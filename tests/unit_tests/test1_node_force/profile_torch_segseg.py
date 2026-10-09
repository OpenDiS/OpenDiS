"""What the torch segment-segment kernel actually launches on a GPU.

    python3 profile_torch_segseg.py --mode eager
    python3 profile_torch_segseg.py --mode dynamic
    python3 profile_torch_segseg.py --mode static --pairs 1e5

One compile mode per run, one pair count per run, tensors already on the
device, a fixed number of timed calls. Reports GPU busy per call, kernel
launches per call, and the split between Triton and eager kernels.

WHY NOT JUST READ THIS OFF bench_segseg_force.py

Because that benchmark times 'torch', 'torch+comp' and 'torch dev' in one
process. Those run different kernels, so a total divided by an assumed
call count describes none of them. Worse, the obvious repair of summing
triton_* for the compiled column and everything else for the eager one is
wrong too: a compiled call still emits plenty of eager kernels. Only
_general_branch goes through torch.compile. The classification above it
(t, tp, c, spindex, the nonzero() gather and the index_copy scatter) and
the whole near-parallel formula body stay eager by design, because they
branch on data; see enable_compile() in
core/pydis/python/pydis/calforce/compute_stress_force_analytic_torch.py.
On a realistic mix 5.3% of pairs take that path, which is a real
contribution rather than a rounding error. The correction recursion, on
the other hand, is compiled, since it re-enters _general_dispatch.

So the split has to come from isolating the mode, not from the names.
The names are still worth reporting inside a mode: they say how much of a
compiled call is actually fused.

THE LAUNCH COUNT IS THE INTERESTING NUMBER

Not the ratio. ExaDiS runs one fused kernel per call, measured at 0.19 ms
for 100k pairs on an A100. Eager torch runs some hundreds of elementwise
kernels over the same data. How far Triton closes that gap is the whole
fusion question, and a per-call launch count answers it directly where a
milliseconds-per-call ratio does not.

FOUR TRAPS, HANDLED HERE

  - A dynamic-shape kernel inherits its tuning from whichever shape
    compiled it. Warming at 1000 pairs and running 100k gave a kernel
    about 1.8x slower on CPU. This warms at the shape it profiles.
  - Inductor caches kernels on disk, outliving the process, so a kernel
    built at some other shape is reused silently. This points
    TORCHINDUCTOR_CACHE_DIR at a fresh directory per run and prints it.
  - First-call autotuning pollutes the window. Warm-up happens before the
    profiler is started, and the window is opened and closed explicitly.
  - This is fp64. An A100 does 9.7 TFLOP/s in double precision and
    Triton's fp64 elementwise codegen is not especially tuned, so a gap
    to hand-written CUDA is plausible on those grounds alone, separately
    from anything about fusion. Read a Triton-vs-ExaDiS gap with that in
    mind before attributing it to global-memory round trips.

READING THE RESULT AGAINST EXADIS

ExaDiS at 100k pairs is ~1e8 flops in 190 us, some 20x off fp64 peak,
which is sane for a register-heavy kernel. If the Triton kernel lands
near that, fusion has closed the gap and the "one fused kernel against
hundreds of round trips" story applies only to eager torch. If it stays
in the milliseconds, ExaDiS is genuinely ahead on kernel time, and which
of fp64 codegen or global-memory traffic explains it is the next
question.

Keep that separate from call time, which has the opposite answer. Around
its 0.19 ms kernel ExaDiS spends 81.6 ms per pyexadis call, against about
5.1 ms for torch on device tensors. Whether kernel time or call time is
the one that matters depends entirely on whether the caller is python per
timestep or C++ inside a loop.

UNDER NSYS

The timed window is bracketed by cudaProfilerApi, so nsys can be told to
record only that:

    nsys profile -t cuda --capture-range=cudaProfilerApi \\
      --capture-range-end=stop -f true -o torch_dynamic \\
      python3 profile_torch_segseg.py --mode dynamic
    nsys stats --report cuda_gpu_kern_sum torch_dynamic.nsys-rep

which needs no division by an assumed call count and no filtering of
warm-up. The numbers this script prints come from torch.profiler and
should agree; nsys is the authority if they do not.
"""

import argparse
import os
import sys
import tempfile


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split('\n')[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', default='dynamic',
                    choices=['eager', 'dynamic', 'static'],
                    help="one per run. eager is no torch.compile at all; "
                         "dynamic and static are the two compile modes.")
    ap.add_argument('--pairs', default='1e5',
                    help="pair count, default 1e5")
    ap.add_argument('--geometry', default='realistic',
                    help="general, realistic or sweep. Default realistic, "
                         "which is the near-parallel mix a real run gives.")
    ap.add_argument('--calls', type=int, default=10,
                    help="timed calls inside the profiled window, default 10")
    ap.add_argument('--warmup', type=int, default=5,
                    help="calls before the window, at the same shape, to "
                         "absorb compilation and autotuning. Default 5.")
    ap.add_argument('--device', default='cuda',
                    help="torch device, default cuda")
    ap.add_argument('--keep-cache', action='store_true', default=False,
                    help="reuse whatever TORCHINDUCTOR_CACHE_DIR is set to "
                         "instead of a fresh one. Off by default: a kernel "
                         "compiled at another shape would be reused silently.")
    args = ap.parse_args(argv)

    n = int(float(args.pairs))

    # Inductor shells out to $CC to build the triton launcher, and a vendor
    # compiler rejects the flags it passes: nvc from the hpc_sdk module,
    # which is loaded by default on mc3, fails on -Wno-psabi. The module
    # catches that and falls back to eager, which turns a compiled run into
    # an eager one wearing its label, so say so before anything is measured
    # rather than leaving it to the fallback notice mid-run.
    if args.mode != 'eager':
        cc = os.environ.get('CC', '')
        if cc and os.path.basename(cc).split('-')[0] not in ('gcc', 'cc',
                                                             'clang'):
            print("warning: CC=%s is not a gcc-like compiler. Inductor "
                  "builds the triton launcher with it and vendor compilers "
                  "reject its flags, in which case this run silently "
                  "measures eager. Re-run as CC=gcc %s ..."
                  % (cc, os.path.basename(sys.argv[0])))

    # before torch is imported, so inductor cannot have read the old value
    if not args.keep_cache:
        cache = tempfile.mkdtemp(prefix='inductor-%s-%d-' % (args.mode, n))
        os.environ['TORCHINDUCTOR_CACHE_DIR'] = cache
    else:
        cache = os.environ.get('TORCHINDUCTOR_CACHE_DIR', '<torch default>')

    # bench_segseg_force puts pydis on the path and owns the pair geometry,
    # so the two scripts cannot drift apart on what they are measuring
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from bench_segseg_force import make_pairs, near_parallel_fraction, \
        MU, NU, A

    import torch
    import pydis.calforce.compute_stress_force_analytic_torch as K

    dev = K.resolve_device(args.device)
    if dev.type != 'cuda':
        print("this profiles GPU kernels; --device resolved to %s" % dev)
        return 1

    K.enable_compile(None if args.mode == 'eager' else args.mode)

    print("torch %s on %s" % (torch.__version__, torch.cuda.get_device_name(dev)))
    print("mode %s, %d pairs, geometry %s, %d warm-up + %d timed calls"
          % (args.mode, n, args.geometry, args.warmup, args.calls))
    print("inductor cache %s" % cache)

    pairs = make_pairs(n, args.geometry)
    print("%.1f%% of pairs take the near-parallel branch, which stays eager "
          "in every mode" % (100.0*near_parallel_fraction(pairs)))

    ts = [torch.as_tensor(x, dtype=torch.float64, device=dev) for x in pairs]
    call = lambda: K.torch_segseg_force_vec(*ts, MU, NU, A)

    # compilation and autotuning happen here, at the shape being profiled,
    # so neither lands in the window below
    for _ in range(args.warmup):
        call()
    torch.cuda.synchronize()
    if args.mode != 'eager' and K._compile_mode is None:
        print("torch.compile fell back to eager; the module said why above")

    report(*measure(torch, call, args.calls), calls=args.calls)
    return 0


def measure(torch, call, calls):
    """measure: CUDA events for `calls` calls, and their wall time"""
    import time
    from torch.profiler import profile, ProfilerActivity

    torch.cuda.synchronize()
    torch.cuda.profiler.start()          # the nsys capture range opens here
    t0 = time.perf_counter()
    with profile(activities=[ProfilerActivity.CUDA]) as prof:
        for _ in range(calls):
            call()
        torch.cuda.synchronize()
    wall = time.perf_counter() - t0
    torch.cuda.profiler.stop()
    return prof, wall


def _device_time(evt):
    """_device_time: self GPU time in us, across torch profiler versions"""
    for attr in ('self_device_time_total', 'self_cuda_time_total'):
        v = getattr(evt, attr, None)
        if v is not None:
            return v
    return 0.0


def report(prof, wall, calls):
    """report: GPU busy, launches, and how much of it Triton fused"""
    rows = [(e.key, e.count, _device_time(e)) for e in prof.key_averages()]
    rows = [r for r in rows if r[2] > 0.0]
    # memory traffic is not a kernel launch and should not be counted as one
    kern = [r for r in rows if not r[0].startswith(('Memcpy', 'Memset'))]
    mem = [r for r in rows if r[0].startswith(('Memcpy', 'Memset'))]

    triton = [r for r in kern if r[0].startswith('triton_')]
    eager = [r for r in kern if not r[0].startswith('triton_')]

    def tot(rs):
        return sum(c for _, c, _ in rs), sum(t for _, _, t in rs)

    nk, tk = tot(kern)
    nt, tt = tot(triton)
    ne, te = tot(eager)
    nm, tm = tot(mem)

    print("")
    print("  per call, mean of %d" % calls)
    print("    wall              %8.2f ms" % (1e3*wall/calls))
    print("    GPU busy          %8.2f ms" % (tk/1e3/calls))
    print("    kernel launches   %8.1f" % (nk/calls))
    print("      triton_*        %8.1f  (%6.2f ms)" % (nt/calls, tt/1e3/calls))
    print("      eager           %8.1f  (%6.2f ms)" % (ne/calls, te/1e3/calls))
    if nm:
        print("    memcpy/memset     %8.1f  (%6.2f ms)"
              % (nm/calls, tm/1e3/calls))

    print("")
    print("  ten costliest kernels")
    print("    %10s %8s %10s  %s" % ('total ms', 'count', 'per call', 'name'))
    for name, count, us in sorted(kern, key=lambda r: -r[2])[:10]:
        print("    %10.3f %8d %10.1f  %s"
              % (us/1e3, count, count/calls, name[:78]))


if __name__ == '__main__':
    sys.exit(main())
