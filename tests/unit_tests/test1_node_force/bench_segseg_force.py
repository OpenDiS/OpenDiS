"""Time the four segment-segment force kernels against each other.

    python3 bench_segseg_force.py                     # everything it can find
    python3 bench_segseg_force.py --device cuda       # on a GPU node
    python3 bench_segseg_force.py --geometry all --sizes 1e4,1e5,1e6
    python3 bench_segseg_force.py --compile static   # fixed-shape ceiling
    python3 bench_segseg_force.py --compile off         # eager only

Four implementations, all computing the same isotropic non-singular
segment-segment interaction:

    ParaDiS C   compute_segseg_force_list, the compiled kernel looping
                over pairs in C. Serial.
    numpy       python_segseg_force_vec, vectorized over pairs.
    torch       torch_segseg_force_vec, the same on a torch device.
    ExaDiS      compute_force_segseglist, C++ with Kokkos threading.

Anything not importable is reported and skipped, so this runs on a node
with no ExaDiS build or no torch without special-casing.

GEOMETRY MATTERS MORE THAN SIZE HERE, WHICH IS EASY TO GET WRONG

The kernel takes one of two branches depending on how close to parallel a
pair is, and they cost very different amounts. Benchmarking a mix that
does not resemble a real simulation gives numbers that are not wrong so
much as answers to a different question. Three mixes are offered:

    general     no near-parallel pairs at all, the cheap branch only
    realistic   8% near-parallel, which is what
                examples/02_frank_read_src/test_frank_read_src_pydis_elast.py
                actually produces, measured over 87522 selected pairs
    sweep       1-c^2 log-uniform over twelve decades, so two thirds are
                near-parallel. This is the accuracy-testing distribution
                from test_segseg_force_torch.py and it is the wrong one
                for timing; it is kept because the contrast is
                informative.

The default is `realistic`. On this machine torch is 6.1x the compiled
kernel under `realistic` and 2.0x under `sweep`, purely from the mix, so
quoting one number without saying which mix produced it is meaningless.

THREADS

Pass --threads N. ExaDiS sizes its Kokkos pool from OMP_NUM_THREADS and
torch has its own default, and on a 4-core mac those came out 8 and 4, a
2x handicap invisible in the output. The compiled ParaDiS kernel is
serial regardless, so its column is a per-core number standing next to
multi-core ones however the others are pinned.

TORCH.COMPILE

A 'torch+comp' column is reported by default, --compile dynamic. It fuses
the general branch, which otherwise materialises around a hundred
full-length temporaries and is where most of the time goes. Measured at
100k pairs, 4 threads: 238 ms to 90 on general geometry, 266 to 111 on a
realistic mix. Fusion moves the answer by ~1.5e-12 relative, orders below
the kernel's own noise floor.

Read the column knowing what it excludes: timing is taken after a warm-up
call, so compilation is not charged to it. Compiling costs about 23 s the
first time in either mode; 'static' then pays roughly 13 s again for every
new pair count, while 'dynamic' serves any count from the one kernel. See
DYNAMIC_NOTE below for how the two compare once warm, and for the two
measurement traps worth controlling for.

WHAT THE GPU COLUMNS MEAN

With --device cuda there are two extra torch timings. "torch" takes numpy
arrays and gives numpy back, so it pays a host-to-device copy in and a
device-to-host copy out, which is what CalForce(force_kernel='torch')
does today. "torch dev" is handed tensors already resident and returns
tensors, which is the ceiling that keeping a whole timestep on the device
would approach.

"torch dev" is compiled or not to match "torch+comp", so the difference
between those two columns is purely the transfer, and on a large problem
that difference is what decides whether moving more of the step onto the
device is worth doing.
"""

import argparse
import os
import sys
import time
from pathlib import Path

# Quieten the OpenMP runtime before anything loads it. Intel's libiomp5,
# which a conda numpy or torch may pull in, prints
#   OMP: Info #277: omp_get_nested routine deprecated ...
# from inside ExaDiS on every call, which lands between the rows of the
# timing table. These have to be set before the runtime initialises, so
# before numpy and torch are imported, not in main().
os.environ.setdefault('KMP_WARNINGS', '0')
os.environ.setdefault('KMP_AFFINITY', 'noverbose')

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import quiet_native_output, kokkos_summary

MU, NU, A = 50.0, 0.3, 0.01


def make_pairs(n, geometry='realistic', seed=99):
    """make_pairs: n segment pairs with a chosen near-parallel fraction"""
    rng = np.random.default_rng(seed)
    if geometry == 'general':
        onemc2 = 10.0**rng.uniform(-2, 0, n)
    elif geometry == 'realistic':
        onemc2 = np.where(rng.random(n) < REALISTIC_NEAR_PARALLEL,
                          0.0, 10.0**rng.uniform(-2, 0, n))
    elif geometry == 'sweep':
        onemc2 = 10.0**rng.uniform(-12, 0, n)
    else:
        raise ValueError("geometry must be general, realistic or sweep")

    sign = np.where(rng.random(n) < 0.15, -1.0, 1.0)
    c = sign * np.sqrt(np.clip(1.0 - onemc2, 0.0, 1.0))

    def unit(m):
        v = rng.normal(size=(m, 3))
        return v / np.linalg.norm(v, axis=1)[:, None]

    t = unit(n)
    perp = unit(n)
    perp = perp - np.einsum('ij,ij->i', perp, t)[:, None] * t
    perp = perp / np.linalg.norm(perp, axis=1)[:, None]
    tp = c[:, None]*t + np.sqrt(np.clip(1.0 - c*c, 0.0, 1.0))[:, None]*perp

    p1 = rng.uniform(-50, 50, (n, 3))
    p2 = p1 + rng.uniform(1, 50, n)[:, None] * t
    p3 = p1 + rng.uniform(-40, 40, (n, 3))
    p4 = p3 + rng.uniform(1, 50, n)[:, None] * tp
    return (p1, p2, p3, p4, unit(n), unit(n))


def near_parallel_fraction(pairs):
    """near_parallel_fraction: how many pairs take the near-parallel branch"""
    p1, p2, p3, p4 = pairs[0], pairs[1], pairs[2], pairs[3]
    t = (p2 - p1) / np.linalg.norm(p2 - p1, axis=1)[:, None]
    tp = (p4 - p3) / np.linalg.norm(p4 - p3, axis=1)[:, None]
    return float((1.0 - np.einsum('ij,ij->i', t, tp)**2 < 1e-4).mean())


def timed(fn, reps, sync=None):
    """timed: mean milliseconds over `reps` runs, after one warm-up"""
    fn()
    if sync:
        sync()
    out = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        if sync:
            sync()
        out.append(time.perf_counter() - t0)
    return 1e3 * sum(out) / len(out)


def collect_kernels(device, dtype, want, compile_mode=None):
    """collect_kernels: the implementations available here, in report order

    Returns a list of (name, factory) where factory(pairs) gives a
    zero-argument callable to time, plus a list of skip messages.
    """
    kernels, skipped = [], []

    if 'c' in want:
        try:
            from pydis.calforce.compute_stress_force_analytic_paradis import (
                compute_segseg_force_list)
            kernels.append(('ParaDiS C', lambda pr: (
                lambda: compute_segseg_force_list(*pr, MU, NU, A))))
        except Exception as exc:
            skipped.append("ParaDiS C: %s" % exc)

    if 'numpy' in want:
        try:
            from pydis.calforce.compute_stress_force_analytic_python import (
                python_segseg_force_vec)
            # the numpy kernel writes through its inputs on the parallel
            # branch, so it is given copies; the copy is charged to it,
            # which is what a caller would pay too
            kernels.append(('numpy', lambda pr: (
                lambda: python_segseg_force_vec(
                    *[x.copy() for x in pr], MU, NU, A))))
        except Exception as exc:
            skipped.append("numpy: %s" % exc)

    if 'torch' in want:
        try:
            import torch
            from pydis.calforce.compute_stress_force_analytic_torch import (
                torch_segseg_force_vec, resolve_device)
            dev = resolve_device(device)
            dt = {'float32': torch.float32,
                  'float64': torch.float64}.get(dtype, torch.float64)
            import pydis.calforce.compute_stress_force_analytic_torch as K

            def torch_runner(pr, mode, dev=dev, dt=dt):
                # the module decides eagerly or compiled from a global, so
                # each runner sets it before calling. It is not reset through
                # enable_compile(), which would discard the compiled kernel
                # and pay for compilation again on every timed call.
                def run():
                    K._compile_mode = mode
                    return torch_segseg_force_vec(
                        *pr, MU, NU, A, device=dev, dtype=dt,
                        allow_float32=(dt is torch.float32))
                return run

            kernels.append(('torch', lambda pr: torch_runner(pr, None)))
            if compile_mode:
                # probe before offering the column: torch.compile raises at
                # wrap time on an unsupported python, and this runs by
                # default now, so an old torch must degrade to a note rather
                # than take the whole benchmark down
                try:
                    torch.compile(lambda x: x)
                    kernels.append(('torch+comp',
                                    lambda pr: torch_runner(pr, compile_mode)))
                except Exception as exc:
                    skipped.append("torch.compile: %s" % exc)
            if dev.type != 'cpu':
                # same kernel, but handed tensors already resident, so the
                # host/device copies are excluded; see the module docstring
                def dev_factory(pr, dev=dev, dt=dt, mode=compile_mode,
                                k=torch_segseg_force_vec):
                    ts = [torch.as_tensor(np.ascontiguousarray(x), dtype=dt,
                                          device=dev) for x in pr]

                    def run():
                        # set the mode explicitly. Without this the column
                        # inherits whatever the previously timed one left
                        # behind, so its meaning depends on execution order
                        K._compile_mode = mode
                        return k(*ts, MU, NU, A,
                                 allow_float32=(dt is torch.float32))
                    return run
                kernels.append(('torch dev', dev_factory))
        except Exception as exc:
            skipped.append("torch: %s" % exc)

    if 'exadis' in want:
        try:
            import pyexadis
            from segseg_tables import exadis_network

            def exadis_factory(pr):
                # the network is built outside the timed call: the other
                # kernels are handed arrays and are not charged for
                # construction either
                G, seg_pairs = exadis_network(pr)
                return lambda: pyexadis.compute_force_segseglist(
                    G.net, MU, NU, A, seg_pairs)
            kernels.append(('ExaDiS', exadis_factory))
        except Exception as exc:
            skipped.append("ExaDiS: %s" % exc)

    return kernels, skipped


def describe_environment(device):
    """describe_environment: what does the work, and with how many threads"""
    print("numpy %s" % np.__version__)
    try:
        import torch
        from pydis.calforce.compute_stress_force_analytic_torch import (
            device_report)
        print("%s, %d threads" % (device_report(device),
                                  torch.get_num_threads()))
    except ImportError:
        print("torch not available, so its columns are skipped. Install it "
              "into the environment that\n  can already import pyexadis, "
              "matching the CUDA the node reports:\n"
              "      pip install torch --index-url "
              "https://download.pytorch.org/whl/cu124\n"
              "  Prefer these wheels over conda-forge pytorch, whose MKL "
              "builds pull a second\n  OpenMP runtime alongside the one "
              "ExaDiS links.")
    except Exception as exc:
        print("torch present but unusable: %s" % exc)
    omp = os.environ.get('OMP_NUM_THREADS', 'unset')
    print("OMP_NUM_THREADS = %s; ExaDiS sizes its Kokkos pool from it, and "
          "prints the pool below" % omp)
    if omp == 'unset':
        print("  no thread count pinned. Kokkos and torch choose their own "
              "defaults and they need not match, which makes the comparison "
              "unfair in whichever direction. Pass --threads N.")
    print("the compiled ParaDiS kernel is serial whatever this says, so read "
          "its column as a per-core number")


def run(sizes, geometry, reps, device, dtype, want, compile_mode=None,
        verbose_init=False):
    """run: one table per geometry mix"""
    kernels, skipped = collect_kernels(device, dtype, want, compile_mode)
    for msg in skipped:
        print("  skipped, %s" % msg)
    if not kernels:
        print("nothing to time")
        return False

    names = [n for n, _ in kernels]
    for mix in geometry:
        print("")
        pairs0 = make_pairs(min(sizes), mix)
        print("  geometry '%s': %.1f%% of pairs take the near-parallel branch"
              % (mix, 100.0*near_parallel_fraction(pairs0)))
        print("  milliseconds, mean of %d runs after one warm-up" % reps)
        print("")
        print("  %9s" % 'pairs' + ''.join("%12s" % n for n in names))

        sync = None
        try:
            import torch
            if device and str(device).startswith('cuda'):
                sync = torch.cuda.synchronize
        except Exception:
            pass

        # ExaDiS prints a "Burgers vector is not conserved" warning when it
        # builds a network of open segment pairs, which is expected for a
        # table of test data. Building every runnable first keeps those
        # warnings above the table instead of between its rows.
        runnables = [(n, make_pairs(n, mix)) for n in sizes]
        with quiet_native_output(not verbose_init):
            runnables = [(n, [factory(pr) for _, factory in kernels])
                         for n, pr in runnables]

        # Compile on the largest size before timing anything. A
        # dynamic-shape kernel is generated once and reused for every
        # later shape, and it inherits its tuning from whichever shape
        # compiled it: warmed at 1000 pairs it then runs 100k in 209 ms,
        # warmed at 100k it runs the same work in 114 ms. Timing sizes in
        # ascending order therefore measures a kernel tuned for the
        # smallest one, which is an artifact of the benchmark rather than
        # anything about torch.
        if compile_mode == 'dynamic':
            for fn in runnables[-1][1]:
                fn()

        rows = []
        for n, fns in runnables:
            row = [timed(fn, reps, sync) for fn in fns]
            rows.append((n, row))
            print("  %9d" % n + ''.join("%12.1f" % v for v in row))

        print("")
        print("  %9s" % 'us/pair' + ''.join("%12s" % n for n in names))
        for n, row in rows:
            print("  %9d" % n + ''.join("%12.2f" % (1e3*v/n) for v in row))

        if compile_mode == 'dynamic' and 'torch+comp' in names:
            print("")
            print(DYNAMIC_NOTE)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split('\n')[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--device', default=None,
                    help="torch device: cuda, cuda:1, cpu, mps. Defaults to "
                         "PYDIS_TORCH_DEVICE, then the best available.")
    ap.add_argument('--dtype', default='float64',
                    choices=['float32', 'float64'],
                    help="working precision for torch, default float64. "
                         "float32 is not accurate enough for production; "
                         "see compute_stress_force_analytic_torch.py")
    ap.add_argument('--sizes', default=None,
                    help="comma separated pair counts, e.g. 1e4,1e5,1e6")
    ap.add_argument('--geometry', default='realistic',
                    help="general, realistic, sweep, or all. Default "
                         "realistic, which is the mix a real run produces.")
    ap.add_argument('--compile', default='dynamic',
                    choices=['off', 'dynamic', 'static'],
                    help="time the torch kernel through torch.compile too, "
                         "as an extra 'torch+comp' column. Default dynamic, "
                         "which compiles once for symbolic shapes and is what "
                         "a real run would get, its pair count changing every "
                         "call. 'static' compiles per pair count: faster once "
                         "warm, but it recompiles for every new size. 'off' "
                         "to skip. Needs torch >= 2.4 on python 3.12, >= 2.6 "
                         "on 3.13; skipped with a note otherwise.")
    ap.add_argument('--verbose-init', action='store_true', default=False,
                    help="show the Kokkos startup banner and the Burgers "
                         "vector warnings from building test networks, "
                         "instead of the two lines summarising them")
    ap.add_argument('--threads', type=int, default=None,
                    help="pin every threaded backend to this many threads, "
                         "so the comparison is like for like. Sets "
                         "OMP_NUM_THREADS, which Kokkos honours, and "
                         "torch.set_num_threads. Without it each backend "
                         "picks its own default and they need not agree: "
                         "on a 4-core mac Kokkos took 8 and torch 4.")
    ap.add_argument('--reps', type=int, default=3,
                    help="timed repetitions per point, default 3")
    ap.add_argument('--skip', default='',
                    help="comma separated kernels to leave out: c, numpy, "
                         "torch, exadis. numpy is the slow one.")
    args = ap.parse_args(argv)

    sizes = DEFAULT_SIZES
    if args.sizes:
        sizes = tuple(int(float(s)) for s in args.sizes.split(','))
    mixes = (['general', 'realistic', 'sweep'] if args.geometry == 'all'
             else [args.geometry])
    want = {'c', 'numpy', 'torch', 'exadis'}
    want -= {s.strip() for s in args.skip.split(',') if s.strip()}

    # must happen before pyexadis.initialize(), which is when Kokkos reads
    # OMP_NUM_THREADS and sizes its pool
    if args.threads:
        os.environ['OMP_NUM_THREADS'] = str(args.threads)
        try:
            import torch
            torch.set_num_threads(args.threads)
        except ImportError:
            pass

    print("segment-segment force kernels, wall time")
    describe_environment(args.device)

    started = False
    try:
        import pyexadis
        if 'exadis' in want:
            with quiet_native_output(not args.verbose_init) as buf:
                pyexadis.initialize()
                started = True
            for line in kokkos_summary(buf.text):
                print("  ExaDiS: %s" % line)
    except Exception:
        pass
    try:
        ok = run(sizes, mixes, args.reps, args.device, args.dtype, want,
                 None if args.compile == 'off' else args.compile,
                 args.verbose_init)
    finally:
        if started:
            pyexadis.finalize()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
