"""Segment-segment interaction forces on a GPU, through torch.

Exercises pydis.calforce.compute_stress_force_analytic_torch, which
evaluates the isotropic non-singular (SBA) segment-segment kernel on a
torch device.

Check 1 runs by default. Checks 2 to 5 need --extra-test: they are the
slower and more searching ones, and they generate their own geometry
rather than reading the stored tables.

Worth running whenever the kernel, torch or numpy changes, and worth
having CI pass it. The stored tables are 160 fixed rows, and they passed
unaltered across a torch 2.2.2 to 2.4.1 upgrade that moved the sweeps by
two orders of magnitude, so on their own they cannot detect that class of
change. --extra-test costs about 3 seconds.

  1. against the two stored reference tables, the same input_data/ and
     ref_data/ that test_segseg_force_pydis_exadis.py uses, held to the
     same per-table tolerances; both tests take them from
     segseg_tables.py so they cannot come to disagree about the data
  2. against the numpy kernel itself on randomized sweeps that reach
     geometries the stored tables do not
  3. degenerate geometry: exactly parallel, antiparallel, and an empty
     pair set. The two parallel cases must match numpy bit for bit, not
     merely come back finite
  4. the device and dtype contract: tensors in gives tensors back on
     their own device, numpy in gives numpy back
  5. float32 is refused unless explicitly overridden, and the
     measurement that justifies refusing it

WHAT THIS TEST IS AND IS NOT PINNING

Check 2 is the real one. The torch kernel is a port of the numpy kernel,
so the numpy kernel is what it must reproduce, and any divergence between
them is a porting bug by definition. Check 1 is weaker than it looks: the
stored references were generated from the compiled ParaDiS kernel, and
the numpy lineage does not agree with that kernel to rounding once
segments come within 1-c^2 < 1e-6 of parallel, so check 1 has to be run
at the numpy tolerances rather than the compiled ones.

That gap is a property of the numpy kernel and not of torch. Measured
over 20000 random pairs, torch-vs-C and numpy-vs-C are identical to every
digit in every band of 1-c^2, while torch-vs-numpy is five orders
tighter. It is also why CalForce still defaults to force_kernel='batch':
a 200-step run of examples/02_frank_read_src/test_frank_read_src_pydis_elast.py
ends 5.0e-05 away from the compiled kernel's trajectory on
force_kernel='torch' and 5.0e-05 away on force_kernel='vec', the two
being 3.6e-07 from each other.

HOW CHECK 2 IS JUDGED, AND WHY NOT BY A FIXED NUMBER

This kernel amplifies rounding by four orders more in some geometries
than others, so check 2 compares torch against numpy band by band in
1-c^2, and in each band against the kernel's own noise floor measured
right there: how far the compiled kernel's answer moves when an input
coordinate is nudged by one ulp. Below that floor two implementations of
the same algebra are indistinguishable however each is written; above it
something is wrong.

The floor is measured at run time rather than written down, because a
written-down number does not survive a version bump. The first version of
this test used a fixed 1e-6 on the general branch, sized from what torch
2.2.2 produced. torch 2.4.1 reduces in a different order and lands at
3.7e-06, which failed that tolerance while sitting comfortably under the
1.8e-05 floor, i.e. it was no less correct. Measured against the floor
both releases pass with the same margin to spare, and a genuine porting
bug would still be caught: it would miss by orders, not by a factor.

float32 is not merely less accurate here, it is wrong: median relative
error 5.4e-03 where the geometry is benign and 6.3e-03 near parallel,
with worst cases of 1.8e+04 and 6.6e+03 and every value finite, so
nothing signals it. The kernel refuses it unless asked twice, and check 5
pins that.

If torch is not installed this test reports that it did not run and
passes. torch is an optional dependency, needed only for the GPU path,
and a machine without it is not a machine with a broken build.

    python3 test_segseg_force_torch.py
    python3 test_segseg_force_torch.py --extra-test
    python3 test_segseg_force_torch.py --extra-test --device cuda --bench
"""

import argparse
import os
import sys
import time
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]   # exadis: only for --bench
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close
from segseg_tables import TABLES, load, MU as mu, NU as nu, A as a

# The tables come from segseg_tables.py, shared with
# test_segseg_force_pydis_exadis.py, so the two tests provably run over
# the same input_data/ and ref_data/ rather than over two lists that
# happen to match today. mu, nu and a come from there too, being the
# constants the stored references were generated with.
#
# The per-table tolerance for a numpy-lineage kernel is carried by the
# table, and the torch kernel is held to it unchanged: it is a port of the
# numpy kernel, so it inherits both the lineage and the distance from the
# compiled reference that comes with it.

# How far torch may sit from numpy, as a multiple of the kernel's own
# noise floor in the same geometry. The floor is measured at run time
# rather than written down here; see check_against_numpy.
#
# A fixed number does not work. The first version of this test used
# 1e-6 relative on the general branch, sized from what torch 2.2.2
# happened to produce, and torch 2.4.1 fails it at 3.7e-06 while being
# no less correct: the two versions reduce in a different order and the
# kernel amplifies that. Anything tight enough to be meaningful on one
# torch release is a false alarm on the next.
TOL_FLOOR_MULTIPLE = 3.0

# Bands of 1-c^2 to judge separately, because the kernel's conditioning
# varies by four orders across them and one tolerance cannot serve all.
# Splitting at the 1e-4 branch threshold alone is not enough: the rows
# just above it are ill-conditioned and the rows well above it are not.
SWEEP_BANDS = [(1e-2, 2.0, '1-c^2 >= 1e-2   '),
               (1e-4, 1e-2, '1-c^2 1e-4..1e-2'),
               (1e-6, 1e-4, '1-c^2 1e-6..1e-4'),
               (0.0, 1e-6, '1-c^2 0 .. 1e-6 '),
               (-1.0, 0.0, '1-c^2 == 0      ')]

# one-ulp probes used to measure the floor; the worst movement over this
# many independent nudges is taken
FLOOR_PROBES = 4

SWEEP_SEEDS = (1, 2, 3)
SWEEP_N = 5000

BENCH_SIZES = (1000, 10000, 100000)

EXTRA_NOTE = "    more thorough tests can be enabled by --extra-test"

SWEEP_EXPLANATION = """\
    Randomly generated segment pairs, not the stored tables above. 1-c^2 is
    spread over twelve decades so every regime is sampled, including
    geometries the 160 stored rows never reach.

    The torch kernel is a port of the numpy one, so numpy is what it has to
    reproduce; the stored tables cannot tell you that, since they only say
    how far each sits from the compiled ParaDiS kernel.

      differ    max |torch - numpy| over the band, relative to |f| there
      noise x3  the pass bar: how far the compiled ParaDiS kernel's own
                answer moves when one input coordinate is changed by one
                ulp, the smallest step float64 can take, times 3

    Below that bar two implementations of the same algebra are
    indistinguishable however each is written, so it is the bar rather
    than any number chosen in advance."""

# Padded to the same width as the labels in
# test_segseg_force_pydis_exadis.py, so a reader moving between the two
# tests sees one column layout. Only torch is run here; the compiled
# library and the numpy kernel are that test's business, and repeating
# them would mean two places reporting the same number.
LABEL_TORCH = "pydis : torch implementation       "

# what the compiled kernel would be held to, kept only so the block below
# can say when the table's tolerance is looser and why
TOL_PARADIS = 1e-9

NEARPAR_NOTE = ("float64 conditioning limit near parallel, not an "
                "implementation defect; see segseg_tables.py")


def random_pairs(seed, n):
    """random_pairs: segment pairs spanning every regime of 1-c^2

    1-c^2 log-uniform over twelve decades, so the near-parallel branch,
    the general branch and the switch between them are all covered
    densely, with a twentieth of the rows exactly parallel and a sixth of
    them antiparallel. The stored tables cannot do this: they are fixed
    lists blessed against two independent implementations, which is the
    right shape for pinning absolute values and the wrong shape for
    covering a continuum.
    """
    rng = np.random.default_rng(seed)
    onemc2 = 10.0**rng.uniform(-12, 0, n)
    onemc2[:n//20] = 0.0
    sign = np.where(rng.random(n) < 0.15, -1.0, 1.0)
    c = sign * np.sqrt(np.clip(1.0 - onemc2, 0.0, 1.0))

    def unit(m):
        v = rng.normal(size=(m, 3))
        return v / np.linalg.norm(v, axis=1)[:, None]

    t = unit(n)
    perp = unit(n)
    perp = perp - np.einsum('ij,ij->i', perp, t)[:, None] * t
    perp = perp / np.linalg.norm(perp, axis=1)[:, None]
    tp = c[:, None] * t + np.sqrt(np.clip(1.0 - c*c, 0.0, 1.0))[:, None] * perp

    p1 = rng.uniform(-50, 50, (n, 3))
    p2 = p1 + rng.uniform(1, 50, n)[:, None] * t
    p3 = p1 + rng.uniform(-40, 40, (n, 3))
    p4 = p3 + rng.uniform(1, 50, n)[:, None] * tp
    return (p1, p2, p3, p4, unit(n), unit(n))


def onemc2_of(pairs):
    """onemc2_of: 1-c^2 as the kernels compute it, from the endpoints"""
    p1, p2, p3, p4 = pairs[0], pairs[1], pairs[2], pairs[3]
    t = (p2 - p1) / np.linalg.norm(p2 - p1, axis=1)[:, None]
    tp = (p4 - p3) / np.linalg.norm(p4 - p3, axis=1)[:, None]
    return 1.0 - np.einsum('ij,ij->i', t, tp)**2


def check_against_tables(torch_kernel, device, dtype):
    """check_against_tables: torch over the stored tables

    Laid out the way test_segseg_force_pydis_exadis.py lays it out, a
    banner naming the table and then one line per implementation, so the
    two tests read the same way over the same data.
    """
    ok = True
    for n, table in enumerate(TABLES, start=1):
        print("--- table %d of %d: %s"
              % (n, len(TABLES), table.description))
        pairs, ref_forces = load(table)
        if pairs is None:
            ok &= report("%s: could not be loaded" % LABEL_TORCH, False)
            continue

        f = np.concatenate(torch_kernel(*pairs, mu, nu, a,
                                        device=device, dtype=dtype), axis=1)
        ok &= report_close(LABEL_TORCH, f, ref_forces, table.tol_python)

        if table.tol_python > TOL_PARADIS:
            print("        the numpy lineage is held to %.0e here rather "
                  "than %.0e:\n        %s"
                  % (table.tol_python, TOL_PARADIS, NEARPAR_NOTE))
        print("")
    return ok


def noise_floor(pairs, forces):
    """noise_floor: how far the answer moves under a one-ulp input change

    Per row, relative to the local force magnitude. Nudging every endpoint
    coordinate by one ulp, the smallest change float64 can express, and
    taking the worst over FLOOR_PROBES independent nudges. Measured with
    the compiled ParaDiS kernel because it is the implementation nothing
    here is a port of, so it gives an independent read on how much the
    geometry itself amplifies rounding.

    This is the yardstick two implementations of the same algebra have to
    be judged against. Below it they are indistinguishable no matter how
    each is written; above it something is actually wrong.
    """
    from pydis.calforce.compute_stress_force_analytic_paradis import (
        compute_segseg_force_list)
    rng = np.random.default_rng(5)
    worst = np.zeros(len(forces))
    for _ in range(FLOOR_PROBES):
        q = [x.copy() for x in pairs]
        for j in range(4):
            q[j] = np.nextafter(q[j], q[j] + rng.choice([-1.0, 1.0],
                                                        size=q[j].shape))
        f = np.concatenate(compute_segseg_force_list(*q, mu, nu, a), axis=1)
        worst = np.maximum(worst, np.abs(f - forces).max(axis=1))
    return worst


def check_against_numpy(torch_kernel, numpy_kernel, device, dtype):
    """check_against_numpy: the port must reproduce what it is a port of

    Judged per band of 1-c^2 against the noise floor measured in that same
    band, rather than against a written-down number, so the test says the
    same thing on any torch release.
    """
    from pydis.calforce.compute_stress_force_analytic_paradis import (
        compute_segseg_force_list)
    ok = True
    for seed in SWEEP_SEEDS:
        pairs = random_pairs(seed, SWEEP_N)
        # the numpy kernel writes through its inputs on the parallel
        # branch, so it gets copies and the torch kernel gets the originals
        fN = np.concatenate(
            numpy_kernel(*[x.copy() for x in pairs], mu, nu, a), axis=1)
        fT = np.concatenate(
            torch_kernel(*pairs, mu, nu, a, device=device, dtype=dtype),
            axis=1)
        fC = np.concatenate(
            compute_segseg_force_list(*pairs, mu, nu, a), axis=1)

        mag = np.maximum(np.abs(fN).max(axis=1), 1e-12)
        rel = np.abs(fT - fN).max(axis=1) / mag
        floor = noise_floor(pairs, fC) / mag
        onemc2 = onemc2_of(pairs)

        print("  seed %d, %d randomly generated pairs" % (seed, SWEEP_N))
        for lo, hi, label in SWEEP_BANDS:
            k = (onemc2 == 0) if lo < 0 else ((onemc2 > lo) & (onemc2 <= hi))
            if not k.any():
                continue
            got, bar = rel[k].max(), floor[k].max()*TOL_FLOOR_MULTIPLE
            ok &= report(
                "    %s n=%5d  differ %.3e  vs noise x%.0f %.3e"
                % (label, int(k.sum()), got, TOL_FLOOR_MULTIPLE, bar),
                bool(got <= bar))
    return ok


def check_degenerate(torch_kernel, numpy_kernel, device, dtype):
    """check_degenerate: parallel, antiparallel and empty inputs

    Exactly parallel is the one geometry where the near-parallel formula
    is exact rather than an approximation, so both implementations follow
    a short, branch-free path through it and there is nothing for a
    different operation order to disturb. They come out bit for bit
    identical, on torch 2.2.2 and 2.4.1 alike, so that is what is
    asserted. Finiteness alone would pass on any two wrong answers.

    The distance from the compiled ParaDiS kernel is reported rather than
    asserted, since it is that test's business, but it is ~1e-15 relative
    here, which says the shared answer is also the right one.
    """
    ok = True

    t = np.array([[1.0, 0.0, 0.0]])
    for name, tp in (("exactly parallel    ", t),
                     ("exactly antiparallel", -t)):
        p1 = np.zeros((1, 3))
        p2 = p1 + 10.0*t
        p3 = p1 + np.array([[0.0, 5.0, 0.0]])
        p4 = p3 + 10.0*tp
        b = np.array([[0.5, -0.5, 0.7]])
        pairs = (p1, p2, p3, p4, b, b)

        fT = np.concatenate(torch_kernel(*pairs, mu, nu, a, device=device,
                                         dtype=dtype), axis=1)
        fN = np.concatenate(
            numpy_kernel(*[x.copy() for x in pairs], mu, nu, a), axis=1)
        diff = float(np.abs(fT - fN).max())
        ok &= report(
            "torch : %s: |f| = %.3e, vs numpy %.3e, finite %s"
            % (name, np.abs(fN).max(), diff, bool(np.all(np.isfinite(fT)))),
            bool(np.all(np.isfinite(fT))) and diff == 0.0)

    e = np.zeros((0, 3))
    f = torch_kernel(e, e, e, e, e, e, mu, nu, a,
                     device=device, dtype=dtype)
    ok &= report("torch : empty pair set returns four (0,3) arrays",
                 all(x.shape == (0, 3) for x in f))
    return ok


def check_device_contract(torch, torch_kernel, device, dtype):
    """check_device_contract: numpy in gives numpy out, tensors stay put"""
    ok = True
    pairs = random_pairs(11, 64)

    f = torch_kernel(*pairs, mu, nu, a, device=device, dtype=dtype)
    ok &= report("torch : numpy in gives numpy out",
                 all(isinstance(x, np.ndarray) for x in f))

    dev = torch.device(device) if device else f_default_device(torch)
    tens = [torch.as_tensor(x, dtype=dtype or torch.float64, device=dev)
            for x in pairs]
    g = torch_kernel(*tens, mu, nu, a)
    ok &= report("torch : tensors in gives tensors back on the same device",
                 all(torch.is_tensor(x) and x.device == tens[0].device
                     for x in g))

    # the two routes must agree bit for bit: same kernel, same device, the
    # only difference being who did the host-to-device copy
    same = max(float(np.abs(np.asarray(gi.cpu()) - fi).max())
               for gi, fi in zip(g, f))
    ok &= report("torch : both entry points agree, max diff %.3e" % same,
                 same == 0.0)
    return ok


def f_default_device(torch):
    """f_default_device: what the kernel would pick if asked for nothing"""
    from pydis.calforce.compute_stress_force_analytic_torch import (
        resolve_device)
    return resolve_device(None)


def check_float32_refused(torch, torch_kernel, numpy_kernel, device):
    """check_float32_refused: single precision must not be usable by accident

    The kernel is asked for float32 twice: once plainly, which must
    raise, and once with allow_float32 set, which must run. The second
    call is what justifies the first, so its error is printed rather than
    asserted: it is the evidence, not the requirement.

    This matters more than a normal precision caveat because float32 here
    fails silently. It produces no non-finite value and no warning, just
    numbers carrying two or three correct digits in the median and none
    at all in the worst rows.
    """
    ok = True
    pairs = random_pairs(7, SWEEP_N)

    try:
        torch_kernel(*pairs, mu, nu, a, device=device, dtype=torch.float32)
        ok &= report("torch : float32 is refused without an explicit "
                     "override", False)
    except ValueError:
        ok &= report("torch : float32 is refused without an explicit "
                     "override", True)

    fN = np.concatenate(
        numpy_kernel(*[x.copy() for x in pairs], mu, nu, a), axis=1)
    f32 = np.concatenate(
        torch_kernel(*pairs, mu, nu, a, device=device,
                     dtype=torch.float32, allow_float32=True), axis=1)
    ok &= report("torch : float32 runs when the override is passed",
                 f32.shape == fN.shape)

    mag = np.maximum(np.abs(fN).max(axis=1), 1e-12)
    rel = np.abs(f32.astype(np.float64) - fN).max(axis=1) / mag
    general = onemc2_of(pairs) >= 1e-4
    finite = np.isfinite(f32).all()
    print("  and this is why: float32 vs float64, relative error")
    print("      benign geometry (1-c^2 >= 1e-4): median %.2e, worst %.2e"
          % (np.median(rel[general]), rel[general].max()))
    print("      near parallel   (1-c^2 <  1e-4): median %.2e, worst %.2e"
          % (np.median(rel[~general]), rel[~general].max()))
    print("      all values finite: %s, so nothing announces the problem"
          % bool(finite))
    return ok


def bench(torch, torch_kernel, numpy_kernel, device, dtype):
    """bench: wall time per pair for all four implementations

    Read it with the parallelism in mind, which is printed above the
    table. The compiled ParaDiS kernel loops over pairs in C on one
    thread. exadis runs Kokkos over its own thread pool. numpy and torch
    are as parallel as their BLAS and thread settings make them. So this
    compares implementations as a user meets them, not algorithms at
    equal resource, and the C column in particular is a per-core number
    standing next to multi-core ones.

    exadis is timed on the force call only. Building the ExaDisNet is
    done once outside the loop, since the other three are handed arrays
    and are not charged for construction either.
    """
    from pydis.calforce.compute_stress_force_analytic_paradis import (
        compute_segseg_force_list)

    try:
        import pyexadis
        from segseg_tables import exadis_network
        have_exadis = True
    except ImportError:
        have_exadis = False

    print("")
    print("  threads: torch %d, exadis %s, ParaDiS C 1 (serial)"
          % (torch.get_num_threads(),
             os.environ.get('OMP_NUM_THREADS', 'its Kokkos pool')))
    print("  mean of 3 runs after one warm-up, milliseconds")
    print("")
    print("  %8s %12s %12s %12s %12s" % ('pairs', 'ParaDiS C', 'numpy',
                                         'torch', 'exadis'))

    if have_exadis:
        pyexadis.initialize()
    for n in BENCH_SIZES:
        pairs = random_pairs(99, n)

        def timeit(fn):
            fn()
            ts = []
            for _ in range(3):
                t0 = time.perf_counter()
                fn()
                if device and str(device).startswith('cuda'):
                    torch.cuda.synchronize()
                ts.append(time.perf_counter() - t0)
            return 1e3*sum(ts)/len(ts)

        tc = timeit(lambda: compute_segseg_force_list(*pairs, mu, nu, a))
        tn = timeit(lambda: numpy_kernel(*[x.copy() for x in pairs],
                                         mu, nu, a))
        tt = timeit(lambda: torch_kernel(*pairs, mu, nu, a,
                                         device=device, dtype=dtype))
        if have_exadis:
            G, seg_pairs = exadis_network(pairs)
            te = timeit(lambda: pyexadis.compute_force_segseglist(
                G.net, mu, nu, a, seg_pairs))
            cell = "%12.1f" % te
        else:
            cell = "%12s" % 'n/a'
        print("  %8d %12.1f %12.1f %12.1f%s" % (n, tc, tn, tt, cell))
    if have_exadis:
        pyexadis.finalize()


def main(device=None, dtype=None, do_bench=False, extra=False):
    print("segment-segment forces through torch")
    try:
        import torch
    except ImportError:
        print("torch is not installed; this test did not run.")
        print("It covers the optional GPU path only, so a machine without "
              "torch is not a machine with a broken build.")
        print("    pip install torch")
        return True

    from pydis.calforce.compute_stress_force_analytic_torch import (
        torch_segseg_force_vec, device_report)
    from pydis.calforce.compute_stress_force_analytic_python import (
        python_segseg_force_vec)

    if isinstance(dtype, str):
        dtype = {'float32': torch.float32, 'float64': torch.float64}[dtype]
    print(device_report(device))
    print("mu = %g, nu = %g, a = %g\n" % (mu, nu, a))

    print("--- against the stored reference tables")
    ok = check_against_tables(torch_segseg_force_vec, device, dtype)

    if not extra:
        print(EXTRA_NOTE)
        if do_bench:
            bench(torch, torch_segseg_force_vec, python_segseg_force_vec,
                  device, dtype)
        return bool(ok)

    print("\n--- torch against the numpy kernel it is a port of")
    print(SWEEP_EXPLANATION)
    ok &= check_against_numpy(torch_segseg_force_vec,
                              python_segseg_force_vec, device, dtype)

    print("\n--- degenerate geometry: cases the sweeps cannot hit reliably")
    ok &= check_degenerate(torch_segseg_force_vec,
                           python_segseg_force_vec, device, dtype)

    print("\n--- the device and dtype contract")
    ok &= check_device_contract(torch, torch_segseg_force_vec, device, dtype)

    print("\n--- single precision is refused")
    ok &= check_float32_refused(torch, torch_segseg_force_vec,
                                python_segseg_force_vec, device)

    if do_bench:
        bench(torch, torch_segseg_force_vec, python_segseg_force_vec,
              device, dtype)
    return bool(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--device', default=None,
                    help="torch device, e.g. cuda, cuda:1, cpu, mps. "
                         "Defaults to PYDIS_TORCH_DEVICE, then the best "
                         "available.")
    ap.add_argument('--dtype', default=None,
                    choices=['float32', 'float64'],
                    help="working precision, default float64")
    ap.add_argument('--extra-test', dest='extra', action='store_true',
                    default=False,
                    help="also run the checks beyond the stored tables: "
                         "the randomized sweeps against the numpy kernel, "
                         "degenerate geometry, the device and dtype "
                         "contract, and the float32 refusal")
    ap.add_argument('--bench', action='store_true', default=False,
                    help="also time the kernel against numpy and the "
                         "compiled ParaDiS library")
    args = ap.parse_args()

    passed = main(device=args.device, dtype=args.dtype,
                  do_bench=args.bench, extra=args.extra)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("\ntest " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
