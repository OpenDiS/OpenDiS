"""Segment-segment interaction forces, PyDiS and ExaDiS.

Two stored tables of segment pairs, each row holding two segments
(p1,p2) and (p3,p4) with Burgers vectors b12 and b34, and the reference
forces f1..f4 on the four end nodes. Geometry and forces live in
separate files, following the split the rest of this folder uses:

  input_data/<stem>.dat   18 columns, the segment pairs
  ref_data/<stem>.npz     the reference forces, plus mu, nu, a and a
                          digest of the geometry they belong to

Three implementations of the same isotropic non-singular kernel are run
over those tables and compared against the stored forces:

  pydis   compute_segseg_force_list  the compiled ParaDiS library (SBA)
  pydis   python_segseg_force_vec    the python/numpy implementation
  exadis  SegSegIso                  reached through
                                     compute_force_segseglist

so the test pins each implementation to the reference and thereby to
the others. ExaDiS has no SBN1 quadrature variant, so the SBN1 case
below stays pydis-only.

compute_force_segseglist runs FORCE_SEGSEG_ISO over an explicit pair
list and nothing else: no core, self, or PK term, which is what makes
it comparable to the pydis kernels. The table becomes one network of
two segments per row with no node shared between rows, so each node
carries exactly one of f1..f4. Those are open lines, so exadis warns on
construction that the Burgers vector is not conserved; that is expected
for a table of test data. PyDiS cannot hold such a network at all,
which is why the exadis side builds an ExaDisNet directly rather than
going through DisNetManager.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import argparse

import numpy as np
from pydis.calforce.compute_stress_force_analytic_paradis import (
    compute_segseg_force_list)
from pydis.calforce.compute_stress_force_analytic_python import (
    python_segseg_force_vec)
from framework.testing import (array_digest, report_close,
                               quiet_native_output, kokkos_summary)
from pydis.build_info import bitrepro_math, build_description
from segseg_tables import (TABLES, load, exadis_network,
                           MU as mu, NU as nu, A as a)

# The tables themselves, which pairs they cover and how tightly a kernel
# of the numpy lineage may be held over each, live in segseg_tables.py so
# that this test and test_segseg_force_torch.py cannot drift apart over
# which data they run. mu, nu and a come from there too: they are the
# constants the stored references were generated with, not a choice this
# test gets to make.

# Why the near-parallel table gets 1e-6 for the python/numpy kernel,
# the tolerance carried by that table in segseg_tables.py.
#
# Not because the python kernel is wrong. Near parallel this kernel is
# ill-conditioned, and no two implementations of it can be expected to agree
# in float64 more tightly than the kernel's own noise floor, however each is
# written.
#
# The measurement, which is what the tolerance is sized against rather than
# the observed disagreement. Nudging all four endpoint coordinates by one
# ulp, the smallest change float64 can express, 1.1e-16 relative, and taking
# the worst over 40 independent nudges, the compiled kernel's own answer
# moves by 5.8e-08 on this table. The python kernel differs from it by
# 2.6e-08, which is inside the distance the answer moves under a change too
# small to represent. 1e-6 clears the floor by about 17x, so the test does
# not go red when a different compiler or libm moves the reference within
# its own uncertainty; it is still 100x tighter than anything a simulation
# would notice. Sizing the tolerance to the observed 2.6e-08 instead would
# leave under 2x, which is a tolerance that passes here and nowhere else.
#
# The amplification comes from the near-parallel decomposition: the
# correction terms inside SpecialRemoteNodeForce are 1.4 to 10.2 in
# magnitude against a final force of ~0.1, so the result is a cancellation
# of quantities two orders larger than itself.
#
# Confirmed independently by evaluating the python kernel in np.longdouble
# and treating that as the accurate answer: the compiled kernel sits 3.4e-08
# from it at 1-c^2 = 1e-05 and the python kernel 4.0e-08. Both float64
# implementations are off by the same amount, in different directions.
#
# So the 0.0 that the compiled kernel and exadis both report is not evidence
# that either is accurate here, and tol_paradis is not really an accuracy
# claim on this table. exadis' SegSegIso is a transcription of this same
# ParaDiS C code and performs the identical sequence of float64 operations,
# so the two round identically and cancel identically; what the 1e-9 checks
# for them is that neither operation sequence has changed. The python kernel
# came from the Matlab lineage and orders its arithmetic differently, which
# is the whole of the difference. Holding it to 1e-9 would be asserting a
# precision the kernel does not have in this regime.
#
# On the sharp edge this table stays away from: at 1-c^2 = 1e-4 exactly, the
# threshold at which every kernel switches from the general formula to the
# near-parallel one, the switch is discontinuous, the two formulas differing
# by ~7.7e-05 where they meet, and one ulp moves the answer by 2.6e-04. The
# generator refuses to bless a row within 10% of the threshold for that
# reason. A row there would fail for all three implementations at once, on
# whichever machine first put it on the other side of the switch, and would
# look like a force bug. That discontinuity is a property of the ParaDiS
# method, present identically in the C kernel and in exadis.
NEARPAR_NOTE = ("float64 conditioning limit near parallel, not an "
                "implementation defect; see the comment above")

force_nint = 3

# One tolerance per table, carried by the table in segseg_tables.py and
# applied to every implementation alike. tolC is for the SBN1 quadrature
# variant, which is commented out below.
tolC = 1e-5


def test_pydis(pairs, ref_forces, tol):
    """test_pydis: both pydis kernels against the reference

    The two are held to different tolerances on a bitwise-reproducible build.
    The compiled kernel is SegSegForce.c, built there against the portable
    log/atan, so against a reference blessed by any such build it must agree
    exactly and anything else is a reproducibility regression. The python/numpy
    kernel reaches libm through numpy, which the build does not touch, so it
    keeps the table's own tolerance either way.
    """
    p1, p2, p3, p4, b12, b34 = pairs
    tol_c = 0.0 if bitrepro_math() else tol

    # use ParaDiS library (SBA)
    fA = compute_segseg_force_list(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok = report_close("pydis : ParaDiS library (SBA)      ",
                      np.concatenate(fA, axis=1), ref_forces, tol_c)

    # use Python code (translated from Matlab)
    fB = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok &= report_close("pydis : python/numpy implementation",
                       np.concatenate(fB, axis=1), ref_forces, tol)


    # use ParaDiS library (SBN1)
    #quad_points, weights = np.polynomial.legendre.leggauss(force_nint)
    #fC = compute_segseg_force_SBN1_vec(p1, p2, p3, p4, b12, b34,
    #                                   mu, nu, a, quad_points, weights)
    #ok &= report_close("pydis : ParaDiS library (SBN1)",
    #                   np.concatenate(fC, axis=1), ref_forces, tolC)

    return bool(ok)


OUT_DIR = Path('output')          # relative to the working directory


def write_ref():
    """write_ref: regenerate every table's reference from pydis

    Written by the ParaDiS library kernel reached through pydis, the same
    implementation the first comparison in test_pydis checks, so that
    comparison reads 0 against a freshly blessed reference and the other
    two read their difference from it.

    The geometry is read straight from the fixture rather than through
    load(), which needs a reference that agrees with the fixture and so
    cannot be used to make one. Each file carries the digest of the
    geometry it came from, which is what lets load() tell a stale
    reference from a kernel fault.

    Files land in output/ for inspection, not in ref_data/.
    """
    OUT_DIR.mkdir(exist_ok=True)
    for table in TABLES:
        geometry = np.loadtxt(table.dat)
        pairs = tuple(geometry[:, i:i+3] for i in range(0, 18, 3))
        p1, p2, p3, p4, b12, b34 = pairs
        forces = np.concatenate(
            compute_segseg_force_list(p1, p2, p3, p4, b12, b34, mu, nu, a),
            axis=1)
        out_file = OUT_DIR / (table.stem + '.npz')
        np.savez(out_file, forces=forces, mu=mu, nu=nu, a=a,
                 geometry_digest=np.array(array_digest(geometry)),
                 n_pairs=geometry.shape[0])
        print("write_ref: wrote %s, %d pairs x %d force columns"
              % (out_file, forces.shape[0], forces.shape[1]))
    print("")
    print("  inspect them, then bless them:")
    for table in TABLES:
        print("      cp %s/%s.npz ref_data/" % (OUT_DIR, table.stem))


def exadis_forces(pairs):
    """exadis_forces: the pairs through exadis' SegSegIso, as (N,12)

    Separated from the reporting so the caller can silence the C++ chatter
    around this part alone: building the network warns about Burgers
    vector conservation, and Intel's OpenMP runtime prints an "OMP: Info"
    line from inside the force call. Neither indicates a problem and both
    would otherwise land between the result lines.
    """
    import pyexadis
    n = pairs[0].shape[0]
    G, seg_pairs = exadis_network(pairs)
    f = pyexadis.compute_force_segseglist(G.net, mu, nu, a, seg_pairs)
    return np.array(f).reshape(n, 12)


def main():
    print("segment-segment forces from pydis and exadis, against the "
          "stored references")
    print(build_description())
    print("%d tables to check\n" % len(TABLES))

    tables = []
    ok = True
    for n, table in enumerate(TABLES, start=1):
        print("--- table %d of %d: %s" % (n, len(TABLES), table.description))
        pairs, ref_forces = load(table)
        if pairs is None:
            ok = False
            print("")
            continue
        tables.append((table.description, pairs, ref_forces, table.tol))
        # note ok &= rather than ok = : one failing table must fail the test,
        # and must not be cleared by a later table passing
        ok &= test_pydis(pairs, ref_forces, table.tol)
        print("")

    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; the exadis half of this test did "
              "not run")
        return False

    # initialize once, outside the loop: it is global setup, not per table
    with quiet_native_output() as buf:
        pyexadis.initialize()
    for line in kokkos_summary(buf.text):
        print("exadis: %s" % line)
    for n, (description, pairs, ref_forces, tol) in enumerate(tables,
                                                              start=1):
        print("--- table %d of %d: %s" % (n, len(tables), description))
        with quiet_native_output():
            f = exadis_forces(pairs)
        ok &= report_close("exadis: SegSegIso                  ",
                           f, ref_forces, tol)
        print("")
    with quiet_native_output():
        pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-ref', dest='write_ref', action='store_true',
                        default=False,
                        help='regenerate the references from pydis into '
                             'output/ instead of running the comparison')
    args = parser.parse_args()

    if args.write_ref:
        write_ref()
        sys.exit(0)

    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
