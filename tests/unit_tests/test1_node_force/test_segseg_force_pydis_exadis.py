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
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
from pydis.calforce.compute_stress_force_analytic_paradis import (
    compute_segseg_force_list)
from pydis.calforce.compute_stress_force_analytic_python import (
    python_segseg_force_vec)
from framework.testing import report_close, array_digest

# The tables, by stem: each is input_data/<stem>.dat holding the segment
# pairs, 18 columns of p1, p2, p3, p4, b12, b34, and ref_data/<stem>.npz
# holding the reference forces f1..f4 alongside the constants they were
# computed with and a digest of the geometry they belong to.
#
# The two were one 30-column .dat until 2026-08-16. Splitting them puts this
# test on the same footing as the rest of the suite, where input_data/ is what
# a test reads in and ref_data/ is what it checks against, and it makes the
# blessing step explicit: a regenerated fixture lands in input_data/ directly,
# while its reference waits in output/ until copied.
#
# Two tables are run. The first is the original random-geometry set; its
# smallest 1-c^2 is 3.4e-03, so it never approaches the eps = 1e-4 threshold
# at which the kernels switch to the parallel formula. The second was added
# because three bugs lived in that region of the python kernel undetected:
# inconsistent [nonpar] subsetting, a parallel threshold of 1e-6 instead of
# 1e-4, and two correction branches passing lists where arrays are required.
# Regenerate it with make_parallel_ref_data.py.
#
# The third entry is the tolerance for the python/numpy kernel. It is looser
# on the second table for a reason that is about the kernel rather than about
# any of the three implementations; see NEARPAR_NOTE below.
REF_FILES = [
    ('segsep_min_2.5_max_32.5_iso_randombvecs_a0.010',
     'random geometry, 1-c^2 from 3.4e-03 to 0.95', 1e-9),
    ('parallel_and_nearparallel_iso_a0.010',
     'parallel and near-parallel, 1-c^2 from 0 to 1e-02', 1e-6),
]

# Why the second table gets 1e-6 for the python/numpy kernel.
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
                "implementation defect; see the comment on REF_FILES")

# the constants the reference was generated with
mu = 50.0
nu = 0.3
a = 0.01
force_nint = 3

# cell edge for the exadis network, as a multiple of the extent of the
# table. Nothing depends on it: the cell is non-periodic, so there are
# no images, and compute_force_segseglist builds no neighbor list
BOX_FACTOR = 10.0

# tol_paradis applies to the compiled library on every table; the tolerance
# for the python/numpy kernel is per table, carried in REF_FILES. tolC is for
# the SBN1 quadrature variant, which is commented out below.
tol_paradis = 1e-9
tolC = 1e-5
tol_exadis = 1e-9


def load_ref(stem):
    """load_ref: the segment pairs, and the reference forces f1..f4

    The geometry comes from input_data/<stem>.dat, 18 columns of
    p1, p2, p3, p4, b12, b34, and the forces from ref_data/<stem>.npz.
    Returns (pairs, ref_forces), or (None, None) with an explanation if
    either is missing or if the two describe different pairs.
    """
    dat, npz = input_dir / (stem + '.dat'), ref_dir / (stem + '.npz')
    for path, what in ((dat, 'fixture'), (npz, 'reference')):
        if path.exists():
            continue
        print("load_ref: no %s at %s/%s" % (what, path.parent.name, path.name))
        # a freshly generated reference sits in output/ until it is blessed,
        # so say so rather than sending anyone to regenerate what already
        # exists
        waiting = Path('output') / path.name
        if waiting.exists():
            print("          it is waiting to be blessed:")
            print("              cp %s %s/" % (waiting, path.parent.name))
        else:
            print("          regenerate with 'make segseg_force_ref'")
        return None, None

    geometry = np.loadtxt(dat)
    ref = np.load(npz, allow_pickle=False)
    forces = ref['forces']

    # A reference and a fixture blessed at different times can end up
    # describing different segment pairs: regenerating the fixture writes
    # input_data/ in place while the new reference waits in output/ to be
    # copied. Comparing forces to the wrong geometry would look like a
    # kernel bug, so refuse rather than report a failure that is not one.
    digest = array_digest(geometry)
    if str(ref['geometry_digest']) != digest:
        print("load_ref: '%s' and '%s' describe different segment pairs"
              % (dat.name, npz.name))
        print("          fixture   %s, %d pairs"
              % (digest, geometry.shape[0]))
        print("          reference %s, %d pairs"
              % (str(ref['geometry_digest']), int(ref['n_pairs'])))
        print("          the fixture has been regenerated and the new "
              "reference not yet blessed:")
        print("              cp output/%s ref_data/" % npz.name)
        return None, None

    print("load_ref: '%s', %d segment pairs, mu = %g, nu = %g, a = %g"
          % (stem, geometry.shape[0], mu, nu, a))
    pairs = (geometry[:, 0:3], geometry[:, 3:6], geometry[:, 6:9],
             geometry[:, 9:12], geometry[:, 12:15], geometry[:, 15:18])
    return pairs, forces


def test_pydis(pairs, ref_forces, tol_python):
    """test_pydis: both pydis kernels against the reference"""
    p1, p2, p3, p4, b12, b34 = pairs

    # use ParaDiS library (SBA)
    fA = compute_segseg_force_list(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok = report_close("pydis : ParaDiS library (SBA)     ",
                      np.concatenate(fA, axis=1), ref_forces, tol_paradis)

    # use Python code (translated from Matlab)
    fB = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok &= report_close("pydis : python/numpy implementation",
                       np.concatenate(fB, axis=1), ref_forces, tol_python)
    if tol_python > tol_paradis:
        print("        python/numpy is held to %.0e here rather than %.0e:\n"
              "        %s" % (tol_python, tol_paradis, NEARPAR_NOTE))

    # use ParaDiS library (SBN1)
    #quad_points, weights = np.polynomial.legendre.leggauss(force_nint)
    #fC = compute_segseg_force_SBN1_vec(p1, p2, p3, p4, b12, b34,
    #                                   mu, nu, a, quad_points, weights)
    #ok &= report_close("pydis : ParaDiS library (SBN1)",
    #                   np.concatenate(fC, axis=1), ref_forces, tolC)

    return bool(ok)


def build_network(pairs):
    """build_network: one ExaDisNet holding every pair in the table

    Two segments per row, four nodes per row, nothing shared between
    rows. Returns the network and the pair list to evaluate.
    """
    from pyexadis_base import ExaDisNet
    import pyexadis

    p1, p2, p3, p4, b12, b34 = pairs
    n = p1.shape[0]
    rn = np.stack([p1, p2, p3, p4], axis=1).reshape(4*n, 3)
    segs = np.zeros((2*n, 8))
    segs[0::2, 0] = np.arange(n)*4 + 0
    segs[0::2, 1] = np.arange(n)*4 + 1
    segs[0::2, 2:5] = b12
    segs[1::2, 0] = np.arange(n)*4 + 2
    segs[1::2, 1] = np.arange(n)*4 + 3
    segs[1::2, 2:5] = b34

    L = BOX_FACTOR * float(np.max(rn.max(axis=0) - rn.min(axis=0)))
    origin = 0.5*(rn.min(axis=0) + rn.max(axis=0)) - 0.5*L
    cell = pyexadis.Cell(h=L*np.eye(3), origin=origin,
                         is_periodic=[False, False, False])

    return ExaDisNet(cell, rn, segs), [[2*i, 2*i+1] for i in range(n)]


def test_exadis(pairs, ref_forces):
    """test_exadis: the exadis kernel against the same reference"""
    import pyexadis

    n = pairs[0].shape[0]
    G, seg_pairs = build_network(pairs)
    f = pyexadis.compute_force_segseglist(G.net, mu, nu, a, seg_pairs)
    return report_close("exadis: SegSegIso                 ",
                        np.array(f).reshape(n, 12), ref_forces,
                        tol_exadis)


def main():
    print("segment-segment forces from pydis and exadis, against the "
          "stored references")
    print("%d tables to check\n" % len(REF_FILES))

    tables = []
    ok = True
    for n, (stem, description, tol_python) in enumerate(REF_FILES, start=1):
        print("--- table %d of %d: %s" % (n, len(REF_FILES), description))
        pairs, ref_forces = load_ref(stem)
        if pairs is None:
            ok = False
            print("")
            continue
        tables.append((description, pairs, ref_forces))
        # note ok &= rather than ok = : one failing table must fail the test,
        # and must not be cleared by a later table passing
        ok &= test_pydis(pairs, ref_forces, tol_python)
        print("")

    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; the exadis half of this test did "
              "not run")
        return False

    # initialize once, outside the loop: it is global setup, not per table
    pyexadis.initialize()
    for n, (description, pairs, ref_forces) in enumerate(tables, start=1):
        print("--- table %d of %d: %s" % (n, len(tables), description))
        ok &= test_exadis(pairs, ref_forces)
        print("")
    pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
