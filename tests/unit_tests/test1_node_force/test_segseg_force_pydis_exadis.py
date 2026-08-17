"""Segment-segment interaction forces, PyDiS and ExaDiS.

A stored table of 100 segment pairs, each row holding two segments
(p1,p2) and (p3,p4) with Burgers vectors b12 and b34, followed by the
reference forces f1..f4 on the four end nodes. The pairs span
separations from 2.5 to 32.5 with random Burgers vectors, at a = 0.01.

Three implementations of the same isotropic non-singular kernel are run
over that table and compared against the stored forces:

  pydis   compute_segseg_force_vec   the compiled ParaDiS library (SBA)
  pydis   python_segseg_force_vec    the python translation from Matlab
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
ref_dir = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
from pydis.calforce.compute_stress_force_analytic_paradis import (
    compute_segseg_force_vec)
from pydis.calforce.compute_stress_force_analytic_python import (
    python_segseg_force_vec)
from framework.testing import report_close

# the segment pairs and their reference forces: columns 0:18 are
# p1, p2, p3, p4, b12, b34 and 18:30 are f1..f4
REF_FILE = 'segsep_min_2.5_max_32.5_iso_randombvecs_a0.010.dat'

# the constants the reference was generated with
mu = 50.0
nu = 0.3
a = 0.01
force_nint = 3

# cell edge for the exadis network, as a multiple of the extent of the
# table. Nothing depends on it: the cell is non-periodic, so there are
# no images, and compute_force_segseglist builds no neighbor list
BOX_FACTOR = 10.0

tolA, tolB, tolC = 1e-9, 1e-9, 1e-5
tol_exadis = 1e-9


def load_ref(ref_file=REF_FILE):
    """load_ref: the segment pairs and the reference forces f1..f4"""
    data = np.loadtxt(ref_dir / ref_file)
    print("load_ref: '%s', %d segment pairs, mu = %g, nu = %g, a = %g"
          % (ref_file, data.shape[0], mu, nu, a))
    pairs = (data[:, 0:3], data[:, 3:6], data[:, 6:9], data[:, 9:12],
             data[:, 12:15], data[:, 15:18])
    return pairs, data[:, 18:30]


def test_pydis(pairs, ref_forces):
    """test_pydis: both pydis kernels against the reference"""
    p1, p2, p3, p4, b12, b34 = pairs

    # use ParaDiS library (SBA)
    fA = compute_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok = report_close("pydis : ParaDiS library (SBA)",
                      np.concatenate(fA, axis=1), ref_forces, tolA)

    # use Python code (translated from Matlab)
    fB = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    ok &= report_close("pydis : python translation  ",
                       np.concatenate(fB, axis=1), ref_forces, tolB)

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
    return report_close("exadis: SegSegIso           ",
                        np.array(f).reshape(n, 12), ref_forces,
                        tol_exadis)


def main():
    print("segment-segment forces from pydis and exadis, against the "
          "stored reference")
    pairs, ref_forces = load_ref()

    ok = test_pydis(pairs, ref_forces)

    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; the exadis half of this test did "
              "not run")
        return False

    pyexadis.initialize()
    ok &= test_exadis(pairs, ref_forces)
    pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
