"""Segment-segment interaction forces, computed by ExaDiS.

The counterpart of test_segseg_force_pydis.py, on the same stored table
and the same reference forces. It corresponds to Test A there, the
analytic SBA kernel; ExaDiS has no SBN1 variant, so Test C stays
pydis-only.

The kernel is reached through pyexadis.compute_force_segseglist, which
runs FORCE_SEGSEG_ISO over an explicit pair list and nothing else: no
core, self, or PK term. The table becomes one network of two segments
per row with no node shared between rows, so each node carries exactly
one of f1..f4. Those are open lines, so ExaDiS warns on construction
that the Burgers vector is not conserved; that is expected for a table
of test data. PyDiS cannot hold such a network at all, which is why
this builds an ExaDisNet directly rather than via DisNetManager.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pyexadis_paths = [str(opendis_root / p) for p in
                  ['python', 'core/exadis/python']]
[sys.path.append(p) for p in pyexadis_paths if not p in sys.path]
ref_dir = Path(__file__).resolve().parent / 'ref_data'

import numpy as np

import pyexadis
from pyexadis_base import ExaDisNet
from framework.testing import report_close

# the segment pairs and their reference forces: columns 0:18 are
# p1, p2, p3, p4, b12, b34 and 18:30 are f1..f4
REF_FILE = 'segsep_min_2.5_max_32.5_iso_randombvecs_a0.010.dat'

# the constants the reference was generated with
mu = 50.0
nu = 0.3
a = 0.01

# cell edge, as a multiple of the extent of the segment table. Nothing
# here depends on it: the cell is non-periodic, so there are no images,
# and compute_force_segseglist builds no neighbor list
BOX_FACTOR = 10.0

tolA = 1e-9


def build_network(segseg_data):
    """build_network: one ExaDisNet holding every pair in the table

    Two segments per row, four nodes per row, nothing shared between
    rows. Returns the network and the pair list to evaluate.
    """
    n = segseg_data.shape[0]
    rn = segseg_data[:, 0:12].reshape(4*n, 3)
    segs = np.zeros((2*n, 8))
    segs[0::2, 0] = np.arange(n)*4 + 0
    segs[0::2, 1] = np.arange(n)*4 + 1
    segs[0::2, 2:5] = segseg_data[:, 12:15]
    segs[1::2, 0] = np.arange(n)*4 + 2
    segs[1::2, 1] = np.arange(n)*4 + 3
    segs[1::2, 2:5] = segseg_data[:, 15:18]

    L = BOX_FACTOR * float(np.max(rn.max(axis=0) - rn.min(axis=0)))
    origin = 0.5*(rn.min(axis=0) + rn.max(axis=0)) - 0.5*L
    cell = pyexadis.Cell(h=L*np.eye(3), origin=origin,
                         is_periodic=[False, False, False])

    seg_pairs = [[2*i, 2*i+1] for i in range(n)]
    return ExaDisNet(cell, rn, segs), seg_pairs


def main():
    print("segment-segment forces from exadis, against the stored "
          "reference '%s'" % REF_FILE)
    segseg_data = np.loadtxt(ref_dir / REF_FILE)
    n = segseg_data.shape[0]
    print("%d segment pairs, mu = %g, nu = %g, a = %g"
          % (n, mu, nu, a))

    f1234_ref = segseg_data[:, 18:]

    # Test A: use exadis SegSegIso, the counterpart of the SBA kernel
    G, seg_pairs = build_network(segseg_data)
    fA = pyexadis.compute_force_segseglist(G.net, mu, nu, a, seg_pairs)
    okA = report_close("TestA, exadis SegSegIso",
                       np.array(fA).reshape(n, 12), f1234_ref, tolA)

    return bool(okA)


if __name__ == "__main__":
    pyexadis.initialize()
    passed = main()
    pyexadis.finalize()

    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
