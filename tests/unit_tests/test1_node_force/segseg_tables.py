"""The stored segment-pair tables, and how to load one.

Shared by test_segseg_force_pydis_exadis.py and test_segseg_force_torch.py
so the two run over exactly the same data. They did already, by both
naming the same two stems, but from two separate declarations, which is
the arrangement that quietly drifts: adding a third table or renaming one
would have updated whichever test the author was looking at.

Deliberately not named test_*.py. The makefiles and the CMake glob
collect test_*.py, and this is support code.

Each table is two files, following the split the rest of this folder
uses:

    input_data/<stem>.dat   18 columns of geometry, p1 p2 p3 p4 b12 b34
    ref_data/<stem>.npz     the reference forces f1..f4, plus mu, nu, a,
                            the row count and a digest of the geometry
                            they were computed from

`tol` travels with the table rather than with the test because it is a
property of the geometry: it is that table's own floating-point noise
floor with margin, so it says how closely *any* implementation can be
expected to reproduce the reference over those pairs. Every kernel is
held to it, the compiled ParaDiS library and ExaDiS included.

It did not start that way. The compiled kernel and ExaDiS were held to
1e-9 while the numpy lineage got the table's looser number, on the
grounds that those two reproduce the reference exactly. On the machine
where the references were generated they do, to 0.0. On another machine
they do not, and cannot: a different compiler, different flags, and CUDA
rather than OpenMP mean a different order of operations. Measured on an
A100 node, the same two tables give

    ParaDiS C   2.1772e-11    2.4387e-11
    ExaDiS      2.4473e-10    5.5384e-08

and the last of those failed a 1e-9 tolerance while sitting essentially
on the near-parallel table's measured noise floor of 5.84e-08. A 1e-9
that only holds on the machine that generated the reference is a
bit-reproducibility check wearing an accuracy check's clothing, and it
fails for everyone else without a defect in sight.

What is lost by widening it is real and worth naming: on the machine that
did generate the references, 1e-9 would have caught a change in the
compiled kernel's operation sequence that these tolerances now let
through. That signal was only ever available on one machine, and a test
that fails elsewhere for no reason is the worse trade.
"""

import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parents[2]
for _p in [str(_root / p) for p in ['python', 'lib', 'core/pydis/python']]:
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np
from framework.testing import array_digest

INPUT_DIR = _here / 'input_data'
REF_DIR = _here / 'ref_data'


class Table:
    """Table: one stored segment-pair table"""

    def __init__(self, stem, description, tol):
        self.stem = stem
        self.description = description
        self.tol = tol

    @property
    def dat(self):
        return INPUT_DIR / (self.stem + '.dat')

    @property
    def npz(self):
        return REF_DIR / (self.stem + '.npz')


# The first is the original random-geometry set. Its smallest 1-c^2 is
# 3.4e-03, so it never approaches the eps = 1e-4 threshold at which the
# kernels switch to the parallel formula. The second was added because
# three bugs lived in that region of the python kernel undetected:
# inconsistent [nonpar] subsetting, a parallel threshold of 1e-6 instead
# of 1e-4, and two correction branches passing lists where arrays are
# required. Regenerate it with make_parallel_ref_data.py.
#
# Each tolerance is that table's own noise floor with roughly 20x margin.
# The floor is how far the compiled kernel's answer moves when one input
# coordinate is changed by one ulp, worst case over 40 independent nudges:
#
#   random geometry        floor 3.71e-10   ->  1e-8
#   parallel/near-parallel floor 5.84e-08   ->  1e-6
#
# The near-parallel floor is two orders looser because the kernel is
# ill-conditioned there, not because any implementation is at fault; the
# measurement is in the block comment in
# test_segseg_force_pydis_exadis.py.
TABLES = [
    Table('segsep_min_2.5_max_32.5_iso_randombvecs_a0.010',
          'random geometry, 1-c^2 from 3.4e-03 to 0.95', 1e-8),
    Table('parallel_and_nearparallel_iso_a0.010',
          'parallel and near-parallel, 1-c^2 from 0 to 1e-02', 1e-6),
]

# the constants every stored reference was generated with
MU, NU, A = 50.0, 0.3, 0.01


def load(table, verbose=True):
    """load: the segment pairs and reference forces of one table

    Returns (pairs, forces), or (None, None) with an explanation printed
    if either file is missing or the two describe different pairs.

    The fixture and the reference are blessed at different moments: a
    regenerated fixture lands in input_data/ immediately while its
    reference waits in output/ to be copied. Comparing forces against the
    wrong geometry would look like a kernel bug, so the digest carried in
    the .npz is checked and a mismatch refuses rather than reports a
    failure that is not one.
    """
    for path, what in ((table.dat, 'fixture'), (table.npz, 'reference')):
        if path.exists():
            continue
        print("segseg_tables: no %s at %s/%s"
              % (what, path.parent.name, path.name))
        waiting = Path('output') / path.name
        if waiting.exists():
            print("               it is waiting to be blessed:")
            print("                   cp %s %s/" % (waiting, path.parent.name))
        else:
            print("               regenerate with 'make segseg_force_ref'")
        return None, None

    geometry = np.loadtxt(table.dat)
    ref = np.load(table.npz, allow_pickle=False)
    digest = array_digest(geometry)
    if str(ref['geometry_digest']) != digest:
        print("segseg_tables: '%s' fixture and reference describe different "
              "segment pairs" % table.stem)
        print("               fixture   %s, %d pairs"
              % (digest, geometry.shape[0]))
        print("               reference %s, %d pairs"
              % (str(ref['geometry_digest']), int(ref['n_pairs'])))
        print("               the fixture has been regenerated and the new "
              "reference not yet blessed:")
        print("                   cp output/%s ref_data/" % table.npz.name)
        return None, None

    if verbose:
        print("load: '%s', %d segment pairs, mu = %g, nu = %g, a = %g"
              % (table.stem, geometry.shape[0], MU, NU, A))
    pairs = tuple(geometry[:, i:i+3] for i in range(0, 18, 3))
    return pairs, ref['forces']


# cell edge for an exadis network built from a table, as a multiple of the
# extent of that table. Nothing depends on it: the cell is non-periodic, so
# there are no images, and compute_force_segseglist builds no neighbor list.
BOX_FACTOR = 10.0


def exadis_network(pairs):
    """exadis_network: one ExaDisNet holding every pair of a table

    Two segments per row, four nodes per row, nothing shared between rows,
    so each node carries exactly one of f1..f4. Returns the network and the
    pair list to evaluate.

    Those are open lines, so exadis warns on construction that the Burgers
    vector is not conserved; that is expected for a table of test data.
    PyDiS cannot hold such a network at all, which is why this builds an
    ExaDisNet directly rather than going through DisNetManager.

    Shared because two callers need it: the comparison against the stored
    references, and the benchmark, which has to time exadis on the same
    pairs it times the other kernels on.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet

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
