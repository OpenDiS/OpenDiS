"""Nodal forces on the 10-loop configuration, computed by ExaDiS.

The counterpart of test_disnet_loop_force_pydis.py, on the same
configuration and constants and against the reference that test writes,
so the two codes are pinned to one set of numbers. The force modes map
across as:

  LineTension   ForceSegLT<CoreDefault, selfforce=false>, i.e. core
                force + PK force, the self force excluded
  CUTOFF_MODEL  self force + PK force + the segment-pair sum, with
                Ec = 0.0 to drop the core term Elasticity_SBA does not
                have

The cell is non-periodic, so the forces do not depend on it, but its
size is not free: ExaDiS clamps its neighbor bins to 3 per direction
once cutoff + maxseg > d/3 and then drops pairs silently, landing about
5e-3 below the all-pairs answer here (plan 9.1). Hence a box of
BOX_BINS*(CUTOFF + maxseg), asserted by check_cutoff_maxseg, and
check_cell_independence as the guard if that ever stops holding.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pyexadis_paths = [str(opendis_root / p) for p in
                  ['python', 'lib', 'core/pydis/python',
                   'core/exadis/python']]
[sys.path.append(p) for p in pyexadis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
np.set_printoptions(threshold=20, edgeitems=5)

import pyexadis
from pyexadis_base import ExaDisNet, CalForce
from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import load_force_ref, report, report_close

# the configuration, as node positions and connectivity
RN_FILE    = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the reference, written by test_disnet_loop_force_pydis.py --write-ref
REF_NAME = 'loop_node_force_ref.npz'

# segment-pair cutoff, shared with the pydis test
CUTOFF = 200.0

# cell edge, as a multiple of CUTOFF + maxseg. 3 is the smallest value
# that keeps the neighbor bins unclamped
BOX_BINS = 3.0

# box scale-up used by check_cell_independence, and what the pair sum is
# allowed to move by when the bin count changes under it
BOX_FACTOR_CHECK = 10.0
TOL_SUM = 1.0e-10

state = {"burgmag": 3e-10, "mu": 50.0, "nu": 0.3, "a": 0.01,
         "maxseg": 2.0, "minseg": 0.5, "rann": 3.0}
Ec_linetension = 1.0e6

atol = 1e-6


def init_loop_from_file(rn_file=RN_FILE, links_file=LINKS_FILE,
                        box_factor=1.0):
    """init_loop_from_file: ExaDisNet from the PyDiS test's input data

    The input files carry no simulation cell and no glide planes, because
    neither enters the quantity being tested. The plane normals are left
    at zero: cross(b, t) is not usable on a loop, since it vanishes
    wherever the line runs screw.

    The box is BOX_BINS*(cutoff + maxseg) so the neighbor bins stay
    unclamped, and box_factor scales it so check_cell_independence can
    vary it.
    """
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'"
          % (rn_file, links_file))
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)

    # exadis wants [n1, n2, bx, by, bz, nx, ny, nz] per segment
    segs = np.zeros((links.shape[0], 8))
    segs[:, :5] = links[:, :5]

    L = box_factor * BOX_BINS * (CUTOFF + state["maxseg"])
    h = L * np.eye(3)
    lo, hi = rn.min(axis=0), rn.max(axis=0)
    if np.any(hi - lo > L):
        raise ValueError("init_loop_from_file: box too small for the "
                         "configuration")
    origin = 0.5*(lo + hi) - 0.5*L
    check_cutoff_maxseg(h, CUTOFF, state["maxseg"])
    cell = pyexadis.Cell(h=h, origin=origin,
                         is_periodic=[False, False, False])

    return cell, ExaDisNet(cell, rn, segs)


def node_force(G, force_mode, **kwargs):
    """node_force: nodal forces from one force mode, ordered by node tag

    Returned in tag order so the rows line up with the reference files,
    which are written in the order PyDiS' all_nodes_tags() yields. Both
    codes tag the nodes (0, i) in input file order, so this is the
    identity here, but relying on that silently is how a reordering
    turns into a wrong-looking physics failure.
    """
    calforce = CalForce(state=state, force_mode=force_mode, **kwargs)
    s = dict(state)
    s["applied_stress"] = np.zeros(6)
    s = calforce.NodeForce(DisNetManager(G), s)
    f = np.array(s["nodeforces"])
    tags = np.array(s["nodeforcetags"])
    return f[np.lexsort((tags[:, 1], tags[:, 0]))]


def force_linetension(G):
    """force_linetension: core + PK, matching NodeForce_LineTension"""
    return node_force(G, 'LineTension', Ec=Ec_linetension)


def force_elasticity(G):
    """force_elasticity: self + PK + pair sum, matching Elasticity_SBA"""
    return node_force(G, 'CUTOFF_MODEL', Ec=0.0, cutoff=CUTOFF)


def node_tags(G):
    """node_tags: the network's node tags, in the order forces come out"""
    tags = np.array(G.get_tags())
    return tags[np.lexsort((tags[:, 1], tags[:, 0]))]


def load_ref():
    """load_ref: read the reference the PyDiS test writes

    Returns (tags, force_linetension, force_elasticity), or None.
    """
    hint = ("generate it with\n"
            "    python3 test_disnet_loop_force_pydis.py --write-ref\n"
            "then copy output/%s into %s" % (REF_NAME, ref_dir))
    expected = {"mu": state["mu"], "nu": state["nu"], "a": state["a"],
                "Ec": Ec_linetension, "cutoff": CUTOFF,
                "maxseg": state["maxseg"]}
    ref = load_force_ref(ref_dir / REF_NAME, expected, hint)
    if ref is None:
        return None
    return ref['tags'], ref['force_linetension'], ref['force_elasticity']


def check_cell_independence(f_lt, f_elast):
    """check_cell_independence: the cell must not affect the forces

    Rebuilds the configuration in a BOX_FACTOR_CHECK times larger box.
    There are no periodic images and the cutoff already covers every
    pair, so nothing is left for the cell to change; this is what fails
    if the cutoff or the box ever stop satisfying check_cutoff_maxseg,
    since the pairs dropped then differ between the two box sizes.

    The line tension is per-segment and must match exactly. The pair
    sum must not: a different bin count changes the order the
    contributions are accumulated in, and the threaded reduction is not
    order-independent, hence TOL_SUM.
    """
    _, G = init_loop_from_file(box_factor=BOX_FACTOR_CHECK)
    label = "%gx box" % BOX_FACTOR_CHECK
    ok = report_close(label + ", LineTension", force_linetension(G),
                      f_lt, atol=0.0)
    ok &= report_close(label + ", Elasticity", force_elasticity(G),
                       f_elast, atol=TOL_SUM)
    return bool(ok)


def main():
    _, G = init_loop_from_file()
    print("nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    ref = load_ref()
    if ref is None:
        return False
    ref_tags, ref_lt, ref_elast = ref

    lt_force_array = force_linetension(G)
    elast_force_array = force_elasticity(G)

    ok = True
    # the reference rows are ordered by the tags PyDiS wrote; comparing
    # forces row by row is only meaningful if this network agrees
    ok &= report("node order matches the reference",
                 np.array_equal(node_tags(G), ref_tags))
    ok &= report_close("LineTension", lt_force_array, ref_lt, atol)
    ok &= report_close("Elasticity", elast_force_array, ref_elast, atol)
    ok &= check_cell_independence(lt_force_array, elast_force_array)
    return bool(ok)


if __name__ == "__main__":
    pyexadis.initialize()
    passed = main()
    pyexadis.finalize()

    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
