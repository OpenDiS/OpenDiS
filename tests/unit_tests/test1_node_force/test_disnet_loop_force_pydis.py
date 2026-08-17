"""Nodal forces on the 10-loop configuration, computed by PyDiS.

Two force modes are checked against one stored reference:

  LineTension     core force + Peach-Koehler force from applied stress
  Elasticity_SBA  self force + PK force + the segment-pair sum

test_disnet_loop_force_exadis.py compares against the same reference
with the same constants and the same cutoff, so the two codes are
pinned to one set of numbers rather than only to each other. A change
here needs the same change there.

CUTOFF keeps every pair in this configuration (largest mid-point
separation 135.3), so it truncates nothing; testing the truncation
belongs with a case built for it. atol = 1e-6 leaves room for the
exadis companion, which lands within 5e-10 of these numbers.

Regenerate the reference with

    python3 test_disnet_loop_force_pydis.py --write-ref
    cp output/loop_node_force_ref.npz ref_data/

The .npz records the constants it was made with, and both tests check
them, so a stale reference is reported as stale rather than as a
physics failure.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in
               ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in pydis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
np.set_printoptions(threshold=20, edgeitems=5)

from pydis.disnet import DisNet, Cell
from pydis.calforce.calforce_disnet import CalForce
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import load_force_ref, report, report_close

# the configuration, as node positions and connectivity
RN_FILE    = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the reference both this test and the exadis companion compare against,
# and where --write-ref puts a newly computed one
REF_NAME = 'loop_node_force_ref.npz'
OUT_DIR  = 'output'

# segment-pair cutoff, and the longest segment the mid-point stage of
# that cutoff has to allow for (the longest one here is 1.395)
CUTOFF = 200.0
MAXSEG = 2.0

# cell edge, as a multiple of CUTOFF + MAXSEG. 3 is the smallest value
# that keeps the exadis neighbor bins unclamped, which is what the
# companion test needs; here it only has to hold the configuration
BOX_BINS = 3.0

state = {"burgmag": 3e-10, "mu": 50.0, "nu": 0.3, "a": 0.01,
         "maxseg": MAXSEG, "rann": 3.0}
Ec_linetension = 1.0e6

atol = 1.0e-6


def init_loop_from_file(rn_file=RN_FILE, links_file=LINKS_FILE):
    """init_loop_from_file: build the loop configuration

    The cell is non-periodic, so closest_image is the identity and the
    forces do not depend on it. It is sized at BOX_BINS*(cutoff+maxseg)
    anyway, matching the exadis test, where that size is what keeps the
    ExaDiS neighbor bins unclamped
    (.plan/2026-07-27/plan_pydis_elast.md 9.1). Giving both tests the
    same cell keeps two descriptions of one setup from drifting apart.
    """
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'"
          % (rn_file, links_file))
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)
    L = BOX_BINS * (CUTOFF + MAXSEG)
    h = L * np.eye(3)
    origin = 0.5*(rn.min(axis=0) + rn.max(axis=0)) - 0.5*L
    check_cutoff_maxseg(h, CUTOFF, MAXSEG)
    G = DisNet(cell=Cell(h=h, origin=origin,
                         is_periodic=[False, False, False]))
    G.add_nodes_segments_from_list(rn, links)
    return G


def compute_forces(G):
    """compute_forces: (tags, line tension force, elastic force)

    One CalForce serves both: Ec is read only by the LineTension path,
    and the Elasticity_SBA path carries no core term at all
    (.plan/2026-07-27/plan_pydis_elast.md 8.1).
    """
    calforce = CalForce(state=state, Ec=Ec_linetension, cutoff=CUTOFF)
    tags = list(G.all_nodes_tags())

    nodeforce_dict, _ = calforce.NodeForce_LineTension(
        G, applied_stress=np.zeros(6))
    f_lt = np.array([nodeforce_dict[tag] for tag in tags])

    nodeforce_dict, _ = calforce.NodeForce_Elasticity_SBA(
        G, applied_stress=np.zeros(6))
    f_elast = np.array([nodeforce_dict[tag] for tag in tags])

    return np.array(tags, dtype=int), f_lt, f_elast


def write_ref(tags, f_lt, f_elast, out_dir=OUT_DIR):
    """write_ref: save the forces as the reference the tests compare to

    The constants are stored alongside the forces so the file says what
    it is a reference for. The node tags are stored for the same reason:
    the rows are only meaningful in a known order, and a reordering
    should be reported as a reordering rather than as wrong physics.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / REF_NAME
    np.savez(out_file,
             tags=tags,
             force_linetension=f_lt,
             force_elasticity=f_elast,
             mu=state["mu"], nu=state["nu"], a=state["a"],
             Ec=Ec_linetension, cutoff=CUTOFF, maxseg=MAXSEG,
             source=np.array('pydis'))
    print("write_ref: wrote %s (%d nodes)" % (out_file, tags.shape[0]))
    print("")
    print("to install it, copy it into ref_data/ next to this script:")
    print("    cp %s %s" % (out_file, ref_dir))
    print("")
    print("this reference encodes the current settings of the test")
    print("(mu = %g, nu = %g, a = %g, Ec = %g, cutoff = %g, maxseg = %g);"
          % (state["mu"], state["nu"], state["a"], Ec_linetension,
             CUTOFF, MAXSEG))
    print("regenerate it whenever those change, in both this test and")
    print("test_disnet_loop_force_exadis.py")
    return out_file


def load_ref():
    """load_ref: read the stored reference, checking what it describes

    Returns (tags, force_linetension, force_elasticity), or None if the
    file is absent or was generated under different constants.
    """
    hint = ("generate it with\n"
            "    python3 %s --write-ref\n"
            "then copy output/%s into %s"
            % (Path(__file__).name, REF_NAME, ref_dir))
    expected = {"mu": state["mu"], "nu": state["nu"], "a": state["a"],
                "Ec": Ec_linetension, "cutoff": CUTOFF,
                "maxseg": MAXSEG}
    ref = load_force_ref(ref_dir / REF_NAME, expected, hint)
    if ref is None:
        return None
    return ref['tags'], ref['force_linetension'], ref['force_elasticity']


def main(write_ref_file=False):
    G = init_loop_from_file()
    print("nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    tags, lt_force_array, elast_force_array = compute_forces(G)

    if write_ref_file:
        write_ref(tags, lt_force_array, elast_force_array)
        return True

    ref = load_ref()
    if ref is None:
        return False
    ref_tags, ref_lt, ref_elast = ref

    ok = True
    # the reference rows are only meaningful in a known node order, so
    # a reordering should be reported as one, not as wrong physics
    ok &= report("node order matches the reference",
                 np.array_equal(tags, ref_tags))
    ok &= report_close("LineTension", lt_force_array, ref_lt, atol)
    ok &= report_close("Elasticity", elast_force_array, ref_elast, atol)
    return bool(ok)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='nodal forces on the 10-loop configuration (pydis)')
    parser.add_argument('--write-ref', dest='write_ref',
                        action='store_true', default=False,
                        help='save the computed forces as a reference '
                             '.npz under output/ instead of comparing '
                             'against the stored one')
    args = parser.parse_args()

    passed = main(write_ref_file=args.write_ref)

    if not args.write_ref:
        tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
        print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
