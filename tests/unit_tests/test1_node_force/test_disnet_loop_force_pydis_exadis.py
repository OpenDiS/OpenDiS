"""Nodal forces on the 10-loop configuration, PyDiS and ExaDiS.

Both codes compute the same two force modes on the same configuration,
with the same constants, cell and cutoff, and both are compared against
one stored reference, so they are pinned to one set of numbers rather
than only to each other. The force modes map across as:

  line tension  pydis  NodeForce_LineTension, core force + PK force
                exadis LINE_TENSION_MODEL, ForceSegLT with
                       selfforce=false, i.e. the same two terms
  elasticity    pydis  NodeForce_Elasticity_SBA with Ec = 0.0, self
                       force + PK force + the segment-pair sum
                exadis CUTOFF_MODEL with Ec = 0.0, the same three terms

Both elasticity calls set Ec = 0.0. Elasticity_SBA gained a core term,
so the term now has to be switched off on both sides to compare these
three, where before it was absent from pydis and only switched off in
exadis. Testing the core term itself wants a case built for it: with
these constants it is 3.4e4 times the elastic force, so folding it in
here would leave the elasticity comparison measuring the core term and
almost nothing else.

The reference is written by this test, from the pydis side:

    python3 test_disnet_loop_force_pydis_exadis.py --write-ref
    cp output/loop_node_force_ref.npz ref_data/

It records the constants it was made with, and they are checked before
comparing, so a stale reference is reported as stale rather than as a
physics failure. Since the reference is the pydis result, the exadis
error against it is also the difference between the two codes.

CUTOFF keeps every pair in this configuration, so it truncates nothing;
testing the truncation belongs with a case built for it. The cell is
non-periodic and the forces do not depend on it, but its size is not
free: exadis clamps its neighbor bins to 3 per direction once
cutoff + maxseg > d/3 and then drops pairs silently, landing about 5e-3
below the all-pairs answer here. Hence a box of BOX_BINS*(CUTOFF +
MAXSEG), which check_cutoff_maxseg asserts on every call.
"""

import os
# Kokkos reads this at initialize(), so it must be set before pyexadis is
# imported, and therefore before any framework import that might pull it in.
#
# Assigned rather than setdefault, as in test3_remesh_rule and
# test4_collision_mode. ExaDiS sums the segment-pair forces with a parallel
# reduction, so the thread count changes the order the additions happen in and
# with it the last bits of the answer: measured here, the exadis elasticity
# residual is 1.4211e-14 on one thread and 1.9540e-14 on eight, each stable.
# Neither is a defect, but leaving it to the machine's core count means the
# number quietly depends on where the test runs.
os.environ['OMP_NUM_THREADS'] = '1'

import re
import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
np.set_printoptions(threshold=20, edgeitems=5)

from pydis.disnet import DisNet, Cell
from pydis.calforce.calforce_disnet import CalForce as PyCalForce
from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import (load_force_ref, report, report_close,
                               quiet_native_output, kokkos_summary)
from pydis.build_info import bitrepro_math, build_description
from pydis.calforce.bitrepro_math import ENABLED as numpy_bitrepro_math

# the configuration, as node positions and connectivity
RN_FILE    = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the reference both codes are compared against, and where --write-ref
# puts a newly computed one
REF_NAME = 'loop_node_force_ref.npz'
OUT_DIR  = 'output'

# segment-pair cutoff, and the longest segment the mid-point stage of
# that cutoff has to allow for (the longest one here is 1.395)
CUTOFF = 200.0
MAXSEG = 2.0

# cell edge, as a multiple of CUTOFF + MAXSEG. 3 is the smallest value
# that keeps the exadis neighbor bins unclamped
BOX_BINS = 3.0

state = {"burgmag": 3e-10, "mu": 50.0, "nu": 0.3, "a": 0.01,
         "maxseg": MAXSEG, "minseg": 0.5, "rann": 3.0}
Ec_linetension = 1.0e6
Ec_elasticity = 0.0

atol = 1.0e-6

# Elasticity's self/pair forces come from SegSegForce.c (its Ec=0 core term
# is exactly zero regardless of numpy, since it is multiplied by Ec), so a
# bitwise-reproducible build must reproduce a reference blessed by another
# such build exactly there. Line tension never touches SegSegForce.c at
# all -- it is selfforcevec_LineTension with Ec=1e6, whose np.dot()/
# np.linalg.norm() calls dispatch to the environment's BLAS (see
# pydis.calforce.bitrepro_math) -- so its own exactness depends on that
# separate, Python-side flag instead. ExaDiS is built against its own
# platform libm regardless, and keeps atol for both.
atol_lt     = 0.0 if numpy_bitrepro_math else atol
atol_elast  = 0.0 if bitrepro_math()      else atol

# ExaDiS is held to zero on the line tension and close to it on the elasticity.
#
# Line tension is exact. It is the core force alone, and on a _repro build
# selfforcevec_LineTension evaluates it in ExaDiS' grouping, so the two agree to
# the last bit against a pydis-blessed reference. Zero rather than a tolerance,
# so that losing that grouping shows up as a failure rather than as drift.
#
# Elasticity cannot be exact, and the reason is not a formula or a rounding
# choice: the pair sum is accumulated in a different *shape* by each code, per
# node in ExaDiS against per segment then per node here, and 900 summands per
# node re-associated differently land one to two ULP apart. Measured 1.4211e-14
# single-threaded, with pydis 7.99e-15 and ExaDiS 1.33e-14 from the correctly
# rounded sum, so neither is at fault. Matching it would mean adopting ExaDiS'
# accumulation shape, which was measured and rejected: see the note on
# OMP_NUM_THREADS above for why the thread count has to be pinned for even this
# number to be stable. 1e-12 is the observed value with about 70x margin.
atol_exadis_lt    = 0.0 if numpy_bitrepro_math else atol
atol_exadis_elast = 1.0e-12


def load_config(rn_file=RN_FILE, links_file=LINKS_FILE):
    """load_config: node positions and links, shared by both codes"""
    print("load_config: rn_file = '%s', links_file = '%s'"
          % (rn_file, links_file))
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)
    return rn, links


def cell_geometry(rn):
    """cell_geometry: the cell both codes are given, as (h, origin)"""
    L = BOX_BINS * (CUTOFF + MAXSEG)
    h = L * np.eye(3)
    lo, hi = rn.min(axis=0), rn.max(axis=0)
    if np.any(hi - lo > L):
        raise ValueError("cell_geometry: box too small for the "
                         "configuration")
    check_cutoff_maxseg(h, CUTOFF, MAXSEG)
    return h, 0.5*(lo + hi) - 0.5*L


def forces_pydis():
    """forces_pydis: (tags, line tension force, elastic force)

    One CalForce per mode, because Ec belongs to the object and the two
    modes want different values: Elasticity_SBA reads Ec now, so a
    single object carrying Ec_linetension would put a core term into the
    elastic result that the exadis call is not asked for.
    """
    rn, links = load_config()
    h, origin = cell_geometry(rn)
    G = DisNet(cell=Cell(h=h, origin=origin,
                         is_periodic=[False, False, False]))
    G.add_nodes_segments_from_list(rn, links)
    print("pydis:  nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    tags = list(G.all_nodes_tags())

    calforce = PyCalForce(state=state, Ec=Ec_linetension, cutoff=CUTOFF)
    nodeforce_dict, _ = calforce.NodeForce_LineTension(
        G, applied_stress=np.zeros(6))
    f_lt = np.array([nodeforce_dict[tag] for tag in tags])

    calforce = PyCalForce(state=state, Ec=Ec_elasticity, cutoff=CUTOFF)
    nodeforce_dict, _ = calforce.NodeForce_Elasticity_SBA(
        G, applied_stress=np.zeros(6))
    f_elast = np.array([nodeforce_dict[tag] for tag in tags])

    return np.array(tags, dtype=int), f_lt, f_elast


def forces_exadis():
    """forces_exadis: (tags, line tension force, elastic force)

    The input files carry no glide planes, and none is needed: the
    normals are left at zero because cross(b, t) is not usable on a
    loop, vanishing wherever the line runs screw, and neither force
    mode reads them.
    """
    from pyexadis_base import ExaDisNet, CalForce
    import pyexadis

    rn, links = load_config()
    h, origin = cell_geometry(rn)
    # exadis wants [n1, n2, bx, by, bz, nx, ny, nz] per segment
    segs = np.zeros((links.shape[0], 8))
    segs[:, :5] = links[:, :5]
    cell = pyexadis.Cell(h=h, origin=origin,
                         is_periodic=[False, False, False])
    G = ExaDisNet(cell, rn, segs)
    print("exadis: nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))
    N = DisNetManager(G)

    def node_force(force_mode, **kwargs):
        calforce = CalForce(state=state, force_mode=force_mode, **kwargs)
        s = dict(state)
        s["applied_stress"] = np.zeros(6)
        s = calforce.NodeForce(N, s)
        f = np.array(s["nodeforces"])
        tags = np.array(s["nodeforcetags"])
        # ordered by tag, so the rows line up with the reference; both
        # codes tag the nodes (0, i) in input file order, but relying on
        # that silently is how a reordering turns into wrong physics
        order = np.lexsort((tags[:, 1], tags[:, 0]))
        return tags[order], f[order]

    tags, f_lt = node_force('LineTension', Ec=Ec_linetension)
    _, f_elast = node_force('CUTOFF_MODEL', Ec=Ec_elasticity, cutoff=CUTOFF)
    return tags, f_lt, f_elast


def write_ref(tags, f_lt, f_elast, out_dir=OUT_DIR):
    """write_ref: save the pydis forces as the stored reference

    The constants are stored alongside the forces so the file says what
    it is a reference for. The node tags are stored for the same
    reason: the rows are only meaningful in a known order.
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
    print("to install it, copy it into the ref_data/ next to this "
          "script:")
    print("    cp %s ref_data/%s" % (out_file, REF_NAME))
    print("")
    print("this reference encodes the current settings of the test")
    print("(mu = %g, nu = %g, a = %g, Ec = %g, cutoff = %g, maxseg = %g);"
          % (state["mu"], state["nu"], state["a"], Ec_linetension,
             CUTOFF, MAXSEG))
    print("regenerate it whenever any of those change")
    return out_file


def load_ref():
    """load_ref: read the stored reference, checking what it describes

    Returns (tags, force_linetension, force_elasticity), or None if the
    file is absent or was generated under different constants.
    """
    hint = ("generate it with\n"
            "    make loop_node_force_ref\n"
            "then copy %s/%s into ref_data/" % (OUT_DIR, REF_NAME))
    expected = {"mu": state["mu"], "nu": state["nu"], "a": state["a"],
                "Ec": Ec_linetension, "cutoff": CUTOFF,
                "maxseg": MAXSEG}
    ref = load_force_ref(ref_dir / REF_NAME, expected, hint)
    if ref is None:
        return None
    return ref['tags'], ref['force_linetension'], ref['force_elasticity']


def compare_to_ref(label, tags, f_lt, f_elast, ref, tol_lt=atol, tol_elast=atol,
                   note_elast=None):
    """compare_to_ref: report one code's three comparisons

    Two tolerances, not one: line tension and elasticity are bitwise
    reproducible under different conditions (see atol_lt/atol_elast above),
    so a caller that cares about that distinction needs to pass them
    separately rather than getting one tol applied to both.
    """
    ref_tags, ref_lt, ref_elast = ref
    ok = report("%s: node order matches the reference" % label,
                np.array_equal(tags, ref_tags))
    ok &= report_close("%s: line tension" % label, f_lt, ref_lt, tol_lt)
    if note_elast:
        print(note_elast)
    ok &= report_close("%s: elasticity " % label, f_elast, ref_elast, tol_elast)
    return bool(ok)


def test_pydis(ref):
    """test_pydis: the pydis forces against the stored reference"""
    tags, f_lt, f_elast = forces_pydis()
    return compare_to_ref("pydis ", tags, f_lt, f_elast, ref,
                          tol_lt=atol_lt, tol_elast=atol_elast)


def test_exadis(ref):
    """test_exadis: the exadis forces against the same reference"""
    tags, f_lt, f_elast = forces_exadis()
    # one line on screen; the reasoning is in the note on atol_exadis_elast
    note = ("note: pydis cannot match exadis' summation order exactly, so the "
            "next line is not held to zero")
    return compare_to_ref("exadis", tags, f_lt, f_elast, ref,
                          tol_lt=atol_exadis_lt, tol_elast=atol_exadis_elast,
                          note_elast=note)


def kokkos_threads(banner):
    """kokkos_threads: the thread pool size Kokkos reports, or None

    From the thread_pool_topology[ N x T x V ] line of the startup banner,
    whose middle number is the thread count.
    """
    m = re.search(r'thread_pool_topology\[\s*\d+\s*x\s*(\d+)', banner)
    return int(m.group(1)) if m else None


def check_threads(banner):
    """check_threads: ExaDiS really is running on the one thread it was pinned to

    Read from the Kokkos banner, not from OMP_NUM_THREADS: this file sets that
    variable unconditionally, so reading it back would only confirm the file
    agrees with itself. The banner is what Kokkos took at initialize(), which
    also catches the assignment being moved after the pyexadis import, where it
    would have no effect while still reading back as '1'.
    """
    threads = kokkos_threads(banner)
    return report("exadis: running on 1 thread (the pair sum's accumulation "
                  "order depends on it), got %s" % threads, threads == 1)


def main(write_ref_file=False):
    if write_ref_file:
        write_ref(*forces_pydis())
        return True

    print(build_description())

    ref = load_ref()
    if ref is None:
        return False

    ok = test_pydis(ref)

    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; the exadis half of this test did "
              "not run")
        return False

    with quiet_native_output() as buf:
        pyexadis.initialize()
    for line in kokkos_summary(buf.text):
        print("exadis: %s" % line)
    ok &= check_threads(buf.text)
    ok &= test_exadis(ref)
    pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='nodal forces on the 10-loop configuration, '
                    'pydis and exadis')
    parser.add_argument('--write-ref', dest='write_ref',
                        action='store_true', default=False,
                        help='save the pydis forces as a reference '
                             '.npz under output/ instead of comparing '
                             'against the stored one')
    args = parser.parse_args()

    passed = main(write_ref_file=args.write_ref)

    if not args.write_ref:
        tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
        print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
