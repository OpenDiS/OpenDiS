"""Topology (multi-node split) handling compared step by step, PyDiS against ExaDiS.

The ExaDiS simulation is authoritative: its TopologySerial is the one whose
result the run carries forward. At every step the network is snapshotted
before topology runs, ExaDiS splits it, and the same snapshot is replayed
through PyDiS's Topology(split_mode='Serial') on a throwaway DisNet. A
disagreement at step n therefore does not change the configuration seen at
step n+1, and the per-step comparisons stay independent of each other.
Structurally the same design as test4_collision_mode's comparing driver,
generalized to a different operator; several of its generic helpers
(canonical_form, network_diff, tag/segment set builders) are reused directly
from there rather than duplicated.

WHAT THIS MEASURES, AND WHY A DISAGREEMENT IS NOT AUTOMATICALLY A BUG

PyDiS's 'Serial' split mode is a transcription of ParaDiS SplitMultiNodes and
is the mode ExaDiS's TopologySerial also ports (see
.plan/2026-08-17/debug_topology.md), so unlike collision this comparison runs
the *same* algorithm on both sides and full agreement is the expectation, not
an aspiration. But the split-direction decision inside that algorithm is, on
any node whose local geometry is exactly or near-exactly symmetric, a genuine
mathematical tie: both mirror-image split directions dissipate the same
power, decided only by which way each code's independent floating-point
roundoff happens to fall. This was measured directly, not assumed (see
.plan/2026-08-21/debug_topology_stage2.md section 2): feeding one code's
exact pre-split geometry into the other's own force and mobility bindings
finds the same tie, at the ~1e-11 to 1e-13 relative level, resolved in
opposite directions.

So a step where both codes split the same node into different results is
classified, not simply failed: classify_split_disagreement (below) replays
the trial through each code's own force/mobility machinery and reports
whether either side's own margin was itself near-tied. A near-tied step is
counted, not asserted on -- it is not a defect in either implementation, and
demanding agreement on a coin flip would make this test measure luck rather
than correctness. A disagreement where *neither* side's margin is a near-tie
is a real finding and this test fails on it.

RUNNING IT

    make                       run at the default step count
    make PYTHON=/opt/anaconda3/bin/python3
    python3 test_topology_mode_pydis_exadis.py --max-step 100
    python3 test_topology_mode_pydis_exadis.py --diagnose

Single-threaded by construction: OMP_NUM_THREADS is pinned to 1 before
pyexadis is imported. This matters more here than for collision:
debug_topology_stage2.md measured that unpinned, multi-threaded ExaDiS can
land on either side of a near-tied split from one run to the next, which
would make this test's pass/fail flicker independently of any PyDiS change.
"""

import os
# must precede the pyexadis import, and therefore any framework import
# that might pull it in
os.environ.setdefault('OMP_NUM_THREADS', '1')

import sys
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if p not in sys.path]

# test4_collision_mode contributes no configuration here (changed 2026-08-21,
# see the plan): only its generic, operator-agnostic comparison helpers
# (canonical_form, network_diff, tag/segment set builders, glide_plane_settings)
# are reused, imported below.
test4_dir = str(Path(__file__).resolve().parents[1] / 'test4_collision_mode')
if test4_dir not in sys.path:
    sys.path.insert(0, test4_dir)

import numpy as np

from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import (report, quiet_native_output, kokkos_summary,
                              max_nearest_distance, GREEN, RED, RESET)

import test_collision_mode_pydis_exadis as t4

# ---------------------------------------------------------------- configuration

# The base is examples/03_binary_junction, current tuned state, decided by
# the user 2026-08-21: two dislocation lines whose coincident free centre
# nodes collide into a 4-arm node that TopologySerial then splits into a
# binary junction. This is the exact case
# .plan/2026-08-21/debug_topology_stage2.md's mirror-tie investigation was
# built around, and the plan's sequencing (finish this test, then return to
# tests/full_runs/03_binary_junction) is why: this test's classifier is
# meant to be the tool that finishes that debugging, so it needs to run on
# that same geometry rather than a different one chosen to avoid the tie.
#
# 'Current tuned state' means Ec=2.8e10 and the EPS_ARM_ASYMMETRY=1e-7
# symmetry break currently in both example scripts -- not the pristine
# geometry, and not Ec=0. Both values are known (debug_topology_stage2.md
# sections 3-6) to be empirical and non-general, not physical thresholds;
# they are carried over here only because this test's base is defined to
# match the example scripts as they currently stand, not because they are
# endorsed as correct.
LBOX = 1000.0
Z0 = 0.125 * LBOX
CUTOFF = 0.25 * LBOX
DT = 1.0e-9
MAX_STEP = 300
EC = 2.8e10
EPS_ARM_ASYMMETRY = 1.0e-7

# Second phase, matching tests/full_runs/03_binary_junction: the junction
# forms unstressed, then this stress turns on and pulls it back apart,
# exercising the topology handler on the way apart too, not just on the way
# in. Same magnitude as that test's UNZIP_STRESS (see
# examples/03_binary_junction/test_binary_junction_pydis_elast.py's comment
# on it for how that value was chosen); the activation step is different by
# request (STRESS_STEP=100 here, a fixed step, vs that test's
# max_step // 2) rather than derived from MAX_STEP.
UNZIP_STRESS = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 4.5e9])
STRESS_STEP = 100

# PyDiS split mode. 'Serial' is the one with an ExaDiS counterpart
# (TopologySerial); 'MaxDiss' is pydis's original, unrelated algorithm and
# has nothing to compare against, the same asymmetry test4_collision_mode's
# PYDIS_MODE default records for collision's 'Proximity' vs 'Retroactive'.
PYDIS_MODE = 'Serial'

# Relative margin below which a disagreement is classified as a known tie
# rather than a real difference. Actual measured ties (debug_topology_stage2.md
# section 2) are ~1e-11 to 1e-13 relative; this is set several orders of
# magnitude above that floor. Not yet calibrated against a real (non-tie)
# disagreement to confirm it does not also swallow one -- an open question
# carried over from the plan (section 4).
TIE_REL_THRESHOLD = 1e-4

# Absolute tolerance for the max nearest-node position difference between
# pydis's replayed and exadis' authoritative network, after topology runs
# each step. Matches tests/full_runs/03_binary_junction's own TOL, same
# length scale (LBOX=1000 in both).
POS_DIFF_TOL = 1.0e-6

OUT_DIR = Path(__file__).resolve().parent / 'output'
REF_DIR = Path(__file__).resolve().parent / 'ref_data'
TAG_STATE = REF_DIR / 'exadis_tag_state.npz'


def base_state():
    """base_state: the material and discretization parameters

    Matches examples/03_binary_junction's current state dict exactly,
    including 'crystal': 'bcc' and 'use_glide_planes': True. Without a
    crystal type, ExaDiS never assigns a glide plane to a segment created by
    a topological split (Crystal::initialize, core/exadis/src/crystal.h:126,
    forces use_glide_planes=0 whenever no crystal type is set at all) --
    resolved to be the actual root cause of the tests/full_runs/03_binary_junction
    divergence (.plan/2026-08-21/debug_topology_stage2.md, section 11), not
    a tie. Checked via t4.glide_plane_settings, same as test4_collision_mode,
    rather than assumed.

    With the crystal type declared, ExaDiS's junction-segment plane happens
    to come out equal to PyDiS's own crystal-agnostic one (see the companion
    exadis example script's comment on this state dict for why: it is a
    property of this test's geometry, not a general equivalence between the
    two models).
    """
    return {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0,
            "maxseg": 0.04 * LBOX, "minseg": 0.01 * LBOX, "rann": 3.0,
            "crystal": "bcc", "use_glide_planes": True}


def init_two_disl_lines(state, pbc=False):
    """init_two_disl_lines: the binary-junction geometry, as the example builds it

    Transcribed from examples/03_binary_junction/test_binary_junction_exadis_elast.py
    rather than imported from it, the same choice test4_collision_mode made
    for the Frank-Read source: that file is a __main__ script with argparse
    and module-level globals set by main(), not a library module meant for
    import. Kept identical, including the EPS_ARM_ASYMMETRY symmetry-breaking
    perturbation and its justification; see that file's own comment for why
    it is needed and .plan/2026-08-21/debug_topology_stage2.md for how the
    value was found.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet, NodeConstraints
    from framework.simulation_setup import remesh_initial_config

    b1 = np.array([-1.0, 1.0, 1.0])
    b2 = np.array([1.0, -1.0, 1.0])
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[pbc, pbc, pbc])
    center = np.array(cell.center())

    PINNED = NodeConstraints.PINNED_NODE
    FREE = NodeConstraints.UNCONSTRAINED
    rn = np.array([[0.0, -Z0, -Z0, PINNED],
                   [0.0, 0.0, 0.0, FREE],
                   [0.0, Z0, Z0, PINNED],
                   [-Z0, 0.0, -Z0, PINNED],
                   [0.0, 0.0, 0.0, FREE],
                   [Z0, 0.0, Z0, PINNED]])
    rn[:, 0:3] += center

    # shift line 1's two pinned endpoints along the line's own direction;
    # see test_binary_junction_pydis_elast.py's comment for why.
    line1_dir = rn[2, :3] - rn[1, :3]
    line1_dir = line1_dir / np.linalg.norm(line1_dir)
    rn[0, :3] += EPS_ARM_ASYMMETRY * Z0 * line1_dir
    rn[2, :3] += EPS_ARM_ASYMMETRY * Z0 * line1_dir

    xi1, xi2 = rn[2, :3] - rn[1, :3], rn[5, :3] - rn[4, :3]
    n1, n2 = np.cross(b1, xi1), np.cross(b2, xi2)
    n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0, :] = np.concatenate(([0, 1], b1, n1))
    links[1, :] = np.concatenate(([1, 2], b1, n1))
    links[2, :] = np.concatenate(([3, 4], b2, n2))
    links[3, :] = np.concatenate(([4, 5], b2, n2))

    rn, links = remesh_initial_config(rn, links, state["maxseg"],
                                      LBOX * np.eye(3), [pbc] * 3)
    return DisNetManager(ExaDisNet(cell, rn, links))


def canonical_form_pbc(data, cell_h, is_periodic, tol=1e-6):
    """canonical_form_pbc: t4.canonical_form, folded into the primary cell first

    Folding only applies under full PBC; this base config is not periodic
    (pbc=False, matching examples/03_binary_junction), so it is a no-op
    here and this stays a pass-through to t4.canonical_form. Kept general
    -- rather than dropped -- because it was needed and verified against a
    periodic configuration while this test's base was still the Frank-Read
    source, and folding is the correct behaviour if a periodic base is used
    again later: a split's new-node position is a computed point that can
    legitimately land at either periodic image of the same physical
    location. t4.canonical_form itself does not fold, so this wraps it
    rather than changing a function test4_collision_mode already certifies.
    """
    if not all(is_periodic):
        return t4.canonical_form(data, tol=tol)
    L = np.diag(np.asarray(cell_h, float))
    data = {**data, 'nodes': {**data['nodes'],
                              'positions': np.mod(np.asarray(data['nodes']['positions'], float), L)}}
    return t4.canonical_form(data, tol=tol)


# ------------------------------------------------------------------- the record

class StepRecord:
    """StepRecord: what one step of the comparison observed

    Plain data, collected for every step and summarized at the end, rather
    than asserted step by step, following test4_collision_mode's StepRecord.
    """

    def __init__(self, istep):
        self.istep = istep
        self.n_nodes_before = 0
        self.n_segs_before = 0
        self.exadis_dn = 0
        self.exadis_dseg = 0
        self.pydis_dn = 0
        self.pydis_dseg = 0
        self.exadis_split = False
        self.pydis_split = False
        self.pos_diff = None      # max nearest-node distance, pydis vs exadis, after topology
        self.tag_state_match = None  # live exadis tag pool vs the stored reference recording
        self.tags_match = None
        self.same_geometry = None
        self.same_counts = None
        self.pydis_error = None
        self.pydis_traceback = None
        self.exadis_diff = None
        self.pydis_diff = None
        # tie classification, only set on a differing split step
        self.classification = None       # 'tie' | 'real_disagreement' | None
        self.pydis_margin = None         # relative (vd1-vd2)/max, or None
        self.exadis_margin = None
        self.split_tag = None            # the multi-arm node tag being split
        self.split_neighbor_order = None  # its arms, in the order the replay saw them
        self.ex_arms = None              # arms exadis moved to its new node
        self.py_arms = None              # arms pydis moved to its new node


def split_target_tag(diff):
    """split_target_tag: which pre-existing node a split's new node came from

    A split adds exactly one node and one segment connecting it back to the
    node it was split from; that connecting segment identifies the parent.
    Returns None if the diff does not look like a single split (e.g. a
    step with more than one simultaneous split, or none).
    """
    if len(diff['added_nodes']) != 1 or len(diff['added_segs']) < 1:
        return None
    new_tag = diff['added_nodes'][0]
    for a, b in diff['added_segs']:
        other = b if a == new_tag else (a if b == new_tag else None)
        if other is not None:
            return other
    return None


def neighbor_order_of(before_data, tag):
    """neighbor_order_of: tag's arms, in the order PyDiS's replay network sees them

    This is exactly the order build_split_list enumerates candidate
    partitions in: it comes from re-importing the network ExaDiS exported
    (DisNet.import_data), which is not guaranteed to match the order PyDiS's
    own native network construction would produce for the same physical
    configuration. Printed on a 'tie' step so that is visible rather than
    asserted.
    """
    from pydis import DisNet
    G = DisNet()
    G.import_data(before_data)
    if not G.has_node(tag):
        return None
    return list(G.neighbors_tags(tag))


# ------------------------------------------------------ tie-vs-bug classifier

def pydis_own_margin(before_data, tag, arms_to_new_node, state):
    """pydis_own_margin: PyDiS's own relative margin for this exact trial

    Replays the specific partition PyDiS (or ExaDiS) used, through PyDiS's
    own evaluate_trial_split machinery, and returns the relative difference
    between the two candidate nodes' squared trial velocities: this is
    exactly the quantity split_direction's '>' comparison decides on, so a
    value near zero means PyDiS's own kernel found this a near-tie,
    independent of what it ultimately chose to do overall.

    Returns None if the partition cannot be replayed (e.g. the node no
    longer has this exact arm set, or it fails PyDiS's own short-arm
    exemption and is never evaluated as a trial at all).
    """
    from pydis import DisNet
    from pydis.topology.topology_serial import (TopologyParams, split_direction,
                                                 split_node_reusing_forces,
                                                 has_short_arm)
    from framework.calforce_base import CalForce_Base

    G = DisNet()
    G.import_data(before_data)
    if not G.has_node(tag):
        return None
    nbrs = list(G.neighbors_tags(tag))
    arms = [n for n in nbrs if n in arms_to_new_node]
    if len(arms) != len(arms_to_new_node) or len(arms) == 0 or len(arms) == len(nbrs):
        return None

    params = state['_t5_topology_params']
    if has_short_arm(G, tag, params.short_seg):
        return None

    force = state['_t5_pydis_force']
    mobility = state['_t5_pydis_mobility']
    state_trial = dict(state)
    try:
        state_trial, node1, node2 = split_node_reusing_forces(
            G, state_trial, tag, arms, mobility)
    except Exception:
        return None
    v1, v2 = state_trial["vel_dict"][node1], state_trial["vel_dict"][node2]
    vd1, vd2 = float(np.dot(v1, v1)), float(np.dot(v2, v2))
    if not (np.isfinite(vd1) and np.isfinite(vd2)):
        return None
    return abs(vd1 - vd2) / max(vd1, vd2, 1e-300)


def exadis_own_margin(before_data, tag, arms_to_new_node, state):
    """exadis_own_margin: ExaDiS's own relative margin for the same trial

    Constructs the pre-split geometry as an ExaDisNet, adds a zero-separation
    trial node carrying arms_to_new_node the way ExaDiS's own split_node
    does internally, and reads the two candidate nodes' velocities straight
    out of ExaDiS's own CalForce/MobilityLaw bindings -- the same code its
    C++ topology routine calls. This is the dual-kernel replication
    technique of debug_topology_stage2.md section 7, generalized into a
    reusable function instead of a throwaway script.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet, CalForce, MobilityLaw

    tags = [tuple(int(v) for v in t) for t in before_data['nodes']['tags']]
    if tag not in tags:
        return None
    tag_idx = tags.index(tag)
    rn = np.asarray(before_data['nodes']['positions'], float)
    constraints = np.asarray(before_data['nodes']['constraints'], float).reshape(-1, 1)
    rn = np.hstack([rn, constraints])

    links = []
    for k, (a, b) in enumerate(before_data['segs']['nodeids']):
        burg = np.asarray(before_data['segs']['burgers'][k], float)
        plane = np.asarray(before_data['segs']['planes'][k], float)
        links.append([int(a), int(b)] + list(burg) + list(plane))
    links = np.asarray(links, float)

    new_idx = rn.shape[0]
    rn = np.vstack([rn, rn[tag_idx:tag_idx + 1, :]])
    links_new = []
    arm_idx = {tags.index(n) for n in arms_to_new_node if n in tags}
    for l in links:
        i, j = int(l[0]), int(l[1])
        if i == tag_idx and j in arm_idx:
            l = l.copy(); l[0] = new_idx
        elif j == tag_idx and i in arm_idx:
            l = l.copy(); l[1] = new_idx
        links_new.append(l)
    links_new = np.asarray(links_new, float)
    links_new = np.vstack([links_new,
                           [tag_idx, new_idx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    cell_h = np.asarray(before_data['cell']['h'], float)
    cell = pyexadis.Cell(h=cell_h,
                         is_periodic=[bool(x) for x in before_data['cell']['is_periodic']])
    net = ExaDisNet(cell, rn, links_new)

    ex_state = dict(state)
    ex_state.setdefault('applied_stress', np.zeros(6))
    calforce = CalForce(force_mode='CUTOFF_MODEL', state=ex_state, Ec=EC, cutoff=CUTOFF)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=ex_state)
    N = DisNetManager(net)
    new_tags = net.get_tags()
    t0, t1 = tuple(new_tags[tag_idx]), tuple(new_tags[new_idx])

    try:
        f0 = calforce.OneNodeForce(N, ex_state, t0, update_state=True)
        f1 = calforce.OneNodeForce(N, ex_state, t1, update_state=True)
        v0 = mobility.OneNodeMobility(N, ex_state, t0, f0, update_state=True)
        v1 = mobility.OneNodeMobility(N, ex_state, t1, f1, update_state=True)
    except Exception:
        return None
    vd0, vd1 = float(np.dot(v0, v0)), float(np.dot(v1, v1))
    if not (np.isfinite(vd0) and np.isfinite(vd1)):
        return None
    return abs(vd0 - vd1) / max(vd0, vd1, 1e-300)


def classify_split_disagreement(before_data, rec, state):
    """classify_split_disagreement: tie, or a real difference

    Only called on a step where the two codes' post-split networks differ
    and at least one of them split a node. Identifies the node each side
    split from its own diff, replays the corresponding trial through both
    codes' own force/mobility kernels, and classifies the step as a known
    tie if either side's own margin is below TIE_REL_THRESHOLD.

    Returns 'tie', 'real_disagreement', or None when the step's shape does
    not match a single-node split on at least one side (multiple
    simultaneous splits, a pure position disagreement with no topology
    change, etc.) -- reported separately rather than forced into one of
    the two buckets.
    """
    ex_tag = split_target_tag(rec.exadis_diff) if rec.exadis_split else None
    py_tag = split_target_tag(rec.pydis_diff) if rec.pydis_split else None
    tag = ex_tag or py_tag
    if tag is None or (ex_tag and py_tag and ex_tag != py_tag):
        return None

    def arms_moved(diff):
        # added_segs includes the connecting segment back to the parent
        # (tag) itself, alongside the arms actually transferred; excluded
        # here, or the replay tries to move the parent onto the new node as
        # one of its own arms.
        new_tag = diff['added_nodes'][0]
        return [b if a == new_tag else a for a, b in diff['added_segs']
               if new_tag in (a, b) and tag not in (a, b)]

    rec.split_tag = tag
    rec.split_neighbor_order = neighbor_order_of(before_data, tag)
    rec.ex_arms = arms_moved(rec.exadis_diff) if ex_tag == tag else None
    rec.py_arms = arms_moved(rec.pydis_diff) if py_tag == tag else None

    # Whichever side actually split is the trial to replay; when both did,
    # they may or may not agree on which arms moved (see diagnose()'s
    # printout when they do not -- a different, coarser kind of tie than the
    # direction-within-the-same-partition case this function was built for).
    arms_to_new_node = rec.ex_arms if rec.ex_arms is not None else rec.py_arms

    rec.pydis_margin = pydis_own_margin(before_data, tag, arms_to_new_node, state)
    rec.exadis_margin = exadis_own_margin(before_data, tag, arms_to_new_node, state)

    margins = [m for m in (rec.pydis_margin, rec.exadis_margin) if m is not None]
    if margins and min(margins) < TIE_REL_THRESHOLD:
        return 'tie'
    if margins:
        return 'real_disagreement'
    return None


# ------------------------------------------------------------------ the harness

def make_driver(SimulateNetwork, ExaDisNet, DisNet, DisNode):
    """make_driver: build the comparing driver subclass

    Built inside a function for the same reason as test4_collision_mode's:
    its base class lives in pyexadis_base, which cannot be imported until
    the paths above are set and pyexadis has initialized.
    """

    class SimulateNetworkComparingTopology(SimulateNetwork):
        """run the exadis TopologySerial handler, and check pydis's against it

        exadis stays authoritative: its topology result is what the
        simulation carries forward.
        """

        def __init__(self, *args, pydis_topology=None, tag_state=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.pydis_topology = pydis_topology
            self.tag_state = tag_state
            self.records = []

        def step_update_response(self, N, state):
            """turn on UNZIP_STRESS at STRESS_STEP, matching
            tests/full_runs/03_binary_junction's two-phase run: unstressed
            junction formation, then stress pulling it back apart. This
            class is built directly on pyexadis_base.SimulateNetwork (there
            is no separate pydis-side SimulateNetwork run here to keep in
            step, unlike that test's cross-process comparison), so exadis'
            own istep convention applies with no offset.
            """
            state = super().step_update_response(N, state)
            if state['istep'] == STRESS_STEP:
                state["applied_stress"] = UNZIP_STRESS
            return state

        def seed_tag_state(self, G, maxindex, recycled):
            """seed_tag_state: give the replayed network exadis' tag pool

            Same logic as test4_collision_mode's method of the same name
            (see its docstring for the general reasoning and for
            live_tag_state()'s history); reproduced here rather than called
            on t4's driver instance, since it is a method of that other
            class. maxindex/recycled (top of stack first) are exadis' own
            live values for this step, from live_tag_state() below, reversed
            here because DisNet.get_new_tag() pops from the end of the list.
            """
            if maxindex is None:
                return
            live = {t[1] for t in G.all_nodes_tags()}
            free = [i for i in recycled if i not in live]
            G._max_tag = (0, max(maxindex, max(live)))
            G._recycled_tags = [(0, i) for i in reversed(free)]

        def live_tag_state(self, N):
            """live_tag_state: exadis' own current (maxindex, recycled) pool

            Straight from exadis' own SerialDisNet via the binding added in
            exadis commit 20ea2e8 ("Added binding to SerialDisNet internal
            tag indexing"). Replaces gen_tag_state.py's file as what
            actually seeds the replay; that file (loaded into
            self.tag_state) is kept only as the reference
            check_tag_state_consistency() compares this live call against,
            since unlike test4_collision_mode's recording (made by
            instrumenting exadis' C++ directly), gen_tag_state.py's file was
            itself only ever a Python-level reconstruction (tracking which
            tags disappear between snapshots), not exadis' actual internal
            state -- worth checking against the real thing now that this
            binding exists, rather than assumed to agree.
            """
            sn = N.get_disnet(ExaDisNet).net._get_serial_network()
            return int(sn._maxindex()), [int(i) for i in sn._recycled_indices()]

        def check_tag_state_consistency(self, rec, live_maxindex, live_recycled):
            """check_tag_state_consistency: live binding vs gen_tag_state.py's recording

            See live_tag_state()'s docstring for why the two might not
            agree: the stored recording is a Python-level reconstruction,
            not a direct reading of exadis' internal state. tag_state[i-1]
            indexing matches gen_tag_state.py's own RecordingDriver, which
            appends one entry per step at the same point in the pipeline
            (after collision, before topology) as this method is called
            from below.
            """
            j = rec.istep - 1
            if self.tag_state is None or not 0 <= j < len(self.tag_state):
                return
            ref_maxindex, ref_recycled = self.tag_state[j]
            rec.tag_state_match = (live_maxindex == ref_maxindex
                                   and live_recycled == ref_recycled)

        def step_topological_operations(self, N, state):
            rec = StepRecord(state.get('istep', len(self.records)))

            if self.cross_slip is not None:
                self.cross_slip.Handle(N, state)
            if self.collision is not None:
                self.collision.HandleCol(N, state)

            before = N.get_disnet(ExaDisNet).export_data()
            rec.n_nodes_before = len(before['nodes']['tags'])
            rec.n_segs_before = len(before['segs']['nodeids'])

            # exadis' own tag pool, at this same point (after collision,
            # before topology) gen_tag_state.py's recording was made -- see
            # seed_tag_state()'s and check_tag_state_consistency()'s
            # docstrings.
            live_maxindex, live_recycled = self.live_tag_state(N)
            self.check_tag_state_consistency(rec, live_maxindex, live_recycled)

            # exadis, authoritative
            self.topology.Handle(N, state)
            after = N.get_disnet(ExaDisNet).export_data()
            rec.exadis_dn = len(after['nodes']['tags']) - rec.n_nodes_before
            rec.exadis_dseg = len(after['segs']['nodeids']) - rec.n_segs_before
            rec.exadis_diff = t4.network_diff(before, after)
            rec.exadis_split = not t4.diff_is_empty(rec.exadis_diff)

            self._replay_pydis(rec, before, after, state, live_maxindex, live_recycled)
            self.records.append(rec)

            if self.remesh is not None:
                self.remesh.Remesh(N, state)

        def _replay_pydis(self, rec, before, after, state, live_maxindex, live_recycled):
            """_replay_pydis: run the pydis topology handler on the same input

            A normal pydis-driven simulation populates state['nodeforce_dict']
            and state['vel_dict'] in step_integrate, before topology ever
            runs. This harness is driven by exadis instead, which never
            touches those keys, so they have to be computed here, on the
            replayed copy, before Topology_Serial can read them.
            """
            G = DisNet()
            G.import_data(before)
            self.seed_tag_state(G, live_maxindex, live_recycled)
            DM = DisNetManager(G)

            try:
                state = state['_t5_pydis_force'].NodeForce(DM, state)
                state = state['_t5_pydis_mobility'].Mobility(DM, state)
                self.pydis_topology.Handle(DM, state)
            except Exception as err:                  # measured, not fatal
                import traceback
                rec.pydis_error = "%s: %s" % (type(err).__name__, err)
                rec.pydis_traceback = traceback.format_exc()
                return

            rec.pydis_dn = G.num_nodes() - rec.n_nodes_before
            rec.pydis_dseg = G.num_segments() - rec.n_segs_before
            after_pydis = G.export_data()
            rec.pydis_diff = t4.network_diff(before, after_pydis)
            rec.pydis_split = not t4.diff_is_empty(rec.pydis_diff)
            rec.tags_match = (t4.node_tag_set(after_pydis) == t4.node_tag_set(after))
            cell_h = before['cell']['h']
            is_periodic = before['cell']['is_periodic']
            rec.pos_diff = max_nearest_distance(
                np.asarray(after_pydis['nodes']['positions'], dtype=float),
                np.asarray(after['nodes']['positions'], dtype=float),
                np.asarray(cell_h, dtype=float))
            rec.same_geometry = (canonical_form_pbc(after_pydis, cell_h, is_periodic)
                                 == canonical_form_pbc(after, cell_h, is_periodic))
            rec.same_counts = (rec.pydis_dn == rec.exadis_dn
                               and rec.pydis_dseg == rec.exadis_dseg)

            if (rec.exadis_split or rec.pydis_split) and not rec.same_geometry:
                rec.classification = classify_split_disagreement(before, rec, state)

    return SimulateNetworkComparingTopology


def run(max_step=MAX_STEP, plot=False, print_freq=20, pydis_mode=PYDIS_MODE):
    """run: the comparing simulation, returning the per-step records"""
    import pyexadis
    from pyexadis_base import ExaDisNet, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
    from pyexadis_base import Collision, Remesh, Topology
    from pydis import DisNet, DisNode
    from pydis import Topology as PydisTopology
    from pydis import CalForce as PydisCalForce, MobilityLaw as PydisMobilityLaw
    from pydis.topology.topology_serial import TopologyParams

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    state = base_state()
    check_cutoff_maxseg(LBOX * np.eye(3), CUTOFF, state["maxseg"])
    net = init_two_disl_lines(state, pbc=False)

    calforce = CalForce(force_mode='CUTOFF_MODEL', state=state, Ec=EC, cutoff=CUTOFF)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint = TimeIntegration(integrator='EulerForward', dt=DT, state=state)
    # 'Proximity' on the exadis side is CollisionRetroactive under an alias
    # (plan_collision_mode.md section 2), and is what
    # test_binary_junction_exadis_elast.py itself asks for; matched here for
    # the same reason the rest of this config matches that script.
    collision = Collision(collision_mode='Proximity', state=state)
    topology = Topology(topology_mode='TopologySerial', state=state,
                        force=calforce, mobility=mobility)
    remesh = Remesh(remesh_rule='LengthBased', state=state)

    use_gp, enforce_gp = t4.glide_plane_settings(state)
    setup = {'use_glide_planes': use_gp, 'enforce_glide_planes': enforce_gp,
             'omp_num_threads': os.environ.get('OMP_NUM_THREADS')}
    # pydis' own use_glide_planes flag, read by the fix in topology_serial.py;
    # matched to whatever exadis actually resolved to rather than assumed,
    # for the same reason glide_plane_settings() itself reads back exadis'
    # Params instead of inferring from the absence of a 'crystal' key.
    state['use_glide_planes'] = bool(use_gp)

    # a pydis force/mobility pair for classify_split_disagreement's replay,
    # separate from the ones driving the authoritative exadis simulation.
    # No pydis Cell/CellList is needed anywhere in this harness: unlike
    # collision, Topology does no neighbour search of its own, it only
    # processes nodes the network already flags as multi-arm.
    pydis_force = PydisCalForce(force_mode='Elasticity_SBA', state=state,
                                Ec=EC, cutoff=CUTOFF)
    # vmax=1e15 is REQUIRED for this test to pass -- without it, step 1
    # reports a spurious mirror-symmetric "tie" that isn't real. PydiS'
    # MobilityLaw defaults to vmax=1e9, and the raw (unclamped) velocities
    # at the split trial are of order 1e9-1e10, well above that default. A
    # velocity cap clamps a vector's *magnitude*, so once both candidate
    # nodes' raw speeds exceed vmax, both get rescaled to exactly vmax:
    # vd1 == vd2 == vmax**2 == 1e18 to the last bit, regardless of what the
    # actual forces were. That exact-looking tie was chased for a long time
    # as if it were physical (.plan/2026-08-21/debug_topology_stage2.md) --
    # it was this clamp. ExaDiS' own GLIDE mobility has no such cap, so
    # leaving pydis at its default silently compared a clamped quantity
    # against an unclamped one. Matches the vmax already used in
    # examples/03_binary_junction/test_binary_junction_pydis_elast.py, for
    # the same reason (see the comment there).
    pydis_mobility = PydisMobilityLaw(mobility_law='SimpleGlide', state=state,
                                      vmax=1.0e15)
    state['_t5_pydis_force'] = pydis_force
    state['_t5_pydis_mobility'] = pydis_mobility
    state['_t5_topology_params'] = TopologyParams(state)

    pydis_topology = PydisTopology(split_mode=pydis_mode, state=state,
                                   force=pydis_force, mobility=pydis_mobility)

    tag_state = t4.load_tag_state(path=TAG_STATE)
    Driver = make_driver(SimulateNetwork, ExaDisNet, DisNet, DisNode)
    sim = Driver(calforce=calforce, mobility=mobility, timeint=timeint,
                collision=collision, topology=topology, remesh=remesh,
                vis=VisualizeNetwork() if plot else None,
                state=state, max_step=max_step, loading_mode='stress',
                applied_stress=np.zeros(6),
                print_freq=print_freq, plot_freq=10 if plot else None,
                plot_pause_seconds=0.01, write_freq=None,
                write_dir=str(OUT_DIR),
                pydis_topology=pydis_topology, tag_state=tag_state)
    sim.run(net, state)
    return sim.records, setup


# ------------------------------------------------------------------- reporting

def summarize(records):
    """summarize: turn the per-step records into the findings that matter"""
    n = len(records)
    if n == 0:
        print("no steps recorded")
        return False

    def count(pred):
        return sum(1 for r in records if pred(r))

    print("")
    print("=" * 70)
    print("MEASUREMENTS over %d steps" % n)
    print("=" * 70)

    ex = count(lambda r: r.exadis_split)
    py = count(lambda r: r.pydis_split)
    print("")
    print("-- what the two codes did")
    print("   exadis split a node        : %d step(s)" % ex)
    print("   pydis  split a node        : %d step(s)" % py)
    print("   both agreed on the counts  : %d/%d" % (count(lambda r: r.same_counts), n))
    print("   resulting tags identical   : %d/%d" % (count(lambda r: r.tags_match), n))
    print("   networks identical up to   : %d/%d" % (count(lambda r: r.same_geometry), n))
    print("     relabelling")
    errs = [r for r in records if r.pydis_error]
    print("   pydis raised               : %d step(s)" % len(errs))
    if errs:
        print("     first: step %d, %s" % (errs[0].istep, errs[0].pydis_error))
        if errs[0].pydis_traceback:
            for line in errs[0].pydis_traceback.strip().splitlines()[-8:]:
                print("     | %s" % line)

    differing = [r for r in records
                if (r.exadis_split or r.pydis_split) and not r.same_geometry]
    ties = [r for r in differing if r.classification == 'tie']
    real = [r for r in differing if r.classification == 'real_disagreement']
    unclassified = [r for r in differing if r.classification is None]
    print("")
    print("-- disagreements, classified")
    print("   steps that differ           : %d" % len(differing))
    print("   classified as a known tie   : %d" % len(ties))
    print("   classified as real          : %d" % len(real))
    print("   not classifiable (shape)    : %d" % len(unclassified))

    if ex:
        print("")
        print("-- steps where exadis split a node")
        for r in records:
            if r.exadis_split:
                flag = "" if (r.tags_match and r.same_geometry) \
                    else "  ***** [%s]" % (r.classification or "unclassified")
                print(("   step %4d: exadis dn=%+d dseg=%+d | "
                      "pydis dn=%+d dseg=%+d | tags: %-5s geometry: %-5s%s"
                      % (r.istep, r.exadis_dn, r.exadis_dseg,
                         r.pydis_dn, r.pydis_dseg,
                         r.tags_match, r.same_geometry, flag)).rstrip())

    print("")
    print("=" * 70)
    print("FINDINGS")
    print("=" * 70)
    if ex == 0:
        print("  ExaDiS never split a node in %d steps, so this configuration does" % n)
        print("  not exercise topology. Raise the step count.")
    else:
        print("  ExaDiS split a node on %d of %d steps." % (ex, n))
    return True


def fmt_tag(tag):
    return "%d,%d" % tag


def fmt_pos(p):
    return "(%.1f %.1f %.1f)" % tuple(p)


def describe(diff, label):
    """describe: one code's decisions for one step, in tags

    Reused verbatim from test4_collision_mode's describe(); topology only
    ever adds a node and a segment (never removes one, unlike collision), but
    the general shape is the same.
    """
    return t4.describe(diff, label)


def diagnose(records):
    """diagnose: print what each code decided on every differing step"""
    differing = [r for r in records
                if (r.exadis_split or r.pydis_split) and not r.same_geometry]
    print("")
    print("=" * 70)
    print("DIAGNOSIS: %d step(s) where the two codes differ" % len(differing))
    print("=" * 70)
    for r in differing:
        print("")
        print("  step %d  (exadis dn=%+d dseg=%+d | pydis dn=%+d dseg=%+d)  [%s]"
              % (r.istep, r.exadis_dn, r.exadis_dseg, r.pydis_dn, r.pydis_dseg,
                 r.classification or "unclassified"))
        for line in describe(r.exadis_diff, 'exadis'):
            print(line)
        for line in describe(r.pydis_diff, 'pydis '):
            print(line)
        if r.pydis_margin is not None or r.exadis_margin is not None:
            print("      pydis  own margin : %s"
                  % ("%.3e" % r.pydis_margin if r.pydis_margin is not None else "n/a"))
            print("      exadis own margin : %s"
                  % ("%.3e" % r.exadis_margin if r.exadis_margin is not None else "n/a"))
        if r.split_tag is not None:
            print("      node %s's arms, in the order this replay enumerated them: %s"
                  % (fmt_tag(r.split_tag),
                     ", ".join(fmt_tag(t) for t in r.split_neighbor_order)))
            print("        exadis moved to its new node : %s"
                  % (", ".join(fmt_tag(t) for t in r.ex_arms) if r.ex_arms else "n/a"))
            print("        pydis  moved to its new node : %s"
                  % (", ".join(fmt_tag(t) for t in r.py_arms) if r.py_arms else "n/a"))
            if r.classification == 'tie':
                same_partition = (r.ex_arms is not None and r.py_arms is not None
                                  and set(r.ex_arms) == set(r.py_arms))
                if same_partition:
                    print("        same arms grouped either way -- codes disagree only on")
                    print("        WHICH of the two resulting nodes moves (mirror-image split")
                    print("        directions, both equally valid by symmetry)")
                else:
                    print("        different arms grouped by each code -- a coarser tie, between")
                    print("        which PAIR of arms splits off, not just which way it moves")
                print("        margin this thin is why: it is a near-exact numerical tie, so")
                print("        even the order arms are enumerated in (which depends on how this")
                print("        network was built/imported, not on physics) can decide it")


def check_setup(setup):
    """check_setup: the conditions the comparison is only meaningful under"""
    ok = report("setup: OMP_NUM_THREADS is 1", setup['omp_num_threads'] == '1')
    return bool(ok)


def check_agreement(records):
    """check_agreement: every splitting step must fully agree, or be a known tie

    Both tags_match and same_geometry are required. A split always mints a
    new node with a new tag each code allocates independently, so this
    depends on seed_tag_state() correctly handing PyDiS's replay ExaDiS's
    own recorded tag pool (TAG_STATE, generated by gen_tag_state.py) --
    without that, a step that only fails on the new node's tag number would
    not be a real bug, but with it in place there is no such excuse left:
    a tag mismatch is either a bug in the seeding or a real divergence in
    which node ends up with which identity, and either way it should fail
    loudly rather than being waved through.

    A step that split and disagreed (in tags, in geometry, or both), and
    was not classified as a known mirror-symmetric tie, is a failure.
    """
    split_steps = [r for r in records if r.exadis_split or r.pydis_split]
    bad = [r for r in split_steps
          if not (r.tags_match and r.same_geometry) and r.classification != 'tie']
    label = "every splitting step fully agrees (tags and geometry), or is a known tie"
    if bad:
        label += "  [differs at %s -- see DIAGNOSIS above for what each code did]" % ", ".join(
            "%d (%s)" % (r.istep, r.classification or "unclassified") for r in bad)
    ok = report(label, not bad and bool(split_steps))
    return bool(ok)


def check_position_agreement(records):
    """check_position_agreement: worst-case nodal position agreement, after topology

    Reports the largest max-nearest-node distance between pydis's replayed
    network and exadis' authoritative one, taken over every step (not just
    the steps where a split happened), against POS_DIFF_TOL. This is a
    coarser, absolute-distance check alongside same_geometry's exact
    (tol=1e-6, up to relabelling) one -- same style as
    tests/full_runs/03_binary_junction's own pydis-vs-exadis reporting.
    """
    diffs = [r for r in records if r.pos_diff is not None]
    worst = max((r.pos_diff for r in diffs), default=float('nan'))
    label = ("max nearest-node distance, pydis vs exadis, after topology = %.4e, "
             "tolerance = %.1e" % (worst, POS_DIFF_TOL))
    ok = report(label, bool(diffs) and worst < POS_DIFF_TOL)
    return bool(ok)


def check_tag_state_consistency(records):
    """check_tag_state_consistency: live exadis tag-pool binding vs gen_tag_state.py's recording

    See SimulateNetworkComparingTopology.check_tag_state_consistency for
    what is being compared and why the two sources might disagree. Steps
    the recording has no entry for (r.tag_state_match is None) are not
    counted either way.
    """
    checked = [r for r in records if r.tag_state_match is not None]
    bad = [r.istep for r in checked if not r.tag_state_match]
    label = ("live exadis tag pool (maxindex, recycled) matches the stored "
             "recording at every step it covers")
    if bad:
        label += "  [differs at %s]" % ", ".join(str(i) for i in bad)
    ok = report(label, not bad and bool(checked))
    return bool(ok)


def check_invariants(records):
    """check_invariants: what must hold whichever split mode runs"""
    n = len(records)
    ok = report("the pydis topology handler ran without raising at every step",
               not any(r.pydis_error for r in records))
    return bool(ok) and n > 0


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--max-step', type=int, default=MAX_STEP)
    parser.add_argument('--plot', action='store_true', default=False,
                        help='show the network while it runs, and let the '
                             'native output through')
    parser.add_argument('--print-freq', type=int, default=50)
    parser.add_argument('--diagnose', action='store_true', default=False,
                        help='print what each code decided on every '
                             'differing step')
    parser.add_argument('--pydis-mode', default=PYDIS_MODE,
                        choices=['MaxDiss', 'Serial'],
                        help='which pydis split mode to compare')
    args = parser.parse_args(argv)

    import pyexadis
    with quiet_native_output(enabled=not args.plot) as captured:
        pyexadis.initialize()
        records, setup = run(max_step=args.max_step, plot=args.plot,
                             print_freq=args.print_freq, pydis_mode=args.pydis_mode)
    for line in kokkos_summary(captured.text):
        print("  %s" % line)

    summarize(records)
    # Printed whenever anything differed, not only with --diagnose: a bare
    # "1 step differs" in the measurements above is not enough to tell a
    # known tie from a real bug without seeing what each code actually did.
    differing = [r for r in records
                if (r.exadis_split or r.pydis_split) and not r.same_geometry]
    if args.diagnose or differing:
        diagnose(records)
    print("")
    ok = check_setup(setup)
    ok &= check_invariants(records)
    ok &= check_agreement(records)
    ok &= check_position_agreement(records)
    ok &= check_tag_state_consistency(records)
    print("")
    tag = GREEN + "PASSED" + RESET if ok else RED + "FAILED" + RESET
    print("test_topology_mode_pydis_exadis: %s" % tag)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
