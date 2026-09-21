"""Collision handling compared step by step, PyDiS against ExaDiS.

The ExaDiS simulation is authoritative: its collision is the one whose
result the run carries forward. At every step the network is snapshotted
before collision, ExaDiS collides it, and the same snapshot is replayed
through a PyDiS collision on a throwaway DisNet. A disagreement at step n
therefore does not change the configuration seen at step n+1, and the
per-step comparisons stay independent of each other.

WHAT THIS MEASURES, AND WHAT IT DOES NOT

PyDiS has no retroactive collision yet, so the PyDiS slot holds the
existing 'Proximity' rule, which is a different and much simpler
algorithm: it tests the current minimum distance between two segments
and ignores both velocity and old positions. ExaDiS 'Retroactive' tests
whether two swept segments came within rann during the last step, or
will during the next.

So the two are NOT expected to agree, and this script does not assert
that they do. It reports how far apart they are, and it measures the
things a retroactive implementation would depend on:

  - are old nodal positions present when collision runs, and do they
    still line up with the network
  - is the timestep the one actually taken
  - does every node have a velocity at that moment
  - does segment order survive export_data/import_data
  - do the two codes allocate the same tags
  - how often ExaDiS collides at all on this configuration

Those measurements are the point. The pass/fail status below covers only
the invariants that must hold regardless of which collision rule runs.

RUNNING IT

    make                       run at the default step count
    make PYTHON=$CONDA_PREFIX/bin/python3
    python3 test_collision_mode_pydis_exadis.py --max-step 50
    python3 test_collision_mode_pydis_exadis.py --plot

Single-threaded by construction: OMP_NUM_THREADS is forced to 1 before
pyexadis is imported, because multi-threaded ExaDiS varies run to run
through the order of floating-point reductions, and over hundreds of
steps that grows into visibly different trajectories. check_setup reports
the thread count Kokkos actually took, rather than the variable this file
sets, so the pin cannot appear to hold while having no effect.
"""

import os
import re
# must precede the pyexadis import, and therefore any framework import
# that might pull it in.
#
# Assigned rather than setdefault: an OMP_NUM_THREADS inherited from the
# environment would silently reintroduce the run-to-run variation this test
# exists to compare against, and over hundreds of steps that grows into
# visibly different trajectories. Overridden here rather than merely
# defaulted, so the comparison cannot be made flaky from outside.
os.environ['OMP_NUM_THREADS'] = '1'

import sys
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if p not in sys.path]

import numpy as np

from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.simulation_setup import remesh_initial_config
# exadis' arm order is not carried by export_data, and it breaks ties
# between near-degenerate glide planes, so the replay has to be given it
# explicitly, by apply_arm_order below. Same argument as seed_tag_state
# below; see framework.arm_order.
from framework.arm_order import export_arm_order, apply_arm_order
from framework.testing import (report, report_close, verdict,
                               quiet_native_output, kokkos_summary)

# Tolerance on the paired node positions. Zero: the two codes agree on every
# colliding step to the last bit, so report_close tags the line [BITWISE].
# It was 6.6159e-12 when the rule was first complete, and five changes closed
# that: the Newton solve through an explicit inverse rather than LAPACK,
# closest_image by box-edge subtraction, a frame fix on the reference side, and
# two grouping switches, SEGMENT_INTERP in collision_retroactive.py and
# PTSEG_RATIO in swept_distance.py, each defaulting to the form the reference
# writes. Flipping either switch reopens a residual of order 1e-14, which this
# tolerance then catches.
#
# The frame fix is a local, uncommitted change inside core/exadis, so a
# `git submodule update` will make this assertion fail. That is deliberate: the
# alternative was a tolerance loose enough to hide it.
#
# Exactness here rests on the single thread forced at the top of this file, as
# in test3_remesh_rule: each step hands pydis exadis' own snapshot, so exadis'
# trajectory does not enter the comparison, but its collision arithmetic is not
# thread-invariant. check_setup below checks the count Kokkos actually took,
# rather than leaving a later failure to be puzzled over.
ATOL_POS = 0.0

# ---------------------------------------------------------------- configuration

# The base configuration is the one test3_remesh_rule uses, so that the
# two tests exercise the same trajectory at different points in the step.
# 600 steps under load, with no zero-stress phase: the expanding loop
# meets its own periodic images and collides. Collisions come in two
# bursts, 260-302 and 529-569, 39 steps in all; nothing happens between
# them, which is why 500 steps is a poor stopping point and 600 is not.
#
# ref_data/exadis_tag_state.npz must cover the whole range, or the steps
# past its end lose their tag comparison; it is keyed by step index.
MAX_STEP = 600
LBOX = 1000.0
CUTOFF = 0.25 * LBOX
DT = 1.0e-8
APPLIED_STRESS = np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0])

# Ec=0.0 disables the core-energy term, which PyDiS Elasticity_SBA has no
# equivalent for. Deliberate, and documented in the example this is built
# from; not to be "cleaned up".
EC = 0.0

# cells per side for the PyDiS neighbour list, matching the pydis example
NDIV = [8, 8, 8]

# Which pydis rule occupies the pydis slot. 'Proximity' is the measuring
# baseline: a different algorithm, so disagreement is expected and nothing
# is asserted on it. 'Retroactive' is the rule under test, which should
# agree with exadis step by step.
PYDIS_MODE = 'Retroactive'

# SimulateNetwork.run writes a stress/strain summary at the end whatever
# write_freq says, so it needs somewhere to put it. Next to this script
# rather than relative to the working directory, so the test behaves the
# same under ctest, which runs it from the build tree.
OUT_DIR = Path(__file__).resolve().parent / 'output'
REF_DIR = Path(__file__).resolve().parent / 'ref_data'
TAG_STATE = REF_DIR / 'exadis_tag_state.npz'


def base_state():
    """base_state: the material and discretization parameters

    No 'crystal' key, which is what leaves use_glide_planes and
    enforce_glide_planes at 0 in ExaDiS. That is asserted at setup rather
    than assumed: setting a crystal type would silently move the
    comparison onto a branch PyDiS has no counterpart for.
    """
    return {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
            "maxseg": 0.04 * LBOX, "minseg": 0.01 * LBOX, "rann": 3.0}


def init_frank_read_src_loop(state, pbc=True):
    """init_frank_read_src_loop: the Frank-Read source, as the example builds it

    Conditioned by remesh_initial_config so no segment starts longer than
    maxseg; the inserted nodes are PINNED, leaving the geometry and the
    boundary condition unchanged.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet, NodeConstraints

    arm_length = 0.125 * LBOX
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[pbc, pbc, pbc])
    burg_vec = np.array([1.0, 0.0, 0.0])

    PINNED = NodeConstraints.PINNED_NODE
    FREE = NodeConstraints.UNCONSTRAINED
    rn = np.array([[0.0, -arm_length / 2.0, 0.0, PINNED],
                   [0.0, 0.0, 0.0, FREE],
                   [0.0, arm_length / 2.0, 0.0, PINNED],
                   [0.0, arm_length / 2.0, -arm_length, PINNED],
                   [0.0, -arm_length / 2.0, -arm_length, PINNED]])
    rn[:, 0:3] += cell.center()

    n = rn.shape[0]
    links = np.zeros((n, 8))
    for i in range(n):
        pn = np.cross(burg_vec, rn[(i + 1) % n, :3] - rn[i, :3])
        pn = pn / np.linalg.norm(pn)
        links[i, :] = np.concatenate(([i, (i + 1) % n], burg_vec, pn))

    rn, links = remesh_initial_config(rn, links, state["maxseg"],
                                      LBOX * np.eye(3), [pbc] * 3)
    return DisNetManager(ExaDisNet(cell, rn, links))


# Registered rather than hard-coded, so that adding the binary junction
# configuration later is an entry here and not a rewrite of the harness.
CONFIGS = {'frank_read_elast': init_frank_read_src_loop}


# ------------------------------------------------------------------- the record

class StepRecord:
    """StepRecord: what one step of the comparison observed

    Plain data. Collected for every step and summarized at the end,
    rather than asserted step by step, so that one run produces the whole
    picture instead of stopping at the first difference.
    """

    def __init__(self, istep):
        self.istep = istep
        self.dt = None
        self.n_nodes_before = 0
        self.n_segs_before = 0
        # old positions
        self.xold_present = False
        self.xold_len = None
        self.xold_aligned = None
        self.xold_tags_cover_network = None
        # velocities
        self.vel_present = False
        self.vel_missing = None
        self.vel_dead = None
        # what each code did
        self.exadis_dn = 0
        self.exadis_dseg = 0
        self.pydis_dn = 0
        self.pydis_dseg = 0
        self.exadis_collided = False
        self.pydis_collided = False
        # tag pool: live exadis binding vs the stored reference recording
        self.tag_state_match = None
        # comparison
        self.tags_match = None
        self.same_geometry = None
        self.max_pos_diff = None
        self.match_one_to_one = None
        self.same_counts = None
        self.pydis_error = None
        self.pydis_traceback = None
        # decision-level diagnosis
        self.exadis_diff = None
        self.pydis_diff = None
        self.pydis_events = []
        # round trip
        self.segorder_preserved = None


def load_tag_state(path=TAG_STATE):
    """load_tag_state: exadis' recorded tag recycling state, one entry per step

    exadis carries maxindex and a stack of freed indices from one step to
    the next. The pydis side of this comparison is rebuilt with
    import_data every step, so it starts each step with an empty pool and
    counts up from the highest tag present instead of reusing a freed one.
    That is a property of how the comparison is built, not of either
    collision rule, and without it the resulting-tags count cannot match
    however correct pydis is.

    Originally the only source for this: a recording made once by printing
    network->maxindex and network->recycled_indices at the top of
    CollisionRetroactive::retroactive_collision_parallel (manual C++
    instrumentation, not reproducible from Python). As of exadis commit
    20ea2e8 ("Added binding to SerialDisNet internal tag indexing"),
    SerialDisNet exposes this live: G.net._get_serial_network()._maxindex()
    and ._recycled_indices() (top of stack first, matching this recording's
    own convention). seed_tag_state() below now seeds pydis' replay from
    that live call instead of this file; this file, and this function, are
    kept only as the reference check_tag_state_consistency() compares the
    live binding against on every step, so a regression in either the
    binding or the historical recording is caught rather than silently
    trusted. Regenerate the file (same manual C++ instrumentation, or now
    also possible by recording the live binding's own output) if the
    reference trajectory changes; it is keyed by step index and silently
    unused if absent.

    Returns a list, one entry per step, of (maxindex, [freed indices]) with
    the stack top first, or None when the file is not there.
    """
    if not TAG_STATE.exists():
        return None
    d = np.load(str(path))
    maxindex, flat, off = d['maxindex'], d['recycled_flat'], d['recycled_offsets']
    return [(int(maxindex[i]), [int(v) for v in flat[off[i]:off[i + 1]]])
            for i in range(len(maxindex))]


def node_tag_set(data):
    """node_tag_set: the set of node tags in an export_data payload"""
    return {tuple(int(v) for v in t) for t in data['nodes']['tags']}


def seg_tag_set(data):
    """seg_tag_set: segments as unordered pairs of node tags"""
    return {frozenset(pair) for pair in seg_tag_pairs(data)}


def node_positions(data):
    """node_positions: tag -> position, from an export_data payload"""
    tags = [tuple(int(v) for v in t) for t in data['nodes']['tags']]
    return dict(zip(tags, np.asarray(data['nodes']['positions'], float)))


def match_nodes(data_a, data_b):
    """match_nodes: pair the two networks' nodes by position, ignoring tags

    Returns (gap, one_to_one). gap is the largest distance between paired
    nodes, taken both ways round so a node present in one network and not the
    other cannot hide. one_to_one is False when two nodes of one network claim
    the same partner, in which case the pairing is not an assignment and gap
    understates the disagreement.

    Keyed on position rather than on tags, because the two codes can give the
    same physical node different tags: the run reports steps that "differ only
    in tag allocation, same geometry", and a tag-keyed comparison would call
    those infinitely far apart. Tag agreement is checked separately by
    tags_match, which is the right place for it.

    Nearest neighbour rather than an optimal assignment. The two are the same
    thing while the networks nearly coincide, which is the regime this test
    runs in, and one_to_one is what reports the case where they stop being.

    Distances are minimum-image; this configuration is periodic in all three
    directions.
    """
    ra = np.asarray(data_a['nodes']['positions'], float)
    rb = np.asarray(data_b['nodes']['positions'], float)
    if ra.shape[0] != rb.shape[0] or ra.size == 0:
        return np.inf, False
    h = np.asarray(data_a['cell']['h'], float)
    d = ra[:, None, :] - rb[None, :, :]
    s = d @ np.linalg.inv(h).T
    s -= np.rint(s)
    dist = np.linalg.norm(s @ h.T, axis=-1)

    fwd, rev = dist.argmin(axis=1), dist.argmin(axis=0)
    one_to_one = (len(set(fwd.tolist())) == fwd.size
                  and len(set(rev.tolist())) == rev.size)
    gap = max(dist[np.arange(fwd.size), fwd].max(),
              dist[rev, np.arange(rev.size)].max())
    return float(gap), bool(one_to_one)


def canonical_form(data, tol=1e-6):
    """canonical_form: the network up to relabelling

    Node positions sorted, and segments as sorted pairs of indices into
    that sorted list, so two networks that describe the same geometry and
    topology compare equal whatever tags they gave their nodes. Positions
    are rounded first, since the two codes reach them by the same
    arithmetic on the same inputs and should agree to a few ULP.

    This is the relabelling-invariant comparison the plan asks for; it
    lives here for now and belongs in framework/testing.py once a second
    caller wants it.
    """
    pos = np.asarray(data['nodes']['positions'], float)
    key = np.round(pos / tol).astype(np.int64)
    order = sorted(range(len(key)), key=lambda i: tuple(key[i]))
    rank = {old_i: new_i for new_i, old_i in enumerate(order)}
    nodes = tuple(tuple(key[i]) for i in order)
    segs = tuple(sorted(tuple(sorted((rank[int(a)], rank[int(b)])))
                        for a, b in data['segs']['nodeids']))
    return nodes, segs


def network_diff(before, after):
    """network_diff: what changed, in tags rather than indices

    ExaDiS reports nothing about what its collision decided, so its
    decisions have to be reconstructed by comparing the network before
    and after. This is the raw material: which nodes and segments went,
    which arrived, and where the new ones are. Collision only ever splits
    a segment or merges two nodes, so this is enough to say what it did.
    """
    pos_before, pos_after = node_positions(before), node_positions(after)
    segs_before, segs_after = seg_tag_set(before), seg_tag_set(after)
    removed = sorted(set(pos_before) - set(pos_after))
    added = sorted(set(pos_after) - set(pos_before))
    # a node that survived but was relocated: collision moves the merge
    # survivor to the collision point, so this is where a disagreement
    # about that point shows up with no change of topology at all
    moved = {t: (pos_before[t], pos_after[t])
             for t in set(pos_before) & set(pos_after)
             if np.linalg.norm(pos_after[t] - pos_before[t]) > 1e-9}
    return {
        'moved_nodes': moved,
        'removed_nodes': removed,
        'added_nodes': added,
        'removed_node_pos': {t: pos_before[t] for t in removed},
        'added_node_pos': {t: pos_after[t] for t in added},
        'removed_segs': sorted(tuple(sorted(s)) for s in segs_before - segs_after),
        'added_segs': sorted(tuple(sorted(s)) for s in segs_after - segs_before),
    }


def diff_is_empty(diff):
    return not (diff['removed_nodes'] or diff['added_nodes']
                or diff['removed_segs'] or diff['added_segs']
                or diff['moved_nodes'])


def seg_tag_pairs(data):
    """seg_tag_pairs: segments as ordered pairs of node tags

    Ordered as stored, so that comparing these lists position by position
    answers whether segment order survived a round trip; compare the sets
    instead to ask only whether the same segments are present.
    """
    tags = [tuple(int(v) for v in t) for t in data['nodes']['tags']]
    return [(tags[int(a)], tags[int(b)]) for a, b in data['segs']['nodeids']]


# ------------------------------------------------------------------ the harness

def make_driver(SimulateNetwork, ExaDisNet, DisNet, DisNode):
    """make_driver: build the comparing driver subclass

    Built inside a function because its base class lives in pyexadis_base,
    which cannot be imported until the paths above are set and pyexadis
    has initialized.
    """

    class SimulateNetworkComparingCollision(SimulateNetwork):
        """run the exadis collision, and check a pydis collision against it

        exadis stays authoritative: its collision is the one whose result
        the simulation carries forward.
        """

        def __init__(self, *args, pydis_collision=None, tag_state=None,
                     **kwargs):
            super().__init__(*args, **kwargs)
            self.pydis_collision = pydis_collision
            self.tag_state = tag_state
            self.records = []

        def seed_tag_state(self, G, maxindex, recycled):
            """seed_tag_state: give the replayed network exadis' tag pool

            import_data leaves a fresh network counting up from the highest
            tag it was handed, with nothing recycled, while exadis reaches
            this step holding indices it freed on earlier ones. Replaying
            without them makes every new tag differ, which shows up as a
            tag mismatch on steps whose geometry agrees exactly.

            maxindex and recycled (top of stack first) are exadis' own live
            values for this step, from live_tag_state() below -- see that
            function for where this used to come from.

            pydis pops from the end of its list and exadis pops the top of
            its stack, so recycled, which arrives top first, is reversed
            here.

            Indices still in use are dropped as a guard only. Correctly
            aligned the pool never names a live node, so this should never
            fire; it is kept because popping such an entry would hand
            get_new_tag a tag the network already has, corrupting the
            replay rather than merely relabelling it.
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
            exadis commit 20ea2e8 -- see load_tag_state()'s docstring for
            the history this replaces. recycled is returned top of stack
            first, the same convention load_tag_state()'s file uses, so the
            two are directly comparable in check_tag_state_consistency
            without reordering either one.
            """
            sn = N.get_disnet(ExaDisNet).net._get_serial_network()
            return int(sn._maxindex()), [int(i) for i in sn._recycled_indices()]

        def check_tag_state_consistency(self, rec, live_maxindex, live_recycled):
            """check_tag_state_consistency: live binding vs the stored recording

            self.tag_state (loaded from ref_data/exadis_tag_state.npz) is no
            longer what seeds the replay -- see seed_tag_state()'s docstring
            -- but it is still exadis' own tag pool state at this same point
            in an earlier, presumably-identical run of this trajectory. If
            live and recorded ever disagree, either the live binding is
            wrong, the recording is stale (the reference trajectory
            changed and load_tag_state()'s file needs regenerating), or the
            run has genuinely stopped being reproducible -- any of which is
            worth surfacing rather than silently trusting one source.
            Leaves rec.tag_state_match as None (not counted either way) on
            a step the recording has no entry for.
            """
            j = rec.istep - 1
            if self.tag_state is None or not 0 <= j < len(self.tag_state):
                return
            ref_maxindex, ref_recycled = self.tag_state[j]
            rec.tag_state_match = (live_maxindex == ref_maxindex
                                   and live_recycled == ref_recycled)

        def step_topological_operations(self, N, state):
            rec = StepRecord(state.get('istep', len(self.records)))
            rec.dt = state.get('dt')

            if self.cross_slip is not None:
                self.cross_slip.Handle(N, state)

            before = N.get_disnet(ExaDisNet).export_data()
            rec.n_nodes_before = len(before['nodes']['tags'])
            rec.n_segs_before = len(before['segs']['nodeids'])

            # exadis' own tag pool, at this same point (top of the collision
            # handler) the old recording was made -- see seed_tag_state()'s
            # and check_tag_state_consistency()'s docstrings.
            live_maxindex, live_recycled = self.live_tag_state(N)
            self.check_tag_state_consistency(rec, live_maxindex, live_recycled)
            # before HandleCol, not after: exadis updates conn in place as it
            # splits and merges, so this has to be the pre-collision order
            live_arms = export_arm_order(N)

            self._measure_inputs(rec, before, state)

            # exadis, authoritative
            self.collision.HandleCol(N, state)
            after = N.get_disnet(ExaDisNet).export_data()
            rec.exadis_dn = len(after['nodes']['tags']) - rec.n_nodes_before
            rec.exadis_dseg = len(after['segs']['nodeids']) - rec.n_segs_before
            rec.exadis_diff = network_diff(before, after)
            rec.exadis_collided = not diff_is_empty(rec.exadis_diff)

            # pydis, on a throwaway copy of the same input
            self._replay_pydis(rec, before, after, state, live_maxindex,
                               live_recycled, live_arms)

            self.records.append(rec)

            if self.topology is not None:
                self.topology.Handle(N, state)
            if self.remesh is not None:
                self.remesh.Remesh(N, state)

        def _measure_inputs(self, rec, before, state):
            """_measure_inputs: are old positions and velocities usable here

            This is the measurement the whole first run exists for. A
            retroactive rule is defined by comparing current positions
            against positions from the start of the step, so if those are
            absent, misaligned, or not mappable onto the network, the
            design has to change before any of it is written.
            """
            n_nodes = rec.n_nodes_before
            tags_now = node_tag_set(before)

            oldnodes = state.get('oldnodes_dict')
            rec.xold_present = oldnodes is not None
            if rec.xold_present:
                xold = oldnodes['positions']
                rec.xold_len = len(xold)
                # exadis passes xold positionally and checks only length
                rec.xold_aligned = (rec.xold_len == n_nodes)
                if 'tags' in oldnodes:
                    old_tags = {tuple(int(v) for v in t)
                                for t in oldnodes['tags']}
                    # every node the rule must handle needs an old position
                    rec.xold_tags_cover_network = tags_now.issubset(old_tags)

            veltags = state.get('nodeveltags')
            rec.vel_present = veltags is not None
            if rec.vel_present:
                vt = {tuple(int(v) for v in t) for t in np.asarray(veltags)}
                rec.vel_missing = len(tags_now - vt)
                rec.vel_dead = len(vt - tags_now)

        def _replay_pydis(self, rec, before, after, state, live_maxindex,
                          live_recycled, live_arms):
            """_replay_pydis: run the pydis collision on the same input"""
            G = DisNet()
            G.import_data(before)
            self.seed_tag_state(G, live_maxindex, live_recycled)
            apply_arm_order(G, live_arms)

            round_trip = G.export_data()
            rec.segorder_preserved = (seg_tag_pairs(round_trip)
                                      == seg_tag_pairs(before))

            # The Proximity rule reads a per-node flag dict that only
            # pydis' Topology creates, and this run is driven by exadis
            # with topology=None. Supplying it here keeps the measurement
            # possible without changing pydis: every node starts clear,
            # which is what init_topology_exemptions does.
            state['nodeflag_dict'] = {tag: DisNode.Flags.CLEAR
                                      for tag in G.all_nodes_tags()}

            try:
                self.pydis_collision.HandleCol(DisNetManager(G), state)
            except Exception as err:                  # measured, not fatal
                import traceback
                rec.pydis_error = "%s: %s" % (type(err).__name__, err)
                rec.pydis_traceback = traceback.format_exc()
                return

            rec.pydis_dn = G.num_nodes() - rec.n_nodes_before
            rec.pydis_dseg = G.num_segments() - rec.n_segs_before
            after_pydis = G.export_data()
            rec.pydis_diff = network_diff(before, after_pydis)
            rec.pydis_collided = not diff_is_empty(rec.pydis_diff)
            records = getattr(self.pydis_collision, 'records', None)
            if records:
                rec.pydis_events = records[-1].events
            rec.tags_match = (node_tag_set(after_pydis) == node_tag_set(after))
            rec.same_geometry = (canonical_form(after_pydis)
                                 == canonical_form(after))
            rec.max_pos_diff, rec.match_one_to_one = match_nodes(after_pydis,
                                                                 after)
            rec.same_counts = (rec.pydis_dn == rec.exadis_dn
                               and rec.pydis_dseg == rec.exadis_dseg)

    return SimulateNetworkComparingCollision


def run(config='frank_read_elast', max_step=MAX_STEP, plot=False,
        print_freq=50, pydis_mode=PYDIS_MODE):
    """run: the comparing simulation, returning the per-step records"""
    import pyexadis
    from pyexadis_base import ExaDisNet, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
    from pyexadis_base import Collision, Remesh
    from pydis import DisNet, DisNode, CellList
    from pydis import Cell as PydisCell
    from pydis import Collision as PydisCollision

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    state = base_state()
    check_cutoff_maxseg(LBOX * np.eye(3), CUTOFF, state["maxseg"])
    net = CONFIGS[config](state)

    calforce = CalForce(force_mode='CUTOFF_MODEL', state=state,
                        Ec=EC, cutoff=CUTOFF)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint = TimeIntegration(integrator='EulerForward', dt=DT, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    remesh = Remesh(remesh_rule='LengthBased', state=state)

    # The setup checks decide which branch of exadis the comparison is
    # actually against, so they are collected as data and reported by the
    # caller: this function runs inside the native-output filter, and
    # anything printed here would be swallowed with the Kokkos banner.
    use_gp, enforce_gp = glide_plane_settings(state)
    setup = {'use_glide_planes': use_gp, 'enforce_glide_planes': enforce_gp}

    # The neighbour list must be given a *pydis* Cell, not the exadis one
    # that net.cell returns here. Both expose is_periodic, but exadis' is a
    # method and pydis' is a list, so passing the wrong one fails inside
    # sort_points_to_list with an opaque TypeError rather than at
    # construction. Building it from the exported cell payload is exactly
    # what DisNet.import_data does, so the list and the replayed network
    # are guaranteed to agree.
    cell_data = net.get_disnet(ExaDisNet).export_data()['cell']
    pydis_cell = PydisCell(h=cell_data['h'], origin=cell_data['origin'],
                           is_periodic=cell_data['is_periodic'])
    nbrlist = CellList(cell=pydis_cell, n_div=NDIV)
    pydis_collision = PydisCollision(collision_mode=pydis_mode, state=state,
                                     nbrlist=nbrlist, collision_record=True)

    Driver = make_driver(SimulateNetwork, ExaDisNet, DisNet, DisNode)
    tag_state = load_tag_state()
    sim = Driver(calforce=calforce, mobility=mobility, timeint=timeint,
                 collision=collision, topology=None, remesh=remesh,
                 vis=VisualizeNetwork() if plot else None,
                 state=state, max_step=max_step, loading_mode='stress',
                 applied_stress=APPLIED_STRESS,
                 print_freq=print_freq, plot_freq=10 if plot else None,
                 plot_pause_seconds=0.01, write_freq=None,
                 write_dir=str(OUT_DIR),
                 pydis_collision=pydis_collision, tag_state=tag_state)
    sim.run(net, state)
    return sim.records, setup


def glide_plane_settings(state):
    """glide_plane_settings: the (use, enforce) exadis will actually apply

    Read back out of the Params object rather than inferred from the
    absence of a 'crystal' key, so that a change in exadis' defaults is
    caught here rather than silently comparing the wrong branch.

    The flags are tri-state on Params: -1 means "not set, resolve from the
    crystal type", and the resolution happens later, in Crystal::initialize,
    not on Params. So reading them raw gives -1 here and -1 is truthy,
    which would make a naive test report glide planes as enforced when
    they are not. The resolution is reproduced rather than assumed:

        type == -1                  -> use = enforce = 0
        use == -1                   -> use = 0 for bcc, else 1
        enforce == -1               -> enforce = use
    """
    from pyexadis_base import get_exadis_params
    cp = get_exadis_params(state).crystalparams
    ctype, use, enforce = cp.type, cp.use_glide_planes, cp.enforce_glide_planes
    if ctype == -1:
        return 0, 0
    if use < 0:
        use = 0 if ctype == 1 else 1          # 1 == BCC_CRYSTAL
    if enforce < 0:
        enforce = use
    return use, enforce


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

    print("")
    print("-- old positions (what a retroactive rule needs)")
    ok_present = count(lambda r: r.xold_present)
    ok_aligned = count(lambda r: r.xold_aligned)
    ok_cover = count(lambda r: r.xold_tags_cover_network)
    print("   present in state           : %d/%d" % (ok_present, n))
    print("   length == node count       : %d/%d" % (ok_aligned, n))
    print("   tags cover the network     : %d/%d" % (ok_cover, n))

    print("")
    print("-- timestep and velocities")
    dts = sorted({r.dt for r in records if r.dt is not None})
    print("   dt present                 : %d/%d" % (
        count(lambda r: r.dt is not None), n))
    if dts:
        print("   dt range                   : %.4e .. %.4e" % (dts[0], dts[-1]))
    print("   velocities present         : %d/%d" % (
        count(lambda r: r.vel_present), n))
    print("   steps with a missing vel   : %d" % count(
        lambda r: r.vel_missing), )
    print("   steps with a dead vel tag  : %d" % count(lambda r: r.vel_dead))

    print("")
    print("-- network transfer")
    print("   segment order preserved    : %d/%d" % (
        count(lambda r: r.segorder_preserved), n))

    print("")
    print("-- what the two codes did")
    ex = count(lambda r: r.exadis_collided)
    py = count(lambda r: r.pydis_collided)
    print("   exadis changed the network : %d step(s)" % ex)
    print("   pydis  changed the network : %d step(s)" % py)
    print("   both agreed on the counts  : %d/%d" % (
        count(lambda r: r.same_counts), n))
    print("   resulting tags identical   : %d/%d" % (
        count(lambda r: r.tags_match), n))
    print("   networks identical up to   : %d/%d" % (
        count(lambda r: r.same_geometry), n))
    print("     relabelling")
    errs = [r for r in records if r.pydis_error]
    print("   pydis raised               : %d step(s)" % len(errs))
    if errs:
        print("     first: step %d, %s" % (errs[0].istep, errs[0].pydis_error))
        if errs[0].pydis_traceback:
            for line in errs[0].pydis_traceback.strip().splitlines()[-8:]:
                print("     | %s" % line)

    if ex:
        print("")
        print("-- steps where exadis collided")
        for r in records:
            if r.exadis_collided:
                # marked so a disagreeing step can be found by eye in a
                # long listing, and by grep
                flag = "" if (r.tags_match and r.same_geometry) else "  *****"
                print(("   step %4d: exadis dn=%+d dseg=%+d | "
                       "pydis dn=%+d dseg=%+d | tags: %-5s geometry: %-5s%s"
                       % (r.istep, r.exadis_dn, r.exadis_dseg,
                          r.pydis_dn, r.pydis_dseg,
                          r.tags_match, r.same_geometry, flag)).rstrip())

    print("")
    print("=" * 70)
    print("FINDINGS")
    print("=" * 70)
    if ok_present == n and ok_aligned == n:
        print("  Old positions are present and aligned at every step, so a")
        print("  retroactive rule can be built on them.")
    else:
        print("  Old positions are NOT reliably available. A retroactive rule")
        print("  cannot be built on state['oldnodes_dict'] as it stands; the")
        print("  harness would have to snapshot positions itself.")
    if ex == 0:
        print("  ExaDiS never collided in %d steps, so this configuration does"
              % n)
        print("  not exercise collision. Raise the step count.")
    else:
        print("  ExaDiS collided on %d of %d steps." % (ex, n))
    return True


def fmt_tag(tag):
    return "%d,%d" % tag


def fmt_pos(p):
    return "(%.1f %.1f %.1f)" % tuple(p)


def describe(diff, label):
    """describe: one code's decisions for one step, in tags"""
    lines = []
    if diff['removed_nodes']:
        lines.append("      %s removed nodes : %s" % (
            label, ", ".join("%s %s" % (fmt_tag(t),
                                        fmt_pos(diff['removed_node_pos'][t]))
                             for t in diff['removed_nodes'])))
    if diff['added_nodes']:
        lines.append("      %s added nodes   : %s" % (
            label, ", ".join("%s %s" % (fmt_tag(t),
                                        fmt_pos(diff['added_node_pos'][t]))
                             for t in diff['added_nodes'])))
    if diff['removed_segs']:
        lines.append("      %s removed segs  : %s" % (
            label, ", ".join("%s-%s" % (fmt_tag(a), fmt_tag(b))
                             for a, b in diff['removed_segs'])))
    if diff['added_segs']:
        lines.append("      %s added segs    : %s" % (
            label, ", ".join("%s-%s" % (fmt_tag(a), fmt_tag(b))
                             for a, b in diff['added_segs'])))
    if diff['moved_nodes']:
        lines.append("      %s moved nodes   : %s" % (
            label, ", ".join("%s %s -> %s" % (fmt_tag(t), fmt_pos(a),
                                              fmt_pos(b))
                             for t, (a, b) in
                             sorted(diff['moved_nodes'].items()))))
    return lines


def diagnose(records):
    """diagnose: attribute each differing step to something specific

    Prints, for every step where the two codes did not produce the same
    network, what each of them decided. ExaDiS' side is reconstructed from
    the before/after diff; pydis' side is its own decision record, which
    additionally says which interval fired and why anything was rejected.
    """
    differing = [r for r in records
                 if (r.exadis_collided or r.pydis_collided)
                 and not r.same_geometry]
    tag_only = [r for r in records
                if (r.exadis_collided or r.pydis_collided)
                and r.same_geometry and not r.tags_match]
    print("")
    print("  %d step(s) differ only in tag allocation, same geometry"
          % len(tag_only))
    print("")
    print("=" * 70)
    print("DIAGNOSIS: %d step(s) where the two codes differ" % len(differing))
    print("=" * 70)

    for r in differing:
        n_ex = len(r.pydis_events)
        print("")
        print("  step %d  (exadis dn=%+d dseg=%+d | pydis dn=%+d dseg=%+d)"
              % (r.istep, r.exadis_dn, r.exadis_dseg, r.pydis_dn, r.pydis_dseg))
        for line in describe(r.exadis_diff, 'exadis'):
            print(line)
        for line in describe(r.pydis_diff, 'pydis '):
            print(line)
        if r.pydis_events:
            print("      pydis  decisions  : %d" % n_ex)
            for e in r.pydis_events:
                seg1 = "%s-%s" % (fmt_tag(e['seg1'][0]), fmt_tag(e['seg1'][1]))
                seg2 = "%s-%s" % (fmt_tag(e['seg2'][0]), fmt_tag(e['seg2'][1]))
                extra = ""
                if e.get('position') is not None:
                    extra = " at %s" % fmt_pos(e['position'])
                detail = ""
                if 'L1' in e:
                    detail = (" L=(%.3f,%.3f) merge %s(%d arms)+%s(%d arms)"
                              " planes=%d"
                              % (e['L1'], e['L2'], fmt_tag(e['node1']),
                                 e['arms1'], fmt_tag(e['node2']), e['arms2'],
                                 e['planes']))
                print("        %-8s %-11s %-14s %s | %s%s%s"
                      % (e['kind'], e['interval'], e['outcome'],
                         seg1, seg2, extra, detail))
        else:
            print("      pydis  decisions  : none")

    classify(differing)
    print("")
    print("  decisions on steps that AGREE     : %s"
          % (interval_split(records, True) or "none"))
    print("  decisions on steps that DIFFER    : %s"
          % (interval_split(records, False) or "none"))


def classify(differing):
    """classify: sort the differences into the shapes the plan predicts"""
    both_acted = [r for r in differing
                  if r.exadis_collided and r.pydis_collided]
    only_exadis = [r for r in differing
                   if r.exadis_collided and not r.pydis_collided]
    only_pydis = [r for r in differing
                  if r.pydis_collided and not r.exadis_collided]
    multi = [r for r in differing if len(r.pydis_events) > 1]

    print("")
    print("  " + "-" * 66)
    print("  both acted but differently        : %d" % len(both_acted))
    print("  exadis acted, pydis did not       : %d" % len(only_exadis))
    print("  pydis acted, exadis did not       : %d" % len(only_pydis))
    print("  steps with >1 pydis decision      : %d  (parallel-detection"
          " exposure)" % len(multi))

    intervals = {}
    outcomes = {}
    for r in differing:
        for e in r.pydis_events:
            intervals[e['interval']] = intervals.get(e['interval'], 0) + 1
            outcomes[e['outcome']] = outcomes.get(e['outcome'], 0) + 1
    if intervals:
        print("  pydis decisions by interval       : %s"
              % ", ".join("%s %d" % kv for kv in sorted(intervals.items())))
        print("  pydis decisions by outcome        : %s"
              % ", ".join("%s %d" % kv for kv in sorted(outcomes.items())))
    return intervals


def interval_split(records, agreeing):
    """interval_split: which interval fires on steps that agree vs differ

    If the retroactive half were under-firing, pairs would fall through to
    the predictive test, which finds the crossing in the next step rather
    than the last one and so reports a different point along each segment.
    That changes the close-to-endpoint decision, and with it whether a
    segment is split or an existing node reused. The counts below are the
    cheap check on that story.
    """
    want = [r for r in records
            if (r.exadis_collided or r.pydis_collided)
            and bool(r.same_geometry) == agreeing]
    counts = {}
    for r in want:
        for e in r.pydis_events:
            counts[e['interval']] = counts.get(e['interval'], 0) + 1
    return counts


def kokkos_threads(banner):
    """kokkos_threads: the thread pool size Kokkos reports, or None

    From the thread_pool_topology[ N x T x V ] line of the startup banner,
    whose middle number is the thread count.
    """
    m = re.search(r'thread_pool_topology\[\s*\d+\s*x\s*(\d+)', banner)
    return int(m.group(1)) if m else None


def check_setup(setup, banner):
    """check_setup: the conditions the comparison is only meaningful under

    The thread count comes from the Kokkos banner, not from OMP_NUM_THREADS.
    This file sets that variable unconditionally, so reading it back would only
    confirm the file agrees with itself; the banner is what Kokkos took at
    initialize(). It also catches the assignment being moved after the pyexadis
    import, where it would have no effect while still reading back as '1'.
    """
    ok = report("setup: enforce_glide_planes resolves to 0",
                setup['enforce_glide_planes'] == 0)
    ok &= report("setup: use_glide_planes resolves to 0",
                 setup['use_glide_planes'] == 0)
    threads = kokkos_threads(banner)
    ok &= report("setup: exadis is running on 1 thread, got %s" % threads,
                 threads == 1)
    return bool(ok)


def check_agreement(records):
    """check_agreement: every step that collided must agree, both ways

    Asserted rather than merely reported. The two codes are expected to
    match now, so a step that collides and does not agree is a failure and
    not an observation. Steps where neither code touched the network are
    not interesting here and are left out.

    Both halves are required. Geometry alone would pass a step whose
    survivor carries the wrong tag, and tags alone would pass one whose
    node sits in the wrong place.
    """
    collided = [r for r in records if r.exadis_collided or r.pydis_collided]
    bad = [r.istep for r in collided
           if not (r.tags_match and r.same_geometry)]
    label = "every colliding step agrees in tags and geometry"
    if bad:
        label += "  [differs at %s]" % ", ".join(str(i) for i in bad)
    ok = report(label, not bad and bool(collided))

    # How far apart the positions are, not just whether canonical_form called
    # them equal. That comparison rounds to its own tol, so it cannot report a
    # margin; a run sitting just inside the rounding would read the same as one
    # agreeing to the last bit.
    measured = [(r.istep, r.max_pos_diff) for r in collided
                if r.max_pos_diff is not None]
    if not measured:
        return bool(report("node positions on colliding steps: nothing to "
                           "compare", False) and ok)

    # a pairing that is not one to one makes its gap a lower bound, so say so
    degenerate = [r.istep for r in collided if r.match_one_to_one is False]
    if degenerate:
        ok &= report("node pairing is one to one on every colliding step "
                     "[not at %s]"
                     % ", ".join(str(i) for i in degenerate), False)

    step, worst = max(measured, key=lambda t: t[1])
    where = "" if worst == 0.0 else ", at step %d" % step
    # "after one collision" because that is the whole of what is compared:
    # _replay_pydis imports exadis' snapshot and runs the pydis rule on it
    # once, so neither side accumulates and a step's number is that step's
    # disagreement rather than a running total
    ok &= report_close("node positions after one collision, worst of %d "
                       "colliding steps%s" % (len(measured), where),
                       np.array([worst]), np.zeros(1), ATOL_POS)
    return bool(ok)


def check_tag_state_consistency(records):
    """check_tag_state_consistency: live exadis tag-pool binding vs the stored recording

    See SimulateNetworkComparingCollision.check_tag_state_consistency for
    what is actually being compared. Steps the recording has no entry for
    (r.tag_state_match is None) are not counted either way -- the recording
    only covers the trajectory it was made against, and a shorter or longer
    run than that is not itself a disagreement.
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
    """check_invariants: what must hold whichever collision rule runs

    Plumbing only: that the inputs a retroactive rule needs were present
    and that nothing raised. Agreement between the two codes is asserted
    separately, in check_agreement.
    """
    n = len(records)
    ok = report("old positions present at every step",
                all(r.xold_present for r in records))
    ok &= report("old positions aligned with the network at every step",
                 all(r.xold_aligned for r in records))
    ok &= report("a timestep is available at every step",
                 all(r.dt for r in records))
    ok &= report("velocities present at every step",
                 all(r.vel_present for r in records))
    ok &= report("no node lacks a velocity",
                 all(not r.vel_missing for r in records))
    ok &= report("segment order survives export/import at every step",
                 all(r.segorder_preserved for r in records))
    ok &= report("the pydis collision ran without raising at every step",
                 not any(r.pydis_error for r in records))
    return bool(ok) and n > 0


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--max-step', type=int, default=MAX_STEP)
    parser.add_argument('--config', default='frank_read_elast',
                        choices=sorted(CONFIGS))
    parser.add_argument('--plot', action='store_true', default=False,
                        help='show the network while it runs, and let the '
                             'native output through')
    parser.add_argument('--print-freq', type=int, default=50)
    parser.add_argument('--diagnose', action='store_true', default=False,
                        help='print what each code decided on every '
                             'differing step')
    parser.add_argument('--pydis-mode', default=PYDIS_MODE,
                        choices=['Proximity', 'Retroactive'],
                        help='which pydis collision rule to compare')
    args = parser.parse_args(argv)

    import pyexadis
    with quiet_native_output(enabled=not args.plot) as captured:
        pyexadis.initialize()
        records, setup = run(config=args.config, max_step=args.max_step,
                             plot=args.plot, print_freq=args.print_freq,
                             pydis_mode=args.pydis_mode)
    for line in kokkos_summary(captured.text):
        print("  %s" % line)

    summarize(records)
    if args.diagnose:
        diagnose(records)
    print("")
    ok = check_setup(setup, captured.text)
    ok &= check_invariants(records)
    ok &= check_agreement(records)
    ok &= check_tag_state_consistency(records)
    print("")
    verdict("test_collision_mode_pydis_exadis", ok)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
