"""LengthBased remesh, pydis against exadis, step by step.

An exadis Frank-Read simulation with elasticity. At every step the network is
snapshotted after collision, the exadis remesh runs and stays authoritative,
and the pydis remesh is replayed on a copy of the same input, so a
disagreement at one step cannot carry into the next.

The run is two-phase: under load the loop only expands and only refinement
fires, so a second phase at zero stress lets it contract and exercises
coarsening.

Both coarsening algorithms are covered, selected by --coarsen-mode: 0
segment-centric, 1 node-centric (the default here, and exadis' and ParaDiS'
own default). Whichever is chosen is passed to both codes from one place, so
the two are never compared across different algorithms. `make` runs both.

Parity is claimed for enforce_glide_planes=0 only.
"""

import os
import re
import sys
from pathlib import Path

# Kokkos reads this at initialize(), so it must be set before pyexadis is
# imported. Multi-threaded exadis is not reproducible run to run.
#
# Assigned rather than setdefault: this test compares the two remesh rules at
# zero tolerance, and exadis' remesh arithmetic is not thread-invariant. Forced
# to 8 threads, one run in two diverged by 3.5527e-15, about one ULP, at step
# 228. An OMP_NUM_THREADS inherited from the environment would therefore turn
# the comparison flaky, so it is overridden here rather than merely defaulted.
# Run the multi-threaded case deliberately if you want it, by editing this
# line; do not expect the exact comparison to survive it.
os.environ['OMP_NUM_THREADS'] = '1'

# this script lives in tests/unit_tests/test3_remesh_rule/, so the repository
# root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python', 'examples/02_frank_read_src']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np

from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import report, report_close, same_network
from framework.testing import quiet_native_output, kokkos_summary
from pydis.disnet import DisNet
from pydis import Remesh as PyDiS_Remesh

from test_frank_read_src_exadis_elast import init_frank_read_src_loop

MAX_STEP = 250          # run twice: under load, then relaxing at zero stress
CUTOFF_FRAC = 0.25      # cutoff as a fraction of the box edge
DETAIL_STEPS = 3        # differing steps reported in full before summarising

# Node positions are compared exactly. same_network already requires that,
# through DisNode.is_equivalent comparing R with ==, so this is the tolerance
# every step has always been held to; naming it puts the number on screen.
#
# Exactness here rests on the single thread forced at the top of this file, not
# on remesh being inherently order-independent. Each step hands pydis exadis' own snapshot, so
# the exadis trajectory's reproducibility does not enter the comparison; but
# exadis' remesh arithmetic itself is not thread-invariant. Forced to 8 threads,
# one run in two diverged at step 228 by 3.5527e-15, about one ULP, and
# same_network failed with it. So the pin is load-bearing, and check_threads
# below reports it rather than leaving a later failure to be puzzled over.
ATOL_POS = 0.0
PRINT_FREQ = 50         # steps between exadis progress lines
STRESS = np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0])


def init_frank_read_src():
    """init_frank_read_src: the elasticity example's configuration

    Its init function conditions the network so no segment starts longer than
    maxseg: the back edges run between pinned nodes and neither code refines
    those.
    """
    Lbox = 1000.0
    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
             "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}
    cutoff = CUTOFF_FRAC*Lbox
    check_cutoff_maxseg(Lbox*np.eye(3), cutoff, state["maxseg"])
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox,
                                   pbc=True, maxseg=state["maxseg"])
    return net, state, cutoff


CONFIGS = {'frank_read': init_frank_read_src}


def as_disnet(data):
    """as_disnet: an exported network as a pydis DisNet"""
    G = DisNet()
    G.import_data(data)
    return G


def max_position_diff(G, G_ref):
    """max_position_diff: greatest node position difference, matched by tag

    Returns inf when the tag sets differ, since there is then no node-to-node
    correspondence to measure and the two are not the same network anyway.
    """
    tags = sorted(G.all_nodes_tags())
    if tags != sorted(G_ref.all_nodes_tags()):
        return np.inf
    if not tags:
        return 0.0
    a = np.array([G.nodes(t).R for t in tags])
    b = np.array([G_ref.nodes(t).R for t in tags])
    return float(np.max(np.abs(a - b)))


def counts(data):
    """counts: nodes and segments in an exported network"""
    return (data["nodes"]["tags"].shape[0],
            data["segs"]["nodeids"].shape[0])


def kokkos_threads(banner):
    """kokkos_threads: the thread pool size Kokkos reports, or None

    From the thread_pool_topology[ N x T x V ] line of the startup banner, whose
    middle number is the thread count.
    """
    m = re.search(r'thread_pool_topology\[\s*\d+\s*x\s*(\d+)', banner)
    return int(m.group(1)) if m else None


def check_threads(banner):
    """check_threads: exadis really is running on the one thread it needs

    Read from the Kokkos banner, not from OMP_NUM_THREADS. This file now sets
    that variable unconditionally, so checking it would only confirm this file
    agrees with itself; the banner is what Kokkos actually took at
    initialize(), which is the thing the exact comparison depends on. It would
    also catch the assignment being moved after the pyexadis import, where it
    would have no effect while still reading back as '1'.
    """
    threads = kokkos_threads(banner)
    return report("setup: exadis is running on 1 thread (the exact comparison "
                  "needs it), got %s" % threads, threads == 1)


def check_scope(state):
    """check_scope: the configuration this comparison is valid for

    enforce_glide_planes resolves to 0 only while no crystal type is set; with
    one, exadis takes a coarsen branch pydis has no counterpart for.
    """
    return report("scope: no crystal type set, so enforce_glide_planes is 0",
                  state.get("crystal") is None
                  and not state.get("enforce_glide_planes", 0)
                  and not state.get("use_glide_planes", 0))


def build_sim(state, cutoff, exadis_rule, pydis_rule, plot, coarsen_mode):
    """build_sim: the comparing driver, wired up

    Built here rather than at module scope so importing this file does not
    require pyexadis.
    """
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
    from pyexadis_base import Collision, Remesh as ExaDiS_Remesh
    from pyexadis_base import SimulateNetwork, ExaDisNet, VisualizeNetwork

    class RemeshCompare(SimulateNetwork):
        """Run the exadis remesh, then replay the pydis remesh on a copy."""

        def __init__(self, *args, **kwargs):
            self.pydis_remesh = kwargs.pop('pydis_remesh')
            super().__init__(*args, **kwargs)
            self.record = []
            self.phase = 'load'

        def live_tag_state(self, N):
            """live_tag_state: exadis' own current (maxindex, recycled) pool

            Straight from exadis' own SerialDisNet via the binding added in
            exadis commit 20ea2e8 ("Added binding to SerialDisNet internal
            tag indexing"); same call as test4_collision_mode's and
            test5_topology_mode's methods of the same name. Unlike those two,
            this test has no stored ref_data/exadis_tag_state.npz recording
            to check against -- LengthBased remesh's own tag allocation
            already agreed without seeding pydis' replay from exadis' pool
            at all, so this is added defensively (matching test4/test5's
            pattern for consistency) rather than because a disagreement was
            observed.
            """
            sn = N.get_disnet(ExaDisNet).net._get_serial_network()
            return int(sn._maxindex()), [int(i) for i in sn._recycled_indices()]

        def seed_tag_state(self, G, maxindex, recycled):
            """seed_tag_state: give the replayed network exadis' tag pool

            Same logic as test4_collision_mode's/test5_topology_mode's
            methods of the same name. maxindex/recycled (top of stack
            first) are exadis' own live values for this step, reversed here
            because DisNet.get_new_tag() pops from the end of the list.
            """
            if maxindex is None:
                return
            live = {t[1] for t in G.all_nodes_tags()}
            free = [i for i in recycled if i not in live]
            G._max_tag = (0, max(maxindex, max(live)))
            G._recycled_tags = [(0, i) for i in reversed(free)]

        def step_topological_operations(self, N, state):
            if self.collision is not None:
                self.collision.HandleCol(N, state)

            before = N.get_disnet(ExaDisNet).export_data()
            live_maxindex, live_recycled = self.live_tag_state(N)
            self.remesh.Remesh(N, state)
            after = N.get_disnet(ExaDisNet).export_data()
            self.record.append(self.compare(before, after, state,
                                            live_maxindex, live_recycled))

        def compare(self, before, after, state, live_maxindex, live_recycled):
            """compare: pydis remesh on the same input, against exadis"""
            entry = {'before': counts(before), 'exadis': counts(after),
                     'pydis': None, 'same': False, 'error': None,
                     'maxdiff': None, 'phase': self.phase}
            try:
                G = as_disnet(before)
                self.seed_tag_state(G, live_maxindex, live_recycled)
                self.pydis_remesh.Remesh(DisNetManager(G), state)
                entry['pydis'] = (G.num_nodes(), G.num_segments())
                G_exadis = as_disnet(after)
                entry['same'] = same_network(G, G_exadis, verbose=False)
                entry['maxdiff'] = max_position_diff(G, G_exadis)
            except Exception as err:
                entry['error'] = "%s: %s" % (type(err).__name__, err)
            return entry

    # Ec=0.0 matches pydis Elasticity_SBA, which has no core-energy term.
    # CUTOFF_MODEL is the mode pydis can match: both truncate the pair sum at
    # the minimum image.
    # coarsen_mode is passed to both sides from one place, so the two cannot be
    # compared across different coarsening algorithms. It has to be passed
    # explicitly either way: exadis' python wrapper defaults it to 1 and pydis
    # defaults it to 0. pydis reads it from state.
    state["coarsen_mode"] = coarsen_mode
    return RemeshCompare(
        state=state,
        calforce=CalForce(force_mode='CUTOFF_MODEL', state=state, Ec=0.0,
                          cutoff=cutoff),
        mobility=MobilityLaw(mobility_law='SimpleGlide', state=state),
        timeint=TimeIntegration(integrator='EulerForward', dt=1.0e-8,
                                state=state),
        collision=Collision(collision_mode='Retroactive', state=state),
        topology=None,
        remesh=ExaDiS_Remesh(remesh_rule=exadis_rule, state=state,
                             coarsen_mode=coarsen_mode),
        pydis_remesh=PyDiS_Remesh(state=state, remesh_rule=pydis_rule),
        vis=VisualizeNetwork() if plot else None,
        max_step=MAX_STEP, loading_mode='stress', applied_stress=STRESS,
        print_freq=PRINT_FREQ, write_freq=None,
        plot_freq=10 if plot else None, plot_pause_seconds=0.01,
        write_dir='output')


def summarise(record):
    """summarise: what the run found, and whether it agreed throughout"""
    bad = [(i, e) for i, e in enumerate(record) if e['error'] or not e['same']]
    acted = [e for e in record if e['before'] != e['exadis']]

    def split(entries):
        """split: how many of these refined, and how many coarsened"""
        return (sum(1 for e in entries if e['exadis'][0] > e['before'][0]),
                sum(1 for e in entries if e['exadis'][0] < e['before'][0]))

    print("")
    print("steps compared              : %d" % len(record))
    print("steps where exadis remeshed : %d (refine %d, coarsen %d)"
          % ((len(acted),) + split(acted)))
    for phase in sorted({e['phase'] for e in record}):
        inphase = [e for e in acted if e['phase'] == phase]
        print("    phase '%-5s'            : %d (refine %d, coarsen %d)"
              % ((phase, len(inphase)) + split(inphase)))
    print("steps where networks differ : %d" % len(bad))
    # a run that only refines proves little: the two codes are known to agree
    # on refinement and to have differed on coarsening
    if not split(acted)[1]:
        print("note: no coarsening occurred, so only the refine path was "
              "compared")

    for i, e in bad[:DETAIL_STEPS]:
        print("")
        print("  step %d: nodes/segs before %s, exadis -> %s"
              % (i, e['before'], e['exadis']))
        print("    pydis -> %s" % (e['error'] or e['pydis'],))
    if len(bad) > DETAIL_STEPS:
        print("")
        print("  ... and %d more" % (len(bad) - DETAIL_STEPS))

    ok = report("pydis and exadis LengthBased agree at every step", not bad)

    # How far apart they are, not just whether they differ. A step that raised
    # has no measurement to contribute and is already counted in `bad` above;
    # if every step raised there is nothing to compare, which is a failure
    # rather than a vacuous pass.
    measured = [(i, e['maxdiff']) for i, e in enumerate(record)
                if e['maxdiff'] is not None]
    if not measured:
        return bool(report("pydis vs exadis node positions: no step produced "
                           "a comparable network", False) and ok)
    step, worst = max(measured, key=lambda t: t[1])
    # naming a step is only informative when one of them stands out; with every
    # step at zero, calling out the first would suggest it differs from the rest
    where = "" if worst == 0.0 else ", first at step %d" % step
    # "after one remesh operation" because that is the whole of what is
    # compared: pydis is handed exadis' snapshot and remeshes it once, so
    # neither side accumulates, and a step's number is that step's
    # disagreement rather than a running total
    ok &= report_close("pydis vs exadis node positions after one remesh "
                       "operation, worst of %d steps%s"
                       % (len(measured), where),
                       np.array([worst]), np.zeros(1), ATOL_POS)
    return bool(ok)


def main(config='frank_read', max_step=MAX_STEP, exadis_rule='LengthBased',
         pydis_rule='LengthBased', plot=False, coarsen_mode=1):
    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; this test did not run")
        return False

    with quiet_native_output() as buf:
        pyexadis.initialize()
    for line in kokkos_summary(buf.text):
        print("exadis: %s" % line)

    net, state, cutoff = CONFIGS[config]()
    # named up front, not just in the banner below: the two coarsening modes are
    # different algorithms, so which one ran is the first thing to know when
    # reading a result, passing or failing
    print("setup: coarsen_mode = %d (%s), passed to both codes"
          % (coarsen_mode,
             'node-centric, the exadis and ParaDiS default'
             if coarsen_mode == 1 else 'segment-centric'))
    # both are setup checks, so && them rather than gating one on the other:
    # a wrong thread count and a wrong glide-plane scope should both be named
    passed = check_threads(buf.text)
    passed = check_scope(state) and passed
    if passed:
        sim = build_sim(state, cutoff, exadis_rule, pydis_rule, plot,
                        coarsen_mode)
        sim.max_step = max_step
        print("comparing '%s' remesh, exadis '%s' against pydis '%s', "
              "coarsen_mode=%d (%s)"
              % (config, exadis_rule, pydis_rule, coarsen_mode,
                 'node-centric' if coarsen_mode == 1 else 'segment-centric'))

        print("phase 'load': %d steps under applied stress" % max_step)
        sim.run(net, state)

        # under load the loop only expands, so segments only ever lengthen and
        # coarsening never fires; removing the stress lets it contract
        sim.phase = 'relax'
        sim.num_steps = max_step
        state["applied_stress"] = np.zeros(6)
        sim.applied_stress = np.zeros(6)
        print("phase 'relax': %d further steps at zero stress" % max_step)
        sim.run(net, state)

        passed = summarise(sim.record)

    with quiet_native_output():
        pyexadis.finalize()
    return bool(passed)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='LengthBased remesh, pydis against exadis, step by step')
    parser.add_argument('--max-step', dest='max_step', type=int,
                        default=MAX_STEP)
    parser.add_argument('--config', dest='config', default='frank_read',
                        choices=sorted(CONFIGS))
    parser.add_argument('--pydis-rule', dest='pydis_rule',
                        default='LengthBased')
    parser.add_argument('--coarsen-mode', dest='coarsen_mode', type=int,
                        default=1, choices=[0, 1],
                        help='coarsening algorithm, passed to both codes: 0 '
                             'segment-centric, 1 node-centric (the exadis and '
                             'ParaDiS default)')
    parser.add_argument('--plot', dest='plot', action='store_true',
                        default=False,
                        help='show the network while it runs; off by '
                             'default so the test can run headless')
    args = parser.parse_args()

    passed = main(config=args.config, max_step=args.max_step,
                  pydis_rule=args.pydis_rule, plot=args.plot,
                  coarsen_mode=args.coarsen_mode)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
