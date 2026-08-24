"""LengthBased remesh, pydis against exadis, step by step.

An exadis Frank-Read simulation with elasticity. At every step the network is
snapshotted after collision, the exadis remesh runs and stays authoritative,
and the pydis remesh is replayed on a copy of the same input, so a
disagreement at one step cannot carry into the next.

The run is two-phase: under load the loop only expands and only refinement
fires, so a second phase at zero stress lets it contract and exercises
coarsening.

Parity is claimed for coarsen_mode=0 and enforce_glide_planes=0 only.
"""

import os
import sys
from pathlib import Path

# Kokkos reads this at initialize(), so it must be set before pyexadis is
# imported. Multi-threaded exadis is not reproducible run to run.
os.environ.setdefault('OMP_NUM_THREADS', '1')

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
from framework.testing import report, same_network
from framework.testing import quiet_native_output, kokkos_summary
from pydis.disnet import DisNet
from pydis import Remesh as PyDiS_Remesh

from test_frank_read_src_exadis_elast import init_frank_read_src_loop

MAX_STEP = 250          # run twice: under load, then relaxing at zero stress
CUTOFF_FRAC = 0.25      # cutoff as a fraction of the box edge
DETAIL_STEPS = 3        # differing steps reported in full before summarising
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


def counts(data):
    """counts: nodes and segments in an exported network"""
    return (data["nodes"]["tags"].shape[0],
            data["segs"]["nodeids"].shape[0])


def check_scope(state):
    """check_scope: the configuration this comparison is valid for

    enforce_glide_planes resolves to 0 only while no crystal type is set; with
    one, exadis takes a coarsen branch pydis has no counterpart for.
    """
    return report("scope: no crystal type set, so enforce_glide_planes is 0",
                  state.get("crystal") is None
                  and not state.get("enforce_glide_planes", 0)
                  and not state.get("use_glide_planes", 0))


def build_sim(state, cutoff, exadis_rule, pydis_rule, plot):
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
                     'phase': self.phase}
            try:
                G = as_disnet(before)
                self.seed_tag_state(G, live_maxindex, live_recycled)
                self.pydis_remesh.Remesh(DisNetManager(G), state)
                entry['pydis'] = (G.num_nodes(), G.num_segments())
                entry['same'] = same_network(G, as_disnet(after),
                                             verbose=False)
            except Exception as err:
                entry['error'] = "%s: %s" % (type(err).__name__, err)
            return entry

    # Ec=0.0 matches pydis Elasticity_SBA, which has no core-energy term.
    # CUTOFF_MODEL is the mode pydis can match: both truncate the pair sum at
    # the minimum image.
    # coarsen_mode=0 is the segment-centric branch, the one pydis implements;
    # the python wrapper defaults it to 1, so it has to be passed.
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
                             coarsen_mode=0),
        pydis_remesh=PyDiS_Remesh(state=state, remesh_rule=pydis_rule),
        vis=VisualizeNetwork() if plot else None,
        max_step=MAX_STEP, loading_mode='stress', applied_stress=STRESS,
        print_freq=PRINT_FREQ, write_freq=None,
        plot_freq=10 if plot else None, plot_pause_seconds=0.01,
        write_dir=str(Path(__file__).resolve().parent / 'output'))


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

    return report("pydis and exadis LengthBased agree at every step", not bad)


def main(config='frank_read', max_step=MAX_STEP, exadis_rule='LengthBased',
         pydis_rule='LengthBased', plot=False):
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
    passed = check_scope(state)
    if passed:
        sim = build_sim(state, cutoff, exadis_rule, pydis_rule, plot)
        sim.max_step = max_step
        print("comparing '%s' remesh, exadis '%s' against pydis '%s'"
              % (config, exadis_rule, pydis_rule))

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
    parser.add_argument('--plot', dest='plot', action='store_true',
                        default=False,
                        help='show the network while it runs; off by '
                             'default so the test can run headless')
    args = parser.parse_args()

    passed = main(config=args.config, max_step=args.max_step,
                  pydis_rule=args.pydis_rule, plot=args.plot)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
