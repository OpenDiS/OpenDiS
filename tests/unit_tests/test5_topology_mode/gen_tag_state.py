"""Record ExaDiS' own tag pool state at every step, for seeding PyDiS' replay.

Analogous to test4_collision_mode's exadis_tag_state.npz
(see load_tag_state()/seed_tag_state() in test_collision_mode_pydis_exadis.py),
but derived differently. That file was made by instrumenting ExaDiS' C++
(printing network->maxindex and network->recycled_indices at the top of the
collision handler) because collision there can recycle a tag on any of many
steps, and a tag freed on one step can still be sitting in ExaDiS' pool many
steps later -- reconstructing that history from Python-level snapshots alone
would mean tracking the pool across the whole 600-step run.

This test's situation is simpler and needs no C++ instrumentation: ExaDiS
splits a node exactly once in the whole 200-step run (see the plan,
.plan/2026-08-20/plan_topology_mode.md section 8), immediately preceded by
the one collision merge that creates the 4-arm node, both in the same step.
So the tag ExaDiS's split reuses is exactly whichever tag the merge freed
*that same step*, computable directly from the tag sets ExaDiS's own network
holds immediately before and after collision.HandleCol() runs, using data
this test's harness already has -- no recording, no instrumentation.

maxindex[i] and recycled_indices[i] are ExaDiS' pool state as it stood
immediately before topology.Handle() runs at step i (the same point
test_topology_mode_pydis_exadis.py's driver calls it 'before'), matching
what a PyDiS-side seed_tag_state() would need to feed DisNet.get_new_tag()
so it reuses the same tag ExaDiS does. Written top-of-stack first, like
test4_collision_mode's file, even though this run only ever has 0 or 1
recycled tags at a time and the order does not matter here.

Run once, like test4_collision_mode's own recording; regenerate if the base
configuration in test_topology_mode_pydis_exadis.py changes.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')

import sys
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if p not in sys.path]

import numpy as np

import test_topology_mode_pydis_exadis as t5

OUT_DIR = Path(__file__).resolve().parent / 'output'
OUT_FILE = OUT_DIR / 'exadis_tag_state.npz'


def record(max_step=t5.MAX_STEP):
    import pyexadis
    from pyexadis_base import ExaDisNet, SimulateNetwork, CalForce, MobilityLaw
    from pyexadis_base import TimeIntegration, Collision, Remesh, Topology
    from framework.simulation_setup import check_cutoff_maxseg

    state = t5.base_state()
    check_cutoff_maxseg(t5.LBOX * np.eye(3), t5.CUTOFF, state["maxseg"])
    net = t5.init_two_disl_lines(state, pbc=False)

    calforce = CalForce(force_mode='CUTOFF_MODEL', state=state, Ec=t5.EC, cutoff=t5.CUTOFF)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint = TimeIntegration(integrator='EulerForward', dt=t5.DT, state=state)
    collision = Collision(collision_mode='Proximity', state=state)
    topology = Topology(topology_mode='TopologySerial', state=state,
                        force=calforce, mobility=mobility)
    remesh = Remesh(remesh_rule='LengthBased', state=state)

    maxindex_list, recycled_list = [], []
    # seeded with the initial (pre-step-1) network's own tags: without this,
    # a tag that the very first step's collision removes before this
    # function's own loop ever takes a snapshot would never be recorded as
    # having existed at all, and step 1's recycled tag would be missed
    # entirely -- exactly what happened on the first attempt at this script.
    ever_seen = {tuple(int(v) for v in t)
                for t in net.get_disnet(ExaDisNet).export_data()['nodes']['tags']}

    class RecordingDriver(SimulateNetwork):
        def step_topological_operations(self, N, state):
            if self.cross_slip is not None:
                self.cross_slip.Handle(N, state)
            if self.collision is not None:
                self.collision.HandleCol(N, state)

            tags_now = {tuple(int(v) for v in t)
                       for t in N.get_disnet(ExaDisNet).export_data()['nodes']['tags']}
            ever_seen.update(tags_now)
            maxindex = max(t[1] for t in ever_seen)
            recycled = sorted(t[1] for t in ever_seen if t not in tags_now)
            maxindex_list.append(maxindex)
            recycled_list.append(list(reversed(recycled)))  # top of stack first

            if self.topology is not None:
                self.topology.Handle(N, state)
            if self.remesh is not None:
                self.remesh.Remesh(N, state)

    sim = RecordingDriver(calforce=calforce, mobility=mobility, timeint=timeint,
                          collision=collision, topology=topology, remesh=remesh,
                          state=state, max_step=max_step, loading_mode='stress',
                          applied_stress=np.zeros(6),
                          print_freq=max_step, write_freq=None,
                          write_dir=str(OUT_DIR))
    sim.run(net, state)
    return maxindex_list, recycled_list


def main():
    import pyexadis
    from framework.testing import quiet_native_output, kokkos_summary

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with quiet_native_output() as captured:
        pyexadis.initialize()
        maxindex_list, recycled_list = record()
    for line in kokkos_summary(captured.text):
        print("  %s" % line)

    maxindex = np.array(maxindex_list, dtype=np.int64)
    offsets = np.zeros(len(recycled_list) + 1, dtype=np.int64)
    for i, r in enumerate(recycled_list):
        offsets[i + 1] = offsets[i] + len(r)
    flat = np.array([v for r in recycled_list for v in r], dtype=np.int64)

    np.savez(str(OUT_FILE), maxindex=maxindex, recycled_flat=flat,
            recycled_offsets=offsets)
    n_nonempty = sum(1 for r in recycled_list if r)
    print("wrote %s: %d steps, %d with a non-empty recycled pool"
          % (OUT_FILE, len(maxindex_list), n_nonempty))
    for i, r in enumerate(recycled_list):
        if r:
            print("  step %d: maxindex=%d recycled=%s" % (i, maxindex_list[i], r))


if __name__ == '__main__':
    main()
