"""Record ExaDiS' own tag pool state at every step, for seeding PyDiS' replay.

Analogous to test4_collision_mode's exadis_tag_state.npz
(see load_tag_state()/seed_tag_state() in test_collision_mode_pydis_exadis.py).

Originally (before exadis commit 20ea2e8, "Added binding to SerialDisNet
internal tag indexing") this file had no direct way to read ExaDiS' actual
internal maxindex/recycled_indices, unlike test4_collision_mode's recording
(made by instrumenting ExaDiS' C++ directly), so it reconstructed an
estimate from Python-level snapshots: tracking every tag ever seen and
inferring which ones must have been freed by noticing they were missing from
the current tag set. That reconstruction is GONE as of this rewrite -- it
turned out to be wrong at many steps once checked against the live binding
(see check_tag_state_consistency() in test_topology_mode_pydis_exadis.py,
which is what caught this). This now just calls the same live binding
test_topology_mode_pydis_exadis.py's own live_tag_state() calls, so what
this file records and what the actual test compares against are, by
construction, reading the same thing -- the file is a snapshot of it at one
point in time, kept as the stored reference that consistency check runs
against on every later run, not a second, independently-derived source.

maxindex[i] and recycled_indices[i] are ExaDiS' pool state as it stood
immediately before topology.Handle() runs at step i (the same point
test_topology_mode_pydis_exadis.py's driver calls it 'before'), matching
what a PyDiS-side seed_tag_state() would need to feed DisNet.get_new_tag()
so it reuses the same tag ExaDiS does. Written top-of-stack first, like
test4_collision_mode's file and like the binding's own _recycled_indices()
return value.

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

OUT_DIR = Path('output')   # cwd-relative, as the tests are
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

    class RecordingDriver(SimulateNetwork):
        def step_update_response(self, N, state):
            # matches SimulateNetworkComparingTopology's own override in
            # test_topology_mode_pydis_exadis.py -- this recording has to see
            # the exact same trajectory the actual test runs, stress switch
            # included, or the tag pool state recorded here will not line up
            # with what the test replays against.
            state = super().step_update_response(N, state)
            if state['istep'] == t5.STRESS_STEP:
                state["applied_stress"] = t5.UNZIP_STRESS
            return state

        def step_topological_operations(self, N, state):
            if self.cross_slip is not None:
                self.cross_slip.Handle(N, state)
            if self.collision is not None:
                self.collision.HandleCol(N, state)

            # straight from ExaDiS' own SerialDisNet, top of stack first --
            # see this file's own docstring and live_tag_state() in
            # test_topology_mode_pydis_exadis.py, which reads the same thing
            # during the actual test.
            sn = N.get_disnet(ExaDisNet).net._get_serial_network()
            maxindex_list.append(int(sn._maxindex()))
            recycled_list.append([int(i) for i in sn._recycled_indices()])

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
