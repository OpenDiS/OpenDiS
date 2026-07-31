import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))

def main(max_step=200):
    global net, sim, state

    Lbox = 1000.0
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}

    # Elasticity_SBA includes the full segment-segment elastic interaction (plus the
    # i==j self term, regularized by the core radius state["a"]), in contrast to the
    # LineTension mode used in test_frank_read_src_pydis.py.
    # Note: this is an O(Nseg^2) double loop in python, so it is much slower than LineTension.
    #
    # The cutoff is stated explicitly, and matched by test_frank_read_src_exadis_elast.py, so
    # that the two runs truncate the segment-segment interaction identically and can be
    # compared. sqrt(3)/2*Lbox is the half diagonal of the cubic cell, which is the largest
    # separation attainable under the minimum image convention: at or above it no pair is
    # dropped, so this is also the cutoff -> infinity limit and gives the same answer as
    # cutoff=None. It is deliberately NOT 0.5*Lbox. That would be a legitimate truncation on
    # the pydis side, but exadis' CUTOFF_MODEL silently discards a further ~75 pairs there --
    # its neighbor list compares segment mid-points using a periodic shift quantized to the
    # bin grid rather than the true minimum image, which under-counts whenever
    # cutoff + maxseg > Lbox/3. See .plan/2026-07-27/plan_pydis_elast.md section 9.1.
    # Revisit once exadis is fixed: comparing at a genuinely truncating cutoff is the more
    # interesting test.
    cutoff    = 0.5*np.sqrt(3.0)*Lbox
    calforce  = CalForce(force_mode='Elasticity_SBA', state=state, cutoff=cutoff)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    # KNOWN LIMITATION: this example only runs to step 269. Topology(split_mode='MaxDiss')
    # calls OneNodeForce for any node with 4 or more arms, and OneNodeForce is not
    # implemented for the Elasticity_* force modes. The bowing source stays free of such
    # nodes until the expanding loop collides with itself; the first 4-arm node appears at
    # step 270, and the run then stops with
    #     NotImplementedError: OneNodeForce_Elasticity_SBA not implemented yet
    # Hence the default max_step=200 below. Topology cannot simply be disabled as a
    # workaround: the Proximity collision handler reads the nodeflag_dict that only
    # Topology.init_topology_exemptions creates, so topology=None fails earlier still with
    # KeyError: 'nodeflag_dict'.
    # To do: implement OneNodeForce_Elasticity_SBA (pydis/calforce/calforce_disnet.py) and
    # raise max_step here. This test case should not be considered closed until then.
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Proximity', state=state, nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=max_step, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output', save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    # max_step > 269 currently fails, see the note on Topology in main()
    parser.add_argument('--max-step', dest='max_step', type=int, default=200)
    args = parser.parse_args()

    main(max_step=args.max_step)

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_pydis_elast_final.json')
