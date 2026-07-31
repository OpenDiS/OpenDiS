import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_circular_loop(radius=100.0, N=20, box_length=1000.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=True):
    '''Generate a circular dislocation loop lying in the z = box_length/2 plane
    '''
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    theta = np.arange(N)*2.0*np.pi/N
    rn = np.vstack([radius*np.cos(theta), radius*np.sin(theta), np.zeros_like(theta)]).T
    rn += cell.center()
    rn = np.hstack([rn, np.full((N,1), DisNode.Constraints.UNCONSTRAINED)])

    # The loop is planar and burg_vec lies in its plane, so every segment glides on the
    # same (0,0,1) plane. Note this cannot be obtained as cross(burg_vec, t) the way it is
    # in 02_frank_read_src: that expression vanishes at the two screw orientations.
    links = np.zeros((N, 8))
    for i in range(N):
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, [0.0, 0.0, 1.0]))

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))

def main(max_step=200, dt=1.0e-9):
    global net, sim, state

    Lbox = 1000.0
    net = init_circular_loop(radius=0.1*Lbox, box_length=Lbox)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    # Geometry and discretization follow 02_frank_read_src (box 1000, core radius a = 1,
    # maxseg/minseg of that order) rather than test_disl_loop_pydis.py, which uses a loop of
    # radius 1 with a = 0.01 in a box of 10. That scale only works with the softened
    # Ec = 1.0e6 line tension used there: with real elasticity the self-force of a radius-1
    # loop is ~2e10 and no attainable applied stress can balance it (the critical stress
    # would be ~0.4*mu), so the loop collapses and annihilates within ~15 steps. At radius
    # 100 the critical stress is ~6.5e8 Pa, which is what the applied stress below is set
    # against. The material constants mu and nu are unchanged from test_disl_loop_pydis.py.
    state = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0, "maxseg": 60.0, "minseg": 20.0, "rann": 3.0}

    # Elasticity_SBA includes the full segment-segment elastic interaction (plus the i==j
    # self term, regularized by the core radius state["a"]), in contrast to the LineTension
    # mode used in test_disl_loop_pydis.py.
    # Note: this is an O(Nseg^2) double loop in python, so it is much slower than LineTension.
    calforce  = CalForce(force_mode='Elasticity_SBA', state=state)
    # SimpleGlide is used here to match 02_frank_read_src, not the 'Relax' law of
    # test_disl_loop_pydis.py. Relax sets v = f with no length normalization and no glide
    # projection, and has no counterpart in exadis, so it cannot be compared across codes.
    # vmax is raised well above its 1e9 default because the nodal velocity grows as the loop
    # expands. Leaving the default would silently clamp the nodes partway through the run and
    # make this disagree with the exadis version, whose GLIDE mobility applies no velocity cap.
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state, vmax=1.0e15)
    timeint   = TimeIntegration(integrator='EulerForward', dt=dt, state=state)
    # Topology is kept enabled only because the Proximity collision handler reads the
    # nodeflag_dict that Topology.init_topology_exemptions creates; with topology=None the
    # run fails with KeyError: 'nodeflag_dict'. Beware that split_mode='MaxDiss' calls
    # OneNodeForce for any node with 4 or more arms, and OneNodeForce is not implemented for
    # the Elasticity_* force modes -- see the same note in
    # 02_frank_read_src/test_frank_read_src_pydis_elast.py. The applied stress below is
    # deliberately super-critical so that the loop expands: a sub-critical stress makes it
    # collapse onto itself, and the resulting collision hits that unimplemented path.
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Proximity', state=state, nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=max_step, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -1.0e9, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output', save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-step', dest='max_step', type=int, default=200)
    parser.add_argument('--dt', dest='dt', type=float, default=1.0e-9)
    args = parser.parse_args()

    main(max_step=args.max_step, dt=args.dt)

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/disl_loop_pydis_elast_final.json')
