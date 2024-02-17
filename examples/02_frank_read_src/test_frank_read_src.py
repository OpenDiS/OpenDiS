import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from pydis.disnet import DisNet, DisNode, Cell
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.collision.collision_disnet import Collision
from pydis.remesh.remesh_disnet import Remesh
from pydis.visualize.vis_disnet import VisualizeNetwork
from pydis.simulate.sim_disnet import SimulateNetwork

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    # To do:
    # add flags (=7 for fixed nodes) (done)
    # add plane_normal to DisEdge
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    G = DisNet(cell=cell)
    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE]])
    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    G.add_nodes_links_from_list(rn, links)
    return G

def main():
    global G, sim
    G = init_frank_read_src_loop(pbc=False)

    bounds = np.array([-0.5*np.diag(G.cell.h), 0.5*np.diag(G.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    calforce = CalForce(mu=160e9, nu=0.31, a=0.01, Ec=1.0e6,
                        applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -2.0e6, 0.0]),
                        force_mode='LineTension')

    mobility  = MobilityLaw(mobility_law='SimpleGlide')
    timeint   = TimeIntegration(integrator='EulerForward')
    collision = Collision(collision_mode='Proximity')
    remesh    = Remesh(remesh_rule='LengthBased', Lmin=0.1, Lmax=0.3)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, collision=collision, remesh=remesh, vis=vis,
                          dt0 = 1.0e-8, max_step=200,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1)
    sim.run(G)

    return G.is_sane()


if __name__ == "__main__":
    main()
