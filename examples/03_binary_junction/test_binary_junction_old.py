import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib'),os.path.abspath('../../core/pydis/python')])

from framework.disnet_manager import DisNetManager
from pydis.disnet import DisNode, Cell, CellList
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.topology.topology_disnet import Topology
from pydis.collision.collision_disnet import Collision
from pydis.remesh.remesh_disnet import Remesh
from pydis.visualize.vis_disnet import VisualizeNetwork
from pydis.simulate.sim_disnet import SimulateNetwork

def init_two_disl_lines(z0=1.0, box_length=8.0, b1=np.array([-1.0,1.0,1.0]), b2=np.array([1.0,-1.0,1.0]), pbc=False):
    # To do:
    print("init_two_disl_lines: z0 = %f" % (z0))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    cell_list = CellList(cell=cell, n_div=[4,4,4])
    rn    = np.array([[0.0, -z0, -z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  z0,  z0,  DisNode.Constraints.PINNED_NODE],
                      [-z0,  0.0,-z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [ z0,  0.0, z0,  DisNode.Constraints.PINNED_NODE]])
    # To do: move calculation of glide plane normals to a separate function
    xi1, xi2 = rn[2,:3] - rn[1,:3], rn[5,:3] - rn[4,:3]
    n1,  n2 = np.cross(b1, xi1),  np.cross(b2, xi2)
    n1,  n2 = n1 / np.linalg.norm(n1),  n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0,:] = np.concatenate(([0, 1], b1, n1))
    links[1,:] = np.concatenate(([1, 2], b1, n1))
    links[2,:] = np.concatenate(([3, 4], b2, n2))
    links[3,:] = np.concatenate(([4, 5], b2, n2))
    net = DisNetManager()
    net.add_disnet(cell=cell, cell_list=cell_list)
    net.add_nodes_links_from_list(rn, links)
    return net

def main():
    global net, sim
    net = init_two_disl_lines(z0=4000, box_length=3.5e4, pbc=False)

    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    calforce = CalForce(mu=160e9, nu=0.31, a=0.1, Ec=1.0e6,
                        applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                        force_mode='LineTension')

    mobility  = MobilityLaw(mobility_law='SimpleGlide', mob=1.0e6)
    topology  = Topology(split_mode='MaxDiss')
    timeint   = TimeIntegration(integrator='EulerForward')
    collision = Collision(collision_mode='Proximity', mindist2=0.1)
    remesh    = Remesh(remesh_rule='LengthBased', Lmin=500, Lmax=2000)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, topology=topology, collision=collision, remesh=remesh, vis=vis,
                          dt0 = 1.0e-8, max_step=200,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1)
    sim.run(net)

    return net.is_sane()


if __name__ == "__main__":
    main()
