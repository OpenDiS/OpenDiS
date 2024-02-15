import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from pydis.disnet import DisNet
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.collision.collision_disnet import Collision
from pydis.remesh.remesh_disnet import Remesh
from pydis.visualize.vis_disnet import VisualizeNetwork
from pydis.simulate.sim_disnet import SimulateNetwork

def init_two_disl_lines(z0=1.0, b1=np.array([-1.0,1.0,1.0]), b2=np.array([1.0,-1.0,1.0])):
    # To do:
    print("init_two_disl_lines: z0 = %f" % (z0))
    G = DisNet()
    rn    = np.array([[0.0, -z0, -z0,  7],
                      [0.0,  0.0, 0.0, 0],
                      [0.0,  z0,  z0,  7],
                      [-z0,  0.0,-z0,  7],
                      [0.0,  0.0, 0.0, 0],
                      [ z0,  0.0, z0,  7]])
    # To do: move calculation of glide plane normals to a separate function
    xi1, xi2 = rn[2,:3] - rn[1,:3], rn[5,:3] - rn[4,:3]
    n1,  n2 = np.cross(b1, xi1),  np.cross(b2, xi2)
    n1,  n2 = n1 / np.linalg.norm(n1),  n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0,:] = np.array([0, 1, b1[0], b1[1], b1[2], n1[0], n1[1], n1[2]])
    links[1,:] = np.array([1, 2, b1[0], b1[1], b1[2], n1[0], n1[1], n1[2]])
    links[2,:] = np.array([3, 4, b2[0], b2[1], b2[2], n2[0], n2[1], n2[2]])
    links[3,:] = np.array([4, 5, b2[0], b2[1], b2[2], n2[0], n2[1], n2[2]])
    G.add_nodes_links_from_list(rn, links)
    if not G.is_sane():
        raise ValueError("sanity check failed")
    return G

def main():
    global G, sim
    G = init_two_disl_lines(z0=4000)

    bounds = np.array([[-1.75e4, -1.75e4, -1.75e4], [1.75e4, 1.75e4, 1.75e4]])
    vis = VisualizeNetwork(bounds=bounds)

    calforce = CalForce(mu=160e9, nu=0.31, a=0.1, Ec=1.0e6,
                        applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -2.0e6, 0.0]),
                        force_mode='LineTension')

    mobility  = MobilityLaw(mobility_law='SimpleGlide')
    timeint   = TimeIntegration(integrator='EulerForward')
    collision = Collision(collision_mode='Proximity', mindist2=0.1)
    remesh    = Remesh(remesh_rule='LengthBased', Lmin=500, Lmax=2000)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, collision=collision, remesh=remesh, vis=vis,
                          dt0 = 1.0e-8, max_step=200,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1)
    sim.run(G)

    return G.is_sane()


if __name__ == "__main__":
    main()
