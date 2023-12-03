from ctypes import *
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

#from pydis.util.paradis_init import paradis_init
from pydis.util.paradis_util import paradis_lib

def init_circular_loop(radius=1.0, N=20, burg_vec=np.array([1.0,0.0,0.0])):
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    G = DisNet()
    theta = np.arange(N)*2.0*np.pi/N
    rn    = np.vstack([radius*np.cos(theta), radius*np.sin(theta), np.zeros_like(theta)]).T
    links = np.zeros((N, 5))
    for i in range(N):
        links[i,:] = np.array([i, (i+1)%N, burg_vec[0], burg_vec[1], burg_vec[2]])
    G.add_nodes_links_from_list(rn, links)
    return G

def main():
    G = init_circular_loop()

    bounds = np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])
    vis = VisualizeNetwork(bounds=bounds)

    calforce = CalForce(mu=160e9, nu=0.31, a=0.01, Ec=1.0e6,
                        applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -1.0e6, 0.0]),
                        force_mode='LineTension')
    mobility  = MobilityLaw(mobility_law='Relax')
    timeint   = TimeIntegration(integrator='EulerForward')
    collision = None #Collision(collision_mode='Proximity')
    remesh    = None #Remesh(remesh_rule='LengthBased')

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, collision=collision, remesh=remesh, vis=vis,
                          dt0 = 1.0e-8, max_step=200, 
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1)

    global home
    #home = paradis_init()
    paradis = paradis_lib()
    home = paradis.paradis_init()

    param = home.contents.param

    param.contents.timestepIntegrator = b'forceBsubcycle'
    param.contents.subInteg0Integ1 = b'RKF-RKF'

    node_data = np.array([3, 
         0,0,      -4000.0000,   500.0000,        6000.0000,   1,   7,
               0,1,       0.5773503000,     0.5773503000,    -0.5773503000,
                          0.,               1.,               1.,
         0,1,       17.0000,     500.0000,        6000.0000,  2,   0,
               0,0,     -0.5773503000,    -0.5773503000,      0.5773503000,
                        0.,              1.,              1.,
               0,2,      0.5773503000,     0.5773503000,     -0.5773503000,
                        0.,              1.,              1.,
         0,2,       4000.0000,    500.0000,        6000.0000,   1,   7,
               0,1,     -0.5773503000,    -0.5773503000,     0.5773503000,
                        0.,              1.,              1.])

    paradis.FreeAllNodes(home)
    paradis.AddNodesFromArray(home, node_data.ctypes.data_as(POINTER(c_double)))

    #sim.run(G)

    #return G.is_sane()


if __name__ == "__main__":
    main()
