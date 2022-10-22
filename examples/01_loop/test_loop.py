import networkx as nx
import numpy as np

from dd_proto import *

ddsim = DDSim()

def set_param():
    print("set_param")
    ddsim.param.mu = 160e9
    ddsim.param.nu = 0.31
    # to do: decide - should bounds be part of param or disnet ?
    ddsim.param.bounds = np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])

    ddsim.param.boundary_cond   = ['inf', 'inf', 'inf']
    ddsim.param.force_mode      = 'LineTension'
    ddsim.param.collision_mode  = 'Proximity'
    ddsim.param.mobility_law    = 'FCC0'
    ddsim.param.time_integrator = 'EulerForward'

def init_loop(radius=1.0, N=20, burg=None, plane_normal=None):
    print("init_loop: create a circular loop")
    G = ddsim.disnet
    bv = np.array([1.0, 0.0, 0.0])
    theta = np.arange(N)*2.0*np.pi/N
    pos   = radius * np.vstack([np.cos(theta), np.sin(theta), np.zeros_like(theta)]).T
    for i in range(N):
        G.add_node((0,i), R = pos[i])
    for i in range(N):
        node = (0,i)
        nbr_node = (0, (i+1)%N)
        G.add_edge(node, nbr_node, burg_vec = bv)
        G.add_edge(nbr_node, node, burg_vec = -bv)

def print_force():
    print("print_force")
    return None

def main():
    set_param()
    init_loop()

    print("calculate node forces")
    ddsim.NodeForce()
    print_force()

    ddsim.plot_disnet()


if __name__ == "__main__":
    main()
