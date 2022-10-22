import networkx as nx
import numpy as np

from dd_proto import *

ddsim = DDSim()

def set_param():
    print("set_param")
    ddsim.param.mu = 160e9
    ddsim.param.nu = 0.31

def init_loop(radius=1.0, N=20, burg=None, plane_normal=None):
    print("init_loop: create a circular loop")
    ddsim.param.bounds = np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])
    G = ddsim.disnet
    bv = np.array([1.0, 0.0, 0.0])
    theta = np.arange(N)*2.0*np.pi/N
    pos   = radius * np.vstack([np.cos(theta), np.sin(theta), np.zeros_like(theta)]).T
    for i in range(N):
        node = (0,i)
        nbr_node = (0, (i+1)%N)
        G.add_node(node, R = pos[i])
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
