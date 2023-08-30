import networkx as nx
import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from dd_proto import *

ddsim = DDSim()

def set_param():
    print("set_param")
    ddsim.param.mu = 160e9
    ddsim.param.nu = 0.31
    ddsim.param.Ec = 1.0e6
    ddsim.param.applied_stress = np.array([0.0, 0.0, 0.0, 0.0, -1.0e6, 0.0])
    # to do: decide - should bounds be part of param or disnet ?
    ddsim.param.bounds = np.array([[-50, -50, -50], [50, 50, 50]])

    ddsim.param.boundary_cond   = ['inf', 'inf', 'inf']
    ddsim.param.force_mode      = 'LineTension'
    ddsim.param.collision_mode  = 'Proximity'
    ddsim.param.mobility_law    = 'FCC0'
    ddsim.param.time_integrator = 'EulerForward'
    ddsim.param.dt0 = 1.0e-8
    ddsim.param.max_step = 200
    ddsim.param.print_freq = 10
    ddsim.param.plot_freq = 10
    ddsim.param.plot_pause_seconds = 0.1

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = ddsim.disnet
    rn = np.loadtxt(rn_file)[:, 1:]
    links = np.loadtxt(links_file)

    N = rn.shape[0]
    num_links = links.shape[0]

    for i in range(N):
        G.add_node((0,i), R = rn[i])

    for j in range(num_links):
        seg = links[j, :2].astype(int)
        bv  = links[j, 2:]
        node = (0,seg[0])
        nbr_node = (0, seg[1])
        G.add_edge(node, nbr_node, burg_vec = bv)
        G.add_edge(nbr_node, node, burg_vec = -bv)

def print_force():
    print("print_force")
    return None

def compare_force(ref_file):
    ref_force = np.loadtxt(ref_file)
    print("ref_force = ", ref_force)

def main():
    set_param()
    init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    print("calculate node forces")
    #ddsim.NodeForce()
    #print_force()

    #compare_force(ref_file = "force_stress_analytic_python.dat")

    ddsim.plot_disnet()
    #ddsim.Run()


if __name__ == "__main__":
    main()
