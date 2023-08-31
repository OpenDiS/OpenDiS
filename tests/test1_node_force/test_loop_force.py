import networkx as nx
import numpy as np
import time
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from dd_proto import *

ddsim = DDSim()

def set_param():
    print("set_param")
    ddsim.param.mu = 50
    ddsim.param.nu = 0.3
    ddsim.param.a  = 0.01
    ddsim.param.Ec = 1.0e6
    ddsim.param.applied_stress = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # to do: decide - should bounds be part of param or disnet ?
    ddsim.param.bounds = np.array([[-50, -50, -50], [50, 50, 50]])

    ddsim.param.boundary_cond   = ['inf', 'inf', 'inf']
    #ddsim.param.force_mode      = 'Elasticity_SBA'
    ddsim.param.force_mode      = 'Elasticity_SBN1_SBA'
    ddsim.param.force_nint      = 3
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

def compare_force(ref_file, rtol = 0.0, atol = 1e-9):
    print("compare_force: ref_file = '%s'" % (ref_file))
    ref_force = np.loadtxt(ref_file)
    node_force = ddsim.Get_NodeForce_Array()
    is_force_close = np.allclose(node_force, ref_force, rtol=rtol, atol=atol)
    print("max error: ", np.max(np.abs(node_force - ref_force)))
    return is_force_close

def main():
    set_param()
    init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    print("calculate node forces")
    tic = time.time()
    ddsim.NodeForce()
    toc = time.time()
    print("\n***********************")
    print("Elapsed time: %g s" %(toc - tic))
    print("***********************\n")

    atol = 1e-4
    is_force_close = compare_force(ref_file = "force_stress_analytic_python.dat", atol=atol)
    if is_force_close:
        print("\033[92mTest Passed!\033[0m")
    else:
        print("\033[91mTest Failed!\033[0m")

    #ddsim.plot_disnet()
    return is_force_close


if __name__ == "__main__":
    is_force_close = main()
    exit(0 if is_force_close else 1)
