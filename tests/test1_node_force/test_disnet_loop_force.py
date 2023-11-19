import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from disnet import DisNet
from calforce_disnet import CalForce

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = DisNet()
    rn = np.loadtxt(rn_file)[:, 1:]
    links = np.loadtxt(links_file)
    G.add_nodes_links_from_list(rn, links)
    return G

def compare_force(ref_file, node_force, rtol = 0.0, atol = 1e-9):
    print("compare_force: ref_file = '%s'" % (ref_file))
    ref_force = np.loadtxt(ref_file)
    is_force_close = np.allclose(node_force, ref_force, rtol=rtol, atol=atol)
    print("max error: ", np.max(np.abs(node_force - ref_force)))
    return is_force_close

def main():
    G = init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    # To do:
    # 1. Add a utility class to plot DisNet
    # vis = VisDisNet
    # vis.plot(G)

    calforce = CalForce(mu=50, nu=0.3, a=0.01, Ec=1.0e6)
    nodeforce_dict = calforce.NodeForce_LineTension(G)
    lt_force_array = np.array([nodeforce_dict[node] for node in G.nodes()])

    atol = 1e-4
    is_lt_force_close = compare_force("force_linetension.dat", lt_force_array, atol=atol)
    if is_lt_force_close:
        print("LineTension Test \033[32mPassed!\033[0m")
    else:
        print("LineTension Test \033[31mFailed!\033[0m")

    nodeforce_dict = calforce.NodeForce_Elasticity_SBA(G)
    elast_force_array = np.array([nodeforce_dict[node] for node in G.nodes()])

    # for debugging purposes:
    # np.savetxt("force.dat", elast_force_array)

    is_elast_force_close = compare_force("force_stress_analytic_python.dat", elast_force_array, atol=atol)
    if is_elast_force_close:
        print("Elasticity Test \033[32mPassed!\033[0m")
    else:
        print("Elasticity Test \033[31mFailed!\033[0m")

    is_force_close = is_lt_force_close and is_elast_force_close
    return is_force_close


if __name__ == "__main__":
    is_force_close = main()
    print("is_force_close = %s" % is_force_close)

    if is_force_close:
        print("test" + '\033[32m' + " PASSED" + '\033[0m')
    else:
        print("test" + '\033[31m' + " FAILED" + '\033[0m')

    exit(0 if is_force_close else 1)
