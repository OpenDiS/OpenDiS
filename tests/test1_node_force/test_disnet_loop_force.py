import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from pydis.disnet import DisNet
from pydis.calforce.calforce_disnet import CalForce

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = DisNet()
    rn = np.loadtxt(rn_file)[:, 1:]
    links = np.loadtxt(links_file)
    G.add_nodes_segments_from_list(rn, links)
    return G

def compare_force(ref_file, node_force, rtol = 0.0, atol = 1e-9):
    print("compare_force: ref_file = '%s'" % (ref_file))
    ref_force = np.loadtxt(ref_file)
    is_force_close = np.allclose(node_force, ref_force, rtol=rtol, atol=atol)
    print("max error: ", np.max(np.abs(node_force - ref_force)))
    return is_force_close

def main():
    G = init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    state = {"burgmag": 3e-10, "mu": 50, "nu": 0.3, "a": 0.01, "rann": 3.0}
    calforce  = CalForce(state=state, Ec=1.0e6)
    #calforce = CalForce(mu=50, nu=0.3, a=0.01, Ec=1.0e6)
    #state = calforce.NodeForce_LineTension(G, state=state)
    nodeforce_dict, segforce_dict = calforce.NodeForce_LineTension(G, applied_stress=np.zeros(6))
    lt_force_array = np.array([nodeforce_dict[tag] for tag in G.all_nodes()])

    atol = 1e-4
    is_lt_force_close = compare_force("force_linetension.dat", lt_force_array, atol=atol)
    if is_lt_force_close:
        print("LineTension Test \033[32mPassed!\033[0m")
    else:
        print("LineTension Test \033[31mFailed!\033[0m")

    state["segforce_dict"] = segforce_dict
    state = calforce.NodeForce_from_SegForce(G, state=state)
    nodeforce_from_segforce_dict = state["nodeforce_dict"]
    lt_force_array = np.array([nodeforce_from_segforce_dict[tag] for tag in G.all_nodes()])
    is_lt_force_from_seg_close = compare_force("force_linetension.dat", lt_force_array, atol=atol)
    if is_lt_force_from_seg_close:
        print("LineTension (from segforce) Test \033[32mPassed!\033[0m")
    else:
        print("LineTension (from segforce) Test \033[31mFailed!\033[0m")

    nodeforce_dict, segforce_dict = calforce.NodeForce_Elasticity_SBA(G, applied_stress=np.zeros(6))
    elast_force_array = np.array([nodeforce_dict[tag] for tag in G.all_nodes()])

    # for debugging purposes:
    # np.savetxt("force.dat", elast_force_array)

    is_elast_force_close = compare_force("force_stress_analytic_python.dat", elast_force_array, atol=atol)
    if is_elast_force_close:
        print("Elasticity Test \033[32mPassed!\033[0m")
    else:
        print("Elasticity Test \033[31mFailed!\033[0m")

    state["segforce_dict"] = segforce_dict
    state = calforce.NodeForce_from_SegForce(G, state=state)
    nodeforce_from_segforce_dict = state["nodeforce_dict"]
    elast_force_array = np.array([nodeforce_from_segforce_dict[tag] for tag in G.all_nodes()])
    is_elast_force_from_seg_close = compare_force("force_stress_analytic_python.dat", elast_force_array, atol=atol)
    if is_elast_force_from_seg_close:
        print("Elasticity (from segforce) Test \033[32mPassed!\033[0m")
    else:
        print("Elasticity (from segforce) Test \033[31mFailed!\033[0m")

    is_force_close = is_lt_force_close and is_lt_force_from_seg_close and is_elast_force_close and is_elast_force_from_seg_close
    return is_force_close


if __name__ == "__main__":
    is_force_close = main()
    print("is_force_close = %s" % is_force_close)

    if is_force_close:
        print("test" + '\033[32m' + " PASSED" + '\033[0m')
    else:
        print("test" + '\033[31m' + " FAILED" + '\033[0m')

    exit(0 if is_force_close else 1)
