import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from disnet import *

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = DisNet()
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
    
    return G

def main():
    G = init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    # To do:
    # 1. Add a utility class to plot DisNet
    # vis = VisDisNet
    # vis.plot(G)

    return G.is_sane()


if __name__ == "__main__":
    sanity_check_passed = main()
    exit(0 if sanity_check_passed else 1)
