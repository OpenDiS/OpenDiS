import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from pydis.disnet import DisNet

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = DisNet()
    rn = np.loadtxt(rn_file)[:, 1:]
    links = np.loadtxt(links_file)
    G.add_nodes_links_from_list(rn, links)
    return G

def main():
    G = init_loop_from_file(rn_file = "loop_rn.dat", links_file = "loop_links.dat")

    return G.is_sane()


if __name__ == "__main__":
    sanity_check_passed = main()
    print("sanity_check_passed = %s" % sanity_check_passed)

    if sanity_check_passed:
        print("test" + '\033[32m' + " PASSED" + '\033[0m')
    else:
        print("test" + '\033[31m' + " FAILED" + '\033[0m')

    exit(0 if sanity_check_passed else 1)
