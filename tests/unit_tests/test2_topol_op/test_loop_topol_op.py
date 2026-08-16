import sys, os
from pathlib import Path
# this script lives in tests/unit_tests/test2_topol_op/, so the repository root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in ['python', 'lib', 'core/pydis/python']]
[sys.path.append(path) for path in pydis_paths if not path in sys.path]
# input data belongs to this test, so it is located relative to the script; any output the
# test may write stays relative to the user's current working directory
input_dir = Path(__file__).resolve().parent / 'input_data'

import numpy as np
from pydis.disnet import DisNet

def init_loop_from_file(rn_file, links_file):
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'" % (rn_file, links_file))
    G = DisNet()
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)
    G.add_nodes_segments_from_list(rn, links)
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
