import sys, os
from pathlib import Path
# this script lives in tests/unit_tests/test1_node_force/, so the repository root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in ['python', 'lib', 'core/pydis/python']]
[sys.path.append(path) for path in pydis_paths if not path in sys.path]
# input and reference data belong to this test, so they are located relative to the script;
# any output the test may write stays relative to the user's current working directory
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import argparse

import numpy as np
from pydis.calforce.compute_stress_analytic_paradis       import compute_seg_stress_coord_dep, compute_seg_stress_coord_indep
from framework.testing import report_close
from pydis.build_info import bitrepro_math, build_description

mu = 1000.0
nu = 0.3
a = 0.01

# Exact on a bitwise-reproducible build, where this result has no reason to
# differ at all from a reference blessed by another such build, and holding it
# to anything looser would let a reproducibility regression through.
#
# StressDueToSeg.c calls no log() or atan(), only arithmetic and sqrt, so it is
# in fact reproducible either way and the _repro build changes nothing about it.
# The tolerance follows the build regardless, so that such a run holds every
# pydis kernel to exact equality with nothing quietly exempt.
atol = 0.0 if bitrepro_math() else 1e-10

seg_data = np.load(input_dir / "seg_data.npy")
p1_list = seg_data[:, 0:3]
p2_list = seg_data[:, 3:6]
b12_list = seg_data[:, 6:9]
x_list = seg_data[:, 9:12]

REF_NAME = "ref_seg_stress.npy"
OUT_DIR = Path('output')          # relative to the working directory

# stacked rather than written into a np.zeros_like(reference), so the shape comes
# from the kernel and a reference is not needed in order to make one
seg_stress = np.array(
    [compute_seg_stress_coord_indep(p1, p2, b12, x, mu, nu, a)
     for p1, p2, b12, x in zip(p1_list, p2_list, b12_list, x_list)])

parser = argparse.ArgumentParser()
parser.add_argument('--write-ref', dest='write_ref', action='store_true',
                    default=False,
                    help='save the computed stresses to output/ as a new '
                         'reference instead of comparing against the stored one')
args = parser.parse_args()

if args.write_ref:
    OUT_DIR.mkdir(exist_ok=True)
    out_file = OUT_DIR / REF_NAME
    np.save(out_file, seg_stress)
    print("write_ref: wrote %s, %d segments" % (out_file, seg_stress.shape[0]))
    print("")
    print("  inspect it, then bless it:")
    print("      cp %s ref_data/" % out_file)
    sys.exit(0)

ref_stress = np.load(ref_dir / REF_NAME)

# report_close prints the error, the tolerance it was judged against and a
# coloured PASSED or FAILED on one line, the same as the other tests here
print(build_description())
passed = report_close("pydis : segment stress at %d field points "
                      "(compute_seg_stress_coord_indep)" % seg_stress.shape[0],
                      seg_stress, ref_stress, atol)

tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
print("test " + tag + '\033[0m')
sys.exit(0 if passed else 1)
