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

import numpy as np
from pydis.calforce.compute_stress_analytic_paradis       import compute_seg_stress_coord_dep, compute_seg_stress_coord_indep

mu = 1000.0
nu = 0.3
a = 0.01

atol = 1e-10

seg_data = np.load(input_dir / "seg_data.npy")
p1_list = seg_data[:, 0:3]
p2_list = seg_data[:, 3:6]
b12_list = seg_data[:, 6:9]
x_list = seg_data[:, 9:12]

ref_stress = np.load(ref_dir / "ref_seg_stress.npy")
seg_stress = np.zeros_like(ref_stress)

for ii, (p1, p2, b12, x) in enumerate(zip(p1_list, p2_list, b12_list, x_list)):
    #seg_stress[ii] = compute_seg_stress_coord_dep(p1, p2, b12, x, mu, nu, a)
    seg_stress[ii] = compute_seg_stress_coord_indep(p1, p2, b12, x, mu, nu, a)

isclose_seg_stress = np.allclose(seg_stress, ref_stress, rtol=0.0, atol=atol)

print(f"max error: {np.max(np.abs(seg_stress - ref_stress)): .6e}")
if isclose_seg_stress:
    print("\033[32mTest Passed!\033[0m")
else:
    print("\033[31mTest Failed!\033[0m")

sys.exit(0 if isclose_seg_stress else 1)
