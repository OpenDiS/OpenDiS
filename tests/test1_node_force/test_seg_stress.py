import networkx as nx
import numpy as np
import sys

from dd_proto import *

ddsim = DDSim()

mu = 1000.0
nu = 0.3
a = 0.01

atol = 1e-10

seg_data = np.load("seg_data.npy")
p1_list = seg_data[:, 0:3]
p2_list = seg_data[:, 3:6]
b12_list = seg_data[:, 6:9]
x_list = seg_data[:, 9:12]

ref_stress = np.load("ref_seg_stress.npy")
seg_stress = np.zeros_like(ref_stress)

for ii, (p1, p2, b12, x) in enumerate(zip(p1_list, p2_list, b12_list, x_list)):
    seg_stress[ii] = compute_seg_stress(p1, p2, b12, x, mu, nu, a)

isclose_seg_stress = np.allclose(seg_stress, ref_stress, rtol=0.0, atol=atol)

print(f"max error: {np.max(np.abs(seg_stress - ref_stress)): .6e}")
if isclose_seg_stress:
    print("\033[92mTest Passed!\033[0m")
else:
    print("\033[91mTest Failed!\033[0m")

sys.exit(0 if isclose_seg_stress else 1)