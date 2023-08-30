import networkx as nx
import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from dd_proto import *

ddsim = DDSim()


tolA, tolB = 1e-9, 1e-9

segseg_data = np.loadtxt("segsep_min_2.5_max_32.5_iso_randombvecs_a0.010.dat")

p1 = segseg_data[:, 0:3]
p2 = segseg_data[:, 3:6]
p3 = segseg_data[:, 6:9]
p4 = segseg_data[:, 9:12]

b12 = segseg_data[:, 12:15]
b34 = segseg_data[:, 15:18]
f1234_ref = segseg_data[:, 18:]

mu = 50.0
nu = 0.3
a = 0.01

# Test A: use ParaDiS library
f1A, f2A, f3A, f4A = compute_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)

max_err_A = np.max(np.abs(np.concatenate((f1A, f2A, f3A, f4A), axis=1) - f1234_ref))
print(f"\nmax error = {max_err_A:.6e}")

testA_PASSED = max_err_A < tolA
print("\033[92mTestA Passed!\033[0m" if testA_PASSED else "\033[91mTestA Failed!\033[0m")

# Test B: use Python code (translated from Matlab)
f1B, f2B, f3B, f4B = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)

max_err_B = np.max(np.abs(np.concatenate((f1B, f2B, f3B, f4B), axis=1) - f1234_ref))
print(f"\nmax error = {max_err_B:.6e}")

testB_PASSED = max_err_B < tolB
print("\033[92mTestB Passed!\033[0m" if testB_PASSED else "\033[91mTestB Failed!\033[0m")

sys.exit(0 if (testA_PASSED and testB_PASSED) else 1)
