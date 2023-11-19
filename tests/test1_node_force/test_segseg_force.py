import networkx as nx
import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from pydis.calforce.compute_stress_force_analytic_paradis import compute_segseg_force_vec
from pydis.calforce.compute_stress_force_analytic_python import python_segseg_force_vec

tolA, tolB, tolC = 1e-9, 1e-9, 1e-5

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
force_nint = 3

# Test A: use ParaDiS library (SBA)
f1A, f2A, f3A, f4A = compute_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)

max_err_A = np.max(np.abs(np.concatenate((f1A, f2A, f3A, f4A), axis=1) - f1234_ref))
print(f"\nmax error = {max_err_A:.6e}")

testA_PASSED = max_err_A < tolA
print("\033[32mTestA Passed!\033[0m" if testA_PASSED else "\033[31mTestA Failed!\033[0m")

# Test B: use Python code (translated from Matlab)
f1B, f2B, f3B, f4B = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)

max_err_B = np.max(np.abs(np.concatenate((f1B, f2B, f3B, f4B), axis=1) - f1234_ref))
print(f"\nmax error = {max_err_B:.6e}")

testB_PASSED = max_err_B < tolB
print("\033[32mTestB Passed!\033[0m" if testB_PASSED else "\033[31mTestB Failed!\033[0m")

sys.exit(0 if (testA_PASSED and testB_PASSED) else 1)

# Test C: use ParaDiS library (SBN1)
#quad_points, weights = np.polynomial.legendre.leggauss(force_nint)
#f1C, f2C, f3C, f4C = compute_segseg_force_SBN1_vec(p1, p2, p3, p4, b12, b34, mu, nu, a, quad_points, weights)
#
#max_err_C = np.max(np.abs(np.concatenate((f1C, f2C, f3C, f4C), axis=1) - f1234_ref))
#print(f"\nmax error = {max_err_C:.6e}")
#
#testC_PASSED = max_err_C < tolC
#print("\033[92mTestC Passed!\033[0m" if testC_PASSED else "\033[91mTestA Failed!\033[0m")
#
#sys.exit(0 if (testA_PASSED and testB_PASSED and testC_PASSED) else 1)
