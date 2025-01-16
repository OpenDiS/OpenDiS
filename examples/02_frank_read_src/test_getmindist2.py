import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from pydis.collision.getmindist2_paradis import GetMinDist2_paradis
from pydis.collision.getmindist2_python  import GetMinDist2_python

test_random_cases = True
#np.random.seed(seed=12)
atol = 1e-12

# test a specific case
p1, v1, p2, v2, p3, v3, p4, v4 = np.array([
    0.55534053, 0.02721592, 0.70769838, 0.45870382, 0.7520356 ,
    0.2865388 , 0.47848469, 0.21513633, 0.89691467, 0.40746523,
    0.21030289, 0.86278309, 0.90734267, 0.37769534, 0.15415964,
    0.0079789 , 0.19020845, 0.81661666, 0.12003181, 0.98307746,
    0.88822869, 0.90810826, 0.602317  , 0.79189523]).reshape(8,3)

dist2, ddist2dt, L1, L2 = GetMinDist2_paradis(p1, v1, p2, v2, p3, v3, p4, v4)
dist2_python, ddist2dt_python, L1_python, L2_python = GetMinDist2_python(p1, v1, p2, v2, p3, v3, p4, v4)
errors_in_this_case = np.array([dist2 - dist2_python, ddist2dt - ddist2dt_python, L1 - L1_python, L2 - L2_python])

if test_random_cases:
    Ntests = 10000
    input_cases = np.zeros((Ntests, 24))
    errors = np.zeros((Ntests, 4))
    for icase in range(Ntests):
        p1, p2 = np.random.rand(3), np.random.rand(3)
        p3, p4 = np.random.rand(3), np.random.rand(3)
        v1, v2 = np.random.rand(3), np.random.rand(3)
        v3, v4 = np.random.rand(3), np.random.rand(3)
    
        input_cases[icase,:] = np.concatenate((p1, v1, p2, v2, p3, v3, p4, v4))
    
        dist2, ddist2dt, L1, L2 = GetMinDist2_paradis(p1, v1, p2, v2, p3, v3, p4, v4)
        dist2_python, ddist2dt_python, L1_python, L2_python = GetMinDist2_python(p1, v1, p2, v2, p3, v3, p4, v4)
    
        if icase < 2 or icase == Ntests - 1:
            print("case %d:" % (icase))
            print("dist2 = %f, dist2_python = %f, diff = %e" % (dist2, dist2_python, dist2 - dist2_python))
            print("ddist2dt = %f, ddist2dt_python = %f, diff = %e" % (ddist2dt, ddist2dt_python, ddist2dt - ddist2dt_python))
            print("L1 = %f, L1_python = %f, diff = %e" % (L1, L1_python, L1 - L1_python))
            print("L2 = %f, L2_python = %f, diff = %e" % (L2, L2_python, L2 - L2_python))
            print("")
    
        errors[icase,:] = np.array([dist2 - dist2_python, ddist2dt - ddist2dt_python, L1 - L1_python, L2 - L2_python])

    print("max error = %e" % (np.max(np.abs(errors))))
    if np.max(np.abs(errors)) < atol:
        print("test" + '\033[32m' + " PASSED" + '\033[0m')
    else:
        print("test" + '\033[31m' + " FAILED" + '\033[0m')
