import numpy as np
from ctypes import c_double
real8 = c_double

try:
    pydis_lib = __import__('pydis_lib')
    found_pydis = True
except ImportError:
    found_pydis = False
    print(' cannot find pydis_lib.py                  ')
    print(' export PYTHONPATH=$HOME/OpenDiS.git/lib:')


def compute_seg_stress_coord_dep(p1, p2, b, x, mu, nu, a):
    """
    dislocation segment from p1 to p2 with Burgers vector b
    field point at x
    """
    sigma=np.ctypeslib.as_ctypes(np.zeros((3,3), dtype=real8))
    pydis_lib.SegmentStress(
        *(mu, nu),
        *(b[0], b[1], b[2]),
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(x[0], x[1], x[2]),
        a,
        sigma
    )

    return np.array(sigma)

def compute_seg_stress_coord_indep(p1, p2, b, x, mu, nu, a):
    """
    dislocation segment from p1 to p2 with Burgers vector b
    field point at x
    """
    sigma_vec=np.ctypeslib.as_ctypes(np.zeros((6), dtype=real8))
    pydis_lib.StressDueToSeg(
        *(x[0], x[1], x[2]),
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(b[0], b[1], b[2]),
        *(a, mu, nu),
        sigma_vec
    )
    sigma = np.array([ [sigma_vec[0], sigma_vec[3], sigma_vec[5]],
                       [sigma_vec[3], sigma_vec[1], sigma_vec[4]],
                       [sigma_vec[5], sigma_vec[4], sigma_vec[2]] ])
    sigma[0,0] = sigma_vec[0]
    sigma[1,1] = sigma_vec[1]
    sigma[2,2] = sigma_vec[2]
    sigma[0,1] = sigma_vec[3]
    sigma[1,0] = sigma_vec[3]
    sigma[1,2] = sigma_vec[4]
    sigma[2,1] = sigma_vec[4]
    sigma[0,2] = sigma_vec[5]
    sigma[2,0] = sigma_vec[5]

    return np.array(sigma)
