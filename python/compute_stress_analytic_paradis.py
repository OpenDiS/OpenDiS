import numpy as np
from ctypes import c_double
real8 = c_double

try:
    opendis_lib = __import__('opendis')
    found_opendis = True
except ImportError:
    found_opendis = False
    print(' cannot find opendis.py                  ')
    print(' export PYTHONPATH=$HOME/OpenDiS.git/lib:')


def compute_seg_stress(p1, p2, b, x, mu, nu, a):
    """
    dislocation segment from p1 to p2 with Burgers vector b
    field point at x
    """
    sigma=np.ctypeslib.as_ctypes(np.zeros((3,3), dtype=real8))
    opendis_lib.SegmentStress(
        *(mu, nu),
        *(b[0], b[1], b[2]),
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(x[0], x[1], x[2]),
        a,
        sigma
    )

    return np.array(sigma)
