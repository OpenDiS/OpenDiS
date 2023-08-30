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


def compute_segseg_force(p1, p2, p3, p4, b1, b2, mu, nu, a, seg12local=1, seg34local=1):
    """
    dislocation segment from p1 to p2 with Burgers vector b1
    dislocation segment from p3 to p4 with Burgers vector b2
    """

    f1x, f1y, f1z = real8(), real8(), real8()
    f2x, f2y, f2z = real8(), real8(), real8()
    f3x, f3y, f3z = real8(), real8(), real8()
    f4x, f4y, f4z = real8(), real8(), real8()
    opendis_lib.SegSegForce(
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(p3[0], p3[1], p3[2]),
        *(p4[0], p4[1], p4[2]),
        *(b1[0], b1[1], b1[2]),
        *(b2[0], b2[1], b2[2]),
        *(a, mu, nu),
        *(seg12local, seg34local),
        *(f1x, f1y, f1z),
        *(f2x, f2y, f2z),
        *(f3x, f3y, f3z),
        *(f4x, f4y, f4z),
    )

    f1 = np.array([f1x.value, f1y.value, f1z.value])
    f2 = np.array([f2x.value, f2y.value, f2z.value])
    f3 = np.array([f3x.value, f3y.value, f3z.value])
    f4 = np.array([f4x.value, f4y.value, f4z.value])

    return f1, f2, f3, f4


# a vectorized version of the above function compute_segseg_force
def compute_segseg_force_vec(
    p1_list, p2_list, p3_list, p4_list, b1_list, b2_list, mu, nu, a
):
    f1_list = np.empty_like(p1)
    f2_list = np.empty_like(p2)
    f3_list = np.empty_like(p3)
    f4_list = np.empty_like(p4)

    for ii, (p1, p2, p3, p4, b1, b2) in enumerate(
        zip(p1_list, p2_list, p3_list, p4_list, b1_list, b2_list)
    ):
        f1_list[ii], f2_list[ii], f3_list[ii], f4_list[ii] = compute_segseg_force(
            p1, p2, p3, p4, b1, b2, mu, nu, a
        )

    return f1_list, f2_list, f3_list, f4_list


def compute_segseg_force_vec(
    p1_list,
    p2_list,
    p3_list,
    p4_list,
    b1_list,
    b2_list,
    mu,
    nu,
    a,
    seg12local=1,
    seg34local=1,
):
    f1 = np.empty_like(p1_list)
    f2 = np.empty_like(p2_list)
    f3 = np.empty_like(p3_list)
    f4 = np.empty_like(p4_list)

    for ii, (p1, p2, p3, p4, b1, b2) in enumerate(
        zip(p1_list, p2_list, p3_list, p4_list, b1_list, b2_list)
    ):
        f1[ii], f2[ii], f3[ii], f4[ii] = compute_segseg_force(
            p1, p2, p3, p4, b1, b2, mu, nu, a, seg12local, seg34local
        )

    return f1, f2, f3, f4
