import numpy as np
from ctypes import c_double, POINTER
real8 = c_double

try:
    pydis_lib = __import__('pydis_lib')
    found_pydis = True
except ImportError:
    found_pydis = False
    raise

def compute_segseg_force_SBN1(p1, p2, p3, p4, b1, b2, mu, nu, a, quad_points, weights, seg12local=1, seg34local=1):
    """
    dislocation segment from p1 to p2 with Burgers vector b1
    dislocation segment from p3 to p4 with Burgers vector b2
    """

    f1x, f1y, f1z = real8(), real8(), real8()
    f2x, f2y, f2z = real8(), real8(), real8()
    f3x, f3y, f3z = real8(), real8(), real8()
    f4x, f4y, f4z = real8(), real8(), real8()
    Nint = quad_points.shape[0]
    quad_points_ptr=np.ctypeslib.as_ctypes(quad_points)
    weights_ptr=np.ctypeslib.as_ctypes(weights)
    pydis_lib.SegSegForce_SBN1(
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(p3[0], p3[1], p3[2]),
        *(p4[0], p4[1], p4[2]),
        *(b1[0], b1[1], b1[2]),
        *(b2[0], b2[1], b2[2]),
        *(a, mu, nu),
        *(Nint, quad_points_ptr, weights_ptr),
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

def compute_segseg_force_SBN1_SBA(p1, p2, p3, p4, b1, b2, mu, nu, a, quad_points, weights, seg12local=1, seg34local=1):
    """
    dislocation segment from p1 to p2 with Burgers vector b1
    dislocation segment from p3 to p4 with Burgers vector b2
    """

    f1x, f1y, f1z = real8(), real8(), real8()
    f2x, f2y, f2z = real8(), real8(), real8()
    f3x, f3y, f3z = real8(), real8(), real8()
    f4x, f4y, f4z = real8(), real8(), real8()
    Nint = quad_points.shape[0]
    quad_points_ptr=np.ctypeslib.as_ctypes(quad_points)
    weights_ptr=np.ctypeslib.as_ctypes(weights)
    pydis_lib.SegSegForce_SBN1_SBA(
        *(p1[0], p1[1], p1[2]),
        *(p2[0], p2[1], p2[2]),
        *(p3[0], p3[1], p3[2]),
        *(p4[0], p4[1], p4[2]),
        *(b1[0], b1[1], b1[2]),
        *(b2[0], b2[1], b2[2]),
        *(a, mu, nu),
        *(Nint, quad_points_ptr, weights_ptr),
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
def compute_segseg_force_SBN1_vec(
    p1_list, p2_list, p3_list, p4_list, b1_list, b2_list, mu, nu, a, quad_points, weights
):
    f1_list = np.empty_like(p1_list)
    f2_list = np.empty_like(p2_list)
    f3_list = np.empty_like(p3_list)
    f4_list = np.empty_like(p4_list)

    for ii, (p1, p2, p3, p4, b1, b2) in enumerate(
        zip(p1_list, p2_list, p3_list, p4_list, b1_list, b2_list)
    ):
        f1_list[ii], f2_list[ii], f3_list[ii], f4_list[ii] = compute_segseg_force_SBN1(
            p1, p2, p3, p4, b1, b2, mu, nu, a, quad_points, weights
        )

    return f1_list, f2_list, f3_list, f4_list

def compute_segseg_force(p1, p2, p3, p4, b1, b2, mu, nu, a, seg12local=1, seg34local=1):
    """
    dislocation segment from p1 to p2 with Burgers vector b1
    dislocation segment from p3 to p4 with Burgers vector b2
    """

    f1x, f1y, f1z = real8(), real8(), real8()
    f2x, f2y, f2z = real8(), real8(), real8()
    f3x, f3y, f3z = real8(), real8(), real8()
    f4x, f4y, f4z = real8(), real8(), real8()
    pydis_lib.SegSegForce(
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


def compute_segseg_force_batch(p1, p2, p3, p4, b1, b2, mu, nu, a,
                               seg12local=1, seg34local=1):
    """Batched segment-segment force: one ctypes call for the whole array.

    Same numbers as compute_segseg_force applied pair by pair, but the loop runs inside
    the compiled library (SegSegForceList) rather than in python, so the ctypes marshalling
    cost is paid once per batch instead of once per pair.

    p1..p4, b1, b2 are (N,3); returns f1..f4 each (N,3).
    """
    p1, p2, p3, p4, b1, b2 = (np.ascontiguousarray(x, dtype=np.float64)
                              for x in (p1, p2, p3, p4, b1, b2))
    n = p1.shape[0]
    f1 = np.empty((n, 3), dtype=np.float64)
    f2 = np.empty((n, 3), dtype=np.float64)
    f3 = np.empty((n, 3), dtype=np.float64)
    f4 = np.empty((n, 3), dtype=np.float64)
    dp = lambda x: x.ctypes.data_as(POINTER(c_double))
    pydis_lib.SegSegForceList(n, dp(p1), dp(p2), dp(p3), dp(p4), dp(b1), dp(b2),
                              a, mu, nu, seg12local, seg34local,
                              dp(f1), dp(f2), dp(f3), dp(f4))
    return f1, f2, f3, f4

# Array-in, array-out wrapper around the scalar compute_segseg_force.
#
# NOT vectorized: it loops in python and calls the compiled scalar routine once
# per pair, so it costs the same as calling compute_segseg_force directly in a
# loop. It was previously named compute_segseg_force_vec, which promised a
# speed-up it does not deliver. For a genuinely vectorized kernel use
# python_segseg_force_vec in compute_stress_force_analytic_python.py, which
# operates on whole batches with numpy.
def compute_segseg_force_list(
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
