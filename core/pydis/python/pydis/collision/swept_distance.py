"""Minimum distance between moving points and segments.

Transcribed from ParaDiS RetroactiveCollision2.c (collisionMethod 4), the
same source ExaDiS CollisionRetroactive ports in
core/exadis/src/collision_types/collision_retroactive.cpp.

Note: this is NOT the algorithm of Sills and Cai (2014). That paper
describes RetroactiveCollisions in RetroactiveCollision.c, ParaDiS
collisionMethod 3, which detects collisions by interval bisection and
handles the earliest one first.

Everything here is pure geometry: numpy arrays in, numbers out, no
dislocation network. Each endpoint is assumed to travel in a straight
line at constant speed from its position at t to its position at tau, so
"during the interval" always means that linear sweep.
"""

import numpy as np

# BLAS-free 3-vector arithmetic on a _repro build, the platform's numpy
# otherwise. Bound once here rather than tested at each call site, because this
# module has too many for that to stay readable. np.dot and np.linalg.norm on a
# 3-vector both dispatch to BLAS, whose summation order is not fixed across
# vendors; see pydis.calforce.bitrepro_math.
from ..calforce.bitrepro_math import (ENABLED as BITREPRO_MATH,
                                      dot_fixed, norm_fixed,
                                      det3, solve3_by_inverse)

if BITREPRO_MATH:
    _dot, _norm = dot_fixed, norm_fixed
else:
    _dot, _norm = np.dot, np.linalg.norm

# Hard-coded small number, as in the C. Distances are squared lengths in
# units of b, so this is not a relative tolerance.
EPS = 1.0e-12

# Newton iteration limits, from the C.
MAX_NEWTON = 20
ERR_TOL = 1.0e-6


def _grow_sphere(center, radius, p):
    """_grow_sphere: smallest sphere containing the old one and p

    Returns (center, radius) unchanged when p is already inside. The
    guard on ds mirrors the C: when p coincides with the centre, the
    division would produce a NaN, which happens in practice when the
    criterion is applied to a hinge, where two of the four points are the
    same node.
    """
    d = p - center
    ds = _norm(d)
    if ds * ds < radius * radius:
        return center, radius
    new_radius = 0.5 * (radius + ds)
    if ds > 1e-20:
        center = p - d / ds * new_radius
    return center, new_radius


def enclosing_sphere(points):
    """enclosing_sphere: a sphere containing every point given

    Not the minimal one: it starts from the sphere on the first two
    points and grows to admit each of the rest in turn, which is what
    ParaDiS FindSphere does for its four points. Only used as a cheap
    conservative filter, so a slightly loose sphere costs a little work
    and never a missed collision.
    """
    points = np.asarray(points, dtype=float)
    center = 0.5 * (points[0] + points[1])
    radius = 0.5 * _norm(points[0] - points[1])
    for p in points[2:]:
        center, radius = _grow_sphere(center, radius, p)
    return center, radius


def point_seg_min_dist(x0, y0, y1):
    """point_seg_min_dist: closest approach of point x0 to segment y0-y1

    Returns (dist2, L), with L in [0, 1] the position along the segment.
    ParaDiS MinPointSegDist.
    """
    diff = x0 - y0
    best = (float(_dot(diff, diff)), 0.0)

    d1 = x0 - y1
    dist2 = float(_dot(d1, d1))
    if dist2 < best[0]:
        best = (dist2, 1.0)

    seg = y0 - y1
    b = float(_dot(seg, seg))
    if b > EPS:
        t = float(-_dot(diff, seg) / b)
        if 0.0 < t < 1.0:
            v = diff + seg * t
            dist2 = float(_dot(v, v))
            if dist2 < best[0]:
                best = (dist2, t)
    return best


def seg_seg_min_dist(x0, x1, y0, y1):
    """seg_seg_min_dist: closest approach of segment x0-x1 to segment y0-y1

    Returns (dist2, L1, L2). ParaDiS MinSegSegDist: the four boundaries of
    the (L1, L2) unit square are checked first, then the interior
    stationary point, which only exists when the segments are not
    parallel.
    """
    diff = x0 - y0
    best = (float(_dot(diff, diff)), 0.0, 0.0)

    # the four edges of the solution square
    for x, fixed_first, fixed in ((x0, True, 0.0), (x1, True, 1.0)):
        dist2, l = point_seg_min_dist(x, y0, y1)
        if dist2 < best[0]:
            best = (dist2, fixed, l)
    for y, fixed in ((y0, 0.0), (y1, 1.0)):
        dist2, l = point_seg_min_dist(y, x0, x1)
        if dist2 < best[0]:
            best = (dist2, l, fixed)

    seg1, seg2 = x1 - x0, y1 - y0
    a = float(_dot(seg1, seg1))
    b = float(_dot(seg2, seg2))
    c = float(_dot(seg1, seg2))
    e = float(_dot(seg1, diff))
    d = float(_dot(seg2, diff))
    g = a * b - c * c

    if abs(g) > EPS:
        l1 = (d * c - e * b) / g
        l2 = (d + c * l1) / b
        if 0.0 < l1 < 1.0 and 0.0 < l2 < 1.0:
            v = diff + seg1 * l1 - seg2 * l2
            dist2 = float(_dot(v, v))
            if dist2 < best[0]:
                best = (dist2, l1, l2)
    return best


def point_point_min_dist_in_time(x1t, x1tau, x2t, x2tau):
    """point_point_min_dist_in_time: closest approach of two moving points

    Returns (dist2, t) with t in [0, 1] the fraction of the interval.
    ParaDiS MindistPtPtInTime.
    """
    diff = x1t - x2t
    best = (float(_dot(diff, diff)), 0.0)

    diff_tau = x1tau - x2tau
    dist2 = float(_dot(diff_tau, diff_tau))
    if dist2 < best[0]:
        best = (dist2, 1.0)

    ddiff = diff_tau - diff
    b = float(_dot(ddiff, ddiff))
    if b > EPS:
        t = float(-_dot(diff, ddiff) / b)
        if 0.0 < t < 1.0:
            v = diff + ddiff * t
            dist2 = float(_dot(v, v))
            if dist2 < best[0]:
                best = (dist2, t)
    return best


def point_seg_min_dist_in_time(x1t, x1tau, x3t, x3tau, x4t, x4tau):
    """point_seg_min_dist_in_time: moving point against moving segment

    Returns (dist2, L2, t). ParaDiS MinDistPtSegInTime.

    The endpoint cases are taken first and used to seed a Newton solve for
    an interior (L2, t). The C does not check the tratio=0 and tratio=1
    faces here, on the stated assumption that the caller has already
    covered them; seg_seg_min_dist_in_time does.
    """
    diff = x1t - x3t
    best = (float(_dot(diff, diff)), 0.0, 0.0)

    for endpoint_t, endpoint_tau, l2 in ((x3t, x3tau, 0.0), (x4t, x4tau, 1.0)):
        dist2, t = point_point_min_dist_in_time(x1t, x1tau, endpoint_t,
                                                endpoint_tau)
        if dist2 < best[0]:
            best = (dist2, l2, t)

    _, l2, t = best
    l34 = x3t - x4t
    dl34 = (x3tau - x3t) + (x4t - x4tau)
    dl13 = (x1tau - x1t) + (x3t - x3tau)

    a = float(_dot(diff, diff))
    b = float(_dot(diff, l34))
    c = float(_dot(diff, dl13))
    d = float(_dot(diff, dl34))
    e = float(_dot(l34, l34))
    f = float(_dot(l34, dl13))
    g = float(_dot(l34, dl34))
    h = float(_dot(dl13, dl13))
    i = float(_dot(dl13, dl34))
    j = float(_dot(dl34, dl34))
    d = d + f

    for _ in range(MAX_NEWTON):
        l2l2 = l2 * l2
        err = np.array([
            c + d * l2 + g * l2l2 + t * (h + 2.0 * i * l2 + j * l2l2),
            b + e * l2 + t * (d + 2.0 * g * l2) + t * t * (i + j * l2),
        ])
        m00 = h + 2.0 * i * l2 + j * l2l2
        m01 = d + 2.0 * g * l2 + 2.0 * t * (i + j * l2)
        m11 = e + 2.0 * g * t + t * t * j
        det = m00 * m11 - m01 * m01
        if abs(det) < EPS or float(_dot(err, err)) < ERR_TOL:
            break
        # Transcribed sign. This adds the correction where the segment
        # against segment solve below subtracts it; the same inconsistency
        # is in ParaDiS and in the ExaDiS port, so it is kept rather than
        # "fixed", which would break agreement with both.
        t = t + (m11 * err[0] - m01 * err[1]) / det
        l2 = l2 + (m00 * err[1] - m01 * err[0]) / det

    if 0.0 < l2 < 1.0 and 0.0 < t < 1.0:
        l2l2 = l2 * l2
        dist2 = ((a + 2.0 * b * l2 + e * l2l2)
                 + 2.0 * t * (c + d * l2 + g * l2l2)
                 + t * t * (h + 2.0 * i * l2 + j * l2l2))
        if dist2 < best[0]:
            best = (dist2, l2, t)
    return best


def seg_seg_min_dist_in_time(x1t, x1tau, x2t, x2tau, x3t, x3tau, x4t, x4tau):
    """seg_seg_min_dist_in_time: two moving segments, closest approach

    Segment 1 runs x1 to x2 and segment 2 runs x3 to x4, each endpoint
    sweeping linearly from its t position to its tau position. Returns
    (dist2, L1, L2, t). ParaDiS MinDistSegSegInTime.

    Six boundary cases are evaluated first: the two segments at each end
    of the interval, and each of the four endpoints against the opposite
    segment throughout it. The best of those seeds a Newton solve for an
    interior (t, L1, L2), which is accepted only if it lands strictly
    inside the unit cube and beats every boundary case.
    """
    dist2, l1, l2 = seg_seg_min_dist(x1t, x2t, x3t, x4t)
    best = (dist2, l1, l2, 0.0)

    dist2, l1, l2 = seg_seg_min_dist(x1tau, x2tau, x3tau, x4tau)
    if dist2 < best[0]:
        best = (dist2, l1, l2, 1.0)

    # each endpoint of one segment against the whole of the other
    for pt, pt_tau, l1_fixed in ((x1t, x1tau, 0.0), (x2t, x2tau, 1.0)):
        dist2, l2, t = point_seg_min_dist_in_time(pt, pt_tau, x3t, x3tau,
                                                  x4t, x4tau)
        if dist2 < best[0]:
            best = (dist2, l1_fixed, l2, t)
    for pt, pt_tau, l2_fixed in ((x3t, x3tau, 0.0), (x4t, x4tau, 1.0)):
        dist2, l1, t = point_seg_min_dist_in_time(pt, pt_tau, x1t, x1tau,
                                                  x2t, x2tau)
        if dist2 < best[0]:
            best = (dist2, l1, l2_fixed, t)

    _, l1, l2, t = best
    l13, l21, l34 = x1t - x3t, x2t - x1t, x3t - x4t
    dl13 = (x1tau - x1t) + (x3t - x3tau)
    dl21 = (x2tau - x2t) + (x1t - x1tau)
    dl34 = (x3tau - x3t) + (x4t - x4tau)

    dot = _dot          # the module alias, not np.dot; see the note at the top
    a, b, c = dot(l13, l13), dot(l13, l21), dot(l13, l34)
    d, e, f = dot(l13, dl13), dot(l13, dl21), dot(l13, dl34)
    g, h = dot(l21, l21), dot(l21, l34)
    i, j, k = dot(l21, dl13), dot(l21, dl21), dot(l21, dl34)
    ll = dot(l34, l34)
    m, n, o = dot(l34, dl13), dot(l34, dl21), dot(l34, dl34)
    pp, q, r = dot(dl13, dl13), dot(dl13, dl21), dot(dl13, dl34)
    ss, tt, u = dot(dl21, dl21), dot(dl21, dl34), dot(dl34, dl34)
    e, f, k = e + i, f + m, k + n

    for _ in range(MAX_NEWTON):
        t2, l1l1, l2l2, l1l2 = t * t, l1 * l1, l2 * l2, l1 * l2
        quad = pp + 2*q*l1 + 2*r*l2 + 2*tt*l1l2 + ss*l1l1 + u*l2l2
        err = np.array([
            d + e*l1 + f*l2 + k*l1l2 + j*l1l1 + o*l2l2 + t*quad,
            b + h*l2 + g*l1 + t*(e + k*l2 + 2*j*l1) + t2*(q + tt*l2 + ss*l1),
            c + h*l1 + ll*l2 + t*(f + k*l1 + 2*o*l2) + t2*(r + tt*l1 + u*l2),
        ])
        m01 = e + k*l2 + 2*j*l1 + 2.0*t*(q + tt*l2 + ss*l1)
        m02 = f + k*l1 + 2*o*l2 + 2.0*t*(r + tt*l1 + u*l2)
        mat = np.array([
            [quad, m01, m02],
            [m01, g + 2*j*t + ss*t2, h + k*t + tt*t2],
            [m02, h + k*t + tt*t2, ll + 2*o*t + u*t2],
        ])
        # det3 and solve3_by_inverse, not np.linalg: the reference forms the
        # determinant as an explicit cofactor sum and the correction through
        # an explicit inverse, where LAPACK factorizes. The two agree to
        # about 1e-16 per iteration, and four iterations of that compound
        # into 1e-13 in the ratio returned here, which is 1e-11 once it is
        # multiplied by a segment length. That was the whole of the residual
        # at step 271 of tests/unit_tests/test4_collision_mode.
        if abs(det3(mat)) < EPS or float(_dot(err, err)) < ERR_TOL:
            break
        corr = solve3_by_inverse(mat, err)
        t, l1, l2 = t - corr[0], l1 - corr[1], l2 - corr[2]

    if 0.0 < t < 1.0 and 0.0 < l1 < 1.0 and 0.0 < l2 < 1.0:
        t2, l1l1, l2l2, l1l2 = t * t, l1 * l1, l2 * l2, l1 * l2
        dist2 = (a + 2*b*l1 + 2*c*l2 + 2*h*l1l2 + g*l1l1 + ll*l2l2
                 + 2*t * (d + e*l1 + f*l2 + k*l1l2 + j*l1l1 + o*l2l2)
                 + t2 * (pp + 2*q*l1 + 2*r*l2 + 2*tt*l1l2 + ss*l1l1 + u*l2l2))
        if dist2 < best[0]:
            best = (float(dist2), float(l1), float(l2), float(t))
    return best


def segs_may_approach(seg1_points, seg2_points, mindist):
    """segs_may_approach: cheap conservative test before the real one

    True when the spheres enclosing each segment's four swept endpoints
    are within mindist of each other. A False is a guarantee that the
    segments stay further apart than mindist for the whole interval; a
    True means nothing on its own.

    Each argument is the four points (start and end of the interval, both
    endpoints) of one segment.
    """
    c1, r1 = enclosing_sphere(seg1_points)
    c2, r2 = enclosing_sphere(seg2_points)
    reach = r1 + r2 + mindist
    return float(_dot(c1 - c2, c1 - c2)) < reach * reach


def swept_seg_seg_collision(mindist, x1t, x1tau, x2t, x2tau,
                            x3t, x3tau, x4t, x4tau):
    """swept_seg_seg_collision: do these two segments meet during the interval

    Returns (collided, dist2, L1, L2). Segment 1 sweeps from x1t-x2t to
    x1tau-x2tau and segment 2 from x3t-x4t to x3tau-x4tau. ParaDiS
    CollisionCriterion.

    The interval is whatever the caller chooses. Retroactive collision
    passes the step just taken with mindist = rann; predictive collision
    passes the step about to be taken with a mindist near zero, so that it
    fires only on segments that will actually cross.

    A collision needs the closest approach to be both within mindist and
    at an interior point of both segments: an L outside [0, 1] means the
    extrapolated lines meet where the segments do not.
    """
    if not segs_may_approach((x1t, x2t, x1tau, x2tau),
                             (x3t, x4t, x3tau, x4tau), mindist):
        return False, np.inf, 0.0, 0.0

    dist2, l1, l2, _ = seg_seg_min_dist_in_time(x1t, x1tau, x2t, x2tau,
                                                x3t, x3tau, x4t, x4tau)
    collided = (dist2 < mindist * mindist
                and 0.0 <= l1 <= 1.0 and 0.0 <= l2 <= 1.0)
    return collided, dist2, l1, l2


def hinge_cos_angle(p1, p3, p4):
    """hinge_cos_angle: cosine of the angle between two arms of a node

    p1 is the shared node, p3 and p4 the far ends. Returns the cosine, and
    leaves the threshold to the caller: how closely two arms must align
    before they are zipped is policy, and this is geometry. ParaDiS
    HingeCollisionCriterion applies 0.98, relaxed to 0.9 when the two far
    nodes are themselves connected.

    Returns 0.0 for a zero-length arm, which no threshold accepts.
    """
    l13, l14 = p1 - p3, p1 - p4
    a = float(_dot(l13, l13))
    b = float(_dot(l14, l14))
    if a < EPS or b < EPS:
        return 0.0
    return float(_dot(l13, l14) / np.sqrt(a) / np.sqrt(b))
