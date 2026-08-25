"""Swept-distance geometry, checked against configurations worked by hand.

These pin the criteria one at a time, which a simulation cannot: a bug in
a branch that a given trajectory never reaches is invisible there. Every
case below has its expected answer written out next to it, derived by
hand rather than captured from a run, so this file is a specification and
not a regression baseline.

No dislocation network is involved, and nothing here imports pyexadis, so
this runs on a machine with no ExaDiS build.

    python3 test_swept_distance.py
"""

import sys
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in ['python', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if p not in sys.path]

import numpy as np

from framework.testing import report, verdict
from pydis.collision.swept_distance import (
    enclosing_sphere, point_seg_min_dist, seg_seg_min_dist,
    point_point_min_dist_in_time, seg_seg_min_dist_in_time,
    segs_may_approach, swept_seg_seg_collision, hinge_cos_angle,
)

# positions come from exact arithmetic on small integers and halves, so
# agreement should be to a few ULP rather than to a physical tolerance
TOL = 1.0e-9

# the thresholds ParaDiS HingeCollisionCriterion applies, which live in
# the caller rather than in the geometry
HINGE_TOL = 0.98
HINGE_TOL_TRI = 0.9


def v(*args):
    return np.array(args, dtype=float)


def close(a, b, tol=TOL):
    return abs(float(a) - float(b)) < tol


# --------------------------------------------------------------- static cases

def test_point_seg():
    """point to segment, the three branches: interior, and each end"""
    seg = (v(-1, 0, 0), v(1, 0, 0))

    # directly above the middle: foot of the perpendicular is interior
    dist2, l = point_seg_min_dist(v(0, 1, 0), *seg)
    ok = report("point-seg: above the middle, dist2=1 at L=0.5",
                close(dist2, 1.0) and close(l, 0.5))

    # beyond the far end: the perpendicular foot is outside, so the
    # nearest point is the endpoint itself
    dist2, l = point_seg_min_dist(v(2, 0, 0), *seg)
    ok &= report("point-seg: beyond the end, dist2=1 at L=1",
                 close(dist2, 1.0) and close(l, 1.0))

    # beyond the near end
    dist2, l = point_seg_min_dist(v(-2, 0, 0), *seg)
    ok &= report("point-seg: before the start, dist2=1 at L=0",
                 close(dist2, 1.0) and close(l, 0.0))

    # degenerate segment: both ends coincide, so only the endpoint
    # branches can fire and the interior one must not divide by zero
    dist2, l = point_seg_min_dist(v(0, 3, 0), v(0, 0, 0), v(0, 0, 0))
    ok &= report("point-seg: zero-length segment, dist2=9",
                 close(dist2, 9.0))
    return ok


def test_seg_seg():
    """segment to segment: skew, parallel, and endpoint-to-endpoint"""
    # perpendicular and offset in z; closest points are both mid-segment
    dist2, l1, l2 = seg_seg_min_dist(v(-1, 0, 0), v(1, 0, 0),
                                     v(0, -1, 1), v(0, 1, 1))
    ok = report("seg-seg: skew pair, dist2=1 at L1=L2=0.5",
                close(dist2, 1.0) and close(l1, 0.5) and close(l2, 0.5))

    # parallel: the interior system is singular, so the answer has to come
    # from the boundary cases
    dist2, l1, l2 = seg_seg_min_dist(v(-1, 0, 0), v(1, 0, 0),
                                     v(-1, 1, 0), v(1, 1, 0))
    ok &= report("seg-seg: parallel pair, dist2=1 (singular interior)",
                 close(dist2, 1.0))

    # collinear and disjoint: nearest approach is end to end
    dist2, l1, l2 = seg_seg_min_dist(v(-1, 0, 0), v(1, 0, 0),
                                     v(2, 0, 0.5), v(4, 0, 0.5))
    ok &= report("seg-seg: end to end, dist2=1.25 at L1=1, L2=0",
                 close(dist2, 1.25) and close(l1, 1.0) and close(l2, 0.0))

    # crossing at a point: zero distance, both interior
    dist2, l1, l2 = seg_seg_min_dist(v(-1, 0, 0), v(1, 0, 0),
                                     v(0, -1, 0), v(0, 1, 0))
    ok &= report("seg-seg: crossing pair, dist2=0 at L1=L2=0.5",
                 close(dist2, 0.0) and close(l1, 0.5) and close(l2, 0.5))
    return ok


def test_enclosing_sphere():
    """the filter sphere contains what it claims to"""
    pts = [v(0, 0, 0), v(2, 0, 0), v(0, 2, 0), v(0, 0, 2)]
    center, radius = enclosing_sphere(pts)
    inside = all(np.linalg.norm(p - center) <= radius + TOL for p in pts)
    ok = report("sphere: contains all four points", inside)

    # coincident points must not produce a NaN centre, which is the case
    # the ParaDiS guard exists for and which a hinge produces in practice
    center, radius = enclosing_sphere([v(1, 1, 1)] * 4)
    ok &= report("sphere: four coincident points give a finite centre",
                 np.all(np.isfinite(center)) and close(radius, 0.0))
    return ok


# ------------------------------------------------------------ moving in time

def test_point_point_in_time():
    """two points, one moving past the other"""
    # x1 sweeps from (0,0,0) to (2,0,0), x2 sits at (1,1,0): closest at
    # the half-way point of the interval, at unit distance
    dist2, t = point_point_min_dist_in_time(v(0, 0, 0), v(2, 0, 0),
                                            v(1, 1, 0), v(1, 1, 0))
    ok = report("pt-pt in time: closest at t=0.5, dist2=1",
                close(dist2, 1.0) and close(t, 0.5))

    # both stationary: the answer must be the static one, at t=0
    dist2, t = point_point_min_dist_in_time(v(0, 0, 0), v(0, 0, 0),
                                            v(3, 0, 0), v(3, 0, 0))
    ok &= report("pt-pt in time: stationary pair, dist2=9 at t=0",
                 close(dist2, 9.0) and close(t, 0.0))
    return ok


def test_static_limit():
    """with tau equal to t, the swept answer is the static answer

    Worth asserting on its own: it is the one case where two independent
    code paths in this module have to produce the same number, so it
    catches a transcription error in the Newton solve that a
    collision-or-not test would not.
    """
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    b0, b1 = v(0, -1, 1), v(0, 1, 1)
    static = seg_seg_min_dist(a0, a1, b0, b1)
    swept = seg_seg_min_dist_in_time(a0, a0, a1, a1, b0, b0, b1, b1)
    return report("swept: tau == t reproduces the static distance",
                  close(swept[0], static[0]) and close(swept[1], static[1])
                  and close(swept[2], static[2]))


def test_crossing_during_interval():
    """segments that pass through each other inside the interval"""
    # segment 1 fixed along x; segment 2 along y, descending through it
    # from z=+0.5 to z=-0.5, so the two are coincident half way through
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    dist2, l1, l2, t = seg_seg_min_dist_in_time(
        a0, a0, a1, a1,
        v(0, -1, 0.5), v(0, -1, -0.5), v(0, 1, 0.5), v(0, 1, -0.5))
    ok = report("swept: crossing pair reaches dist2=0 mid-interval",
                close(dist2, 0.0, 1e-6) and close(t, 0.5, 1e-6))

    collided, dist2, l1, l2 = swept_seg_seg_collision(
        0.1, a0, a0, a1, a1,
        v(0, -1, 0.5), v(0, -1, -0.5), v(0, 1, 0.5), v(0, 1, -0.5))
    ok &= report("criterion: crossing pair collides at mindist=0.1",
                 collided and 0.0 <= l1 <= 1.0 and 0.0 <= l2 <= 1.0)
    return ok


def test_approaching_but_clear():
    """segments that approach and stay further apart than mindist"""
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    # descends from z=2 to z=1, so the closest approach is 1
    args = (a0, a0, a1, a1,
            v(0, -1, 2), v(0, -1, 1), v(0, 1, 2), v(0, 1, 1))
    dist2, _, _, t = seg_seg_min_dist_in_time(*args)
    ok = report("swept: approaching pair bottoms out at dist2=1",
                close(dist2, 1.0, 1e-6) and close(t, 1.0, 1e-6))

    collided, _, _, _ = swept_seg_seg_collision(0.1, *args)
    ok &= report("criterion: approaching pair does not collide", not collided)

    # and it does collide once mindist is opened up past the gap
    collided, _, _, _ = swept_seg_seg_collision(1.5, *args)
    ok &= report("criterion: same pair collides at mindist=1.5", collided)
    return ok


def test_retroactive_versus_predictive():
    """the case the two intervals exist to separate

    A pair that did not meet during the step just taken, but will during
    the next one if the velocities hold. The retroactive interval must
    stay silent and the predictive one must fire; a rule that ran only the
    first would miss this pair entirely.
    """
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    b0, b1 = v(0, -1, 0), v(0, 1, 0)      # positions of segment 2 at t
    # segment 2 descends 2 per step: was at z=3, is at z=1, will be at z=-1
    prev = [b + v(0, 0, 2) for b in (b0, b1)]
    nxt = [b - v(0, 0, 2) for b in (b0, b1)]
    now = [b + v(0, 0, 1) for b in (b0, b1)]

    retro, _, _, _ = swept_seg_seg_collision(
        0.1, a0, a0, a1, a1, prev[0], now[0], prev[1], now[1])
    ok = report("criterion: retroactive interval stays silent", not retro)

    pred, _, l1, l2 = swept_seg_seg_collision(
        1.0e-6, a0, a0, a1, a1, now[0], nxt[0], now[1], nxt[1])
    ok &= report("criterion: predictive interval fires on the same pair",
                 pred and 0.0 <= l1 <= 1.0 and 0.0 <= l2 <= 1.0)
    return ok


def test_lines_meet_but_segments_do_not():
    """the L-outside-[0,1] rejection

    Two segments whose infinite lines intersect, positioned so the meeting
    point lies beyond the end of one of them. Distance alone would accept
    this pair; the L test is what rejects it.
    """
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    # along y at x=5, far outside segment 1's span, and sweeping through
    # the z=0 plane so the lines do meet at some instant
    b0t, b0tau = v(5, -1, 0.5), v(5, -1, -0.5)
    b1t, b1tau = v(5, 1, 0.5), v(5, 1, -0.5)
    collided, dist2, _, _ = swept_seg_seg_collision(
        0.1, a0, a0, a1, a1, b0t, b0tau, b1t, b1tau)
    return report("criterion: intersecting lines, disjoint segments, no "
                  "collision", not collided and dist2 > 0.1 * 0.1)


def test_filter_rejects_distant_pairs():
    """the sphere filter, and that it never rejects a real collision"""
    far = segs_may_approach((v(0, 0, 0), v(1, 0, 0), v(0, 0, 0), v(1, 0, 0)),
                            (v(0, 0, 50), v(1, 0, 50), v(0, 0, 50),
                             v(1, 0, 50)), mindist=1.0)
    ok = report("filter: rejects a pair 50 apart", not far)

    near = segs_may_approach((v(0, 0, 0), v(1, 0, 0), v(0, 0, 0), v(1, 0, 0)),
                             (v(0, 0, 0.5), v(1, 0, 0.5), v(0, 0, 0.5),
                              v(1, 0, 0.5)), mindist=1.0)
    ok &= report("filter: admits a pair 0.5 apart", near)

    # the filter must be conservative: anything it rejects must genuinely
    # be further apart than mindist for the whole interval, so a pair that
    # the full criterion accepts must never be filtered out
    a0, a1 = v(-1, 0, 0), v(1, 0, 0)
    b = ((0, -1, 0.5), (0, -1, -0.5), (0, 1, 0.5), (0, 1, -0.5))
    b0t, b0tau, b1t, b1tau = (v(*x) for x in b)
    admitted = segs_may_approach((a0, a1, a0, a1),
                                 (b0t, b1t, b0tau, b1tau), mindist=0.1)
    ok &= report("filter: admits the crossing pair it must not reject",
                 admitted)
    return ok


def test_hinge_angle():
    """the hinge cosine, and the thresholds the caller applies to it"""
    def arms(theta):
        return (v(0, 0, 0), v(1, 0, 0),
                v(np.cos(theta), np.sin(theta), 0))

    # by construction the cosine of the angle between the two arms
    for theta in (0.05, 0.3, 1.0):
        c = hinge_cos_angle(*arms(theta))
        if not close(c, np.cos(theta)):
            return report("hinge: cosine matches cos(theta)", False)
    ok = report("hinge: cosine matches cos(theta) for several angles", True)

    tight = hinge_cos_angle(*arms(0.1))      # cos = 0.995
    loose = hinge_cos_angle(*arms(0.3))      # cos = 0.955
    wide = hinge_cos_angle(*arms(0.8))       # cos = 0.697

    ok &= report("hinge: a tight hinge passes the 0.98 threshold",
                 tight > HINGE_TOL)
    ok &= report("hinge: a 0.3 rad hinge fails 0.98 but passes the "
                 "triangle-relaxed 0.9",
                 loose < HINGE_TOL and loose > HINGE_TOL_TRI)
    ok &= report("hinge: a wide hinge fails both thresholds",
                 wide < HINGE_TOL_TRI)

    # a zero-length arm has no direction, so no threshold may accept it
    ok &= report("hinge: zero-length arm yields 0.0",
                 close(hinge_cos_angle(v(0, 0, 0), v(0, 0, 0), v(1, 0, 0)),
                       0.0))
    return ok


def main():
    tests = [test_point_seg, test_seg_seg, test_enclosing_sphere,
             test_point_point_in_time, test_static_limit,
             test_crossing_during_interval, test_approaching_but_clear,
             test_retroactive_versus_predictive,
             test_lines_meet_but_segments_do_not,
             test_filter_rejects_distant_pairs, test_hinge_angle]
    ok = True
    for t in tests:
        ok &= bool(t())
    print("")
    verdict("test_swept_distance", ok)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
