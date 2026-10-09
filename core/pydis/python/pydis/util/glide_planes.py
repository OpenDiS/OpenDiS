"""Glide-plane algebra: picking independent planes and projecting onto them.

Transcribed from AdjustCollisionPoint in ParaDiS RetroactiveCollision2.c
and the ExaDiS port in
core/exadis/src/collision_types/collision_retroactive.cpp.

Both functions take plain arrays and know nothing about dislocation
networks, so they serve any caller that has to move a node while
respecting the planes its arms glide on. Collision is the first; merging
and cross-slip want the same thing.
"""

import numpy as np

# Bitwise-reproducible replacements for the numpy calls below, and the flag
# saying whether to use them. Imported from calforce because that is where the
# other half of this machinery already lives; it belongs somewhere neutral once
# a third caller wants it.
from ..calforce.bitrepro_math import ENABLED as BITREPRO_MATH, dot3, solve_repro

# Independence threshold from the C, where it is called newplanecond. Two
# planes count as distinct when |n1.n2| < 0.875, and a third is accepted
# only when the determinant of the three normals is comfortably non-zero.
PLANE_COND = 0.875


class PlaneSet:
    """PlaneSet: a well-conditioned set of glide-plane constraints

    Accepts candidate planes one at a time and keeps only those
    sufficiently independent of the ones already held, up to max_planes.
    Greedy and order-dependent by construction, matching the C: an earlier
    plane can exclude a later one.

    It is an accumulator rather than a filter over a list because the
    constraints come from two nodes in turn and each candidate must be
    tested against everything accepted so far, from either node. Testing
    the two nodes independently lets in contradictory parallel planes,
    which makes the constrained solve singular and silently leaves the
    point where it started.
    """

    def __init__(self, max_planes=3, tol=PLANE_COND):
        self.max_planes = max_planes
        self.npc2 = tol * tol
        self.onemnpc4 = (1.0 - tol) ** 4
        self._normals = np.zeros((3, 3))
        self._offsets = np.zeros(3)
        self.n = 0

    @property
    def normals(self):
        return self._normals[:self.n]

    @property
    def offsets(self):
        return self._offsets[:self.n]

    def accepts(self, normal, allow_first=True):
        """accepts: would this plane add an independent constraint

        allow_first=False rejects a candidate that would be the first
        accepted. That reproduces a quirk of the C, whose plane-selection
        switch omits the zero-planes case when scanning the second node of
        a collision, so a plane offered there is rejected outright if the
        first node offered none. Preserved rather than corrected because
        agreement with the code being compared against is the point; see
        the plan.
        """
        if self.n >= self.max_planes:
            return False
        if self.n == 0:
            return allow_first
        if self.n == 1:
            cosine = float(dot3(self._normals[0], normal) if BITREPRO_MATH
                           else np.dot(self._normals[0], normal))
            return cosine * cosine < self.npc2
        trial = self._normals.copy()
        trial[2] = normal
        det = float(np.linalg.det(trial))
        return det * det > self.onemnpc4

    def add(self, normal, offset, allow_first=True):
        """add: keep this plane if it is independent enough; True if kept"""
        if not self.accepts(normal, allow_first=allow_first):
            return False
        self._normals[self.n] = normal
        self._offsets[self.n] = offset
        self.n += 1
        return True

    def extend(self, candidates, allow_first=True):
        """extend: offer a sequence of (normal, offset) pairs in order"""
        for normal, offset in candidates:
            self.add(normal, offset, allow_first=allow_first)
        return self


def select_independent_planes(candidates, max_planes=3, tol=PLANE_COND,
                              allow_first=True):
    """select_independent_planes: PlaneSet over one sequence, as arrays

    Convenience for the single-group case. Where constraints come from
    more than one source, use PlaneSet directly so that every candidate is
    tested against every plane already accepted.
    """
    planes = PlaneSet(max_planes=max_planes, tol=tol)
    planes.extend(candidates, allow_first=allow_first)
    return planes.normals, planes.offsets


def constrained_plane_point(point, normals, offsets,
                            use_exadis_projection=False):
    """constrained_plane_point: nearest point satisfying every plane

    Minimizes |x - point| subject to normals @ x == offsets, by solving
    the Lagrange system

        [ I    N^T ] [ x      ]   [ point   ]
        [ N    0   ] [ lambda ] = [ offsets ]

    With no planes the point is returned unchanged. A singular system
    means the constraints cannot be met, and the point is likewise
    returned unchanged, which is what the C does when its inversion
    fails.

    use_exadis_projection routes the same system through
    exadis_projected_point instead, which forms the inverse rather than
    solving. See the note there; it is not the default, and no caller
    passes it today.
    """
    if use_exadis_projection:
        return exadis_projected_point(point, normals, offsets)

    point = np.asarray(point, dtype=float)
    k = len(offsets)
    if k == 0:
        return point.copy()

    size = 3 + k
    mat = np.zeros((size, size))
    mat[:3, :3] = np.eye(3)
    mat[3:, :3] = normals
    mat[:3, 3:] = np.asarray(normals).T
    rhs = np.concatenate([point, np.asarray(offsets, dtype=float)])

    try:
        if BITREPRO_MATH:
            return solve_repro(mat, rhs)[:3]
        return np.linalg.solve(mat, rhs)[:3]
    except np.linalg.LinAlgError:
        return point.copy()


# ExaDiS' version of the same projection, selectable by the flag below.
#
# AdjustCollisionPoint (core/exadis/src/collision_types/collision_retroactive.cpp:978)
# builds the same Lagrange system constrained_plane_point does, on a 6x6 padded
# to a fixed lda whatever the plane count, then solves it by forming the
# explicit inverse with Gauss-Jordan (MatrixInvert, :873) and multiplying,
# rather than solving directly. Its pivot rule is unusual: it swaps rows only
# when the diagonal falls below 1e-12, and then takes the first row with a
# larger magnitude rather than the largest.
#
# Off by default, for two measured reasons.
#
#   1. It does not bring the two codes together. On the systems
#      tests/unit_tests/test4_collision_mode actually produces, the two
#      projections agree to the last bit: with the flag on, that test's whole
#      output is byte-identical, residual against ExaDiS still 6.6159e-12. So
#      whatever separates the two codes, it is not this. (The projections are
#      not identical in general: on random single-plane systems they do differ,
#      by up to a few times 1e-13.)
#
#   2. Accuracy is a wash, so it decides nothing either way. Worst error against
#      a long-double reference, coordinates ~500, over thousands of random
#      systems:
#
#          planes   direct solve   ExaDiS' inverse
#          k=1      5.7310e-13     5.6555e-13
#          k=2      3.9380e-12     4.6678e-12
#          k=3      1.5098e-01     1.5098e-01
#
#      k=1 ties, k=2 slightly favours the direct solve, and at k=3 both are
#      swamped by the conditioning of three near-coplanar normals. Only k=1
#      arises in practice: all 117 merges of that test keep exactly one plane.
#
# Reached by passing use_exadis_projection=True to constrained_plane_point.


def exadis_projected_point(point, normals, offsets):
    """exadis_projected_point: AdjustCollisionPoint's projection

    Same result as constrained_plane_point up to rounding, by forming the
    inverse rather than solving. See the note above for why it is not the
    default.
    """
    point = np.asarray(point, dtype=float)
    k = len(offsets)
    if k == 0:
        return point.copy()
    mat = [[0.0]*6 for _ in range(6)]
    mat[0][0] = mat[1][1] = mat[2][2] = 1.0
    n = np.asarray(normals, dtype=float).reshape(k, 3)
    for i in range(k):
        for j in range(3):
            mat[3+i][j] = float(n[i][j])
            mat[j][3+i] = float(n[i][j])
    rhs = [float(point[0]), float(point[1]), float(point[2]), 0.0, 0.0, 0.0]
    for i in range(k):
        rhs[3+i] = float(offsets[i])
    order = k + 3
    inv = _gauss_jordan_inverse(mat, order)
    if inv is None:
        return point.copy()
    out = []
    for i in range(3):
        acc = 0.0
        for j in range(order):
            acc = acc + inv[i][j]*rhs[j]
        out.append(acc)
    return np.array(out)


def _gauss_jordan_inverse(mat, order, lda=6):
    """_gauss_jordan_inverse: exadis MatrixInvert, operation for operation

    Returns None where the C returns 0, so the caller falls back the same way.
    """
    eps = 1.0e-12
    t = [[float(mat[i][j]) for j in range(lda)] for i in range(lda)]
    inv = [[float(i == j) for j in range(lda)] for i in range(lda)]
    for i in range(order):
        fmax = abs(t[i][i])
        if fmax < eps:            # swap only on a tiny pivot, first larger row
            for j in range(i+1, order):
                if abs(t[j][i]) > fmax:
                    fmax = abs(t[j][i])
                    for k in range(order):
                        t[i][k], t[j][k] = t[j][k], t[i][k]
                        inv[i][k], inv[j][k] = inv[j][k], inv[i][k]
                    break
        if fmax < eps:
            return None
        fval = 1.0 / t[i][i]
        for j in range(order):
            t[i][j] = t[i][j] * fval
            inv[i][j] = inv[i][j] * fval
        for k in range(order):
            if k != i:
                fval = t[k][i]
                for j in range(order):
                    t[k][j] = t[k][j] - fval*t[i][j]
                    inv[k][j] = inv[k][j] - fval*inv[i][j]
    return inv


def plane_violation(normal, direction, length):
    """plane_violation: how far a segment departs from its glide plane

    |n . t| * L, the out-of-plane extent of a segment of the given length
    running along direction (a unit vector). Zero when the segment lies in
    the plane.
    """
    d = dot3(normal, direction) if BITREPRO_MATH else np.dot(normal, direction)
    return abs(float(d)) * length
