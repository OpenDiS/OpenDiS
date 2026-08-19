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
            cosine = float(np.dot(self._normals[0], normal))
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


def constrained_plane_point(point, normals, offsets):
    """constrained_plane_point: nearest point satisfying every plane

    Minimizes |x - point| subject to normals @ x == offsets, by solving
    the Lagrange system

        [ I    N^T ] [ x      ]   [ point   ]
        [ N    0   ] [ lambda ] = [ offsets ]

    With no planes the point is returned unchanged. A singular system
    means the constraints cannot be met, and the point is likewise
    returned unchanged, which is what the C does when its inversion
    fails.
    """
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
        return np.linalg.solve(mat, rhs)[:3]
    except np.linalg.LinAlgError:
        return point.copy()


def plane_violation(normal, direction, length):
    """plane_violation: how far a segment departs from its glide plane

    |n . t| * L, the out-of-plane extent of a segment of the given length
    running along direction (a unit vector). Zero when the segment lies in
    the plane.
    """
    return abs(float(np.dot(normal, direction))) * length
