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


def select_independent_planes(candidates, max_planes=3, tol=PLANE_COND,
                              allow_first=True):
    """select_independent_planes: keep a well-conditioned subset

    candidates is an iterable of (normal, offset) with normal a unit
    vector and offset the value of normal.x on the plane. Returns
    (normals, offsets) as a (k, 3) and a (k,) array, k <= max_planes.

    Greedy and order-dependent by construction, matching the C: each
    candidate is accepted only if it is sufficiently independent of those
    already held, so an earlier plane can exclude a later one.

    allow_first=False rejects a candidate that would be the first
    accepted. That reproduces a quirk of the C, which omits the
    zero-planes case when scanning the second node of a collision, so a
    plane offered there is rejected when the first node offered none. It
    is preserved rather than corrected because the point of this
    transcription is agreement with the code being compared against; see
    the plan for the reasoning.
    """
    npc2 = tol * tol
    onemnpc4 = (1.0 - tol) ** 4

    normals = np.zeros((3, 3))
    offsets = np.zeros(3)
    n = 0

    for normal, offset in candidates:
        if n >= max_planes:
            break
        if n == 0:
            accept = allow_first
        elif n == 1:
            cosine = float(np.dot(normals[0], normal))
            accept = cosine * cosine < npc2
        else:
            trial = normals.copy()
            trial[2] = normal
            det = float(np.linalg.det(trial))
            accept = det * det > onemnpc4
        if accept:
            normals[n] = normal
            offsets[n] = offset
            n += 1

    return normals[:n], offsets[:n]


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
