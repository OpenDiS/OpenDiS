"""bitrepro_math: BLAS-free replacements for small-vector numpy calls.

np.dot() and np.linalg.norm() on a float64 array dispatch straight to
whatever BLAS the environment's numpy is linked against -- MKL, Accelerate,
OpenBLAS, ... -- even for a 3-element vector (confirmed via
MKL_VERBOSE=1: a plain np.dot(a, b) on length-3 arrays shows up as an MKL
DDOT(3, ...) call). Different BLAS vendors round their dot/norm
differently by about a ULP, which is invisible on its own but shows up
once a caller multiplies the result by something large, e.g. the
Ec ~ 1e6 core-energy term in selfforcevec_LineTension.

dot3()/norm3()/matvec3()/cross3() below use only +, -, * and math.sqrt in a
fixed, left-to-right order, so no BLAS library sits between the call and the
answer -- the same reasoning as portable_math.c on the C side, one layer up the
stack.

matvec3() and cross3() are here for pkforcevec(), which reaches BLAS through
the @ operator rather than np.dot(). np.cross() does not call BLAS, but it is
written out too: it decides its output dtype and loop order from the shapes it
is handed, and pinning three subtractions costs nothing next to leaving that to
a library.

solve_repro() is the same idea for np.linalg.solve(), which goes to LAPACK.
Used by util/glide_planes.constrained_plane_point for the small KKT system that
projects a collision point onto the intersection of its glide planes.

Defaults to whichever way the compiled library went, via
pydis.build_info.bitrepro_math(): every SYS ending in _repro sets
PYDIS_BITREPRO_MATH ON for the C kernel (core/pydis/c/CMakeLists.txt), and
_repro is meant to turn all of this reproducibility machinery on together,
not just the C half -- there is no separate build-time setting for numpy's
BLAS choice, since that is fixed by the conda environment rather than by
anything pydis configures, so this is the closest a *_repro build gets to
switching it on by itself. The environment variable PYDIS_BITREPRO_MATH
still overrides that default either way, e.g. to exercise one half without
the other.
"""

import os
import math

import numpy as np

from ..build_info import bitrepro_math as _c_bitrepro_math

ENABLED = (os.environ["PYDIS_BITREPRO_MATH"] == "1"
          if "PYDIS_BITREPRO_MATH" in os.environ
          else _c_bitrepro_math())


def dot3(a, b):
    """dot3: a.b for length-3 vectors, without going through BLAS."""
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def norm3(a):
    """norm3: ||a|| for a length-3 vector, without going through BLAS."""
    return math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])


def dot_fixed(a, b):
    """dot_fixed: dot product summed in index order, for any length

    dot3's generalisation. The collision path applies np.dot to 2-vectors as
    well as 3-vectors, the residuals of the small systems it solves for the
    closest-approach parameters, so the drop-in replacement there cannot assume
    a length. Length 3 keeps dot3's written-out form; anything else accumulates
    left to right, which is equally fixed.
    """
    if len(a) == 3:
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    acc = 0.0
    for i in range(len(a)):
        acc = acc + a[i]*b[i]
    return acc


def norm_fixed(a):
    """norm_fixed: ||a|| for any length, without going through BLAS"""
    return math.sqrt(dot_fixed(a, a))


def matvec3(m, v):
    """matvec3: m.v for a 3x3 matrix and a length-3 vector, one dot3 per row

    Replaces the @ operator, which dispatches to BLAS' gemv the same way
    np.dot() dispatches to ddot. Returns a tuple; every caller either indexes
    it or assigns it into a row of an array, and both work as they would with
    an ndarray.
    """
    return (dot3(m[0], v), dot3(m[1], v), dot3(m[2], v))


def det3(m):
    """det3: determinant of a 3x3, by the explicit cofactor sum

    The same expression, in the same order, as exadis' M33_DET macro in
    core/exadis/src/collision_types/collision_retroactive.cpp. np.linalg.det
    factorizes instead, through LAPACK, and the two disagree in the last bits;
    on the Newton iteration of swept_distance.seg_seg_min_dist_in_time that
    disagreement compounds into 1e-13 in the returned ratio. Fixed order, no
    BLAS, no LAPACK.
    """
    return (m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0]
            + m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0]
            - m[0][1]*m[1][0]*m[2][2] - m[0][0]*m[1][2]*m[2][1])


def solve3_by_inverse(m, v):
    """solve3_by_inverse: m^-1 . v, forming the adjugate explicitly

    Follows exadis' M33_INV followed by V3_M33_V3_MUL, operation for
    operation. A singular matrix gives a zero inverse rather than raising,
    which is what that macro does: it replaces 1/det by 0 when det is zero.
    Callers guard on det3 before calling, as the C does.

    Not a general-purpose linear solver, and not the accurate way to solve a
    3x3: forming the inverse loses a little against factorizing. It is here
    because agreeing with the reference to the last bit matters more than the
    last bit of accuracy, and because it avoids LAPACK, whose ordering is not
    fixed across platforms.
    """
    det = det3(m)
    det = 1.0/det if abs(det) > 0.0 else 0.0
    inv = ((det*(m[1][1]*m[2][2] - m[1][2]*m[2][1]),
            det*(m[0][2]*m[2][1] - m[0][1]*m[2][2]),
            det*(m[0][1]*m[1][2] - m[0][2]*m[1][1])),
           (det*(m[1][2]*m[2][0] - m[1][0]*m[2][2]),
            det*(m[0][0]*m[2][2] - m[0][2]*m[2][0]),
            det*(m[0][2]*m[1][0] - m[0][0]*m[1][2])),
           (det*(m[1][0]*m[2][1] - m[1][1]*m[2][0]),
            det*(m[0][1]*m[2][0] - m[0][0]*m[2][1]),
            det*(m[0][0]*m[1][1] - m[0][1]*m[1][0])))
    return matvec3(inv, v)


def cross3(a, b):
    """cross3: a x b for length-3 vectors, component by component."""
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def solve_repro(mat, rhs):
    """solve_repro: dense linear solve with a fixed evaluation order

    Gaussian elimination with partial pivoting, unblocked, in Python floats, so
    every operation is one IEEE-754 add, subtract, multiply or divide in a
    written-down order. The same algorithm class LAPACK's dgesv uses; what is
    removed is the blocking, the vendor-specific ordering and any FMA
    contraction, none of which LAPACK promises to keep constant across
    platforms.

    Not bit-compatible with np.linalg.solve, except that the two agree exactly
    on the 4x4 case this is wanted for. On the 5x5 and 6x6 cases they differ by
    a few ULP while the system is well conditioned, and by more as it
    approaches singular, which is a property of the system rather than of
    either solver.

    Pivoting picks the largest magnitude in the column and the lowest row index
    on a tie, so the pivot sequence is a function of the values alone.

    Raises np.linalg.LinAlgError on a singular matrix, matching what callers
    already catch from np.linalg.solve.
    """
    a = [[float(v) for v in row] for row in np.asarray(mat, dtype=float)]
    b = [float(v) for v in np.asarray(rhs, dtype=float)]
    n = len(b)

    for col in range(n):
        piv, best = col, abs(a[col][col])
        for r in range(col + 1, n):
            if abs(a[r][col]) > best:
                piv, best = r, abs(a[r][col])
        if best == 0.0:
            raise np.linalg.LinAlgError("solve_repro: singular matrix")
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            b[col], b[piv] = b[piv], b[col]
        # one reciprocal per column, so the row loop below multiplies rather
        # than divides, as the blocked libraries also do
        inv = 1.0 / a[col][col]
        for r in range(col + 1, n):
            f = a[r][col] * inv
            if f != 0.0:
                for c in range(col, n):
                    a[r][c] = a[r][c] - f * a[col][c]
                b[r] = b[r] - f * b[col]

    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        acc = b[r]
        for c in range(r + 1, n):
            acc = acc - a[r][c] * x[c]
        x[r] = acc / a[r][r]
    return np.array(x)
