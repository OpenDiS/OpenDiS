"""bitrepro_math: BLAS-free replacements for small-vector numpy calls.

np.dot() and np.linalg.norm() on a float64 array dispatch straight to
whatever BLAS the environment's numpy is linked against -- MKL, Accelerate,
OpenBLAS, ... -- even for a 3-element vector (confirmed via
MKL_VERBOSE=1: a plain np.dot(a, b) on length-3 arrays shows up as an MKL
DDOT(3, ...) call). Different BLAS vendors round their dot/norm
differently by about a ULP, which is invisible on its own but shows up
once a caller multiplies the result by something large, e.g. the
Ec ~ 1e6 core-energy term in selfforcevec_LineTension.

dot3()/norm3() below use only +, * and math.sqrt in a fixed, left-to-right
order, so no BLAS library sits between the call and the answer -- the same
reasoning as portable_math.c on the C side, one layer up the stack.

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
