"""Segment-segment interaction forces, computed by PyDiS.

A stored table of 100 segment pairs, each row holding two segments
(p1,p2) and (p3,p4) with Burgers vectors b12 and b34, followed by the
reference forces f1..f4 on the four end nodes. The pairs span
separations from 2.5 to 32.5 with random Burgers vectors, at a = 0.01.

Two implementations of the same isotropic non-singular kernel are run
over that table and compared against the stored forces:

  A  compute_segseg_force_vec     the compiled ParaDiS library (SBA)
  B  python_segseg_force_vec      the python translation from Matlab

so the test pins the compiled path to the reference and the two
implementations to each other. test_segseg_force_exadis.py checks the
ExaDiS kernel against the same table; A is the case it corresponds to.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in
               ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in pydis_paths if not p in sys.path]
ref_dir = Path(__file__).resolve().parent / 'ref_data'

import numpy as np
from pydis.calforce.compute_stress_force_analytic_paradis import (
    compute_segseg_force_vec)
from pydis.calforce.compute_stress_force_analytic_python import (
    python_segseg_force_vec)
from framework.testing import report_close

# the segment pairs and their reference forces: columns 0:18 are
# p1, p2, p3, p4, b12, b34 and 18:30 are f1..f4
REF_FILE = 'segsep_min_2.5_max_32.5_iso_randombvecs_a0.010.dat'

# the constants the reference was generated with
mu = 50.0
nu = 0.3
a = 0.01
force_nint = 3

tolA, tolB, tolC = 1e-9, 1e-9, 1e-5


def main():
    print("segment-segment forces from two pydis implementations, "
          "against the stored reference '%s'" % REF_FILE)
    segseg_data = np.loadtxt(ref_dir / REF_FILE)
    print("%d segment pairs, mu = %g, nu = %g, a = %g"
          % (segseg_data.shape[0], mu, nu, a))

    p1 = segseg_data[:, 0:3]
    p2 = segseg_data[:, 3:6]
    p3 = segseg_data[:, 6:9]
    p4 = segseg_data[:, 9:12]

    b12 = segseg_data[:, 12:15]
    b34 = segseg_data[:, 15:18]
    f1234_ref = segseg_data[:, 18:]

    # Test A: use ParaDiS library (SBA)
    fA = compute_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    okA = report_close("TestA, ParaDiS library (SBA)",
                       np.concatenate(fA, axis=1), f1234_ref, tolA)

    # Test B: use Python code (translated from Matlab)
    fB = python_segseg_force_vec(p1, p2, p3, p4, b12, b34, mu, nu, a)
    okB = report_close("TestB, python translation  ",
                       np.concatenate(fB, axis=1), f1234_ref, tolB)

    # Test C: use ParaDiS library (SBN1)
    #quad_points, weights = np.polynomial.legendre.leggauss(force_nint)
    #fC = compute_segseg_force_SBN1_vec(p1, p2, p3, p4, b12, b34,
    #                                   mu, nu, a, quad_points, weights)
    #okC = report_close("TestC, ParaDiS library (SBN1)",
    #                   np.concatenate(fC, axis=1), f1234_ref, tolC)

    return bool(okA and okB)


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
