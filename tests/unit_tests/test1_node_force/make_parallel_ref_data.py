"""Generate the parallel / near-parallel segment-pair reference table.

The original table (segsep_min_2.5_max_32.5_iso_randombvecs_a0.010.dat) has
100 rows whose smallest 1-c^2 is 3.4e-03, so nothing in it comes within an
order of magnitude of the eps = 1e-4 threshold at which the kernels switch
from the general formula to the parallel one. Three bugs lived in that
region of pydis' python kernel undetected because of it.

This table straddles that threshold deliberately:

  1-c^2 = 0            exactly parallel and antiparallel
  1e-6, 1e-5           below eps: parallel branch
  5e-5, 9e-5           still below, approaching eps
  1.1e-4, 2e-4         just above, general branch
  1e-3, 1e-2           further above, where the general formula is still
                       ill-conditioned

and for each, one pair with the second segment offset along its own
direction, which is what reaches the near-parallel correction branches in
SpecialRemoteNodeForce.

Note what is deliberately absent: no row sits *on* 1-c^2 = 1e-4. The switch
between the two formulas is discontinuous, the two disagreeing by ~7.7e-05
where they meet, so a pair landing exactly on the threshold has an answer
that moves by 2.6e-04 under a one-ulp change of an input coordinate. Such a
row tests nothing: whether it passes depends on which side of the switch a
given compiler and libm happen to put it, and when it flips it fails for all
three implementations at once and looks like a force bug. Straddling the
threshold means rows either side of it, not on it. 9e-5 and 1.1e-4 are close
enough to exercise both formulas near their common edge and far enough that
float64 keeps them on the side they were generated for.

Reference values are taken from the compiled ParaDiS kernel, but only after
checking they agree with exadis, an independent implementation. Any row
where the two disagree is reported and not blessed: that is a finding, not
something to average over.

Two files come out, because the table is two things. The geometry is a
fixture and goes straight to input_data/ where the test reads it. The forces
are a reference and go to output/, to be inspected and copied into ref_data/
by hand, like every other reference in this suite.

    python3 make_parallel_ref_data.py
    cp output/parallel_and_nearparallel_iso_a0.010.npz ref_data/
"""

import sys
from pathlib import Path

opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]
here = Path(__file__).resolve().parent
input_dir = here / 'input_data'
out_dir = Path('output')          # relative to the working directory

import numpy as np
from pydis.calforce.compute_stress_force_analytic_paradis import (
    compute_segseg_force_list)
from framework.testing import array_digest

# The geometry and the forces are written separately, following the split the
# rest of this folder uses: input_data/ holds what a test reads in, ref_data/
# holds what it checks against.
#
# The .dat is the fixture. It defines which segment pairs the test covers and
# nothing about it needs blessing, so it goes straight to input_data/.
#
# The .npz is the reference, and it goes to output/ for you to inspect and
# copy into ref_data/ by hand, like every other reference here. Regenerating
# the fixture without copying the new reference would leave the two
# describing different pairs, so the .npz records the row count and a hash of
# the geometry it was computed from, and the test refuses to run if they
# disagree.
STEM = 'parallel_and_nearparallel_iso_a0.010'
INPUT_FILE = STEM + '.dat'
REF_FILE = STEM + '.npz'

mu, nu, a = 50.0, 0.3, 0.01
BLESS_TOL = 1e-10          # C and exadis must agree this well to bless a row

# The threshold at which every kernel switches formula, and how close to it a
# row is allowed to sit. The switch is discontinuous, so a row inside this
# band has an answer that depends on which side of it a given compiler puts
# the pair; the check below refuses to write a table containing one. Exactly
# parallel is exempt: it is on the parallel side by construction, not by a
# rounding accident, and it is the one geometry where the formula is exact.
SWITCH_EPS = 1e-4
SWITCH_GUARD = 0.1         # keep |1-c^2 - eps| > 10% of eps


def build_table():
    """build_table: the segment pairs, as (N,18) of p1,p2,p3,p4,b12,b34"""
    rng = np.random.default_rng(20260816)
    # 1e-4 itself is excluded on purpose; see the module docstring
    onemc2 = [0.0, 0.0, 1e-6, 1e-5, 5e-5, 9e-5, 1.1e-4, 2e-4, 1e-3, 1e-2]
    rows = []
    for k, o in enumerate(onemc2):
        anti = (k == 1)                      # the second 0.0 is antiparallel
        c = -1.0 if anti else np.sqrt(max(0.0, 1.0 - o))
        for sep in (2.5, 12.0, 32.5):
            # 7.0 is what reaches the correction branch
            for along in (0.0, 7.0):
                t = np.array([1.0, 0.0, 0.0])
                n = np.array([0.0, 1.0, 0.0])
                tp = c*t + np.sqrt(o)*n
                tp = tp/np.linalg.norm(tp)
                L1, L2 = 10.0, 10.0
                p1 = np.zeros(3)
                p2 = p1 + L1*t
                p3 = p1 + along*t + sep*n + 0.3*np.array([0.0, 0.0, 1.0])
                p4 = p3 + L2*tp
                b12 = rng.normal(size=3); b12 /= np.linalg.norm(b12)
                b34 = rng.normal(size=3); b34 /= np.linalg.norm(b34)
                rows.append(np.concatenate([p1, p2, p3, p4, b12, b34]))
    return np.array(rows)


def exadis_forces(pairs):
    """exadis_forces: the same pairs through exadis' SegSegIso, as (N,12)"""
    import pyexadis
    from pyexadis_base import ExaDisNet
    p1, p2, p3, p4, b12, b34 = pairs
    n = p1.shape[0]
    rn = np.stack([p1, p2, p3, p4], axis=1).reshape(4*n, 3)
    segs = np.zeros((2*n, 8))
    segs[0::2, 0] = np.arange(n)*4 + 0
    segs[0::2, 1] = np.arange(n)*4 + 1
    segs[0::2, 2:5] = b12
    segs[1::2, 0] = np.arange(n)*4 + 2
    segs[1::2, 1] = np.arange(n)*4 + 3
    segs[1::2, 2:5] = b34
    L = 10.0 * float(np.max(rn.max(axis=0) - rn.min(axis=0)))
    origin = 0.5*(rn.min(axis=0) + rn.max(axis=0)) - 0.5*L
    cell = pyexadis.Cell(h=L*np.eye(3), origin=origin,
                         is_periodic=[False, False, False])
    G = ExaDisNet(cell, rn, segs)
    f = pyexadis.compute_force_segseglist(G.net, mu, nu, a,
                                          [[2*i, 2*i+1] for i in range(n)])
    return np.array(f).reshape(n, 12)


def main():
    table = build_table()
    pairs = tuple(table[:, i:i+3] for i in range(0, 18, 3))
    p1, p2, p3, p4 = pairs[0], pairs[1], pairs[2], pairs[3]
    t = (p2-p1)/np.linalg.norm(p2-p1, axis=1)[:, None]
    tp = (p4-p3)/np.linalg.norm(p4-p3, axis=1)[:, None]
    onemc2 = 1.0 - np.einsum('ij,ij->i', t, tp)**2
    print("built %d pairs, 1-c^2 from %.1e to %.1e"
          % (table.shape[0], onemc2.min(), onemc2.max()))

    # refuse to bless a row sitting on the formula switch, however it got
    # there: as generated, or because float64 moved it onto the threshold
    on_switch = np.where((onemc2 > 0.0) &
                         (np.abs(onemc2 - SWITCH_EPS)
                          < SWITCH_GUARD*SWITCH_EPS))[0]
    if on_switch.size:
        print("NOT BLESSED, %d row(s) within %.0f%% of the 1-c^2 = %.0e "
              "formula switch, where the answer is not well defined in "
              "float64:" % (on_switch.size, 100*SWITCH_GUARD, SWITCH_EPS))
        for i in on_switch[:10]:
            print("   row %3d  1-c^2 = %.6e" % (i, onemc2[i]))
        return 1

    fC = np.concatenate(compute_segseg_force_list(*pairs, mu, nu, a), axis=1)

    import pyexadis
    pyexadis.initialize()
    fE = exadis_forces(pairs)
    pyexadis.finalize()

    scale = max(np.abs(fC).max(), 1e-30)
    rel = np.abs(fC - fE).max(axis=1)/scale
    bad = np.where(rel > BLESS_TOL)[0]
    print("C vs exadis: worst relative difference %.3e over %d rows"
          % (rel.max(), len(rel)))
    if bad.size:
        print("NOT BLESSED, %d row(s) disagree beyond %.0e:"
              % (bad.size, BLESS_TOL))
        for i in bad[:10]:
            print("   row %3d  1-c^2 = %.3e  rel diff = %.3e"
                  % (i, onemc2[i], rel[i]))
        return 1

    # the fixture, straight to where the test reads it
    input_dir.mkdir(exist_ok=True)
    np.savetxt(input_dir / INPUT_FILE, table, fmt='%24.16e')

    # hash what the test will actually read, not what is in memory: %24.16e is
    # a round trip, and a digest of the pre-write values could disagree with
    # the post-read ones in the last place
    geometry = np.loadtxt(input_dir / INPUT_FILE)
    digest = array_digest(geometry)

    # the reference, to output/ for a human to bless
    out_dir.mkdir(exist_ok=True)
    np.savez(out_dir / REF_FILE, forces=fC, mu=mu, nu=nu, a=a,
             n_pairs=table.shape[0], geometry_digest=digest,
             source='compute_segseg_force_list, cross-checked against '
                    'exadis SegSegIso to %.0e' % BLESS_TOL)

    print("wrote %s/%s: %d rows x %d geometry columns"
          % (input_dir.name, INPUT_FILE, table.shape[0], table.shape[1]))
    print("wrote %s/%s: %d rows x %d force columns"
          % (out_dir, REF_FILE, fC.shape[0], fC.shape[1]))
    print("")
    print("  the fixture is in place; the reference is not. To bless it:")
    print("      cp %s/%s ref_data/" % (out_dir, REF_FILE))
    print("  until you do, the test will report the fixture and the reference")
    print("  describing different pairs and refuse to compare them.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
