"""Full-run comparison: binary junction with elasticity, PyDiS vs ExaDiS.

A wrapper, like the one in tests/full_runs/02_frank_read_src/. It does not
reimplement the simulation; it runs the two example scripts in
examples/03_binary_junction/ as subprocesses and checks what they produce, so
the examples stay the single source of truth for how the run is set up.

WHAT THIS TEST CAN AND CANNOT DO TODAY

The cross-code comparison is not possible yet, and the reason is structural
rather than a tolerance being too tight.

Two dislocation lines start with a node each at the box centre, so the first
collision merges them into a node with four arms. Splitting that node into two
three-arm nodes joined by a junction segment IS the process this case exists to
model. ExaDiS does it in C++ without difficulty. PyDiS cannot: its
Topology(split_mode='MaxDiss') calls OneNodeForce on any node with four or more
arms, and OneNodeForce is unimplemented for the Elasticity_* force modes, so
test_binary_junction_pydis_elast.py stops at step 0 with NotImplementedError.
Disabling Topology does not help, since the Proximity collision handler reads
the nodeflag_dict that only Topology.init_topology_exemptions creates.

So this test does what can be done now:

  1. both examples must run to completion
  2. each must form a junction by the midpoint of the run (checked from the
     intermediate checkpoint the examples write there): segments carrying
     b1 + b2, bounded by two three-arm nodes. That is a physics assertion,
     not an "it ran" one
  3. the two must agree at that midpoint
  4. UNZIP_STRESS, which turns on for the second half of the run (see the
     examples' own comments on it), must destroy the junction it formed --
     checked from the final state
  5. the two must agree at the end too

Nothing here is optional. An earlier version of this file treated PyDiS as
allowed to fail, because at the time OneNodeForce was unimplemented for the
Elasticity_* modes and PyDiS genuinely could not run the case. That was a
mistake even then: it reported PASSED while PyDiS was crashing, so a real
breakage in PyDiS looked exactly like the known limitation. It has since been
implemented, and the allowance is gone.

Output is written under the working directory this test is run from, not next
to the sources.
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# this script lives in tests/full_runs/03_binary_junction/, so the repository
# root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
sys.path.append(str(opendis_root / 'python'))

import numpy as np
from framework.testing import report, run_script, compare_configs

examples_dir = opendis_root / 'examples' / '03_binary_junction'

MAX_STEP = 300
PRINT_FREQ = 100
WRITE_FREQ = MAX_STEP

# Matches each example's own `stress_step = max_step // 2`: the step at
# which it writes its intermediate checkpoint and UNZIP_STRESS turns on.
STRESS_STEP = MAX_STEP // 2

PYDIS_SCRIPT = 'test_binary_junction_pydis_elast.py'
EXADIS_SCRIPT = 'test_binary_junction_exadis_elast.py'

PYDIS_JSON = Path('output') / 'binary_junction_pydis_elast_final.json'
EXADIS_JSON = Path('output') / 'binary_junction_exadis_elast_final.json'
PYDIS_MID_JSON = (Path('output')
                 / ('binary_junction_pydis_elast_step%d.json' % STRESS_STEP))
EXADIS_MID_JSON = (Path('output')
                  / ('binary_junction_exadis_elast_step%d.json' % STRESS_STEP))

# Burgers vectors of the two lines, and of the junction they form. b1 + b2 is
# what a binary junction carries, and finding it is the check that a junction
# formed rather than the two lines simply passing through or annihilating.
B1 = np.array([-1.0, 1.0, 1.0])
B2 = np.array([1.0, -1.0, 1.0])
B_JUNCTION = B1 + B2

# A junction is bounded by exactly two three-arm nodes, one at each end of the
# junction segment. Fewer means it did not form; more would mean something
# else happened.
N_JUNCTION_ENDS = 2

TOL = 1.0e-6


def junction_summary(json_file):
    """junction_summary: what the run produced, as (n_junction_segs, length,
    n_three_arm)

    A junction segment is one whose Burgers vector is +-(b1+b2); the sign
    depends on which way round the segment was written, which carries no
    physical meaning.
    """
    with open(json_file) as f:
        data = json.load(f)
    pos = np.array(data['nodes']['positions'], dtype=float)
    ids = np.array(data['segs']['nodeids'], dtype=int)
    burg = np.array(data['segs']['burgers'], dtype=float)

    is_junction = np.all(np.isclose(np.abs(burg), np.abs(B_JUNCTION),
                                    atol=1e-6), axis=1)
    length = sum(float(np.linalg.norm(pos[ids[k, 1]] - pos[ids[k, 0]]))
                 for k in np.where(is_junction)[0])
    degree = Counter(ids.flatten())
    n_three_arm = sum(1 for d in degree.values() if d == 3)
    return int(is_junction.sum()), length, n_three_arm


def check_junction_formed(label, json_file):
    """check_junction_formed: require that this run has a junction by now"""
    n_seg, length, n_ends = junction_summary(json_file)
    print("%s: %d junction segments (b = %s), total length %.2f, "
          "%d three-arm nodes"
          % (label, n_seg, np.array2string(B_JUNCTION, precision=0),
             length, n_ends))
    ok = report("%s: a junction formed" % label, n_seg > 0 and length > 0.0)
    ok &= report("%s: bounded by %d three-arm nodes"
                 % (label, N_JUNCTION_ENDS), n_ends == N_JUNCTION_ENDS)
    return ok


def check_junction_destroyed(label, json_file):
    """check_junction_destroyed: require UNZIP_STRESS has undone the junction"""
    n_seg, length, n_ends = junction_summary(json_file)
    print("%s: %d junction segments (b = %s), total length %.2f, "
          "%d three-arm nodes"
          % (label, n_seg, np.array2string(B_JUNCTION, precision=0),
             length, n_ends))
    return report("%s: UNZIP_STRESS destroyed the junction" % label,
                  n_seg == 0 and length == 0.0)


def main(plot=False):
    # --force-mode=CUTOFF_MODEL is what makes the exadis run comparable to
    # pydis' Elasticity_SBA: DDD_FFT_MODEL adds a long-range FFT contribution
    # that pydis has no counterpart for. Both example scripts default to
    # plot=True and expose --no-plot to turn it off, so plot=True here means
    # passing neither flag rather than a --plot one.
    common_args = (['--max-step', MAX_STEP,
                    '--print-freq', PRINT_FREQ,
                    '--write-freq', WRITE_FREQ]
                   + ([] if plot else ['--no-plot']))
    pydis_args = list(common_args)
    exadis_args = ['--force-mode=CUTOFF_MODEL'] + common_args

    ok = True

    ok &= report("run exadis example",
                 run_script(examples_dir / EXADIS_SCRIPT, exadis_args,
                           headless=not plot))
    if not ok:
        print("the exadis example failed; nothing else can be checked")
        return False
    if not EXADIS_JSON.is_file() or not EXADIS_MID_JSON.is_file():
        print("expected output not found under %s" % EXADIS_JSON.parent)
        return False

    ok &= check_junction_formed("exadis", EXADIS_MID_JSON)
    ok &= check_junction_destroyed("exadis", EXADIS_JSON)

    ok &= report("run pydis  example",
                 run_script(examples_dir / PYDIS_SCRIPT, pydis_args,
                           headless=not plot))
    print("")
    if not ok:
        print("the pydis example failed, so nothing below can be checked")
        return False
    if not PYDIS_JSON.is_file() or not PYDIS_MID_JSON.is_file():
        print("expected output not found under %s" % PYDIS_JSON.parent)
        return False

    ok &= check_junction_formed("pydis ", PYDIS_MID_JSON)
    ok &= check_junction_destroyed("pydis ", PYDIS_JSON)

    n_pydis, n_exadis, d_mid = compare_configs(PYDIS_MID_JSON, EXADIS_MID_JSON)
    print("at step %d: nodes: pydis = %d, exadis = %d"
          % (STRESS_STEP, n_pydis, n_exadis))
    print("at step %d: max nearest-node distance, pydis vs exadis = %.4e"
          % (STRESS_STEP, d_mid))
    ok &= report("node counts agree at step %d" % STRESS_STEP,
                 n_pydis == n_exadis)
    ok &= report("configurations agree within %.1e at step %d"
                 % (TOL, STRESS_STEP), d_mid < TOL)

    n_pydis, n_exadis, d_pair = compare_configs(PYDIS_JSON, EXADIS_JSON)
    print("at step %d: nodes: pydis = %d, exadis = %d"
          % (MAX_STEP, n_pydis, n_exadis))
    print("at step %d: max nearest-node distance, pydis vs exadis = %.4e"
          % (MAX_STEP, d_pair))
    ok &= report("node counts agree at step %d" % MAX_STEP,
                 n_pydis == n_exadis)
    ok &= report("configurations agree within %.1e at step %d"
                 % (TOL, MAX_STEP), d_pair < TOL)

    # No stored reference yet. Blessing one now would freeze the exadis-only
    # behaviour of a case whose point is the comparison, and it would have to
    # be regenerated the moment pydis can run. Add it, and the ref_data
    # machinery of 02_frank_read_src, once both codes complete the run.
    print("")
    print("note: no stored reference for this case yet, so both codes are "
          "checked against")
    print("      each other and against the junction assertions, not against "
          "a blessed run")
    return bool(ok)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', dest='plot', action='store_true',
                        help='let the two example scripts show their plots '
                             '(off by default here, unlike running them '
                             'directly)')
    args = parser.parse_args()

    passed = main(plot=args.plot)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
