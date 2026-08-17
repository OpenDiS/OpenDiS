"""Full-run comparison: dislocation loop with elasticity, PyDiS vs ExaDiS.

This is a wrapper. It does not reimplement the simulation. It runs the two
example scripts in examples/01_loop/ as subprocesses and compares the
configurations they write, so the examples stay the single source of truth for
how the run is set up.

The two must agree because they are configured to compute the same thing:
  - PyDiS  Elasticity_SBA  with cutoff = 0.25*Lbox
  - ExaDiS CUTOFF_MODEL    with the same cutoff and Ec = 0.0
ExaDiS' default force mode is DDD_FFT_MODEL, which legitimately does NOT
agree, because it adds the long-range FFT contribution PyDiS does not have.
Hence --force-mode=CUTOFF_MODEL below.

Output is written under the working directory this test is run from, not next
to the sources.
"""

import sys
from pathlib import Path

# this script lives in tests/full_runs/02_frank_read_src/, so the repository
# root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
sys.path.append(str(opendis_root / 'python'))
from framework.testing import report, run_script
from framework.testing import compare_configs, compare_to_ref

examples_dir = opendis_root / 'examples' / '01_loop'

# Tolerance on the node-to-node agreement, in units of b, in a box of 1000.
# Set well above the measured round-off: summation order in the threaded ExaDiS
# kernels is not deterministic, so an exact figure would be flaky, while
# anything approaching this bound would mean a real divergence.
TOL = 1.0e-6

MAX_STEP = 200

# the examples print a progress line every PRINT_FREQ steps; keep the test
# output short without silencing it entirely, so a hung run is still visible
PRINT_FREQ = 100

# steps between intermediate configuration dumps. The comparison uses only the
# final configuration, written after the run regardless of this setting, so the
# intermediate dumps are pure noise here. MAX_STEP writes one at the last step.
WRITE_FREQ = MAX_STEP

PYDIS_SCRIPT = 'test_disl_loop_pydis_elast.py'
EXADIS_SCRIPT = 'test_disl_loop_exadis_elast.py'

# how each example is invoked. --force-mode=CUTOFF_MODEL is what makes the
# exadis run comparable; see the note at the top of this file.
PYDIS_ARGS = ['--no-plot',
              '--max-step', MAX_STEP,
              '--print-freq', PRINT_FREQ,
              '--write-freq', WRITE_FREQ]
EXADIS_ARGS = ['--no-plot',
               '--force-mode=CUTOFF_MODEL',
               '--max-step', MAX_STEP,
               '--print-freq', PRINT_FREQ,
               '--write-freq', WRITE_FREQ]

PYDIS_JSON = Path('output') / 'disl_loop_pydis_elast_final.json'
EXADIS_JSON = Path('output') / 'disl_loop_exadis_elast_final.json'

# Both codes are checked against a stored reference rather than only against
# each other. Agreeing with each other is a weaker statement: a change that
# shifted both equally would pass unnoticed, which is exactly how the maxseg
# bound defect behaved, both codes dropping the same pairs and agreeing while
# both were wrong. Regenerate with 'make disl_loop_elast_ref'.
# anchored to the script, not the working directory: the reference is input
# data belonging to this test, whereas output/ belongs to whoever ran it
REF_NPZ = (Path(__file__).resolve().parent / 'ref_data'
           / 'disl_loop_elast_ref.npz')


def main():
    ok = True

    ok &= report("run pydis  example",
                 run_script(examples_dir / PYDIS_SCRIPT, PYDIS_ARGS))
    ok &= report("run exadis example",
                 run_script(examples_dir / EXADIS_SCRIPT, EXADIS_ARGS))
    if not ok:
        print("an example script failed; skipping the comparison")
        return False

    for f in (PYDIS_JSON, EXADIS_JSON):
        if not f.is_file():
            print("expected output not found: %s" % f)
            return False

    n_pydis, n_exadis, d_pair = compare_configs(PYDIS_JSON, EXADIS_JSON)
    print("nodes: pydis = %d, exadis = %d" % (n_pydis, n_exadis))
    print("max nearest-node distance, pydis vs exadis = %.4e" % d_pair)

    ok &= report("node counts agree", n_pydis == n_exadis)

    if not REF_NPZ.is_file():
        print("reference not found: %s" % REF_NPZ)
        print("generate it with 'make disl_loop_elast_ref' and copy it "
              "into ref_data/")
        return False

    n_run, n_ref, d_pydis = compare_to_ref(PYDIS_JSON, REF_NPZ)
    print("max nearest-node distance, pydis  vs ref = %.4e (nodes %d vs %d)"
          % (d_pydis, n_run, n_ref))
    n_run, n_ref, d_exadis = compare_to_ref(EXADIS_JSON, REF_NPZ)
    print("max nearest-node distance, exadis vs ref = %.4e (nodes %d vs %d)"
          % (d_exadis, n_run, n_ref))

    ok &= report("pydis  matches ref within %.1e" % TOL, d_pydis < TOL)
    ok &= report("exadis matches ref within %.1e" % TOL, d_exadis < TOL)
    return bool(ok)


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
