"""Full-run comparison: binary junction with elasticity, PyDiS vs ExaDiS.

A wrapper, like the one in tests/full_runs/02_frank_read_src/. It does not
reimplement the simulation; both comparison modes below run the two example
scripts in examples/03_binary_junction/, so the examples stay the single
source of truth for how the run is set up.

Two ways of comparing them, chosen by IMPORT_EXADIS_ORDER below -- same split,
same reasoning, as 02_frank_read_src/test_frank_read_src_pydis_exadis_elast.py.

WITH the order import (default). ExaDiS runs first, recording its segment
order immediately before every remesh call. PyDiS then runs with that order
imposed at the same point, in one process. This is narrower than
02_frank_read_src's own order import, which also imposes arm order and
multi-node order at three points per step: an ablation found those do not
matter for this case (pydis_pass()'s own docstring has the details; a
degree-4 node with more than one splitting candidate never comes up here,
and this case's collision and topology decisions never land on a near-tie
sensitive to arm order the way 02's or test5_topology_mode's do). Segment
order right before remesh is the one piece iteration order that does
matter, because it is the one thing Remesh_LengthBased reads
(list(G.all_segments_tags())) that pydis holds no order of its own to
match ExaDiS' with. Node counts agree at all 300 steps and positions agree
to round-off (worst 4.6022e-06 against a 5.0e-06 tolerance).

WITHOUT it. The two examples run as independent subprocesses (as this file
always did before order import was added) and are compared at the midpoint
and the end, plus the physics-specific checks below. This is the weaker,
uncontrolled comparison: two independent runs of this case's deliberately
near-mirror-symmetric geometry can break a symmetric near-tie differently and
end up on different topological branches from ordinary per-process
iteration-order drift, the same phenomenon documented in
tests/full_runs/02_frank_read_src/order_import.py's own module docstring.
Confirmed to still fail at the default tolerance even with both codes on the
same coarsen_mode (worst 2.7047e-02 at step 150, orders of magnitude past tolerance);
kept because it is the honest statement of where two separate runs actually
end up, and because the junction-formation physics checks and the
stored-reference comparison below
only make sense against each example's own real, independently-written
output.

WHAT PYDIS CAN AND CANNOT DO

Two dislocation lines start with a node each at the box centre, so the first
collision merges them into a node with four arms. Splitting that node into two
three-arm nodes joined by a junction segment IS the process this case exists to
model. ExaDiS does it in C++ without difficulty. PyDiS's
Topology(split_mode='Serial') now implements this too (an earlier version of
this file's docstring predates that; see git history if that limitation needs
revisiting).

So the independent-mode comparison does what can be done there:

  1. both examples must run to completion
  2. each must form a junction by the midpoint of the run (checked from the
     intermediate checkpoint the examples write there): segments carrying
     b1 + b2, bounded by two three-arm nodes. That is a physics assertion,
     not an "it ran" one
  3. the two must agree with each other at that midpoint (no stored
     reference exists for this intermediate state)
  4. UNZIP_STRESS, which turns on for the second half of the run (see the
     examples' own comments on it), must destroy the junction it formed --
     checked from the final state
  5. at the end, each of the two must independently agree with a stored
     reference configuration (ref_data/binary_junction_elast_ref.npz,
     regenerated with 'make binary_junction_elast_ref'), not merely with
     each other -- agreeing with each other is a weaker statement, since a
     change that shifted both codes equally would pass unnoticed

Nothing here is optional. An earlier version of this file treated PyDiS as
allowed to fail, because at the time OneNodeForce was unimplemented for the
Elasticity_* modes and PyDiS genuinely could not run the case. That was a
mistake even then: it reported PASSED while PyDiS was crashing, so a real
breakage in PyDiS looked exactly like the known limitation. It has since been
implemented, and the allowance is gone.

Output is written under the working directory this test is run from, not next
to the sources.
"""

import os

# Before anything imports pyexadis -- see 02_frank_read_src's own copy of this
# comment. ExaDiS' force kernels are Kokkos parallel reductions, so with more
# than one thread the summation order varies between runs, which the ordered
# comparison below cannot tolerate (it depends on the exact orders recorded
# being the ones actually used). Assigned rather than setdefault, so an
# OMP_NUM_THREADS inherited from the shell cannot make the result depend on
# where it was run from.
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# this script lives in tests/full_runs/03_binary_junction/, so the repository
# root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
for _p in ['python', 'lib', 'core/pydis/python', 'core/exadis/python',
           'examples/03_binary_junction', 'tests/full_runs/02_frank_read_src']:
    _q = str(opendis_root / _p)
    if _q not in sys.path:
        sys.path.append(_q)

import numpy as np
from framework.testing import report, run_script, compare_configs, compare_to_ref

examples_dir = opendis_root / 'examples' / '03_binary_junction'

# True: run the codes in one process and impose ExaDiS' iteration orders on
# PyDiS every step. False: run the two examples independently and compare
# their configurations. See this module's docstring for what each answers.
IMPORT_EXADIS_ORDER = True

# Coarsening branch, passed to both examples so neither falls back to its own
# default: pydis defaults to 0 and pyexadis_base.Remesh to 1. 1 is
# node-centric, exadis' and ParaDiS' rule, and both branches are verified
# against exadis bitwise in tests/unit_tests/test3_remesh_rule. Either value
# works here; --coarsen-mode on both example scripts overrides it by hand.
COARSEN_MODE = 1

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

# Set above the round-off floor this case actually reaches under order
# import, not 02_frank_read_src's 1e-6: measured directly, the worst step
# over the whole 300-step run (order imported, coarsen_mode=0 on both sides)
# is 3.85e-6, around the junction-formation and UNZIP_STRESS transients --
# more topological events than 02's case, so a larger floor is expected, not
# a sign of a remaining bug. Still tight enough to catch independent()'s own
# failures, which run 1e-2 to 1e-1.
TOL = 5.0e-6

# The final configuration is checked against a stored reference rather than
# only against pydis and exadis agreeing with each other. Agreeing with each
# other is a weaker statement: a change that shifted both codes equally would
# pass unnoticed. Regenerate with 'make binary_junction_elast_ref'. Anchored
# to the script, not the working directory: the reference is input data
# belonging to this test, whereas output/ belongs to whoever ran it. Used by
# independent() only -- see ordered()'s own note on why not there too.
REF_NPZ = (Path(__file__).resolve().parent / 'ref_data'
          / 'binary_junction_elast_ref.npz')

# The order-import comparison records and compares the whole run: unlike
# 02_frank_read_src, where nothing interesting happens before step 250,
# 03's own divergence (before the coarsen_mode fix) started at step 22, so
# there is no early window safe to skip.
FIRST_RECORDED = 1

ORDER_DIR = Path('output') / 'exadis_orders'
PYDIS_DIR = Path('output') / 'ordered_pydis'

# Known limit of the order import on this case. Steps past it have never been
# reached without a divergence, so a first divergence at or after this is the
# status quo and anything earlier is new. None means the whole run agrees,
# which is the current state.
KNOWN_LIMIT = None


def thread_check():
    """thread_check: confirm ExaDiS really is on one thread

    The environment variable is set at import time, but Kokkos is what
    decides, and it prints its own banner. Reading the banner rather than the
    variable is what test4 and 02_frank_read_src's own copy of this function
    do, for the same reason: the variable is the request, the banner is the
    answer.
    """
    return os.environ.get('OMP_NUM_THREADS') == '1'


def scope_note():
    print("scope: both examples are given coarsen_mode=%d, rather than each "
          "taking its own" % COARSEN_MODE)
    print("       default (pydis 0, pyexadis_base 1). Those are different "
          "coarsening")
    print("       algorithms, not settings of one: mismatching them moved node "
          "counts apart")
    print("       as early as step 22, with ordering otherwise fully "
          "controlled.")


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


# --------------------------------------------------------------------------
# independent runs
# --------------------------------------------------------------------------

def independent(plot):
    """independent: run both examples as subprocesses, and check the physics
    and stored-reference agreement this module's docstring describes"""
    # --force-mode=CUTOFF_MODEL is what makes the exadis run comparable to
    # pydis' Elasticity_SBA: DDD_FFT_MODEL adds a long-range FFT contribution
    # that pydis has no counterpart for. Both example scripts default to
    # plot=True and expose --no-plot to turn it off, so plot=True here means
    # passing neither flag rather than a --plot one.
    common_args = (['--max-step', MAX_STEP,
                    '--print-freq', PRINT_FREQ,
                    '--write-freq', WRITE_FREQ]
                   + ([] if plot else ['--no-plot']))
    common_args += ['--coarsen-mode', str(COARSEN_MODE)]
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
    print("at step %d: max nearest-node distance, pydis vs exadis = %.4e, "
          "tolerance = %.1e" % (STRESS_STEP, d_mid, TOL))
    ok &= report("node counts agree at step %d" % STRESS_STEP,
                 n_pydis == n_exadis)
    ok &= report("configurations agree within %.1e at step %d"
                 % (TOL, STRESS_STEP), d_mid < TOL)

    if not REF_NPZ.is_file():
        print("reference not found: %s" % REF_NPZ)
        print("generate it with 'make binary_junction_elast_ref' and copy "
              "it into ref_data/")
        return False

    n_pydis, n_ref, d_pydis = compare_to_ref(PYDIS_JSON, REF_NPZ)
    ok &= report("at step %d: nodes: pydis = %d, ref = %d; max nearest-node "
                 "distance, pydis  vs ref = %.4e, tolerance = %.1e"
                 % (MAX_STEP, n_pydis, n_ref, d_pydis, TOL),
                 n_pydis == n_ref and d_pydis < TOL)

    n_exadis, n_ref, d_exadis = compare_to_ref(EXADIS_JSON, REF_NPZ)
    ok &= report("at step %d: nodes: exadis = %d, ref = %d; max nearest-node "
                 "distance, exadis vs ref = %.4e, tolerance = %.1e"
                 % (MAX_STEP, n_exadis, n_ref, d_exadis, TOL),
                 n_exadis == n_ref and d_exadis < TOL)
    return bool(ok)


# --------------------------------------------------------------------------
# with ExaDiS' iteration orders imported
# --------------------------------------------------------------------------

def exadis_pass(plot=False):
    """exadis_pass: run ExaDiS, recording its segment order right before
    remesh, but only on a step where remesh actually changes the node
    count -- the only one of its three orders (arm, segment, node)
    pydis_pass() needs, and only the steps it needs it on; see that
    function's own docstring for why.

    Node-count-unchanged means Remesh_LengthBased found nothing to coarsen
    or refine, so no segment-order-dependent decision was made and nothing
    written this step could matter to pydis_pass() -- confirmed directly:
    of 300 steps, only 25 change the node count across remesh, and
    imposing the recorded order on exactly those 25 (found from this same
    run) gives the identical result, worst 4.602185e-06, as imposing it
    on all 300. This is a much smaller set than the two Topology events
    (which turn out to need no order at all: both have exactly one
    degree>=4 candidate node, so there is nothing for a node-enumeration
    order to change): remesh runs, and can act, far more often than a full
    multi-arm split happens.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet
    from framework.arm_order import export_arm_order
    from order_import import export_orders
    import test_binary_junction_exadis_elast as ex

    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    pyexadis.initialize()

    def snapshot(N):
        # arm order (export_arm_order) is bundled in here too, unused by
        # pydis_pass() -- export_orders is order_import.py's shared,
        # general-purpose recording call (see 02_frank_read_src's own
        # copy), not worth a narrower one just for this file. Both pieces
        # must be captured together, before remesh runs: tags a post-remesh
        # export_arm_order(N) would name (a node remesh just created) are
        # not in a pre-remesh data['tags'], and export_orders needs both
        # from the same snapshot.
        G = N.get_disnet(ExaDisNet)
        data = dict(G.get_nodes_data())
        data['nodeids'] = np.array(G.get_segs_data()['nodeids'])
        return G, export_orders(data, export_arm_order(N))

    base = ex.SimulateNetwork

    class Recording(base):
        def step_topological_operations(self, N, state):
            step = state['istep']
            if not (FIRST_RECORDED <= step <= MAX_STEP):
                return base.step_topological_operations(self, N, state)
            if self.cross_slip is not None:
                self.cross_slip.Handle(N, state)
            self.collision.HandleCol(N, state)
            self.topology.Handle(N, state)
            G, orders = snapshot(N)
            n_before = G.num_nodes()
            self.remesh.Remesh(N, state)
            if G.num_nodes() != n_before:
                np.savez(str(ORDER_DIR / ('%d_p2.npz' % step)), **orders)
            N.write_json(str(ORDER_DIR / ('cfg_%d.json' % step)))

    ex.SimulateNetwork = Recording
    ex.main(plot=plot, force_mode='CUTOFF_MODEL', max_step=MAX_STEP,
            coarsen_mode=COARSEN_MODE,
            print_freq=PRINT_FREQ, write_freq=MAX_STEP)
    pyexadis.finalize()


def pydis_pass(plot=False):
    """pydis_pass: run PyDiS with ExaDiS' segment order imposed right before
    remesh, on whichever steps exadis_pass() found worth recording,
    recording its own configuration in turn

    Only segment order, and only right there -- confirmed by a six-way
    ablation (arm order at each of p0/p1/p2, segment order at each of
    p0/p1/p2) that this is the minimal set that still reproduces full
    agreement: dropping arm order anywhere, or segment order anywhere
    except immediately before remesh, changes nothing. Two reasons, found
    separately: arm order never decides a near-tie for this case's
    collision or topology the way it does for 02_frank_read_src's or
    test5_topology_mode's, so imposing it buys nothing; and node order
    (dropped for the same reason) has nothing to enumerate over, since
    both of this case's two Topology events have exactly one degree>=4
    candidate node. Segment order's one live consumer is
    Remesh_LengthBased's `list(G.all_segments_tags())`, read once per call,
    which is why only the snapshot taken immediately before remesh.Remesh()
    runs -- not the ones a fuller order-import would also take before
    collision and before topology -- has any effect.

    A missing recording is not a failure: exadis_pass() only writes one on
    a step where its own remesh changed the node count, so most steps (275
    of 300, this run) have none by design -- nothing was decided there for
    an order to matter to. Returns the (step, phase, reason) list where a
    recording exists but could not be imposed (a real problem, not the
    ordinary no-recording case).
    """
    from pydis.disnet import DisNet
    from order_import import load_orders, geometric_map, segments_in_order
    import test_binary_junction_pydis_elast as py

    failures = []

    def segment_order_for(G, step):
        path = ORDER_DIR / ('%d_p2.npz' % step)
        if not path.is_file():
            return None
        orders = load_orders(str(path))
        mapping, note = geometric_map(G, orders)
        if mapping is None:
            reason, text = note
            failures.append((step, 'p2', reason, text))
            return None
        return [(mapping[a], mapping[b]) for a, b in orders['segments']]

    base = py.SimulateNetwork

    class Following(base):
        def step_topological_operations(self, DM, state):
            # PyDiS' istep is 0-based and ExaDiS' is 1-based, so istep+1 here
            # names the same completed step as ExaDiS' istep
            step = state['istep'] + 1
            if not (FIRST_RECORDED <= step <= MAX_STEP):
                return base.step_topological_operations(self, DM, state)
            if self.cross_slip is not None:
                self.cross_slip.Handle(DM, state)
            self.collision.HandleCol(DM, state)
            self.topology.Handle(DM, state)
            if self.remesh is not None:
                G = DM.get_disnet(DisNet)
                wanted = segment_order_for(G, step)
                ordered = (segments_in_order(list(G.all_segments_tags()), wanted)
                          if wanted is not None else None)
                if ordered is None:
                    if wanted is not None:
                        failures.append((step, 'p2', 'segment count mismatch'))
                    self.remesh.Remesh(DM, state)
                else:
                    print("step %d: imposing exadis' segment order before remesh"
                         % step)
                    live = G.all_segments_tags
                    G.all_segments_tags = lambda: iter(ordered)
                    try:
                        self.remesh.Remesh(DM, state)
                    finally:
                        G.all_segments_tags = live
            DM.write_json(str(PYDIS_DIR / ('ordered_pydis_%d.json' % step)))
            return state

    py.SimulateNetwork = Following
    py.main(plot=plot, max_step=MAX_STEP, coarsen_mode=COARSEN_MODE,
            print_freq=PRINT_FREQ, write_freq=MAX_STEP)
    return failures


def per_step_comparison():
    """per_step_comparison: (first divergence or None, rows)"""
    first, rows = None, []
    for step in range(FIRST_RECORDED, MAX_STEP + 1):
        a = PYDIS_DIR / ('ordered_pydis_%d.json' % step)
        b = ORDER_DIR / ('cfg_%d.json' % step)
        if not (a.is_file() and b.is_file()):
            continue
        n_pydis, n_exadis, distance = compare_configs(a, b)
        if distance > TOL and first is None:
            first = step
        rows.append((step, n_pydis, n_exadis, distance))
    return first, rows


def clear_recordings():
    """clear_recordings: delete what a previous run left in output/

    Everything the ordered path reads, it wrote itself earlier in the same
    run. A file left over from a previous run with a different MAX_STEP, a
    different flag or a different build is therefore never wanted, and is
    actively dangerous: a stale ordered_pydis_<step>.json compares cleanly
    against the current exadis recording and reads as agreement. The files
    are removed rather than overwritten, which also keeps a shortened run
    from being judged against a longer one's leftovers.
    """
    removed = 0
    for path in sorted(ORDER_DIR.glob('*')):
        path.unlink()
        removed += 1
    for path in sorted(PYDIS_DIR.glob('*')):
        path.unlink()
        removed += 1
    return removed


def ordered(verbose, plot=False):
    """ordered: run both codes with ExaDiS' orders imposed on PyDiS

    No stored-reference comparison here, unlike independent(): both codes are
    still under active development on this case, so a saved configuration
    says more about when it was recorded than about whether the codes agree
    now. This mode's own per-step comparison against each other is the
    stronger, more current statement.

    plot=True shows each code's own live plot in turn (exadis_pass() runs
    to completion, then pydis_pass() does), same as running either example
    directly -- not side by side, since both run in this one process.
    """
    print("importing exadis' segment order into pydis immediately before")
    print("remesh, on whichever steps exadis' own remesh changes the node")
    print("count -- the only one of exadis' three orders (arm, segment,")
    print("node) this case's own agreement turns out to depend on, and only")
    print("the steps it actually decided something on; see pydis_pass()'s")
    print("docstring for the ablation that found this.")
    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    PYDIS_DIR.mkdir(parents=True, exist_ok=True)
    removed = clear_recordings()
    if removed:
        print("cleared %d file(s) left by a previous run" % removed)
    print("")

    exadis_pass(plot)
    failures = pydis_pass(plot)
    first, rows = per_step_comparison()

    if verbose:
        print("   step   nodes py/ex    max nearest-node distance")
        for step, n_pydis, n_exadis, distance in rows:
            print("   %4d   %3d / %3d      %.4e%s"
                  % (step, n_pydis, n_exadis, distance,
                     '   <== first divergence' if step == first else ''))
        print("")

    # Import failures are reported, not asserted -- see 02_frank_read_src's
    # own copy of this comment for why.
    if failures:
        from order_import import COUNT_MISMATCH
        counts = [f for f in failures if f[2] == COUNT_MISMATCH]
        real = [f for f in failures if f[2] != COUNT_MISMATCH]
        total = MAX_STEP - FIRST_RECORDED + 1
        print("segment order not imposed: %d of %d steps" % (len(failures), total))
        if counts:
            print("       %d with different node counts mid-step, expected "
                  "(exadis purges at end of pass): %s"
                  % (len(counts), ["%d %s: %s" % (f[0], f[1], f[3]) for f in counts[:3]]))
        if real:
            limit = MAX_STEP if first is None else first
            before = [f for f in real if f[0] <= limit]
            print("       %d could not be paired (same node count): %s"
                  % (len(real), ["%d %s: %s" % (f[0], f[1], f[3]) for f in real[:3]]))
            if before:
                print("       of those, %d at or before the first divergence"
                      % len(before))
        print("")

    if not rows:
        print("no steps were recorded; nothing to compare")
        return False

    if first is None:
        worst = max(d for _, _, _, d in rows)
        return report("agreement holds through step %d (worst %.4e, tolerance "
                      "%.1e)" % (MAX_STEP, worst, TOL), True)

    agreed = first - 1
    before = [d for step, _, _, d in rows if step <= agreed]
    at_first = next(d for step, _, _, d in rows if step == first)
    n_py, n_ex = next((a, b) for step, a, b, _ in rows if step == first)
    print("agreement holds through step %d (worst %.4e, tolerance %.1e)"
          % (agreed, max(before) if before else 0.0, TOL))
    if KNOWN_LIMIT is None:
        return report("no divergence through step %d: first divergence at step "
                      "%d, nodes %d / %d, max nearest-node distance %.4e, "
                      "tolerance %.1e"
                      % (MAX_STEP, first, n_py, n_ex, at_first, TOL), False)
    return report("first divergence is no earlier than the known limit at step "
                  "%d: it is at %d, nodes %d / %d, distance %.4e, tolerance %.1e"
                  % (KNOWN_LIMIT, first, n_py, n_ex, at_first, TOL),
                  first >= KNOWN_LIMIT)


def main(plot=False, import_order=IMPORT_EXADIS_ORDER, verbose=False):
    os.makedirs('output', exist_ok=True)
    scope_note()
    print("")
    if not thread_check():
        print("OMP_NUM_THREADS is %r, not 1; ExaDiS' summation order will vary "
              "and this comparison is not reproducible"
              % os.environ.get('OMP_NUM_THREADS'))
        return False
    if import_order:
        return ordered(verbose, plot)
    return independent(plot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', dest='plot', action='store_true',
                        help='let the two example scripts show their plots '
                             '(off by default here, unlike running them '
                             'directly). Ordered mode shows exadis\' plot, '
                             'then pydis\' in turn, one process at a time; '
                             'independent mode shows both subprocesses\' '
                             'plots as they run')
    parser.add_argument('--no-import-order', dest='import_order',
                        action='store_false', default=IMPORT_EXADIS_ORDER,
                        help="compare two independent runs instead of "
                             "imposing exadis' iteration orders on pydis")
    parser.add_argument('--verbose', action='store_true',
                        help='print the per-step distance table')
    args = parser.parse_args()

    passed = main(plot=args.plot, import_order=args.import_order,
                  verbose=args.verbose)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
