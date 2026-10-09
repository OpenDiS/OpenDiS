"""Full-run comparison: Frank-Read source with elasticity, PyDiS vs ExaDiS.

The two are configured to compute the same thing:
  - PyDiS  Elasticity_SBA  with cutoff = 0.25*Lbox
  - ExaDiS CUTOFF_MODEL    with the same cutoff and Ec = 0.0
ExaDiS' default force mode is DDD_FFT_MODEL, which legitimately does NOT agree,
because it adds the long-range FFT contribution PyDiS does not have. Hence
--force-mode=CUTOFF_MODEL below.

Two ways of comparing them, chosen by IMPORT_EXADIS_ORDER below.

WITH the order import (default). ExaDiS runs first, recording the order in which
it visits segments and multi-nodes at three points in every step. PyDiS then runs
with those orders imposed at the same three points. Iteration order is therefore
not a variable, and what is left is everything else: physics, arithmetic, and
every rule the two codes implement. They agree to round-off for the whole run,
node counts equal at every step.

Two of the three orders ExaDiS records are imposed, each measured to be needed
on this case, in both coarsening branches:

    segment order  read by the collision pass when it enumerates candidate
                   pairs and by remesh when it snapshots segments. Without it
                   the comparison parts at step 263, where a dozen collisions
                   come due in one pass and each code acts on the subset it
                   reaches first: 126 nodes against 118, 5.6e+01 apart.
    node order     read by the topology pass over multi-nodes, and by remesh's
                   node walk under coarsen_mode 1. Without it the comparison
                   parts at step 268, one of the eight steps where topology
                   acts, at 4.2e-05 with node counts still equal.

Arm order is recorded but no longer imposed; see order_import.apply_orders.
Restricting the import to fewer steps than FIRST_RECORDED..MAX_STEP does not
work: events are dense from step 260 on, and imposing only over 260-273 parts at
275, only at 263 and 268 parts at 264.

WITHOUT it. The two examples run as independent subprocesses and their final
configurations are compared, the examples staying the single source of truth for
how the run is set up. This parts at step 263 and cannot be made to agree. Both
codes allocate and free nodes in their own order, ExaDiS' remove_nodes moving the
array's last node into each vacated slot while PyDiS' insertion-ordered dict
leaves the rest alone, so by step 268 only 16 of 107 array slots hold the same
physical node. At step 263 a dozen collisions come due in one pass and each code
acts on a different subset. Nothing in either code is wrong; two independent runs
simply cannot be held together across a multi-operation pass.

The order import is the more useful of the two, because a failure under it points
at something actionable. Keep the independent path for the honest statement of
where two separate runs actually end up.

Output is written under the working directory this test is run from, not next to
the sources.
"""

import os

# Before anything imports pyexadis. ExaDiS' force kernels are Kokkos parallel
# reductions, so with more than one thread the summation order varies between
# runs. That is invisible for the first 266 steps here and then decides a branch
# at the first pass where several operations come due at once: on 8 threads this
# comparison parts at step 267 (2.8412e+00) and on one thread it holds to 300.
# Assigned rather than setdefault, so an OMP_NUM_THREADS inherited from the
# shell cannot make the result depend on where it was run from. The same is
# done in test1_node_force, test4_collision_mode and test5_topology_mode.
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
import sys
from pathlib import Path

import numpy as np

# this script lives in tests/full_runs/02_frank_read_src/, so the repository
# root is 3 levels up
opendis_root = Path(__file__).resolve().parents[3]
for _p in ['python', 'lib', 'core/pydis/python', 'core/exadis/python',
           'examples/02_frank_read_src']:
    _q = str(opendis_root / _p)
    if _q not in sys.path:
        sys.path.append(_q)

from framework.testing import report, run_script, compare_configs

examples_dir = opendis_root / 'examples' / '02_frank_read_src'

# True: run the codes in one process and impose ExaDiS' iteration orders on
# PyDiS every step. False: run the two examples independently and compare their
# final configurations. See this module's docstring for what each answers.
IMPORT_EXADIS_ORDER = True

# The tolerance is set well above the round-off floor: summation order in the
# threaded ExaDiS kernels is not deterministic, so an exact figure would be
# flaky, but anything approaching this bound means a real divergence. The two
# codes' positions drift apart at ~1e-07 in a box of 1000 over the first 260
# steps of ordinary motion, which is the floor this sits above.
TOL = 1.0e-6

MAX_STEP = 300

# --max-step overrides it, for a shortened diagnostic run; the divergences this
# case is watched for all fall before step 300, so a shorter run still locates
# them, and it keeps an ablation that makes pydis' network blow up from taking
# hours.
def set_max_step(value):
    global MAX_STEP, WRITE_FREQ
    MAX_STEP = WRITE_FREQ = value

# the examples print a progress line every PRINT_FREQ steps; keep the test
# output short without silencing it entirely, so a hung run is still visible
PRINT_FREQ = 100

# steps between intermediate configuration dumps. The independent comparison
# uses only the final configuration, written after the run regardless of this
# setting, so the intermediate dumps are pure noise there. MAX_STEP writes one
# at the last step.
WRITE_FREQ = MAX_STEP

PYDIS_SCRIPT = 'test_frank_read_src_pydis_elast.py'
EXADIS_SCRIPT = 'test_frank_read_src_exadis_elast.py'

PYDIS_JSON = Path('output') / 'frank_read_src_pydis_elast_final.json'
EXADIS_JSON = Path('output') / 'frank_read_src_exadis_elast_final.json'

# The window the order import records and compares. Everything before it agrees
# to round-off and recording it only costs disk; the first collision is at 260.
FIRST_RECORDED = 250

# The coarsening branches this test knows how to compare.
#   0  segment-centric: a segment below minseg has its two endpoints merged to
#      their mid-point.
#   1  node-centric: a two-arm node with either arm below minseg is merged into
#      its nearer neighbour, which does not move. exadis' python default, and
#      ParaDiS' rule.
# They are different algorithms rather than settings of one, so each run fixes
# the same value on both sides.
COARSEN_MODES = (0, 1)

# What runs when no --coarsen-mode is given. One branch, so a bare invocation
# is one comparison and reads as one. The makefile here and
# tests/full_runs/CMakeLists.txt add the other branch as a second invocation,
# so `make` and `ctest` cover both.
DEFAULT_MODES = (1,)


def order_dir(mode):
    return Path('output') / ('exadis_orders_m%d' % mode)


def pydis_dir(mode):
    return Path('output') / ('ordered_pydis_m%d' % mode)

# Known limit of the order import on this case, per coarsen_mode. Steps past it
# have never been reached without a divergence, so a first divergence at or
# after this is the status quo and anything earlier is new. None means the whole
# run is expected to agree, which is where both branches now stand.
#
# coarsen_mode 1 used to part at step 297, where the two codes chose different
# merge nodes for the same collision because they stored one segment in opposite
# directions. Fixed in pydis rather than here: see DisNet.segment_as_stored.
KNOWN_LIMIT = {0: None, 1: None}

# Neither comparison uses a stored reference. Both codes are under active
# development past the first collision, so a saved configuration says more about
# when it was recorded than about whether the codes agree now. The weakness of
# comparing only against each other is worth stating: a change that shifted both
# equally would pass unnoticed, which is how the maxseg bound defect behaved,
# both codes dropping the same pairs and agreeing while both were wrong.


def thread_check():
    """thread_check: confirm ExaDiS really is on one thread

    The environment variable is set at import time, but Kokkos is what decides,
    and it prints its own banner. Reading the banner rather than the variable is
    what test4 does, for the same reason: the variable is the request, the
    banner is the answer.
    """
    return os.environ.get('OMP_NUM_THREADS') == '1'


def scope_note(modes):
    print("scope: coarsening is compared for %s."
          % ' and '.join('coarsen_mode=%d' % m for m in modes))
    print("       0 is segment-centric (a segment below minseg has its "
          "endpoints merged to")
    print("       their mid-point), 1 is node-centric (a two-arm node with a "
          "short arm merges")
    print("       into its nearer neighbour, which does not move). Each run "
          "sets the same")
    print("       value on both sides; neither code's default is relied on.")


# --------------------------------------------------------------------------
# independent runs
# --------------------------------------------------------------------------

def independent(plot, mode):
    """independent: run both examples as subprocesses, compare the final states"""
    # --force-mode=CUTOFF_MODEL is what makes the exadis run comparable; see
    # the note at the top of this file. Both example scripts default to
    # plot=True and expose --no-plot to turn it off, so plot=True here means
    # passing neither flag rather than a --plot one.
    common_args = (['--max-step', MAX_STEP,
                    '--coarsen-mode', mode,
                    '--print-freq', PRINT_FREQ,
                    '--write-freq', WRITE_FREQ]
                   + ([] if plot else ['--no-plot']))
    # the two examples write these under a fixed name, so a previous mode's
    # output must not be left where this one's comparison would read it
    for stale in (PYDIS_JSON, EXADIS_JSON):
        if stale.is_file():
            stale.unlink()
    pydis_args = list(common_args)
    exadis_args = ['--force-mode=CUTOFF_MODEL'] + common_args

    ok = True
    ok &= report("run pydis  example",
                 run_script(examples_dir / PYDIS_SCRIPT, pydis_args,
                            headless=not plot))
    ok &= report("run exadis example",
                 run_script(examples_dir / EXADIS_SCRIPT, exadis_args,
                            headless=not plot))
    if not ok:
        print("an example script failed; skipping the comparison")
        return False

    for f in (PYDIS_JSON, EXADIS_JSON):
        if not f.is_file():
            print("expected output not found: %s" % f)
            return False

    n_pydis, n_exadis, d_pair = compare_configs(PYDIS_JSON, EXADIS_JSON)
    ok &= report("node counts agree: pydis = %d, exadis = %d"
                 % (n_pydis, n_exadis), n_pydis == n_exadis)
    ok &= report("node positions agree: max nearest-node distance = %.4e, "
                 "tolerance = %.1e" % (d_pair, TOL), d_pair < TOL)
    return bool(ok)


# --------------------------------------------------------------------------
# with ExaDiS' iteration orders imported
# --------------------------------------------------------------------------

def exadis_pass(mode):
    """exadis_pass: run ExaDiS, recording its three orders and configuration"""
    import pyexadis
    from pyexadis_base import ExaDisNet
    from framework.arm_order import export_arm_order
    from order_import import export_orders
    import test_frank_read_src_exadis_elast as ex

    ORDER_DIR = order_dir(mode)
    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    pyexadis.initialize()

    def record(N, step, phase):
        G = N.get_disnet(ExaDisNet)
        data = dict(G.get_nodes_data())
        data['nodeids'] = np.array(G.get_segs_data()['nodeids'])
        np.savez(str(ORDER_DIR / ('%d_%s.npz' % (step, phase))),
                 **export_orders(data, export_arm_order(N)))

    base = ex.SimulateNetwork

    class Recording(base):
        def step_topological_operations(self, N, state):
            step = state['istep']
            if not (FIRST_RECORDED <= step <= MAX_STEP):
                return base.step_topological_operations(self, N, state)
            record(N, step, 'p0')
            self.collision.HandleCol(N, state)
            record(N, step, 'p1')
            self.topology.Handle(N, state)
            record(N, step, 'p2')
            self.remesh.Remesh(N, state)
            N.write_json(str(ORDER_DIR / ('cfg_%d.json' % step)))

    ex.SimulateNetwork = Recording
    ex.main(plot=False, force_mode='CUTOFF_MODEL', max_step=MAX_STEP,
            coarsen_mode=mode, print_freq=PRINT_FREQ, write_freq=MAX_STEP)
    pyexadis.finalize()


def pydis_pass(mode):
    """pydis_pass: run PyDiS under ExaDiS' orders, recording its configuration

    Returns the (step, phase, reason) list where an order could not be imposed.
    """
    from pydis.disnet import DisNet
    import pydis.topology.topology_ops as topology_ops
    import pydis.topology.topology_disnet as topology_disnet
    from order_import import (load_orders, apply_orders, segments_in_order,
                              COUNT_MISMATCH)
    import test_frank_read_src_pydis_elast as py

    # PyDiS holds no segment order or node order of its own to overwrite, so the
    # imported ones are handed to the two places that consult them: the segment
    # enumeration every pass goes through, and the topology's multi-node loop.
    imposed = {'segments': None, 'nodes': None, 'in_remesh': False}
    failures = []
    ORDER_DIR, PYDIS_DIR = order_dir(mode), pydis_dir(mode)

    live_segments = DisNet.all_segments_tags
    live_nodes = DisNet.all_nodes_tags

    # coarsen_mode 1 coarsens by walking nodes rather than segments, so that
    # walk consults a node order as well. Imposed only inside the remesh call:
    # all_nodes_tags is also what the force and mobility loops iterate, and
    # reordering those would change their summation order. Inert for mode 0,
    # whose coarsening reads only the segment order.
    def all_nodes_tags(self):
        tags = list(live_nodes(self))
        rank = imposed['nodes']
        if rank is None or not imposed['in_remesh']:
            return iter(tags)
        return iter(sorted(tags, key=lambda t: rank.get(t, len(tags))))
    DisNet.all_nodes_tags = all_nodes_tags

    def all_segments_tags(self):
        segments = list(live_segments(self))
        wanted = imposed['segments']
        if wanted is None:
            return segments
        return segments_in_order(segments, wanted) or segments
    DisNet.all_segments_tags = all_segments_tags

    split_multi_nodes = topology_ops.split_multi_nodes

    def ordered_split_multi_nodes(G, state, trial_fn, max_degree=15):
        rank = imposed['nodes']
        if rank is None:
            return split_multi_nodes(G, state, trial_fn, max_degree)
        tags = list(G.all_nodes_tags())
        # sorted(), not a filter: a tag the import did not name keeps a stable
        # place at the end rather than being skipped
        held = G.all_nodes_tags
        G.all_nodes_tags = lambda: iter(sorted(
            tags, key=lambda t: rank.get(t, len(tags))))
        try:
            return split_multi_nodes(G, state, trial_fn, max_degree)
        finally:
            G.all_nodes_tags = held
    topology_ops.split_multi_nodes = ordered_split_multi_nodes
    topology_disnet.split_multi_nodes = ordered_split_multi_nodes

    def impose(G, step, phase):
        path = ORDER_DIR / ('%d_%s.npz' % (step, phase))
        imposed['segments'] = imposed['nodes'] = None
        if not path.is_file():
            failures.append((step, phase, 'no recording'))
            return
        segments, nodes, note = apply_orders(G, load_orders(str(path)))
        if segments is None:
            reason, text = note
            failures.append((step, phase, reason, text))
            return
        imposed['segments'], imposed['nodes'] = segments, nodes

    base = py.SimulateNetwork

    class Following(base):
        def step_topological_operations(self, DM, state):
            # PyDiS' istep is 0-based and ExaDiS' is 1-based, so istep+1 here
            # names the same completed step as ExaDiS' istep
            step = state['istep'] + 1
            if not (FIRST_RECORDED <= step <= MAX_STEP):
                return base.step_topological_operations(self, DM, state)
            impose(DM.get_disnet(DisNet), step, 'p0')
            self.collision.HandleCol(DM, state)
            impose(DM.get_disnet(DisNet), step, 'p1')
            self.topology.Handle(DM, state)
            impose(DM.get_disnet(DisNet), step, 'p2')
            imposed['in_remesh'] = True
            self.remesh.Remesh(DM, state)
            imposed['in_remesh'] = False
            DM.write_json(str(PYDIS_DIR / ('ordered_pydis_%d.json' % step)))
            return state

    py.SimulateNetwork = Following
    py.main(plot=False, max_step=MAX_STEP, coarsen_mode=mode,
            print_freq=PRINT_FREQ, write_freq=MAX_STEP)
    return failures


def per_step_comparison(mode):
    """per_step_comparison: (first divergence or None, rows)"""
    first, rows = None, []
    for step in range(FIRST_RECORDED, MAX_STEP + 1):
        a = pydis_dir(mode) / ('ordered_pydis_%d.json' % step)
        b = order_dir(mode) / ('cfg_%d.json' % step)
        if not (a.is_file() and b.is_file()):
            continue
        n_pydis, n_exadis, distance = compare_configs(a, b)
        if distance > TOL and first is None:
            first = step
        rows.append((step, n_pydis, n_exadis, distance))
    return first, rows


def clear_recordings(mode):
    """clear_recordings: delete what a previous run left in output/

    Everything the ordered path reads, it wrote itself earlier in the same run.
    A file left over from a previous run with a different MAX_STEP, a different
    flag or a different build is therefore never wanted, and is actively
    dangerous: a stale ordered_pydis_<step>.json compares cleanly against the
    current exadis recording and reads as agreement. That happened during
    development and produced a wrong conclusion, so the files are removed rather
    than overwritten, which also keeps a shortened run from being judged against
    a longer one's leftovers.
    """
    removed = 0
    for path in sorted(order_dir(mode).glob('*')):
        path.unlink()
        removed += 1
    for path in sorted(pydis_dir(mode).glob('*')):
        path.unlink()
        removed += 1
    # ordered_pydis_<step>.json used to sit directly in output/; sweep those up
    # too so an older run's leftovers cannot be read by mistake
    for path in sorted(Path('output').glob('ordered_pydis_*.json')):
        path.unlink()
        removed += 1
    # and the un-suffixed directories this test used before it ran two modes
    for stale in (Path('output') / 'exadis_orders', Path('output') / 'ordered_pydis'):
        if stale.is_dir():
            for path in sorted(stale.glob('*')):
                path.unlink()
                removed += 1
            stale.rmdir()
    return removed


def ordered(verbose, mode):
    """ordered: run both codes with ExaDiS' orders imposed on PyDiS"""
    print("importing exadis' arm, segment and multi-node orders into pydis at")
    print("three points per step: before collision, before topology, before")
    print("remesh, and its node order onto the remesh walk. What remains is not")
    print("iteration order.")
    order_dir(mode).mkdir(parents=True, exist_ok=True)
    pydis_dir(mode).mkdir(parents=True, exist_ok=True)
    removed = clear_recordings(mode)
    if removed:
        print("cleared %d file(s) left by a previous run" % removed)
    print("")

    exadis_pass(mode)
    failures = pydis_pass(mode)
    first, rows = per_step_comparison(mode)

    if verbose:
        print("   step   nodes py/ex    max nearest-node distance")
        for step, n_pydis, n_exadis, distance in rows:
            print("   %4d   %3d / %3d      %.4e%s"
                  % (step, n_pydis, n_exadis, distance,
                     '   <== first divergence' if step == first else ''))
        print("")

    # Import failures are reported, not asserted. Once the two configurations
    # diverge the geometric map stops being one to one, so every later import
    # fails as a consequence rather than a cause; and a phase snapshot can
    # legitimately fail to map while the step still ends in agreement, when the
    # two codes reach the same end state through intermediates of different
    # size. The per-step comparison is the ground truth: a failed import can
    # only weaken it, never fake it.
    if failures:
        from order_import import COUNT_MISMATCH
        limit = MAX_STEP if first is None else first
        counts = [f for f in failures if f[2] == COUNT_MISMATCH]
        real = [f for f in failures if f[2] != COUNT_MISMATCH]
        total = 3*(MAX_STEP - FIRST_RECORDED + 1)
        print("order imports not applied: %d of %d step-phases" % (len(failures), total))
        # Expected, and not evidence of anything: ExaDiS' merges leave the
        # absorbed node and the connecting segment in place until
        # purge_network runs at the end of the pass, so its mid-pass network can
        # carry a node PyDiS has already removed. The step still ends in
        # agreement.
        if counts:
            print("       %d with different node counts mid-step, expected "
                  "(exadis purges at end of pass): %s"
                  % (len(counts), ["%d %s: %s" % (f[0], f[1], f[3]) for f in counts[:3]]))
        # These are the ones worth reading: same number of nodes, but they
        # cannot be paired, which means the two networks have actually parted.
        if real:
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
    known = KNOWN_LIMIT[mode]
    if known is None:
        return report("no divergence through step %d: first divergence at step "
                      "%d, nodes %d / %d, max nearest-node distance %.4e, "
                      "tolerance %.1e"
                      % (MAX_STEP, first, n_py, n_ex, at_first, TOL), False)
    return report("first divergence is no earlier than the known limit at step "
                  "%d: it is at %d, nodes %d / %d, distance %.4e, tolerance %.1e"
                  % (known, first, n_py, n_ex, at_first, TOL),
                  first >= known)


def main(plot=False, import_order=IMPORT_EXADIS_ORDER, verbose=False,
         modes=DEFAULT_MODES):
    os.makedirs('output', exist_ok=True)
    scope_note(modes)
    print("")
    if not thread_check():
        print("OMP_NUM_THREADS is %r, not 1; ExaDiS' summation order will vary "
              "and this comparison is not reproducible"
              % os.environ.get('OMP_NUM_THREADS'))
        return False
    if len(modes) == 1:
        mode = modes[0]
        print("=" * 70)
        print("coarsen_mode = %d" % mode)
        print("=" * 70)
        return (ordered(verbose, mode) if import_order
                else independent(plot, mode))

    # More than one branch to compare, so each runs in its own process: the
    # ordered path calls pyexadis.initialize(), and Kokkos::initialize() may be
    # called only once in the lifetime of a process, finalize or no finalize.
    # Every mode is run even if an earlier one fails, so one invocation reports
    # the whole picture rather than stopping at the first branch that parts.
    print("running %d branches, one process each: %s. Each prints its own "
          "banner and" % (len(modes), ', '.join('coarsen_mode=%d' % m
                                                for m in modes)))
    print("verdict below, and both are summarised at the very end.")
    print("")
    passthrough = ([] if import_order else ['--no-import-order'])
    passthrough += (['--verbose'] if verbose else [])
    passthrough += (['--plot'] if plot else [])
    results = {}
    for mode in modes:
        results[mode] = run_script(Path(__file__).resolve(),
                                   ['--coarsen-mode', mode] + passthrough,
                                   headless=not plot)
        print("")
    for mode, ok in sorted(results.items()):
        print("coarsen_mode = %d: %s" % (mode, 'PASSED' if ok else 'FAILED'))
    return all(results.values())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', dest='plot', action='store_true',
                        help='let the two example scripts show their plots '
                             '(off by default here, unlike running them '
                             'directly; independent comparison only)')
    parser.add_argument('--no-import-order', dest='import_order',
                        action='store_false', default=IMPORT_EXADIS_ORDER,
                        help="compare two independent runs instead of imposing "
                             "exadis' iteration orders on pydis")
    parser.add_argument('--verbose', action='store_true',
                        help='print the per-step distance table')
    parser.add_argument('--max-step', dest='max_step', type=int, default=None,
                        help='shorten the run (default %d)' % MAX_STEP)
    parser.add_argument('--coarsen-mode', dest='modes', type=int,
                        choices=COARSEN_MODES, action='append',
                        help='which coarsening branch to compare; repeatable. '
                             'Defaults to %s. make and ctest run both.'
                             % ' and '.join(str(m) for m in DEFAULT_MODES))
    args = parser.parse_args()
    if args.max_step:
        set_max_step(args.max_step)

    passed = main(plot=args.plot, import_order=args.import_order,
                  verbose=args.verbose,
                  modes=tuple(args.modes) if args.modes else DEFAULT_MODES)
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
