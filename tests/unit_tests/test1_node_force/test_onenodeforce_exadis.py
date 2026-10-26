"""OneNodeForce against NodeForce, ExaDiS.

The ExaDiS counterpart of test_onenodeforce_pydis.py, pinning the same
invariant on the other code:

    OneNodeForce(tag) == the NodeForce() entry for tag, for every node

CalForce offers two routes to a nodal force. NodeForce evaluates the whole
network and returns a force per node; OneNodeForce evaluates one node, through
compute_node_force, which is what TopologySerial calls when it scores a trial
split. Two routes to one quantity is the arrangement that drifts, and here it
did: the serial node_force queried the neighbor bins once at the node position,
reused that list for every arm, and summed every returned pair with no distance
test, while the full path filters on get_min_dist2_segseg < cutoff2.

MATCH_GLOBAL

This test passes match_global=True, which applies the whole-network path's
cutoff test per arm. The default keeps the single node-position query on
purpose, to agree with the team node_force that Topology uses, so it is only
reported here, not asserted.

THE CUTOFF IS THE POINT

That defect is invisible while the whole configuration fits inside the cutoff,
since then no pair the bins return should be rejected anyway. It only appears
once a pair lands inside the +-1 bin neighborhood and outside the cutoff, which
is the ordinary case in any run large enough to need a cutoff. So this test runs
CUTOFF_MODEL twice: once with a cutoff that keeps every pair, and once with one
that keeps a minority of them. Only the second separates the two branches: the
default parts by 3.7e-02 relative there, and agrees on every all-pairs line.

Only the loop exercises that filtering. Every pair of segments in a star shares
the central node, so a star has no non-adjacent pairs for a cutoff to reject and
its tight-cutoff lines pass either way. The stars are here for the arm count
instead: they are what puts the per-node route on a node Topology would split.

Equality is required to rounding, not to a physical tolerance. Both routes sum
the same pair contributions in the same units, so a correct OneNodeForce differs
from its NodeForce row only in the order the terms are added. Anything larger is
a missing or duplicated term. The tolerance scales with the largest force in the
network, these being absolute quantities whose size depends on the
configuration.

Multi-arm nodes are here for the reason they are in the pydis test: the stored
loop is 450 two-arm nodes, and OneNodeForce exists for the four-and-more-arm
case that Topology acts on, which the loop never presents.

DDD_FFT_MODEL is out of scope. Its node force reads a long-range grid built by
pre_compute for the whole network, so per-node and whole-network agreement means
something different there and wants a case built for it, with a periodic cell.

    python3 test_onenodeforce_exadis.py
"""

import os
# Kokkos reads this at initialize(), so it must be set before pyexadis is
# imported, and therefore before any framework import that might pull it in.
# ExaDiS sums the pair forces with a parallel reduction, so the thread count
# decides the order of the additions and with it the last bits of the answer.
os.environ['OMP_NUM_THREADS'] = '1'

import contextlib
import io
import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'

import numpy as np
from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.testing import (report, report_close, quiet_native_output,
                               kokkos_summary)

RN_FILE = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the same material constants the other tests in this folder use
STATE = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0,
         "maxseg": 2.0, "minseg": 0.5, "rann": 0.5}

# A non-zero applied stress, because the Peach-Koehler term is part of what
# OneNodeForce has to reproduce and a zero stress would let an implementation
# omit it and still pass. Voigt order: xx, yy, zz, yz, xz, xy.
APPLIED_STRESS = np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0])

# Ec is left at the ExaDiS default, which is the pydis default too: both read
# mu/(4*pi)*log(a/0.1) when none is given, ParaDiS Initialize.c:625 and ExaDiS
# force_core.h CoreDefault being the same expression. So the core term is in the
# comparison here exactly as it is in the pydis test.

# The cells the pydis test uses, so that the two run on the same geometry: ten
# times the configuration's largest extent for the loop, which is 1027.2, and a
# fixed 2000 for the stars. Non-periodic in both, so no image enters.
LOOP_BOX_FACTOR = 10.0
STAR_BOX = 2000.0

# CUTOFF_MODEL needs a finite cutoff where the pydis test simply passes None, so
# 300 stands in for "no truncation": the loop measures 168 across and a star
# 200, so every pair is inside it and both routes see the set pydis sees.
CUTOFF_ALL_PAIRS = 300.0

# BEYOND THE PYDIS TEST. Nothing on the pydis side corresponds to this, and it
# is the case with teeth, so it is kept as an extra: 20 keeps about 12% of the
# loop's pairs, which is what puts a pair inside the neighbor bins and outside
# the cutoff. See the module docstring.
CUTOFF_TIGHT = 20.0

# cell edge, as a multiple of cutoff + maxseg. 3 is the smallest value that
# keeps the ExaDiS neighbor bins unclamped; below it ExaDiS drops pairs silently
# and the two routes would then disagree for a reason that is not the one under
# test. Both cells above clear it, and check_cutoff_maxseg asserts it per call.
BOX_BINS = 3.0

# Modes where the two routes agree bitwise, and are held to it. LineTension
# reaches a node's force from its own two arms by the same additions in the same
# order on both routes, with no pair sum to re-associate, so it is exact on every
# configuration here including the multi-arm stars. Held to zero so that stops
# being an accident nobody would notice losing.
EXACT_MODES = ('LineTension',)

# Relative to the largest nodal force in the network. Not a rounding
# tolerance: the CUTOFF_MODEL rows carry about 1.1e-08, and the cause is
# measured. The pair kernel is not symmetric in its arguments: passing the same
# two segments the other way round moves a node's force by up to 9.0e-09 in
# pydis' port of it. NodeForce evaluates each pair once with the lower index
# first, while OneNodeForce always puts the node's own arm first, so about half
# the pairs are evaluated the other way round. 1e-7 leaves about 8x margin and
# sits far below the 3.7e-02 the default route gives at the tight cutoff.
TOL_REL = 1e-7

# Multi-arm configurations, the same stars the pydis test uses so that the two
# tests speak about the same shapes. Each entry is (arm directions, Burgers
# vectors) for one free node joined to that many pinned ends. The Burgers
# vectors sum to zero, as conservation at the node requires, and none is
# parallel to its own arm, the glide plane b x xi being undefined for a screw
# arm. Three arms is what a junction leaves behind, four is what a collision of
# two lines creates and what Topology acts on, five is one past that.
STARS = {
    3: ([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[1, 1, 0], [0, -1, 1], [-1, 0, -1]]),
    4: ([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
        [[0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]]),
    5: ([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0]],
        [[0, 1, 1], [0, -1, 1], [1, 0, 1], [0, 1, -1], [-1, -1, -2]]),
}

STAR_ARM_LENGTH = 100.0

PINNED, FREE = 7, 0


def default_route_residual(calforce, N, state, tags, many, scale):
    """default_route_residual: match_global=False's residual, relative to scale

    Reported, not asserted; see MATCH_GLOBAL in the module docstring.
    """
    one = np.array([calforce.OneNodeForce(N, state, t, update_state=False)
                    for t in tags])
    return float(np.max(np.abs(one - many)))/scale


def cutoff_note(cutoff):
    """cutoff_note: what this block's cutoff does, for the header line

    So a reader meets "cutoff = 20, which truncates" before the
    CUTOFF_MODEL(20) rows it labels, rather than after them.
    """
    if not cutoff:
        return ', no cutoff'
    if cutoff >= CUTOFF_ALL_PAIRS:
        return ', cutoff = %g, which keeps every pair' % cutoff
    return ', cutoff = %g, which truncates most pairs' % cutoff


def quiet_check_cutoff_maxseg(h, cutoff, maxseg):
    """quiet_check_cutoff_maxseg: the neighbor-bin check without its report

    The check still raises when the cell is too small for the bins; only its
    four lines of commentary are dropped, so this test's output stays the shape
    of test_onenodeforce_pydis.py's.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        check_cutoff_maxseg(h, cutoff, maxseg)


def loop_config():
    """loop_config: the stored loop, as (nodes, segs, longest segment)"""
    rn = np.loadtxt(input_dir / RN_FILE)[:, 1:4]
    links = np.loadtxt(input_dir / LINKS_FILE)
    nodes = np.hstack((rn, np.full((rn.shape[0], 1), FREE, dtype=float)))
    span = float(np.max(rn.max(axis=0) - rn.min(axis=0)))
    # exadis wants [n1, n2, bx, by, bz, nx, ny, nz]. The glide planes are left
    # at zero: cross(b, t) is not usable on a loop, vanishing wherever the line
    # runs screw, and neither force mode reads them.
    segs = np.zeros((links.shape[0], 8))
    segs[:, :5] = links[:, :5]
    ends = rn[links[:, :2].astype(int)]
    maxseg = float(np.max(np.linalg.norm(ends[:, 1] - ends[:, 0], axis=1)))
    return nodes, segs, maxseg, LOOP_BOX_FACTOR*span


def star_config(arm_dirs, burgs):
    """star_config: one free node of len(arm_dirs) arms, pinned at the ends"""
    k = len(arm_dirs)
    nodes = np.zeros((k+1, 4))
    nodes[0, 3] = FREE
    nodes[1:, 0:3] = STAR_ARM_LENGTH*np.array(arm_dirs, dtype=float)
    nodes[1:, 3] = PINNED

    segs = np.zeros((k, 8))
    for i, (d, b) in enumerate(zip(arm_dirs, burgs)):
        n = np.cross(np.array(b, dtype=float), np.array(d, dtype=float))
        segs[i, :] = np.concatenate(([0, i+1], b, n/np.linalg.norm(n)))
    return nodes, segs, STAR_ARM_LENGTH, STAR_BOX


def network(nodes, segs, maxseg, box, cutoff):
    """network: the configuration centred in the cell the pydis test gives it

    Non-periodic, so no image enters and both routes see the same pairs. The
    cell is the configuration's own, ten times its largest extent for the loop
    and STAR_BOX for a star, so the two tests run the same geometry. The only
    extra demand here is that it leave ExaDiS' neighbor bins unclamped, which
    both of those clear.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet

    L = box
    h = L*np.eye(3)
    if cutoff:
        quiet_check_cutoff_maxseg(h, cutoff, maxseg)
    nodes = nodes.copy()
    nodes[:, :3] += 0.5*L - 0.5*(nodes[:, :3].max(axis=0)
                                 + nodes[:, :3].min(axis=0))
    cell = pyexadis.Cell(h=h, is_periodic=[False]*3)
    # quiet: a star's pinned ends carry one arm each, so ExaDiS reports Burgers
    # vector non-conservation there. That is the configuration being what it is,
    # not a fault, and it says nothing about the two routes.
    with quiet_native_output():
        return DisNetManager(ExaDisNet(cell, nodes, segs))


def check_case(mode, shown, config, label, cutoff=None):
    """check_case: OneNodeForce against NodeForce on one (mode, configuration)

    Returns (ran, ok): ran is False when a route did not run at all, in which
    case nothing was checked and ok is False. Neither is skipped, for the reason
    the pydis test gives: a green suite while Topology has no usable
    OneNodeForce answers "is this working" with yes.
    """
    from pyexadis_base import CalForce

    nodes, segs, maxseg, box = config
    N = network(nodes, segs, maxseg, box, cutoff)
    if label == 'loop':
        print("load_network: %d nodes, %d links, box = %.1f%s"
              % (nodes.shape[0], segs.shape[0], box, cutoff_note(cutoff)))
    # the extra cutoff rides in the mode column, so the columns stay the two
    # the pydis test prints
    name = "  %-22s %-13s OneNodeForce == NodeForce," % (shown, label)

    state = dict(STATE, applied_stress=APPLIED_STRESS)
    kwargs = {'cutoff': cutoff} if cutoff else {}
    # quiet: constructing a force module prints ExaDiS' "Setting rtol" line
    with quiet_native_output():
        calforce = CalForce(state=state, force_mode=mode, **kwargs)
    try:
        # NodeForce first, as a run does it: the pair list CUTOFF_MODEL's
        # per-node route reads is built by the pre_compute inside this call
        state = calforce.NodeForce(N, state)
    except Exception as exc:
        print("     %s: %s" % (type(exc).__name__, str(exc).splitlines()[0][:100]))
        return False, report("%s NodeForce runs" % name, False)

    many = np.array(state["nodeforces"])
    tags = [tuple(int(x) for x in t) for t in state["nodeforcetags"]]
    try:
        # update_state=False, or each call would overwrite the very row it is
        # about to be compared against
        one = np.array([calforce.OneNodeForce(N, state, t, update_state=False,
                                              match_global=True)
                        for t in tags])
    except Exception as exc:
        print("     %s: %s" % (type(exc).__name__, str(exc).splitlines()[0][:100]))
        return False, report("%s OneNodeForce runs" % name, False)

    scale = max(float(np.max(np.abs(many))), 1e-30)
    tol = 0.0 if mode in EXACT_MODES else TOL_REL*scale
    ok = report_close("%s %d nodes" % (name, len(tags)), one, many, tol)
    if cutoff and label == 'loop':
        print("     match_global=False: %.2e, not asserted"
              % default_route_residual(calforce, N, state, tags, many, scale))
    if not ok:
        worst = int(np.argmax(np.abs(one - many).max(axis=1)))
        print("     worst node %s: OneNodeForce %s"
              % (str(tags[worst]), np.array2string(one[worst], precision=6)))
        print("     %22s NodeForce    %s"
              % ('', np.array2string(many[worst], precision=6)))
        # rounding or a missing term: the two are orders apart, see TOL_REL
        print("     %22s relative to the largest force: %.4e"
              % ('', float(np.max(np.abs(one - many)))/scale))
    return True, ok


def cases():
    """cases: the (mode, column label, cutoff) rows and the configurations

    The first two rows are the pydis test's LineTension and elasticity rows,
    CUTOFF_MODEL being what Elasticity_SBA maps to. The third has no
    counterpart there; see CUTOFF_TIGHT.
    """
    modes = [('LineTension', 'LineTension', None),
             ('CUTOFF_MODEL', 'CUTOFF_MODEL', CUTOFF_ALL_PAIRS),
             ('CUTOFF_MODEL', 'CUTOFF_MODEL(%g)' % CUTOFF_TIGHT, CUTOFF_TIGHT)]
    configs = [(loop_config(), 'loop')]
    configs += [(star_config(*STARS[k]), '%d-arm star' % k)
                for k in sorted(STARS)]
    return modes, configs


def main():
    import pyexadis
    print("OneNodeForce against NodeForce, exadis")
    print("applied stress (Voigt) = %s"
          % np.array2string(APPLIED_STRESS, precision=1))

    with quiet_native_output() as buf:
        pyexadis.initialize()
    for line in kokkos_summary(buf.text):
        print(line)
    print("")

    modes, configs = cases()
    ok, checked = True, []
    for mode, shown, cutoff in modes:
        for config, label in configs:
            ran, case_ok = check_case(mode, shown, config, label, cutoff=cutoff)
            ok &= case_ok
            if ran:
                checked.append((shown, label))
        print("")

    pyexadis.finalize()

    want = len(modes)*len(configs)
    print("checked %d of %d (mode, configuration) combinations"
          % (len(checked), want))
    if len(checked) != want:
        # counted as failures, not skipped: a green result while a mode has no
        # usable OneNodeForce would answer "is this working" with yes while
        # Topology still cannot score a trial split
        done = set(checked)
        missing = [(shown, l) for _, shown, _ in modes for _, l in configs
                   if (shown, l) not in done]
        print("NOT checked, counted as failures: %s"
              % ', '.join('%s/%s' % ml for ml in missing))
    return bool(ok) and len(checked) == want


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("\ntest " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
