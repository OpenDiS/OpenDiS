"""OneNodeForce against NodeForce, for every force mode that has one.

CalForce offers two ways to get the force on a node. NodeForce evaluates
the whole network and returns a force per node; OneNodeForce evaluates a
single node. Topology needs the second one: split_multi_nodes tries a
trial split of every node with four or more arms and asks for the force
on each of the two candidate nodes, which would be absurdly expensive if
it had to re-evaluate the network twice per candidate.

Two functions computing the same quantity by different routes is exactly
the arrangement that drifts, so this pins the invariant:

    OneNodeForce(tag) == the NodeForce() entry for tag, for every node

Equality is required to rounding rather than to a physical tolerance.
Both routes sum the same segment contributions in the same units, and a
correct OneNodeForce differs from its NodeForce row only by the order the
terms are added, so anything larger is a missing or duplicated term
rather than arithmetic. The tolerance is scaled by the largest force in
the network, since these are absolute quantities whose size depends on
the configuration.

EVERY MODE IS REQUIRED, NOT OPTIONAL

A mode whose OneNodeForce raises NotImplementedError fails this test, and
so does one whose NodeForce cannot run. Neither is skipped. Skipping
would let the suite go green while Topology still has no usable
OneNodeForce, which is the condition that stops
examples/03_binary_junction/test_binary_junction_pydis_elast.py at step 0
and 02_frank_read_src at step 270. A green test that coexists with a
missing function is worse than no test, because it answers the question
"is this working" with yes.

    python3 test_onenodeforce_pydis.py
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'

import numpy as np
from pydis.disnet import DisNet, DisNode, Cell
from pydis.calforce.calforce_disnet import CalForce
from framework.disnet_manager import DisNetManager
from framework.testing import report, report_close

RN_FILE = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the same material constants the other tests in this folder use
STATE = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0,
         "maxseg": 2.0, "minseg": 0.5, "rann": 0.5}

# A non-zero applied stress, because the Peach-Koehler term is part of
# what OneNodeForce has to reproduce and a zero stress would let an
# implementation omit it and still pass. Voigt order: xx, yy, zz, yz, xz,
# xy.
APPLIED_STRESS = np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0])

# No cutoff, so OneNodeForce and NodeForce sum over the same set of pairs
# and nothing about the truncation can influence the comparison. Whether
# the cutoff selects the right pairs is a separate question with its own
# test; this one is only about the two routes agreeing.
CUTOFF = None

# relative to the largest nodal force in the network; see the module
# docstring for why this is a rounding tolerance and not a physical one
TOL_REL = 1e-12

# Modes where the two routes agree exactly, and are held to it.
#
# LineTension is pure numpy on both routes and reaches each node's force by the
# same additions in the same order, so it is bitwise, on every configuration
# here including the multi-arm stars. Held to zero so that stops being an
# accident nobody would notice losing.
#
# The elasticity modes cannot be: OneNodeForce sums only the pairs touching the
# node while NodeForce sums all of them and accumulates, so the same
# contributions arrive re-associated. That shows up on the 450-node loop, at
# 1.1e-04 for Elasticity_SBA against a scale of 1.7e+13, and it is a rounding
# difference rather than a disagreement.
EXACT_MODES = ('LineTension',)

MODES = ['LineTension', 'Elasticity_SBA', 'Elasticity_SBN1_SBA']

# Multi-arm configurations, and the reason they are here.
#
# The stored loop below is 450 nodes and every one of them has two arms, so
# on its own this test never exercises the case OneNodeForce exists for:
# Topology calls it only on nodes with four or more arms, when it tries a
# split. A two-arm node reaches a different, shorter path through the
# accumulation, so 450 of them is 450 samples of the easy case.
#
# Each entry is (arm directions, Burgers vectors) for a star: one free node
# at the centre joined to that many pinned ends. The Burgers vectors sum to
# zero, as conservation at the node requires, and none is parallel to its own
# arm, since the glide plane is b x xi and would be undefined for a screw
# arm. Three arms is what a junction leaves behind, four is what a collision
# of two lines creates and what Topology acts on, five is one past that.
STARS = {
    3: ([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[1, 1, 0], [0, -1, 1], [-1, 0, -1]]),
    4: ([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
        [[0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]]),
    5: ([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0]],
        [[0, 1, 1], [0, -1, 1], [1, 0, 1], [0, 1, -1], [-1, -1, -2]]),
}

STAR_ARM_LENGTH = 100.0
STAR_BOX = 2000.0


def load_network():
    """load_network: the stored loop, as a DisNetManager"""
    rn = np.loadtxt(input_dir / RN_FILE)[:, 1:]
    links = np.loadtxt(input_dir / LINKS_FILE)
    # a box comfortably larger than the loop, non-periodic, so that no
    # image or cutoff effect enters and the two routes see the same pairs
    extent = float(np.max(rn[:, :3].max(axis=0) - rn[:, :3].min(axis=0)))
    cell = Cell(h=10.0*extent*np.eye(3), is_periodic=[False]*3)
    rn = rn.copy()
    rn[:, 0:3] += cell.center()
    print("load_network: %d nodes, %d links, box = %.1f"
          % (rn.shape[0], links.shape[0], 10.0*extent))
    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))


def star_network(arm_dirs, burgs):
    """star_network: one free node of len(arm_dirs) arms, pinned at the ends"""
    PINNED = DisNode.Constraints.PINNED_NODE
    FREE = DisNode.Constraints.UNCONSTRAINED
    cell = Cell(h=STAR_BOX*np.eye(3), is_periodic=[False]*3)
    k = len(arm_dirs)
    rn = np.zeros((k+1, 4))
    rn[0, 3] = FREE
    for i, d in enumerate(arm_dirs):
        rn[i+1, 0:3] = STAR_ARM_LENGTH*np.array(d, dtype=float)
        rn[i+1, 3] = PINNED
    rn[:, 0:3] += cell.center()

    links = np.zeros((k, 8))
    for i, (d, b) in enumerate(zip(arm_dirs, burgs)):
        xi = np.array(d, dtype=float)
        b = np.array(b, dtype=float)
        n = np.cross(b, xi)
        n = n / np.linalg.norm(n)
        links[i, :] = np.concatenate(([0, i+1], b, n))
    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))


def check_mode(mode, network=None, label=None):
    """check_mode: one force mode, or a report that it has no OneNodeForce

    Returns (ran, ok): ran is False when the mode has no OneNodeForce
    implementation, in which case ok is True and nothing was checked.
    """
    net = network if network is not None else load_network()
    state = dict(STATE)
    state["applied_stress"] = APPLIED_STRESS
    calforce = CalForce(force_mode=mode, state=state, cutoff=CUTOFF)

    try:
        state = calforce.NodeForce(net, state)
    except Exception as exc:
        # A mode whose NodeForce does not run has no force to compare
        # against, so nothing here can be verified. That is a failure, not a
        # reason to pass quietly.
        print("     %s: %s" % (type(exc).__name__,
                               str(exc).splitlines()[0][:100]))
        return False, report("  %-22s %-13s NodeForce runs"
                             % (mode, label or 'loop'), False)
    nodeforce_dict = state["nodeforce_dict"]
    tags = list(nodeforce_dict.keys())
    scale = max(float(np.max(np.abs(np.array(list(
        nodeforce_dict.values()))))), 1e-30)

    # every node, not a sample: the cheap case to get wrong is a node with
    # an unusual arm count, and on this loop they are all two-arm, so the
    # cost of being thorough is one pass either way
    one = np.zeros((len(tags), 3))
    for i, tag in enumerate(tags):
        try:
            one[i] = calforce.OneNodeForce(net, state, tag,
                                           update_state=False)
        except NotImplementedError as exc:
            return False, report("  %-22s %-13s OneNodeForce is implemented"
                                 % (mode, label or 'loop'), False)

    many = np.array([nodeforce_dict[t] for t in tags])
    tol = 0.0 if mode in EXACT_MODES else TOL_REL*scale
    ok = report_close("  %-22s %-13s OneNodeForce == NodeForce, %d nodes"
                      % (mode, label or 'loop', len(tags)),
                      one, many, tol)
    if not ok:
        worst = int(np.argmax(np.abs(one - many).max(axis=1)))
        print("     worst node %s: OneNodeForce %s"
              % (str(tags[worst]), np.array2string(one[worst], precision=6)))
        print("     %22s NodeForce    %s"
              % ('', np.array2string(many[worst], precision=6)))
    return True, ok


def main():
    print("OneNodeForce against NodeForce, pydis")
    print("applied stress (Voigt) = %s\n"
          % np.array2string(APPLIED_STRESS, precision=1))

    configs = [(None, 'loop')]
    configs += [(star_network(*STARS[k]), '%d-arm star' % k)
                for k in sorted(STARS)]

    ok = True
    checked = []
    for mode in MODES:
        for net, label in configs:
            ran, mode_ok = check_mode(mode, network=net, label=label)
            ok &= mode_ok
            if ran:
                checked.append((mode, label))
        print("")

    want = len(MODES) * len(configs)
    print("checked %d of %d (mode, configuration) combinations"
          % (len(checked), want))
    if len(checked) != want:
        # counted as failures, not skipped: a green result while a mode has
        # no usable OneNodeForce would answer "is this working" with yes
        # while Topology still cannot split a multi-arm node
        done = set(checked)
        missing = [(m, l) for m in MODES for _, l in configs
                   if (m, l) not in done]
        print("NOT checked, counted as failures: %s"
              % ', '.join('%s/%s' % ml for ml in missing))
    return bool(ok) and len(checked) == want

if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("\ntest " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
