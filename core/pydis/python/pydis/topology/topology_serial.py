"""@package docstring
Topology_Serial: split multi nodes by maximum power dissipation, ParaDiS criterion

Matches the behavior of the 'TopologySerial' topology model of ExaDiS, and
follows ParaDiS SplitMultiNodes in Topology.c, which is where the criterion
comes from.

Tries every way of splitting a multi-arm node and takes the one that releases
the most power, as Topology_MaxDiss does, but measures that power with the two
new nodes moved apart rather than on top of each other. The zero-separation
split serves only to read a separation direction off the velocities the
mobility law then gives them; they are moved apart by split_dist along that
direction and judged there. This is the difference that matters, because at
zero separation the trial cannot see the energy released by pulling two
dislocation lines apart.
"""

import numpy as np
from copy import deepcopy

from ..disnet import Tag
from framework.disnet_manager import DisNetManager
from .topology_ops import (build_split_list, exempt_from_collisions,
                           split_node_and_update_forces)

# ParaDiS eps in SplitMultiNodes, used both to break the power comparison and
# to push the separation just past split_dist.
SPLIT_EPS = 1.0e-12


class TopologyParams:
    """TopologyParams: lengths, tolerances and noise floors for splitting multi nodes

    Each attribute is named for the ParaDiS parameter it stands for, so the
    criterion below reads against Topology.c without a translation table. Used
    by 'Serial' only; 'MaxDiss' needs nothing from state.
    """

    def __init__(self, state: dict) -> None:
        self.rann = state.get("rann", None)
        self.minseg = state.get("minseg", None)
        self.a = state.get("a", None)
        self.mob = state.get("mob", 1.0)
        self.dt = state.get("dt", None)

        # ParaDiS Initialize.c sets rTol from the core radius and derives rann
        # from rTol, not the other way round. Taking rtol from rann instead
        # gives a much larger value wherever the two are set independently,
        # and vnoise below is proportional to it.
        self.rtol = state.get("rtol", None)
        if self.rtol is None and self.a is not None:
            self.rtol = 0.25*self.a

        # Both noise floors are otherwise set by the expressions below. They
        # are left settable because those expressions carry ParaDiS units.
        self._epsvel = state.get("epsvel", None)
        self._vnoise = state.get("vnoise", None)

        missing = [name for name in ("rann", "minseg", "rtol")
                   if getattr(self, name) is None]
        if missing:
            raise ValueError("TopologyParams: state has no %s"
                             % ", ".join(missing))

    @property
    def split_dist(self) -> float:
        """how far apart a split moves the two nodes, ParaDiS minDist"""
        return 2.0*self.rann

    @property
    def short_seg(self) -> float:
        """a node with an arm shorter than this is left alone"""
        return min(5.0, 0.1*self.minseg)

    @property
    def vnoise(self) -> float:
        """force noise floor subtracted in the power test, ParaDiS rTol/deltaTT

        Zero until the first time step has been taken, since dt is not known
        before then.
        """
        if self._vnoise is not None:
            return self._vnoise
        return 0.0 if not self.dt else self.rtol/self.dt

    @property
    def epsvel(self) -> float:
        """a trial split whose nodes move slower than this is not a candidate

        ParaDiS reads this as the velocity reached under an applied stress of
        0.001 MPa.
        """
        if self._epsvel is not None:
            return self._epsvel
        return 0.001e6*self.mob


def store_node_force(state: dict, tag: Tag, f: np.ndarray) -> None:
    """store_node_force: record a nodal force the way OneNodeForce does

    The mobility law rebuilds its dictionary from the nodeforces array, so a
    force written only to the dictionary does not survive the next mobility
    call.
    """
    state["nodeforce_dict"][tag] = f
    if "nodeforces" in state and "nodeforcetags" in state:
        tags = state["nodeforcetags"]
        ind = np.where((tags[:,0]==tag[0]) & (tags[:,1]==tag[1]))[0]
        if ind.size == 1:
            state["nodeforces"][ind[0]] = f
        else:
            state["nodeforces"] = np.vstack((state["nodeforces"], f))
            state["nodeforcetags"] = np.vstack((state["nodeforcetags"], tag))
    else:
        state["nodeforces"] = np.array([f])
        state["nodeforcetags"] = np.array([tag])


def node_force_from_arms(G, segforce_dict: dict, node: Tag, orig_tag: Tag) -> np.ndarray:
    """node_force_from_arms: sum the stored segment forces of a node's arms

    The arms were attached to orig_tag before the split, which is how they are
    keyed in segforce_dict. The segment joining the two split nodes has no
    entry, and contributes nothing while it has zero length.
    """
    f = np.zeros(3)
    for nbr in G.neighbors_tags(node):
        if (orig_tag, nbr) in segforce_dict:
            f += segforce_dict[(orig_tag, nbr)][0:3]
        elif (nbr, orig_tag) in segforce_dict:
            f += segforce_dict[(nbr, orig_tag)][3:6]
    return f


def split_node_reusing_forces(G, state, tag, nbrs_to_split, mobility):
    """split_node_reusing_forces: split a node in place without calling the force kernel

    Both nodes stay where the original was, so every segment force is unchanged
    and each node's force is just the sum over its own arms. ParaDiS reads these
    from its backup copy for the same reason. Asking the kernel instead would
    mean evaluating the force on the zero-length segment between the two nodes,
    which is not defined and comes back as nan.
    """
    pos0 = G.nodes(tag).R.copy()
    node1, node2 = G.split_node(tag, pos0.copy(), pos0.copy(), nbrs_to_split)
    segforce_dict = state["segforce_dict"]
    for node in (node1, node2):
        store_node_force(state, node,
                         node_force_from_arms(G, segforce_dict, node, tag))
    state = mobility.Mobility(DisNetManager(G), state)
    return state, node1, node2


def has_short_arm(G, tag: Tag, short_seg: float) -> bool:
    """has_short_arm: whether any arm of a node is shorter than short_seg

    A node with a very short arm is left unsplit. Splitting it would leave an
    even shorter segment, which oscillates at high velocity and drags the
    timestep down with it. ParaDiS calls this a temporary kludge and compares
    the squared arm length against short_seg rather than the length; the same
    comparison is kept here.
    """
    for nbr in G.neighbors_tags(tag):
        d = G.seg_vector(tag, nbr)
        if np.dot(d, d) < short_seg:
            return True
    return False


def split_direction(vel1: np.ndarray, vel2: np.ndarray, epsvel: float):
    """split_direction: unit vector along which to move two split nodes apart

    Taken from whichever of the two moves faster, pointing the way that node is
    already going, so the separation follows the motion the mobility law asked
    for rather than an arbitrary direction.

    Returns (dirvec, move_first), where move_first says which of the two nodes
    is the one displaced, or None when neither node moves. A partition whose
    nodes both sit still is not a candidate.
    """
    vd1, vd2 = np.dot(vel1, vel1), np.dot(vel2, vel2)
    # a comparison against nan is False, so non-finite velocities would
    # otherwise pass the test below rather than failing it
    if not (np.isfinite(vd1) and np.isfinite(vd2)):
        return None
    if max(vd1, vd2) <= epsvel*epsvel:
        return None
    if vd1 > vd2:
        return -vel1/np.sqrt(vd1), True
    return vel2/np.sqrt(vd2), False


def set_connecting_plane(G, tag1: Tag, tag2: Tag, dirvec: np.ndarray) -> None:
    """set_connecting_plane: give the segment joining two split nodes the glide
    plane implied by a separation direction

    A trial split leaves this segment without a plane, because the two nodes
    coincide and there is no line direction to derive one from. Once a direction
    is known the plane follows, and the velocities are worth recomputing under
    the glide constraint it imposes.
    """
    if not G.has_segment(tag1, tag2):
        return
    seg = G.segments((tag1, tag2))
    seg.plane_normal = G.find_precise_glide_plane(seg.burg_vec_from(tag1),
                                                  dirvec)


def separated_positions(pos0: np.ndarray, dirvec: np.ndarray,
                        move_first: bool, split_dist: float):
    """separated_positions: where the two split nodes end up

    One node moves and the other stays, which is what ParaDiS does, so the pair
    straddles the original position rather than being centred on it. Either way
    the second node ends up split_dist along dirvec from the first.
    """
    step = split_dist*(1.0 + SPLIT_EPS)*dirvec
    if move_first:
        return pos0 - step, pos0.copy()
    return pos0.copy(), pos0 + step


def evaluate_trial_split(G, state, tag: Tag, nbrs_to_split: list,
                         force, mobility, params: TopologyParams):
    """evaluate_trial_split: the power released by one way of splitting a node

    Returns (power, pos1, pos2), or None when this partition is not a
    candidate. Runs on a copy, so the network passed in is left alone.

    Three stages, following ParaDiS: split with both nodes at the original
    position, read a separation direction off the velocities that follow, then
    split again with the nodes moved apart and judge the result there.
    """
    pos0 = G.nodes(tag).R.copy()

    G_trial = G.copy()
    state_trial = deepcopy(state)
    state_trial, node1, node2 = split_node_reusing_forces(
        G_trial, state_trial, tag, nbrs_to_split, mobility)

    chosen = split_direction(state_trial["vel_dict"][node1],
                             state_trial["vel_dict"][node2], params.epsvel)
    if chosen is None:
        return None

    # the direction fixes the glide plane of the connecting segment, and the
    # velocities are re-read with that constraint in force
    set_connecting_plane(G_trial, node1, node2, chosen[0])
    state_trial = mobility.Mobility(DisNetManager(G_trial), state_trial)

    chosen = split_direction(state_trial["vel_dict"][node1],
                             state_trial["vel_dict"][node2], params.epsvel)
    if chosen is None:
        return None
    dirvec, move_first = chosen

    pos1, pos2 = separated_positions(pos0, dirvec, move_first, params.split_dist)

    G_split = G.copy()
    state_split = deepcopy(state)
    state_split, node1, node2 = split_node_and_update_forces(
        G_split, state_split, tag, pos1, pos2, nbrs_to_split, force, mobility)

    f1, f2 = state_split["nodeforce_dict"][node1], state_split["nodeforce_dict"][node2]
    v1, v2 = state_split["vel_dict"][node1], state_split["vel_dict"][node2]

    # keep only splits that go on separating; a pair moving back together would
    # be undone on the next step
    if np.dot(v2 - v1, dirvec) <= 0.0:
        return None

    power = (np.dot(f1, v1) + np.dot(f2, v2)
             - params.vnoise*(np.linalg.norm(f1) + np.linalg.norm(f2)))
    if not np.isfinite(power):
        return None
    return power, pos1, pos2


def Topology_Serial(G, tag: Tag, state: dict, force, mobility,
                    params: TopologyParams) -> dict:
    """Topology_Serial: try to split multi-arm node in different ways
        and select the way that maximizes the power dissipation
    """
    if has_short_arm(G, tag, params.short_seg):
        return state

    nbrs = list(G.neighbors_tags(tag))
    nbr_idx_list = build_split_list(G.out_degree(tag))

    # the unsplit node is the baseline the trials have to beat
    power0 = np.dot(state["nodeforce_dict"][tag], state["vel_dict"][tag])

    best = None
    for nbr_idx in nbr_idx_list:
        nbrs_to_split = [nbrs[i] for i in nbr_idx]
        trial = evaluate_trial_split(G, state, tag, nbrs_to_split,
                                     force, mobility, params)
        if trial is None:
            continue
        if best is None or trial[0] > best[0]:
            best = trial + (nbrs_to_split,)

    if best is None or best[0] - power0 <= SPLIT_EPS:
        return state

    power, pos1, pos2, nbrs_to_split = best
    state, split_node1, split_node2 = split_node_and_update_forces(
        G, state, tag, pos1, pos2, nbrs_to_split, force, mobility)
    exempt_from_collisions(state, split_node1, split_node2)

    return state
