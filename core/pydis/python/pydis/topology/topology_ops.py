"""@package docstring
Topology_Ops: low-level operations shared by the topology split modes

Kept apart from topology_disnet so that each split mode lives in its own module
and can import these without importing the mode dispatch that selects it.
"""

import itertools

from ..disnet import DisNode, Tag
from ..util.arm_order import split_node_matching_exadis
from framework.disnet_manager import DisNetManager


def build_split_list(n: int) -> list:
    """build_split_list: build a list of length n with True's and False's
       there must be at least two True's and at least 2 False's
    """
    indices = list(range(n))
    all_bool_list = list(itertools.product([True, False], repeat=n))
    # the first element must be selected to avoid double counting
    selected_list = [list(itertools.compress(indices,item)) for item in all_bool_list if item[0] and sum(item) >= 2 and sum(item) <= n-2]
    return selected_list


def init_topology_exemptions(G, state) -> None:
    """init_topology_exemptions: initialize the topology exemptions
    """
    nodeflag_dict = {}
    for tag in G.all_nodes_tags():
        nodeflag_dict[tag]  = DisNode.Flags.CLEAR
        nodeflag_dict[tag] &= ~(DisNode.Flags.NO_COLLISIONS | DisNode.Flags.NO_MESH_COARSEN)
    state["nodeflag_dict"] = nodeflag_dict
    return state


def exempt_from_collisions(state, *tags) -> None:
    """exempt_from_collisions: mark nodes as exempt from collisions this time step

    A split is otherwise liable to be undone by the collision handler that runs
    after it.
    """
    for tag in tags:
        if not tag in state['nodeflag_dict']:
            state['nodeflag_dict'][tag] = DisNode.Flags.CLEAR
        state['nodeflag_dict'][tag] |= DisNode.Flags.NO_COLLISIONS


def ensure_node_forces(G, state, force, mobility):
    """ensure_node_forces: fill in forces and velocities for nodes that have none

    The collision handler runs before the split when collide_before_split is
    set, and merges nodes without touching nodeforce_dict, so a multi-arm node
    created by a collision reaches the split with no force recorded. Both split
    modes need one to score the unsplit node.

    The fill-in is network-wide rather than just the node about to be split,
    because Mobility rebuilds every velocity from nodeforce_dict and raises on
    any node missing from it. Nothing happens when the entries are already
    there, so a run whose collisions never produce a multi-node is unaffected.
    """
    forces, vels = state["nodeforce_dict"], state.get("vel_dict", {})
    stale = [tag for tag in G.all_nodes_tags()
             if tag not in forces or tag not in vels]
    if not stale:
        return state
    for tag in stale:
        if tag not in forces:
            force.OneNodeForce(DisNetManager(G), state, tag, update_state=True)
    return mobility.Mobility(DisNetManager(G), state)


def split_node_and_update_forces(G, state, tag, pos1, pos2, nbrs_to_split, force, mobility):
    """split_node_and_update_forces: split a node and refresh both new nodes' forces
    """
    split_node1, split_node2 = split_node_matching_exadis(G, tag, pos1, pos2, nbrs_to_split)
    # calculate nodal forces and velocities for the trial split
    # To do: pass DisNetManager instead of DisNet in these static methods
    f1 = force.OneNodeForce(DisNetManager(G), state, split_node1, update_state=True)
    f2 = force.OneNodeForce(DisNetManager(G), state, split_node2, update_state=True)
    state = mobility.Mobility(DisNetManager(G), state)
    return state, split_node1, split_node2


def split_multi_nodes(G, state: dict, trial_fn, max_degree=15) -> None:
    """split_multi_nodes: examines all nodes with at least four arms and decides
       if the node should be split and some of the node's arms moved to a new node.
       guarantees sanity after operation

       trial_fn(G, tag, state) decides and performs the split for one node, and
       is what distinguishes the split modes.

       This function calls the lower level split_node() function.
    """
    nodes = list(G.all_nodes_tags())
    for tag in nodes:
        n_degree = G.out_degree(tag)

        if n_degree < 4:
            continue
        elif n_degree > max_degree:
            raise ValueError("split_multi_node: Node %s has more than %d arms" % (str(tag), n_degree))

        state = trial_fn(G, tag, state)

    return state
