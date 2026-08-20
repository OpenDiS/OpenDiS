"""@package docstring
Topology_MaxDiss: split multi nodes by maximum power dissipation

Tries every way of splitting a multi-arm node and takes the one that releases
the most power, measured with the two new nodes left on top of each other.

Kept as it was so that results produced with it stay reproducible. Two
consequences of measuring at zero separation are worth knowing. The trial
cannot see the energy released by pulling two dislocation lines apart, so a
junction node tends to be left alone; and the force kernel is asked for the
force on the zero-length segment between the two nodes, which is not defined
and comes back as nan, whereupon np.max propagates it and the comparison below
is False. Topology_Serial is the mode that separates the nodes first and avoids
both.
"""

import numpy as np
from copy import deepcopy

from ..disnet import Tag
from .topology_ops import (build_split_list, exempt_from_collisions,
                           split_node_and_update_forces)


def Topology_MaxDiss(G, tag: Tag, state: dict, force, mobility, power_th=1e-3) -> dict:
    """Topology_MaxDiss: try to split multi-arm node in different ways
        and select the way that maximizes the power dissipation
    """
    n_degree = G.out_degree(tag)
    nbrs = list(G.neighbors_tags(tag))
    nbr_idx_list = build_split_list(n_degree)

    power0 = np.dot(state["nodeforce_dict"][tag], state["vel_dict"][tag])

    pos0 = G.nodes(tag).R
    n_splits = len(nbr_idx_list)
    power_diss = np.zeros(n_splits)
    for k in range(n_splits):
        nbrs_to_split = [nbrs[i] for i in nbr_idx_list[k]]

        # make a copy of the network G to make trial splits
        G_trial = G.copy()
        state_trial = deepcopy(state)
        state_trial, split_node1, split_node2 = split_node_and_update_forces(G_trial, state_trial, tag, pos0.copy(), pos0.copy(), nbrs_to_split, force, mobility)

        power_diss[k] = np.dot(state_trial["nodeforce_dict"][split_node1], state_trial["vel_dict"][split_node1]) \
                      + np.dot(state_trial["nodeforce_dict"][split_node2], state_trial["vel_dict"][split_node2])

    if np.max(power_diss) - power0 > power_th:
        # select the split that leads to the maximum power dissipation
        k_sel = np.argmax(power_diss)
        do_split = True
    else:
        do_split = False

    if do_split:
        nbrs_to_split = [nbrs[i] for i in nbr_idx_list[k_sel]]
        state, split_node1, split_node2 = split_node_and_update_forces(G, state, tag, pos0.copy(), pos0.copy(), nbrs_to_split, force, mobility)
        exempt_from_collisions(state, split_node1, split_node2)

    return state
