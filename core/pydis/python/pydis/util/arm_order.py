"""Splitting a segment the way ExaDiS orders the resulting arms.

The order in which a node's arms are visited is not physics, but it is not
free either: a collision walks a node's arms and keeps the first glide plane
sufficiently independent of those already kept, so when two candidate planes
are near-degenerate the arm order decides which one survives, and the
collision point moves by an ulp.

DisNet.insert_node_between and ExaDiS' SerialDisNet::split_seg build the same
topology in different arm orders. This module holds the wrapper that makes
the two agree, for callers comparing against ExaDiS.

The neutral half of the same problem, reading an arm order out of an ExaDiS
network and imposing it on a DisNet built from exported data, is in
framework.arm_order, whose set_node_arm_order this uses.
"""

from framework.arm_order import set_node_arm_order


def insert_node_matching_exadis(G, tag1, tag2, new_tag, R):
    """insert_node_matching_exadis: split a segment, ExaDiS' arm order kept

    DisNet.insert_node_between, followed by the three reorderings that make
    the result match ExaDiS' SerialDisNet::split_seg. ExaDiS rewires both old
    nodes' connections in place, keeping their slots, and builds the new
    node's connections far node first:

        cnew.add_connection(n2, snew,  1);
        cnew.add_connection(n1, i,   -1);

    PyDiS removes and re-adds instead, which leaves each old node's new arm
    at the end of its list, and orders the new node's arms near node first.
    Both differences are arm order only, and arm order breaks ties between
    near-degenerate glide planes, so both move node positions by an ulp.

    A free function rather than a change to insert_node_between, so that
    disnet.py stays as it is; callers that want ExaDiS' ordering ask for it
    here.
    """
    order1 = [new_tag if nbr == tag2 else nbr
              for nbr in G.neighbors_tags(tag1)]
    order2 = [new_tag if nbr == tag1 else nbr
              for nbr in G.neighbors_tags(tag2)]

    G.insert_node_between(tag1, tag2, new_tag, R)

    set_node_arm_order(G, tag1, order1)
    set_node_arm_order(G, tag2, order2)
    set_node_arm_order(G, new_tag, [tag2, tag1])
