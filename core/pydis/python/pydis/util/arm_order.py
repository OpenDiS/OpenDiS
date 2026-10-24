"""Splitting a node the way ExaDiS orders the resulting arms.

The order in which a node's arms are visited is not physics, but it is not
free either: a collision walks a node's arms and keeps the first glide plane
sufficiently independent of those already kept, so when two candidate planes
are near-degenerate the arm order decides which one survives, and the
collision point moves by an ulp. For a multi-arm split, arm order feeds a
node force's own summation, which the split-direction decision is built on
directly, so the same sensitivity reaches the split position too.

DisNet.insert_node_between and ExaDiS' SerialDisNet::split_seg build the same
topology in different arm orders, and so do DisNet.split_node and
SerialDisNet::split_node. This module holds the wrappers that make the two
agree, for callers comparing against ExaDiS.

The neutral half of the same problem, reading an arm order out of an ExaDiS
network and imposing it on a DisNet built from exported data, is in
framework.arm_order, whose set_node_arm_order this uses.
"""

from framework.arm_order import set_node_arm_order


def insert_node_matching_exadis(G, tag1, tag2, new_tag, R):
    """insert_node_matching_exadis: split a segment, ExaDiS' arm order kept

    GUARANTEE, load-bearing for keeping this equivalent to plain PyDiS
    rather than a second, divergent implementation: this function calls
    DisNet.insert_node_between UNMODIFIED first, for the actual topology
    change, and only then calls set_node_arm_order, which deregisters and
    re-registers the SAME Edge objects already present -- it cannot create,
    drop, or mutate an edge, only change which key a dict yields first. So
    whatever insert_node_between guarantees (Burgers balance, sane
    structure) still holds after this call; the only observable effect is
    on floating-point summation order in code that later walks a node's
    arms (force, mobility, glide-plane selection), which is the point.
    Keep new callers to this shape -- call the real PyDiS operation, then
    reorder -- rather than adding independent logic here, or this stops
    being a safe drop-in and becomes a second topology implementation to
    keep in sync by hand.

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


def split_node_matching_exadis(G, tag, pos1, pos2, nbrs_to_split, eps_b=1e-3):
    """split_node_matching_exadis: split a multi-arm node, ExaDiS' arm order kept

    GUARANTEE, load-bearing for keeping this equivalent to plain PyDiS
    rather than a second, divergent implementation: this function calls
    DisNet.split_node UNMODIFIED first, for the actual topology change, and
    only then calls set_node_arm_order, which deregisters and re-registers
    the SAME Edge objects already present -- it cannot create, drop, or
    mutate an edge, only change which key a dict yields first. So whatever
    split_node guarantees (Burgers balance, sane structure) still holds
    after this call; the only observable effect is on floating-point
    summation order in code that later walks a node's arms (force,
    mobility, split direction), which is the point. Keep new callers to
    this shape -- call the real PyDiS operation, then reorder -- rather
    than adding independent logic here, or this stops being a safe drop-in
    and becomes a second topology implementation to keep in sync by hand.

    DisNet.split_node, followed by the one reordering that makes the result
    match ExaDiS' SerialDisNet::split_node (core/exadis/src/network.cpp).
    ExaDiS visits the arms being moved in *descending* original-index order
    (`for (int k = arms.size()-1; k >= 0; k--)`), having sorted them
    ascending first, so the new node's connections are built far-index
    first. PyDiS' split_node visits nbrs_to_split in the order given, which
    Topology_Serial builds ascending (build_split_list's itertools.compress
    preserves the order of the neighbour list it draws from). The old
    node's remaining arms need no fix: both codes drop the moved arms and
    keep the rest in their existing relative order (Conn::remove_connection
    shifts left; a dict pop does the same), and both append a new
    connecting segment last, when the split is Burgers-imbalanced.

    Reordering only the new node's arms is therefore the entire fix, same
    principle as insert_node_matching_exadis for split_seg: a free function
    rather than a change to split_node, so callers that want ExaDiS'
    ordering ask for it here. Assumes nbrs_to_split already reflects the
    node's ExaDiS-seeded arm order (framework.arm_order.apply_arm_order),
    since reversing an arbitrarily-ordered list would not reproduce
    anything ExaDiS does.
    """
    split_node1, split_node2 = G.split_node(tag, pos1, pos2, nbrs_to_split, eps_b)
    set_node_arm_order(G, split_node2, list(reversed(nbrs_to_split)))
    return split_node1, split_node2
