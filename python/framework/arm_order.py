"""@package docstring
arm_order: carry a node's arm ordering from an ExaDiS network to a PyDiS one

The order in which a node's arms are visited is not physics, but it is not
free either. Both codes' collision rules walk a node's arms and keep the
first glide plane sufficiently independent of those already kept, so when
two candidate planes are near-degenerate the arm order decides which one
survives, and the collision point moves by an ulp.

That ordering does not survive a network export. ExaDisNet.export_data()
carries cell, nodes and segments and no connectivity, so a PyDiS DisNet
rebuilt from it takes its arm order from the segment array. ExaDiS' own
order is the accumulated history of every operation that touched the node:
Conn::add_connection appends and Conn::remove_connection shifts left
(core/exadis/src/network.h), starting from generate_connectivity(), which
sweeps segments. The two therefore agree on a freshly built network and
drift apart from the first topological change. Measured on
tests/unit_tests/test4_collision_mode at step 547, node 0,52 is
[0,104 0,73] in ExaDiS against [0,73 0,104] in PyDiS, segments [28, 1]
against [1, 28], and the collision point differs by 1.1369e-13 in z.

Note that this is a different thing from the segment array's own order,
which export/import does preserve.

Two functions, plus the primitive both need:

    order = export_arm_order(N)   # read what ExaDiS holds
    apply_arm_order(G, order)     # impose it on a freshly imported DisNet

That is only half of matching the two codes. Seeding the imported network
fixes the state a step starts from; a step that then splits a segment
diverges again, because PyDiS and ExaDiS update a node's arm list
differently. Seeding alone leaves 5.6843e-14 on that test. The other half is
pydis.util.arm_order.insert_node_matching_exadis, which uses
set_node_arm_order below.

Everything here works through DisNet's and the graph's own methods, so that
neither core/pydis/python/pydis/disnet.py nor core/exadis/ has to change.
Both are more awkward to modify than they look: the graph is a submodule of
a separate project, and core/exadis is a submodule too.
"""

import numpy as np


def export_arm_order(N):
    """export_arm_order: each node's arm order, as ExaDiS holds it

    Returns {node tag: [neighbour tags]}, the neighbours in ExaDiS' own conn
    order. Accepts a DisNetManager, an ExaDisNet, or the pyexadis network
    itself.

    Read it at the point the order matters. ExaDiS updates conn in place as
    it splits and merges, so an order taken after an operation describes the
    network after it.
    """
    net = _exadis_net(N)
    serial = net._get_serial_network()
    tags = [(int(row[0]), int(row[1]))
            for row in np.array(net.get_nodes_array())]
    order = {}
    for i, tag in enumerate(tags):
        conn = serial.conn(i)
        order[tag] = [tags[conn.node(j)] for j in range(conn.num)]
    return order


def apply_arm_order(G, order):
    """apply_arm_order: put a PyDiS DisNet's arms in a given order

    order is {node tag: [neighbour tags]} as export_arm_order() returns it.
    Returns the number of nodes reordered.

    A tag the network does not have is skipped, and a neighbour named for a
    node that no longer has that arm is ignored, so an order taken before a
    topological change stays usable after one.
    """
    aligned = 0
    for tag, nbr_tags in order.items():
        if G.has_node(tag):
            set_node_arm_order(G, tag, nbr_tags)
            aligned += 1
    return aligned


def set_node_arm_order(G, tag, nbr_tags):
    """set_node_arm_order: put one node's arms in the given neighbour order

    Arms whose neighbour is not named keep their relative order, after the
    named ones.

    A node holds its arms in an insertion-ordered dict, so reordering means
    re-inserting: deregister every arm, then register them again in the
    wanted order.
    """
    rank = {nbr: k for k, nbr in enumerate(nbr_tags)}
    node = G.tags_to_nodes[tag]

    def key(edge):
        other = (edge.source.tag if edge.source.tag != tag
                 else edge.target.tag)
        return rank.get(other, len(rank))

    # list() first: the loops below mutate what edges() iterates over
    ordered = sorted(list(node.edges()), key=key)
    for edge in ordered:
        node._deregister_edge(edge)
    for edge in ordered:
        node._register_edge(edge)


def _exadis_net(N):
    """_exadis_net: the pyexadis network behind a manager, a wrapper, or itself

    Imported here rather than at module scope so that this module stays
    importable without ExaDiS present, as the rest of the framework does.
    """
    if hasattr(N, 'get_disnet'):
        from pyexadis_base import ExaDisNet
        N = N.get_disnet(ExaDisNet)
    return getattr(N, 'net', N)
