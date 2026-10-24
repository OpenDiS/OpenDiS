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


def export_data_with_arms_first(G, tags):
    """export_data_with_arms_first: G.export_data(), reordered for ExaDiS'
    generate_connectivity() to reproduce named tags' own arm order

    The companion of export_arm_order()/apply_arm_order(), aimed the other
    way: feeding an ExaDisNet from PyDiS, not reading one. ExaDiS'
    generate_connectivity() (network.cpp), run inside
    SerialDisNet::import_data() and the ExaDisNet(cell, nodes, segs)
    constructor alike, builds each node's conn by scanning the segment
    array once and appending a connection to both endpoints wherever they
    occur -- so a node's rebuilt conn order is exactly the order its own
    segments occur in that array, not any node's own view of it.
    G.export_data()'s own segment order is the DisNet graph's global
    edge-insertion history, which need not agree with any one node's own
    per-node order (the same drift export_arm_order()'s docstring
    describes, met from the export side instead of the import side).

    Places each tag in tags' own arms (G.neighbors_tags(tag)'s order) at
    the front of the segment array, processing tags in the order given and
    each segment on its first mention, so a tag processed before its arms
    are otherwise claimed sees generate_connectivity() reproduce its exact
    order. A segment two named tags both claim (the connecting segment a
    topology trial's two candidate nodes share) goes to whichever is
    processed first; callers that need both should list the one whose
    order matters more first, or avoid relying on that one shared slot's
    order at all.

    Only the segment array's row order changes; nodeids/burgers/planes
    move together, and cell/nodes are untouched.
    """
    data = G.export_data()
    tags_list = [tuple(int(v) for v in t) for t in data['nodes']['tags']]
    tag_to_idx = {t: i for i, t in enumerate(tags_list)}
    nodeids = data['segs']['nodeids']

    from collections import deque
    by_pair = {}
    for k, (a, b) in enumerate(nodeids):
        by_pair.setdefault(frozenset((int(a), int(b))), deque()).append(k)

    order, placed = [], set()
    for tag in tags:
        tag_idx = tag_to_idx.get(tag)
        if tag_idx is None:
            continue
        for nbr in G.neighbors_tags(tag):
            nbr_idx = tag_to_idx.get(nbr)
            if nbr_idx is None:
                continue
            q = by_pair.get(frozenset((tag_idx, nbr_idx)))
            if not q:
                continue
            for _ in range(len(q)):
                k = q.popleft()
                if k not in placed:
                    order.append(k)
                    placed.add(k)
                    break
                q.append(k)
    order.extend(k for k in range(len(nodeids)) if k not in placed)

    data['segs'] = {key: np.asarray(val)[order] for key, val in data['segs'].items()}
    return data


def rebase_export_data(G, before):
    """rebase_export_data: G.export_data(), with every segment before also
    has kept at before's own original array row

    export_data_with_arms_first fixes which of a node's own arms sees
    generate_connectivity() first, and that turned out not to be the whole
    story: measured directly on test5_topology_mode step 194, comparing
    PyDiS' resulting split direction against ExaDiS' own true one (backed
    out of its recorded final position), export_data_with_arms_first left a
    2.16e-14 relative gap despite every node's own arm order already being
    right. The remaining gap is the *other* nodes' arms: ExaDiS'
    generate_connectivity() rebuilds every node's conn by a single pass
    over the segment array, so a segment two nodes both still hold moving
    position -- which export_data_with_arms_first's per-tag placement does
    to segments it has not yet claimed -- also reorders whichever of those
    two nodes was not the one being placed. A whole-network pair sum (the
    force this feeds) sums over all of that, not just one node's arms, so
    it is sensitive to the array's global order, not only to each node's
    local view of it.

    before is an ExaDisNet.export_data() snapshot of the network as it
    stood before whatever produced G (a topology trial's split, here).
    Segments unmodified by that change keep before's own row position
    exactly; only the ones actually touched (an endpoint tag rewired to a
    new node) are updated in place, matched to their original row by the
    endpoint tag that did not move, restricted to candidates whose other
    endpoint is a tag before did not have at all -- a split only ever
    rewires an arm onto a node it just created, never onto one that already
    existed, so that restriction is exact, unlike matching on the Burgers
    vector alone (tried first; rejected because a crystal's segments reuse
    a small set of glide directions, so two unrelated segments sharing a
    node can carry the identical vector, and the wrong one matched first in
    a direct test on step 194's real trial). This is exadis_own_margin's
    own construction technique (test_topology_mode_pydis_exadis.py),
    generalized: that function hand-writes the same reasoning for one
    specific trial; this one derives it from a before/after diff so any
    caller can use it. Segments new since before (a topology trial's
    connecting segment, with both endpoints new) are appended at the end,
    in G's own order, since before has no position for them to preserve.

    Matching a row is not enough on its own: DisNet stores each edge with
    a source/target that need not agree with which end before called n1,
    so a segment picked up by tag alone can come back with its two ends
    swapped and its Burgers vector negated to match -- physically the same
    segment, but not the bit-identical row before had, since a sum that
    later walks the array in n1/n2 order is not guaranteed invariant under
    swapping both together (found directly on step 194's real trial: one
    of the two rewired segments came back swapped, and that alone fully
    accounted for the gap between this function's first version and
    exadis_own_margin's hand-built construction). Every matched row is
    therefore re-oriented to start from before's own n1 tag (or, for a
    rewired row, whichever of its two tags before also used), negating the
    Burgers vector to match; the plane normal is a property of the glide
    plane, not of n1 vs n2, and is left alone. Confirmed directly against
    exadis_own_margin's construction for step 194's real trial, after this
    fix: bit-identical result.
    """
    data = G.export_data()
    before_tags = [tuple(int(v) for v in t) for t in before['nodes']['tags']]
    cur_tags = [tuple(int(v) for v in t) for t in data['nodes']['tags']]
    new_tags = set(cur_tags) - set(before_tags)

    before_nodeids = np.asarray(before['segs']['nodeids'])
    cur_nodeids = np.asarray(data['segs']['nodeids'])

    cur_by_pair = {}
    for k, (a, b) in enumerate(cur_nodeids):
        pair = frozenset((cur_tags[int(a)], cur_tags[int(b)]))
        cur_by_pair.setdefault(pair, []).append(k)

    matched = set()
    order = []
    flip = {}
    for a, b in before_nodeids:
        ta, tb = before_tags[int(a)], before_tags[int(b)]
        j = next((c for c in cur_by_pair.get(frozenset((ta, tb)), ())
                  if c not in matched), None)
        if j is None:
            # One endpoint was rewired to a new node; find the row that
            # took its place by the endpoint that did not move, requiring
            # the other endpoint to be a tag before had no node for at all.
            for anchor in (ta, tb):
                for c, (aa, bb) in enumerate(cur_nodeids):
                    if c in matched:
                        continue
                    pair = (cur_tags[int(aa)], cur_tags[int(bb)])
                    if anchor in pair and any(t in new_tags for t in pair):
                        j = c
                        break
                if j is not None:
                    break
        if j is not None:
            # Orient the matched row to start from before's own n1 tag
            # (ta), or, if ta itself was the one rewired away, from tb --
            # the endpoint before actually put first either way.
            cj, cjb = cur_tags[int(cur_nodeids[j][0])], cur_tags[int(cur_nodeids[j][1])]
            want_first = ta if ta in (cj, cjb) else tb
            flip[j] = (cj != want_first)
            order.append(j)
            matched.add(j)
    order.extend(k for k in range(len(cur_nodeids)) if k not in matched)

    segs = {key: np.array(val)[order] for key, val in data['segs'].items()}
    for row, j in enumerate(order):
        if flip.get(j):
            segs['nodeids'][row] = segs['nodeids'][row][::-1]
            segs['burgers'][row] = -segs['burgers'][row]
    data['segs'] = segs
    return data


def _exadis_net(N):
    """_exadis_net: the pyexadis network behind a manager, a wrapper, or itself

    Imported here rather than at module scope so that this module stays
    importable without ExaDiS present, as the rest of the framework does.
    """
    if hasattr(N, 'get_disnet'):
        from pyexadis_base import ExaDisNet
        N = N.get_disnet(ExaDisNet)
    return getattr(N, 'net', N)
