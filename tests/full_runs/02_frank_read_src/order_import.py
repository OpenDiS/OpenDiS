"""Carrying ExaDiS' iteration orders into PyDiS, step by step.

PyDiS and ExaDiS agree on the physics of this case to round-off, but they
disagree on the order in which they visit things, and several operations are
order-dependent when a single pass performs more than one of them. Three orders
matter, and each is consulted at a different point in the step:

    arm order      per node, read by the collision rule when it picks a glide
                   plane and by a trial split's force summation
    segment order  read by the collision pass when it enumerates candidate
                   pairs, and by remesh when it snapshots segments
    node order     read by the topology pass when it enumerates multi-nodes

The orders drift because ExaDiS' SerialDisNet::remove_nodes removes a node by
moving the array's last node into the vacated slot ("Remove the nodes while not
preserving the order (faster)"), while PyDiS holds nodes in an insertion-ordered
dict where a removal leaves the others alone. They part at the first removal.
Measured on this case at step 268: only 16 of 107 slots hold the same physical
node in both codes, though the networks agree to 1.5e-07. remove_nodes also
frees tags, so the recycled tag pools drift too and there is no common labelling
to match on. The map below is therefore geometric, not by tag.

This is an instrument, not a fix. It answers "is anything other than ordering
still different?", and on this case the answer is no: with all three orders
imported the two codes agree to ~1.5e-07 through step 268, eight steps past the
first collision, where without it they part at 263. It cannot repair drift that
happens WITHIN a pass, because ExaDiS permutes its own arrays as it goes.
"""

import numpy as np

from framework.arm_order import apply_arm_order

BOX = 1000.0

# Largest offset, in b, between a PyDiS node and the ExaDiS node taken to be the
# same one. Generous on purpose: the map only has to be unambiguous, and the
# networks are typically within 1e-07. A failure to build the map is reported
# rather than silently skipped, because a silently unaligned step reads as an
# ordering-independent difference and is the easiest way to draw a wrong
# conclusion here.
MAP_TOLERANCE = 15.0


def phases():
    """phases: the three points in a step at which an order is imported

    p0 before collision, p1 after collision and before topology, p2 after
    topology and before remesh.
    """
    return ('p0', 'p1', 'p2')


def export_orders(exadis_net_data, arms_by_tag):
    """export_orders: ExaDiS' three orders, as plain arrays ready for np.savez

    exadis_net_data is get_nodes_data() plus a nodeids array; arms_by_tag is
    framework.arm_order.export_arm_order's output. The node order is implicit
    in the row order of tags and positions, and the segment order in that of
    nodeids.
    """
    tags = np.array([[int(t[0]), int(t[1])] for t in exadis_net_data['tags']])
    index = {(int(t[0]), int(t[1])): i for i, t in enumerate(tags)}
    counts, flat = [], []
    for t in tags:
        neighbours = arms_by_tag[(int(t[0]), int(t[1]))]
        counts.append(len(neighbours))
        flat.extend(index[n] for n in neighbours)
    return dict(tags=tags,
                positions=np.array(exadis_net_data['positions'], dtype=float),
                nodeids=np.array(exadis_net_data['nodeids'])[:, :2].astype(int),
                arm_counts=np.array(counts, dtype=int),
                arm_neighbours=np.array(flat, dtype=int))


def load_orders(path):
    """load_orders: read back what export_orders wrote, keyed by tag"""
    z = np.load(path)
    tags = [tuple(int(x) for x in t) for t in z['tags']]
    arms, k = {}, 0
    for i, tag in enumerate(tags):
        n = int(z['arm_counts'][i])
        arms[tag] = [tags[j] for j in z['arm_neighbours'][k:k+n]]
        k += n
    return dict(tags=tags,
                positions=z['positions'][:, :3],
                segments=[(tags[a], tags[b]) for a, b in z['nodeids']],
                arms=arms)


def _minimum_image(delta):
    return delta - BOX*np.rint(delta/BOX)


# Why a map could not be built. COUNT_MISMATCH is a different kind of event
# from the other two and worth separating: ExaDiS' merges do not delete the
# absorbed node or the connecting segment, leaving that to purge_network at the
# end of the pass, so its mid-pass network legitimately carries nodes PyDiS has
# already removed. Seeing it at a phase boundary says nothing about whether the
# two codes agree. AMBIGUOUS and TOO_FAR are the interesting failures: the
# networks are the same size but their nodes cannot be paired, which means they
# have actually diverged.
COUNT_MISMATCH = 'count'
AMBIGUOUS = 'ambiguous'
TOO_FAR = 'far'


def geometric_map(G, orders):
    """geometric_map: {ExaDiS tag: PyDiS tag}, matched by position

    Returns (mapping, note). mapping is None when no map could be built and note
    then says why, as (reason, text) with reason one of COUNT_MISMATCH,
    AMBIGUOUS or TOO_FAR. On success note is (None, text) carrying the worst
    offset, so a caller can report how close the match was. Keyed geometrically
    because the two codes' tag pools have drifted; see this module's docstring.
    """
    pydis_tags = list(G.all_nodes_tags())
    if len(pydis_tags) != len(orders['tags']):
        return None, (COUNT_MISMATCH, 'node counts differ: pydis %d, exadis %d'
                      % (len(pydis_tags), len(orders['tags'])))
    pydis_pos = np.array([G.nodes(t).R for t in pydis_tags])
    mapping, worst, claimed = {}, 0.0, set()
    for exadis_tag, exadis_pos in zip(orders['tags'], orders['positions']):
        distance = np.linalg.norm(_minimum_image(pydis_pos - exadis_pos), axis=1)
        j = int(np.argmin(distance))
        mapping[exadis_tag] = pydis_tags[j]
        claimed.add(j)
        worst = max(worst, float(distance[j]))
    if len(claimed) != len(pydis_tags):
        return None, (AMBIGUOUS, 'not one to one: %d of %d pydis nodes claimed, '
                      'worst offset %.2e'
                      % (len(claimed), len(pydis_tags), worst))
    if worst > MAP_TOLERANCE:
        return None, (TOO_FAR, 'worst offset %.2e exceeds %.1f'
                      % (worst, MAP_TOLERANCE))
    return mapping, (None, 'worst offset %.2e' % worst)


def apply_orders(G, orders):
    """apply_orders: impose ExaDiS' three orders on a PyDiS network

    Returns (segment_order, node_rank, note), or (None, None, note) if the map
    could not be built; see geometric_map for what note carries. segment_order
    is a list of PyDiS tag pairs in ExaDiS' segment order, and node_rank is
    {PyDiS tag: ExaDiS node index}; both are for the callers that consult them,
    since PyDiS holds no such order of its own to overwrite. The arm order is
    applied here, in place, because a PyDiS node does hold one.
    """
    mapping, note = geometric_map(G, orders)
    if mapping is None:
        return None, None, note
    apply_arm_order(G, {mapping[t]: [mapping[n] for n in neighbours]
                        for t, neighbours in orders['arms'].items()})
    segment_order = [(mapping[a], mapping[b]) for a, b in orders['segments']]
    node_rank = {mapping[t]: i for i, t in enumerate(orders['tags'])}
    return segment_order, node_rank, note


def segments_in_order(live_segments, wanted):
    """segments_in_order: live segments, reordered to follow wanted

    Segments named in wanted come first, in that order; any the network has but
    wanted does not are appended, keeping their own relative order, so no
    segment is ever dropped. Returns None if the two do not describe the same
    number of segments, so the caller can leave the order alone rather than
    impose a partial one.
    """
    live = [tuple(s) for s in live_segments]
    remaining = {frozenset(s): s for s in live}
    ordered = [remaining.pop(frozenset(w)) for w in wanted
               if frozenset(w) in remaining]
    ordered += list(remaining.values())
    return ordered if len(ordered) == len(live) else None
