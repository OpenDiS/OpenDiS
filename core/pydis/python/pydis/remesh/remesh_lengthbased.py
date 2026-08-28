"""@package docstring
Remesh_LengthBased: remesh by segment length

Coarsen segments shorter than minseg, refine segments longer than maxseg.
Matches the exadis rule of the same name, step for step, on a network with no
junctions.
"""

import numpy as np
from ..disnet import DisNet, DisNode


def _removable(G: DisNet, tag) -> bool:
    """_removable: whether coarsening may delete this node

    Mirrors exadis SerialDisNet::constrained_node: a node is left alone unless
    it has exactly two arms and is not pinned.
    """
    return (G.out_degree(tag) == 2
            and G.nodes(tag).constraint != DisNode.Constraints.PINNED_NODE)


def _movable(G: DisNet, tag, merged_already) -> bool:
    """_movable: whether coarsening may delete this node or move it

    _removable, and additionally not already the survivor of a merge earlier in
    this same pass. ExaDiS' merge_nodes_position (core/exadis/src/network.cpp)
    does not delete the absorbed node or the segment joining the two: for that
    self-connection it only zeroes the segment's Burgers vector and leaves the
    connection in place, with the removal deferred to a later purge. The
    survivor therefore carries a dead zero-Burgers arm for the rest of the
    pass, so its arm count reads 3 rather than 2 and constrained_node() calls
    it constrained.

    SIDE EFFECT, reproduced here deliberately: that freezes the node for the
    remainder of the pass. A further merge on one of its segments keeps it
    exactly where it is instead of moving it to the new mid-point, and a
    segment whose BOTH endpoints have already been merged is not coarsened at
    all, since merge_nodes returns without doing anything when both nodes are
    constrained. Any check that skipped arms with a zero Burgers vector when
    counting a node's connections would remove the freezing.

    Measured on tests/full_runs/02_frank_read_src at step 266, coarsening the
    chain A-B-C where B has both arms below minseg: ExaDiS merges B into C at
    their mid-point and then absorbs A without moving that survivor, landing
    at midpoint(B,C), while an unfrozen mid-point rule lands at
    midpoint(A, midpoint(B,C)), 4.69 away.
    """
    return (_removable(G, tag) and tag not in merged_already)


def Remesh_LengthBased(G: DisNet, params) -> None:
    """Remesh_LengthBased: coarsen below minseg, refine above maxseg

    Two shapes, chosen by params.interleave_coarsen_refine.

    True (default), matching exadis: ONE pass over the segment snapshot, each
    segment bisected or merged according to its length at the moment it is
    reached. This is the shape of exadis' refine_coarsen
    (core/exadis/src/remesh.h), and it is what makes a segment lengthened by an
    earlier merge in the same pass get bisected, with no rule saying so.

    False, matching ParaDiS: coarsening over the whole snapshot, then
    refinement, as ParaDiS' RemeshRule_2 does with MeshCoarsen(home) then
    MeshRefine(home) (external/paradis/src/RemeshRule_2.c). Refinement then
    takes its candidates per refine_from_fresh_scan, which only applies here.

    The two references genuinely differ, so pydis has to pick one, and the pick
    is measurable. On tests/full_runs/02_frank_read_src, with exadis' iteration
    orders imported so that ordering is not a variable, interleaving tracks
    exadis to round-off for the whole run (~1e-07, node counts equal at every
    step, 24/24 nodes and 2.8483e-08 at step 300). Two passes part from it at
    step 269: exadis bisects a 48.24 long segment (maxseg 40) that coarsening
    had just lengthened, and no choice of refinement candidates reproduces
    that, because at step 268 exadis leaves an equivalent 42.12 one alone. The
    difference is the pass structure, not the refinement rule.

    Lengths and merge points are read live in both shapes, at the moment each
    segment is visited. Only the candidate list is fixed in advance.

    Both shapes above are coarsen_mode 0. coarsen_mode 1 is a different
    algorithm and takes its own path: refine every over-long segment, then
    coarsen node by node. See _coarsen_nodes.
    """
    if getattr(params, 'coarsen_mode', 0) == 1:
        # recycle=True: this refinement runs before any coarsening, so the
        # tag pool holds only what earlier passes freed, which exadis reuses
        _refine(G, params, list(G.all_segments_tags()), recycle=True)
        _coarsen_nodes(G, params)
        return

    all_segments_list = list(G.all_segments_tags())
    if getattr(params, 'interleave_coarsen_refine', True):
        _interleave(G, params, all_segments_list)
        return
    _coarsen(G, params, all_segments_list)
    if getattr(params, 'refine_from_fresh_scan', False):
        all_segments_list = list(G.all_segments_tags())
    _refine(G, params, all_segments_list)


def _follow_merges(merged_into, tag):
    """_follow_merges: where a tag ended up after earlier merges in this pass

    merged_into maps a removed node's tag to the tag that absorbed it. A
    survivor can itself be absorbed later, so the chain is followed to its end.
    """
    seen = set()
    while tag in merged_into and tag not in seen:
        seen.add(tag)
        tag = merged_into[tag]
    return tag


def _live_segment(G: DisNet, merged_into, tag1, tag2):
    """_live_segment: the segment a snapshot entry now names, or None

    Follows both endpoints through the merges made so far in this pass and
    checks the segment still exists. Returns (tag1, tag2, r1, r2, length) with
    r2 taken as the image of tag2 nearest tag1, all read live.
    """
    tag1 = _follow_merges(merged_into, tag1)
    tag2 = _follow_merges(merged_into, tag2)
    if tag1 == tag2:
        return None
    if not (G.has_node(tag1) and G.has_node(tag2) and G.has_segment(tag1, tag2)):
        return None
    r1 = G.nodes(tag1).R.copy()
    r2 = G.cell.closest_image(Rref=r1, R=G.nodes(tag2).R.copy())
    return tag1, tag2, r1, r2, float(np.linalg.norm(r2 - r1))


def _coarsen_segment(G: DisNet, tag1, tag2, r1, r2,
                     merged_into, merged_already) -> None:
    """_coarsen_segment: merge one under-length segment's endpoints

    RULE CHANGE (LengthBased): the second endpoint is removed and the first
    moved to the mid-point, where before the first was removed and the second
    left in place. Matches exadis SerialDisNet::merge_nodes
    (core/exadis/src/network.cpp), which merges n2 into n1 at the mid-point.

    _movable rather than _removable, so that a node already merged once in this
    pass keeps its position instead of moving again; see _movable for the exadis
    behavior this reproduces and its side effect.
    """
    if _movable(G, tag2, merged_already):
        tag, survivor, R = tag2, tag1, (0.5*(r1+r2)
                                        if _movable(G, tag1, merged_already)
                                        else r1)
    elif _movable(G, tag1, merged_already):
        tag, survivor, R = tag1, tag2, r2
    else:
        return
    if G.out_degree(tag) != 2:
        return
    G.remove_two_arm_node(tag)
    merged_into[tag] = survivor
    merged_already.add(survivor)
    # the removal can orphan the survivor and take it with it
    if G.has_node(survivor):
        G.nodes(survivor).R = R.copy()


def _refine_segment(G: DisNet, tag1, tag2, r1, r2, recycle=False) -> None:
    """_refine_segment: bisect one over-length segment

    A segment between two pinned nodes is left alone.

    RULE CHANGE (LengthBased): recycle=False by default, so a node inserted
    here cannot take a tag freed by coarsening earlier in the same pass. exadis
    frees tags only in SerialDisNet::purge_network, once the whole pass is
    over, so a tag freed during a pass is not available within it.

    recycle=True is for a refinement pass that runs before any coarsening, as
    coarsen_mode 1's does. Nothing has been freed within the pass yet there, so
    the pool holds only tags freed by earlier passes, which exadis' split_seg
    does draw on. Refusing them would be wrong rather than conservative:
    measured on tests/unit_tests/test3_remesh_rule at steps 302 and 380, where
    exadis reused a tag and pydis took maxindex+1.
    """
    if (G.nodes(tag1).constraint == DisNode.Constraints.PINNED_NODE and
            G.nodes(tag2).constraint == DisNode.Constraints.PINNED_NODE):
        return
    G.insert_node_between(tag1, tag2, G.get_new_tag(recycle=recycle),
                          (r1 + r2)/2.0)
    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed after a bisection")


def _interleave(G: DisNet, params, all_segments_list) -> None:
    """_interleave: one pass, each segment refined or coarsened as it is reached

    exadis' loop shape, and the reason it matters. ExaDiS' refine_coarsen
    (core/exadis/src/remesh.h) is a single loop over segment indices that tests
    each segment's length at the moment it reaches it and either bisects or
    merges. So a segment that an earlier merge in the same pass lengthened is
    bisected when the cursor gets to it, with no rule saying so; and one
    lengthened by a later merge is not. Running coarsening and refinement as two
    separate passes cannot express either case, which is what
    tests/full_runs/02_frank_read_src showed: at step 268 exadis leaves a 42.12
    long segment (maxseg 40) that coarsening had just lengthened, and at step
    269 it bisects a 48.24 one, in otherwise identical circumstances.

    ParaDiS is two passes, MeshCoarsen(home) then MeshRefine(home)
    (external/paradis/src/RemeshRule_2.c), so this is where the two references
    part company and pydis has to pick one. See Remesh_LengthBased.
    """
    merged_into = {}
    merged_already = set()
    for tag1, tag2 in all_segments_list:
        live = _live_segment(G, merged_into, tag1, tag2)
        if live is None:
            continue
        tag1, tag2, r1, r2, length = live
        if length > params.maxseg:
            _refine_segment(G, tag1, tag2, r1, r2)
        elif length < params.minseg:
            _coarsen_segment(G, tag1, tag2, r1, r2, merged_into, merged_already)
    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed after interleave")


def _coarsen_nodes(G: DisNet, params) -> None:
    """_coarsen_nodes: coarsen_mode 1, node-centric coarsening

    Walks 2-arm unconstrained nodes rather than segments. A node whose shorter
    arm is under minseg is merged into its NEARER neighbour, and that neighbour
    stays where it is. Mirrors the coarsen_mode == 1 block of exadis'
    refine_coarsen (core/exadis/src/remesh.h), which is also the shape of
    ParaDiS' MeshCoarsen (external/paradis/src/RemeshRule_2.c), and is exadis'
    python-side default.

    Three ways it differs from mode 0, all of them consequences of walking
    nodes instead of segments:

    - NOTHING MOVES. The survivor is placed at r0, its own position taken as
      the image nearest the node being removed, so mode 1 only ever deletes
      nodes. Mode 0 puts the survivor at the mid-point of the merged segment.
    - The test is on the node, so a node needs only ONE arm under minseg to go;
      mode 0 needs the one segment it is looking at to be under minseg.
    - Constrained means `constraint != UNCONSTRAINED`, which also excludes
      surface and corner nodes; mode 0's constrained_node only excludes pinned
      ones.

    Two exadis details reproduced deliberately:

    - `frozen` mirrors the dead-arm side effect documented in _movable. A
      survivor keeps a zeroed-Burgers arm to the node it absorbed until the
      end-of-pass purge, so its arm count reads 3 and the `conn[i].num != 2`
      test skips it for the rest of the pass.
    - exadis also skips a candidate whose neighbour has no connections left.
      That cannot arise here: pydis deletes an absorbed node outright, so it
      cannot still be someone's neighbour.

    Neighbour order does not matter, unlike elsewhere in this file: it only
    decides which of l0, l1 is which, and the choice is made on the lengths.
    Only an exact tie would depend on the order.
    """
    frozen = set()
    for tag in list(G.all_nodes_tags()):
        if not G.has_node(tag) or tag in frozen:
            continue
        if G.out_degree(tag) != 2:
            continue
        if G.nodes(tag).constraint != DisNode.Constraints.UNCONSTRAINED:
            continue

        ri = G.nodes(tag).R.copy()
        nbrs = list(G.neighbors_tags(tag))
        r = [G.cell.closest_image(Rref=ri, R=G.nodes(n).R.copy()) for n in nbrs]
        length = [float(np.linalg.norm(ri_n - ri)) for ri_n in r]
        if min(length) > params.minseg:
            continue

        near = 0 if length[0] < length[1] else 1
        survivor, R = nbrs[near], r[near]
        G.remove_two_arm_node(tag)
        frozen.add(survivor)
        # the removal can orphan the survivor and take it with it
        if G.has_node(survivor):
            G.nodes(survivor).R = R.copy()

    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed after "
                         "node-centric coarsening")


def _coarsen(G: DisNet, params, all_segments_list) -> None:
    """_coarsen: merge the endpoints of segments shorter than minseg

    all_segments_list fixes *which* segments are candidates (the pre-pass
    snapshot Remesh_LengthBased captured), but each one's length and merge
    point are read live, at the moment that segment is visited, not from a
    frozen snapshot. Matches exadis SerialDisNet::refine_coarsen
    (core/exadis/src/remesh.h:64-137): its single loop re-reads
    network->nodes[n].pos fresh on every iteration, so when two coarsening
    segments share an endpoint, the one visited second sees the first
    merge's result. Collecting every candidate from a frozen snapshot and
    applying them afterward -- the previous approach here -- computes both
    merge points from pre-coarsen positions instead, which disagrees with
    exadis whenever that sharing happens. Confirmed as a real cross-code
    divergence in tests/full_runs/03_binary_junction: two segments on
    either side of the same node were both under minseg at once.
    """
    # RULE CHANGE, to match ParaDiS: a snapshot entry naming a node that an
    # earlier merge in this same pass removed is followed to the node that
    # absorbed it, instead of being skipped. ParaDiS' MeshCoarsen
    # (RemeshRule_2.c) walks node keys by index over the live array and
    # re-reads each node's neighbours, so a merge earlier in the pass never
    # stops a later one; ExaDiS' refine_coarsen does the same walking segments
    # by index, whose endpoints a merge rewires to the survivor. Only this
    # tag-keyed snapshot lost the entry, and the segment it named survives the
    # merge reconnected to the survivor and is often still under minseg, so a
    # coarsening cascade stopped partway. Measured consequence: a closed
    # four-node loop of perimeter ~60 with two arms under minseg collapsed
    # completely in ParaDiS and ExaDiS and was left as a three-node loop here.
    #
    # Old version:
    #   for tag1, tag2 in all_segments_list:
    #       if not (G.has_node(tag1) and G.has_node(tag2)
    #              and G.has_segment(tag1, tag2)):
    #           continue
    merged_into = {}
    merged_already = set()
    for tag1, tag2 in all_segments_list:
        live = _live_segment(G, merged_into, tag1, tag2)
        if live is None:
            continue
        tag1, tag2, r1, r2, length = live
        if length < params.minseg:
            _coarsen_segment(G, tag1, tag2, r1, r2, merged_into, merged_already)

    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed 1")


def _refine(G: DisNet, params, all_segments_list, recycle=False) -> None:
    """_refine: bisect segments longer than maxseg

    all_segments_list is the pre-coarsen snapshot from Remesh_LengthBased,
    not a fresh scan; see that function's docstring for why. A segment from
    that snapshot may no longer exist (one of its endpoints could have been
    the node coarsening just removed) and is skipped rather than treated as
    an error: exadis's equivalent loop simply never revisits a segment
    coarsening has already consumed either.
    """
    for tag1, tag2 in all_segments_list:
        live = _live_segment(G, {}, tag1, tag2)
        if live is None:
            continue
        tag1, tag2, r1, r2, length = live
        if length > params.maxseg:
            _refine_segment(G, tag1, tag2, r1, r2, recycle=recycle)

    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed 2")
