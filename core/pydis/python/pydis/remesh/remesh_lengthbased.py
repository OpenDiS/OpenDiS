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


def Remesh_LengthBased(G: DisNet, params) -> None:
    """Remesh_LengthBased: coarsen below minseg, refine above maxseg

    Segments are snapshotted before coarsening, and that snapshot -- not a
    fresh scan -- is what refine works from. Matches exadis
    (core/exadis/src/remesh.h:63-64, refine_coarsen), which fixes
    nsegs = network->number_of_segs() before its single loop over both
    operations, so a segment a coarsen step lengthens by merging away a
    neighbor is never revisited for refinement within the same pass; that
    catches up on the next remesh call instead. Without this, pydis refines
    such a segment immediately, one step ahead of exadis. Confirmed as the
    cause of a real cross-code divergence in tests/full_runs/03_binary_junction
    (a node coarsened away left its neighbor's segment longer than maxseg).
    """
    all_segments_list = list(G.all_segments_tags())
    _coarsen(G, params, all_segments_list)
    _refine(G, params, all_segments_list)


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
    for tag1, tag2 in all_segments_list:
        if not (G.has_node(tag1) and G.has_node(tag2)
               and G.has_segment(tag1, tag2)):
            continue
        r1, r2 = G.nodes(tag1).R.copy(), G.nodes(tag2).R.copy()
        # apply PBC
        r2 = G.cell.closest_image(Rref=r1, R=r2)
        L = np.linalg.norm(r2-r1)
        if not (L < params.minseg):
            continue
        # RULE CHANGE (LengthBased): the second endpoint of the segment is
        # removed and the first is moved to the mid-point, where before the
        # first was removed and the second left in place. Matches exadis
        # SerialDisNet::merge_nodes in core/exadis/src/network.cpp, which
        # merges n2 into n1 at the mid-point of the two.
        if _removable(G, tag2):
            tag, survivor, R = tag2, tag1, (0.5*(r1+r2)
                                            if _removable(G, tag1) else r1)
        elif _removable(G, tag1):
            tag, survivor, R = tag1, tag2, r2
        else:
            continue
        if G.out_degree(tag) != 2:
            continue
        G.remove_two_arm_node(tag)
        # the removal can orphan the survivor and take it with it
        if G.has_node(survivor):
            G.nodes(survivor).R = R.copy()

    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed 1")


def _refine(G: DisNet, params, all_segments_list) -> None:
    """_refine: bisect segments longer than maxseg

    all_segments_list is the pre-coarsen snapshot from Remesh_LengthBased,
    not a fresh scan; see that function's docstring for why. A segment from
    that snapshot may no longer exist (one of its endpoints could have been
    the node coarsening just removed) and is skipped rather than treated as
    an error: exadis's equivalent loop simply never revisits a segment
    coarsening has already consumed either.
    """
    for tag1, tag2 in all_segments_list:
        if not (G.has_node(tag1) and G.has_node(tag2)
               and G.has_segment(tag1, tag2)):
            continue
        node1, node2 = G.nodes(tag1), G.nodes(tag2)
        r1, r2 = node1.R.copy(), node2.R.copy()
        # apply PBC
        r2 = G.cell.closest_image(Rref=r1, R=r2)
        L = np.linalg.norm(r2-r1)
        if (L > params.maxseg) and ((node1.constraint != DisNode.Constraints.PINNED_NODE) or (node2.constraint != DisNode.Constraints.PINNED_NODE)):
            # RULE CHANGE (LengthBased): recycle=False, so a node inserted
            # here cannot take a tag freed by the coarsening above. exadis
            # frees tags only in SerialDisNet::purge_network, once the whole
            # pass is over, so a tag freed during a pass is not available
            # within it.
            new_tag = G.get_new_tag(recycle=False)
            r = (r1 + r2)/2.0
            G.insert_node_between(tag1, tag2, new_tag, r)
            if not G.is_sane():
                raise ValueError("Remesh_LengthBased: sanity check failed 1a")

    if not G.is_sane():
        raise ValueError("Remesh_LengthBased: sanity check failed 2")
