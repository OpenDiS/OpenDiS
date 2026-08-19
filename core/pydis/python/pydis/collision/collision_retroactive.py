"""Retroactive collision detection and handling.

Follows ParaDiS RetroactiveCollisions2 in
external/ParaDiS.git/src/RetroactiveCollision2.c, which is ParaDiS
collisionMethod 4 and the algorithm ExaDiS CollisionRetroactive ports in
core/exadis/src/collision_types/collision_retroactive.cpp.

Note: this is NOT the algorithm of Sills and Cai (2014). That paper
describes RetroactiveCollisions in RetroactiveCollision.c, ParaDiS
collisionMethod 3, which detects collisions by interval bisection and
handles the earliest one first. Method 4 instead tests whether two swept
segments come within rann during an interval, and has no hit time.

Two passes in a fixed order: segment pairs, then hinges. Nodes are
addressed by tag throughout; no array index outlives the candidate
generator that produced it.

Where ExaDiS and ParaDiS differ, this follows ExaDiS, because ExaDiS is
what it is compared against step by step. Each such place says so in a
comment naming both behaviours. The hinge pass and the collision latch
are one; merging two constrained nodes is another; when merged nodes and
annihilated segments are removed is a third, and the only one carrying a
switch, the purge_at_end argument of RetroactiveCollision.
"""

from collections import Counter

import numpy as np

from ..disnet import DisNet, DisNode
from ..util.glide_planes import PlaneSet, constrained_plane_point
from ..util.state_access import velocities_by_tag, old_positions_by_tag
from .swept_distance import swept_seg_seg_collision, hinge_cos_angle

# mindist for the predictive interval: near zero, so that the second pass
# fires only on segments that will genuinely cross rather than approach
PREDICTIVE_MINDIST = 1.0e-6

# hinge thresholds on the cosine between two arms of a node, relaxed when
# the two far nodes are themselves connected (a triangle)
HINGE_TOL = 0.98
HINGE_TOL_TRI = 0.9

# a node is thrown too far if the adjusted collision point is further than
# this multiple of rann from it
THROW_FACTOR = 4.0

ZERO_LENGTH2 = 1.0e-20


class CollisionRecord:
    """CollisionRecord: what the rule decided, for tests and diagnostics

    Off by default: the production path should not pay for a diagnostic
    only the comparison uses.
    """

    def __init__(self):
        self.events = []
        self.reasons = Counter()

    def add(self, kind, interval, outcome, **details):
        self.reasons["%s:%s:%s" % (kind, interval, outcome)] += 1
        self.events.append(dict(kind=kind, interval=interval,
                                outcome=outcome, **details))

    def summary(self):
        return dict(self.reasons)


class RetroactiveCollision:
    """RetroactiveCollision: the rule, as a small orchestration of helpers

    Holds only what the pass needs: the parameters, and the per-pass
    latch. Constructed per call rather than kept on the Collision object,
    so nothing leaks between steps.

    purge_at_end defers removing merged nodes and annihilated segments to
    a single purge_network at the end of the pass, which is what exadis
    does, so the hinge pass sees the connectivity a merge leaves behind.
    False deletes them at the merge, which is what ParaDiS does. The two
    give different results; see DisNet.merge_node.
    """

    def __init__(self, G, state, record=None, purge_at_end=True):
        self.G = G
        self.rann = state["rann"]
        self.mindist2 = self.rann ** 2
        self.dt = state.get("dt")
        self.use_glide_planes = bool(state.get("use_glide_planes", False))
        # exadis defers node and segment removal to the end of the pass, so
        # its later passes see annihilated arms as live connectivity. False
        # restores the ParaDiS behaviour of deleting them at the merge.
        self.purge_at_end = bool(purge_at_end)
        self.record = record

        self.velocities = velocities_by_tag(state)
        self.old_positions = old_positions_by_tag(state)
        self.missing_velocity = 0
        self.missing_old_position = 0
        self.last_planes = 0

        # segments already involved in a collision this pass, as frozensets
        # of their two tags. Per-pass state, so it is local and not a flag
        # on the network.
        self.done = set()

        # nodes eligible for the hinge pass, filled in by run()
        self.pre_pass_nodes = frozenset()

    # ------------------------------------------------------------ lookups

    def velocity(self, tag):
        """velocity: nodal velocity, zero if the mobility law never saw it"""
        v = self.velocities.get(tag) if self.velocities else None
        if v is None:
            self.missing_velocity += 1
            return np.zeros(3)
        return v

    def old_position(self, tag):
        """old_position: where a node was at the start of the step

        A node absent from the set is one created after the snapshot was
        taken, so it has not moved and its current position is its old
        one. Counted, because with the standard driver it should not
        happen and a non-zero count means an invariant changed.
        """
        r = self.old_positions.get(tag)
        if r is None:
            self.missing_old_position += 1
            return self.G.nodes(tag).R
        return r

    def swept(self, tag, ref):
        """swept: (start, end) positions of a node over the past interval

        Both reduced to the periodic image nearest ref, so that a segment
        pair is compared in one consistent frame.
        """
        now = self.G.cell.closest_image(Rref=ref, R=self.G.nodes(tag).R)
        was = self.G.cell.closest_image(Rref=now, R=self.old_position(tag))
        return was, now

    def predicted(self, tag, now):
        """predicted: where a node will be after another step at its speed"""
        return now + self.dt * self.velocity(tag)

    # -------------------------------------------------------- the criterion

    def detect(self, seg1, seg2):
        """detect: does this pair of segments collide, over either interval

        Returns (interval, L1, L2) with interval 'retroactive' or
        'predictive', or None. The retroactive interval is tried first and
        the predictive one only if it stays silent, which is the order the
        C uses and which matters: the two can both fire and report
        different points along the segments.
        """
        (t1, t2), (t3, t4) = seg1, seg2
        p1, p2 = self.G.nodes(t1).R, None
        p1_old, p1_now = self.swept(t1, p1)
        p2_old, p2_now = self.swept(t2, p1_now)
        p3_old, p3_now = self.swept(t3, p1_now)
        p4_old, p4_now = self.swept(t4, p3_now)

        if (np.dot(p1_now - p2_now, p1_now - p2_now) < ZERO_LENGTH2
                or np.dot(p3_now - p4_now, p3_now - p4_now) < ZERO_LENGTH2):
            return None

        hit, _, l1, l2 = swept_seg_seg_collision(
            self.rann, p1_old, p1_now, p2_old, p2_now,
            p3_old, p3_now, p4_old, p4_now)
        if hit:
            return 'retroactive', l1, l2

        if not self.dt:
            return None
        n1, n2 = self.predicted(t1, p1_now), self.predicted(t2, p2_now)
        n3, n4 = self.predicted(t3, p3_now), self.predicted(t4, p4_now)
        hit, _, l1, l2 = swept_seg_seg_collision(
            PREDICTIVE_MINDIST, p1_now, n1, p2_now, n2, p3_now, n3, p4_now, n4)
        if hit:
            return 'predictive', l1, l2
        return None

    # --------------------------------------------------------- the operators

    def merge_node_for_segment(self, tag1, tag2, ratio):
        """merge_node_for_segment: which node this segment contributes

        Reuses an endpoint when the collision point falls within rann of
        it, and otherwise splits the segment there and returns the new
        node. Returns (tag, velocity) or (None, None) if the segment has
        gone.

        Called once per segment of a pair, so four times per collision in
        the worst case, which is why it is a function and not inline.
        """
        if not self.G.has_segment(tag1, tag2):
            return None, None

        r1 = self.G.nodes(tag1).R
        vec = self.G.seg_vector(tag1, tag2)
        length2 = float(np.dot(vec, vec))

        if ratio * ratio * length2 < self.mindist2:
            return tag1, self.velocity(tag1)
        if (1.0 - ratio) ** 2 * length2 < self.mindist2:
            return tag2, self.velocity(tag2)

        new_tag = self.G.get_new_tag()
        # folded into the primary cell, as ParaDiS folds every split and
        # merge position with FoldBox. Without it a node placed across a
        # periodic face keeps an out-of-cell coordinate: the same point,
        # written differently, which compares unequal.
        pnew = self.G.cell.fold(r1 + ratio * vec)
        self.G.insert_node_between(tag1, tag2, new_tag, pnew)
        v = ((1.0 - ratio) * self.velocity(tag1)
             + ratio * self.velocity(tag2))
        return new_tag, v

    def arm_planes(self, tag, velocity):
        """arm_planes: (normal, offset) candidates from one node's arms

        One per arm, in the C's order of preference: the segment's own
        glide plane when glide planes are in use, else the plane spanned by
        the line and the velocity, else by the line and the Burgers vector.
        An arm that defines none is skipped.
        """
        r = self.G.nodes(tag).R
        for _, vec, burg, plane in self.G.arm_vectors(tag):
            length = float(np.linalg.norm(vec))
            if length < 1.0e-20:
                continue
            direction = -vec / length          # from neighbour to this node
            options = []
            if self.use_glide_planes and plane is not None:
                options.append(np.asarray(plane, dtype=float))
            options.append(np.cross(direction, velocity))
            options.append(np.cross(direction, burg))
            for normal in options:
                norm2 = float(np.dot(normal, normal))
                if norm2 > 1.0e-12:
                    unit = normal / np.sqrt(norm2)
                    yield unit, float(np.dot(unit, r))
                    break

    def constrained_node(self, tag):
        """constrained_node: is this node ineligible to be relocated

        True for a node that does not have exactly two arms, or that
        carries any constraint. A junction or a pinned end is held in
        place; an ordinary interior node is free to move.
        """
        return (self.G.out_degree(tag) != 2
                or self.G.nodes(tag).constraint
                != DisNode.Constraints.UNCONSTRAINED)

    def collision_point(self, tag1, v1, tag2, v2):
        """collision_point: where the two merging nodes should meet

        The base point, moved onto the intersection of the glide planes of
        both nodes' arms. A node with a single arm cannot be relocated, so
        its position is used unchanged.
        """
        r1 = self.G.nodes(tag1).R
        r2 = self.G.cell.closest_image(Rref=r1, R=self.G.nodes(tag2).R)
        if self.G.out_degree(tag1) == 1:
            return r1
        if self.G.out_degree(tag2) == 1:
            return r2

        # Change new position to be one of the node position so that
        # constrained node does not need to move
        c1, c2 = self.constrained_node(tag1), self.constrained_node(tag2)
        if c1 and not c2:
            final_pos = r1
        elif c2 and not c1:
            final_pos = r2
        else:
            final_pos = 0.5 * (r1 + r2)
        # one accumulator across both nodes, so that a plane offered by
        # the second is tested against those already taken from the first.
        # Two contradictory parallel planes would otherwise both be kept
        # and the solve would go singular, leaving the point unmoved.
        planes = PlaneSet()
        planes.extend(self.arm_planes(tag1, v1))
        planes.extend(self.arm_planes(tag2, v2), allow_first=False)
        self.last_planes = planes.n
        return self.G.cell.fold(
            constrained_plane_point(final_pos, planes.normals, planes.offsets))

    def thrown_too_far(self, tag1, tag2, newpos):
        """thrown_too_far: would this collision fling both nodes

        The guard is on both nodes together, as in ExaDiS: a collision
        that moves one node a long way is allowed provided the other stays
        put, which is the pinned case.
        """
        limit = THROW_FACTOR ** 2 * self.mindist2
        d1 = self.G.nodes(tag1).R - self.G.cell.closest_image(
            Rref=self.G.nodes(tag1).R, R=newpos)
        d2 = self.G.nodes(tag2).R - self.G.cell.closest_image(
            Rref=self.G.nodes(tag2).R, R=newpos)
        return float(np.dot(d1, d1)) > limit and float(np.dot(d2, d2)) > limit

    def merge(self, tag1, tag2, newpos):
        """merge: collapse the two nodes at newpos

        The first node survives unless the second is constrained, in which
        case the two exchange roles so that the constrained one is kept.
        This is the reference's rule as written, which holds here because
        the pairs are now ordered the same way (see segment_pairs). Note
        DisNet.merge_node deletes its *first* argument, the opposite
        convention, hence the call order below.

        The merge goes ahead even when both nodes are constrained, to match
        exadis, whose merge never inspects constraints. DIFFERENT FROM
        ParaDiS, which unpins one of them first so that its own merge,
        which refuses to delete a pinned node, can proceed:

            if (mergenode1->constraint == PINNED_NODE &&
                mergenode2->constraint == PINNED_NODE)
                mergenode1->constraint &= ~PINNED_NODE;

        Unpinning silently changes a boundary condition the caller set, so
        the exadis behaviour is taken instead: keep both constraints and
        merge anyway. Refusing outright, which is what PyDiS did before, is
        the one option neither code takes: it drops a collision that should
        have formed a junction, and says nothing.
        """
        if self.G.nodes(tag2).constraint != DisNode.Constraints.UNCONSTRAINED:
            survivor, dead = tag2, tag1
        else:
            survivor, dead = tag1, tag2
        return self.G.merge_node(dead, survivor, position=newpos,
                                 ignore_constraints=True,
                                 defer_purge=self.purge_at_end)

    def latch(self, *tags):
        """latch: bar every segment on these nodes from colliding again"""
        for tag in tags:
            if self.G.has_node(tag):
                for nbr in self.G.neighbors_tags(tag):
                    self.done.add(frozenset((tag, nbr)))

    # ------------------------------------------------------------- the passes

    def collide(self, seg1, seg2, interval, l1, l2, kind='segment'):
        """collide: carry out one detected collision

        Shared by both passes, since they differ only in how they choose
        the pair and the ratios along it.
        """
        node1, v1 = self.merge_node_for_segment(seg1[0], seg1[1], l1)
        if node1 is None:
            return
        node2, v2 = self.merge_node_for_segment(seg2[0], seg2[1], l2)
        if node2 is None or node2 == node1:
            return

        arms1 = self.G.out_degree(node1)
        arms2 = self.G.out_degree(node2)
        newpos = self.collision_point(node1, v1, node2, v2)
        if self.thrown_too_far(node1, node2, newpos):
            self.report(kind, interval, 'rejected_throw', seg1, seg2,
                        L1=l1, L2=l2, node1=node1, node2=node2,
                        arms1=arms1, arms2=arms2, planes=self.last_planes)
            return

        self.latch(node1, node2)
        merged, status = self.merge(node1, node2, newpos)
        outcome = 'merged' if merged is not None else 'merge_failed'
        self.report(kind, interval, outcome, seg1, seg2,
                    merged=merged, status=status, position=newpos,
                    L1=l1, L2=l2, node1=node1, node2=node2,
                    arms1=arms1, arms2=arms2, planes=self.last_planes)

    def report(self, kind, interval, outcome, seg1, seg2, **details):
        if self.record is not None:
            self.record.add(kind, interval, outcome,
                            seg1=seg1, seg2=seg2, **details)

    def candidate_cutoff(self, segments):
        """candidate_cutoff: how far apart two segments can be and still meet

        sqrt(4*dr2max + l2max/2) + rann, the bound exadis uses to size its
        neighbour list, where dr2max is the largest squared nodal
        displacement over the step and l2max the largest squared segment
        length. Two segments whose mid-points are further apart than this
        cannot have come within rann at any point of the interval.

        Reproducing the bound rather than inventing one matters for more
        than speed: it makes the two codes consider the same candidate
        pairs, so a disagreement is a disagreement about the criterion and
        not about which pairs were looked at.
        """
        dr2max = 0.0
        for tag in self.G.all_nodes_tags():
            now = self.G.nodes(tag).R
            was = self.G.cell.closest_image(Rref=now, R=self.old_position(tag))
            dr2max = max(dr2max, float(np.dot(now - was, now - was)))
        l2max = 0.0
        for t1, t2 in segments:
            vec = self.G.seg_vector(t1, t2)
            l2max = max(l2max, float(np.dot(vec, vec)))
        return np.sqrt(4.0 * dr2max + 0.5 * l2max) + self.rann

    def segment_pairs(self):
        """segment_pairs: candidate pairs, as tag pairs and never indices

        Pairs sharing a node are hinges and belong to the other pass. The
        rest are filtered on mid-point separation before the expensive
        swept-distance criterion sees them; the filter is vectorized
        because it is the only part of this rule that is O(N^2).
        """
        segments = [tuple(sorted(pair))
                    for pair in self.G.all_segments_tags()]
        if len(segments) < 2:
            return

        cutoff = self.candidate_cutoff(segments)
        vecs = np.array([self.G.seg_vector(t1, t2) for t1, t2 in segments])
        mid = np.array([self.G.nodes(t1).R for t1, t2 in segments]) + 0.5 * vecs
        half = 0.5 * np.linalg.norm(vecs, axis=1)
        # minimum image on the mid-point separations, so the filter is not
        # fooled by a pair either side of a periodic face
        deltas = mid[:, None, :] - mid[None, :, :]
        shape = deltas.shape
        reduced = np.array([self.G.cell.map(d)
                            for d in deltas.reshape(-1, 3)]).reshape(shape)
        # The cutoff bounds how far apart two segments can be, not how far
        # apart their mid-points are. A segment reaches half its length
        # beyond its own mid-point, so both half-lengths have to be added
        # before the comparison; without them a long pair that genuinely
        # collides near its ends is dropped before the criterion sees it.
        reach = cutoff + half[:, None] + half[None, :]
        near = np.einsum('ijk,ijk->ij', reduced, reduced) < reach * reach

        # The higher-indexed segment is yielded first, matching exadis,
        # whose pair loop keeps only k < i and so treats the higher index
        # as the pair's first segment. The criterion itself is symmetric
        # under the exchange, but the order decides which segment is split
        # first and which node survives the merge, so it has to agree.
        for i, seg_i in enumerate(segments):
            for k in range(i):
                if not near[i, k]:
                    continue
                seg_k = segments[k]
                if set(seg_i) & set(seg_k):
                    continue                  # shares a node: a hinge
                yield seg_i, seg_k

    def detect_all(self):
        """detect_all: every colliding pair, on the untouched network

        Detection is completed before any collision is carried out, which
        is how the reference implementation is organized: it finds all
        pairs in parallel and then executes them in sequence. The
        difference is observable. Detecting as we go would let an earlier
        collision move the nodes a later detection reads, and would offer
        the criterion segments split moments before, whose new nodes carry
        no velocity because the mobility law never saw them.
        """
        return [(seg1, seg2) + hit
                for seg1, seg2 in self.segment_pairs()
                for hit in [self.detect(seg1, seg2)] if hit is not None]

    def run_segment_pass(self):
        """run_segment_pass: collide segments that do not share a node

        Detect first, then execute, skipping any pair whose segments an
        earlier collision has already consumed.
        """
        for seg1, seg2, interval, l1, l2 in self.detect_all():
            if frozenset(seg1) in self.done or frozenset(seg2) in self.done:
                continue
            if not (self.G.has_segment(*seg1) and self.G.has_segment(*seg2)):
                continue
            self.collide(seg1, seg2, interval, l1, l2)

    def hinge_pairs(self):
        """hinge_pairs: pairs of arms of the same node, as tags

        Only nodes that existed before the segment pass ran are eligible,
        as hinge node or as either far node. Matches exadis, whose hinge
        loop is bounded by the pre-pass node count and drops an arm whose
        far node is beyond it:

            for (int i = 0; i < nnodes; i++) {
                ...
                if (n3 >= nnodes) continue;
                    ...
                    if (n4 >= nnodes) continue;

        so a node a split created this pass cannot take part in a zip.
        """
        for tag in self.pre_pass_nodes:
            if not self.G.has_node(tag):
                continue
            nbrs = [n for n in self.G.neighbors_tags(tag)
                    if n in self.pre_pass_nodes]
            for i, n3 in enumerate(nbrs):
                for n4 in nbrs[i + 1:]:
                    yield tag, n3, n4

    def detect_hinge(self, tag, n3, n4):
        """detect_hinge: are these two arms closed enough to zip

        The threshold is relaxed when the two far nodes are connected, so
        that a triangle about to collapse is caught earlier.
        """
        p1 = self.G.nodes(tag).R
        p3 = self.G.cell.closest_image(Rref=p1, R=self.G.nodes(n3).R)
        p4 = self.G.cell.closest_image(Rref=p1, R=self.G.nodes(n4).R)
        tol = HINGE_TOL_TRI if self.G.has_segment(n3, n4) else HINGE_TOL

        if hinge_cos_angle(p1, p3, p4) > tol:
            return 'retroactive'
        if not self.dt:
            return None
        q1 = self.predicted(tag, p1)
        q3, q4 = self.predicted(n3, p3), self.predicted(n4, p4)
        if hinge_cos_angle(q1, q3, q4) > tol:
            return 'predictive'
        return None

    def run_hinge_pass(self):
        """run_hinge_pass: zip arms of the same node that have closed up

        The shorter arm's far node is one merge node; the other comes from
        the longer arm at the ratio of the two lengths, which is where the
        shorter one reaches to.
        """
        for tag, n3, n4 in self.hinge_pairs():
            if not (self.G.has_segment(tag, n3)
                    and self.G.has_segment(tag, n4)):
                continue
            # No latch check here: changed to match exadis, whose hinge
            # loop carries the equivalent skipseg tests commented out and
            # so zips arms that a collision just latched.
            # DIFFERENT FROM ParaDiS, which sets NO_COLLISIONS on the node
            # a merge produced and skips it in the hinge loop, so it would
            # decline the zip that follows a merge.

            l13 = self.G.seg_length(tag, n3)
            l14 = self.G.seg_length(tag, n4)
            if min(l13, l14) < 1.0e-10:
                continue
            interval = self.detect_hinge(tag, n3, n4)
            if interval is None:
                continue

            short, long_ = (n4, n3) if l14 < l13 else (n3, n4)
            ratio = min(l14, l13) / max(l14, l13)
            node1, v1 = short, self.velocity(short)
            node2, v2 = self.merge_node_for_segment(tag, long_, ratio)
            if node2 is None or node2 == node1:
                continue

            newpos = self.collision_point(node1, v1, node2, v2)
            if self.thrown_too_far(node1, node2, newpos):
                self.report('hinge', interval, 'rejected_throw',
                            (tag, n3), (tag, n4))
                continue
            self.latch(node1, node2)
            merged, status = self.merge(node1, node2, newpos)
            self.report('hinge', interval,
                        'merged' if merged is not None else 'merge_failed',
                        (tag, n3), (tag, n4), merged=merged, status=status,
                        position=newpos)

    def run(self):
        # snapshot before the segment pass, because the hinge pass may not
        # act on nodes that pass creates. See hinge_pairs.
        self.pre_pass_nodes = frozenset(self.G.all_nodes_tags())
        self.run_segment_pass()
        self.run_hinge_pass()
        if self.purge_at_end:
            self.G.purge_network()


def handle_collision_retroactive(G: DisNet, state: dict, record=None,
                                 purge_at_end=True):
    """handle_collision_retroactive: one retroactive collision pass

    Raises when the inputs a retroactive rule is defined by are absent.
    ExaDiS instead substitutes the current positions for the old ones,
    which silently reduces the rule to a proximity test; PyDiS is
    deliberately stricter here, because the degradation is invisible in
    the result.
    """
    if "rann" not in state:
        raise KeyError("handle_collision_retroactive: state has no 'rann'")
    if old_positions_by_tag(state) is None:
        raise ValueError(
            "handle_collision_retroactive: no old nodal positions in state. "
            "Retroactive collision compares current positions against those "
            "at the start of the step; without them it would silently become "
            "a proximity test. Populate state['oldnodes_dict'] (see "
            "SimulateNetwork.save_old_nodes) before calling.")
    if not state.get("dt"):
        raise ValueError(
            "handle_collision_retroactive: state has no non-zero 'dt'. "
            "The predictive half of the criterion sweeps one timestep "
            "ahead, and a zero step would switch it off silently.")

    rule = RetroactiveCollision(G, state, record=record,
                                purge_at_end=purge_at_end)
    rule.run()

    if not G.is_sane():
        raise ValueError("handle_collision_retroactive: sanity check failed")
    return rule
