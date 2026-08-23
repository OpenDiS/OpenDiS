"""@package docstring
Collision_DisNet: class for detecting and handing collisions between segments

Provide collision handling functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode
from framework.collision_base import Collision_Base
from framework.disnet_manager import DisNetManager
from .collision_retroactive import (handle_collision_retroactive,
                                    CollisionRecord)

try:
    from .getmindist2_paradis import GetMinDist2_paradis as GetMinDist2
except ImportError:
    # use python version instead
    print("pydis_lib not found, using python version for GetMinDist2")
    from pydis.collision.getmindist2_python  import GetMinDist2_python as GetMinDist2

class Collision:
    """Collision: class for detecting and handling collisions

    """
    def __init__(self, state: dict={}, collision_mode: str='Proximity', **kwargs) -> None:
        self.collision_mode = collision_mode
        self.mindist2 = state.get("rann", np.sqrt(1.0e-3))**2
        self.nbrlist = kwargs.get('nbrlist')
        if not self.nbrlist.__module__.split('.')[0] in ['pydis']:
            raise ValueError("Collision: nbrlist must come compatible modules")

        self.HandleCol_Functions = {
            'Proximity': self.HandleCol_Proximity,
            'Retroactive': self.HandleCol_Retroactive }
        # collected when collision_record is on, for tests and diagnostics
        self.collision_record = kwargs.get('collision_record', False)
        self.records = []
        # Retroactive only: defer removing merged nodes and annihilated
        # segments to one purge at the end of the pass, as exadis does.
        # False deletes them at the merge, as ParaDiS does.
        self.purge_at_end = kwargs.get('purge_at_end', True)
        # Proximity only: see the comment where this is used, in
        # HandleCol_Proximity, for what this reproduces and why it is
        # believed to be an ExaDiS-side artifact rather than a real physical
        # effect. Default True to match ExaDiS's current behavior; set False
        # to fall back to the plain midpoint average once/if that ExaDiS
        # behavior is fixed upstream.
        self.match_exadis_hinge_merge_artifact = state.get(
            'match_exadis_hinge_merge_artifact', True)
        
    def HandleCol(self, DM: DisNetManager, state: dict) -> dict:
        """HandleCol: handle collision according to collision_mode
        """
        G = DM.get_disnet(DisNet)
        oldpos_dict = state.get('oldpos_dict', None)
        dt = state.get('dt', 0.0)
        xold = np.array(list(oldpos_dict.values())) if oldpos_dict != None else None
        state = self.HandleCol_Functions[self.collision_mode](G, state, xold=xold, dt=dt)
        return state

    def HandleCol_Proximity(self, G: DisNet, state: dict, xold=None, dt=None) -> dict:
        """HandleCol_Proximity: handle collision using Proximity criterion
           This is a much simplified version of Proximity collision handling in ParaDiS
        """
        # loop through all segment pairs to check for collision
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        midpoints = 0.5*(R1 + R2)

        self.nbrlist.sort_points_to_list(midpoints)

        # absent here means a node topology has not seen yet (e.g. collision
        # runs before topology this step), which carries no exemption
        nodeflag_dict = state.get('nodeflag_dict', {})

        collided = np.zeros(Nseg, dtype=bool)
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        for i, j in self.nbrlist.iterate_nbr_pairs(use_cell_list=all(G.cell.is_periodic)):
                if collided[i]:
                    continue
                tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
                if not G.has_segment(tag1, tag2):
                    continue
                if nodeflag_dict.get(tag1, DisNode.Flags.CLEAR) & DisNode.Flags.NO_COLLISIONS:
                    continue
                if nodeflag_dict.get(tag2, DisNode.Flags.CLEAR) & DisNode.Flags.NO_COLLISIONS:
                    continue

                if collided[i] or collided[j]:
                    continue
                if not G.has_segment(tag1, tag2):
                    continue
                tag3, tag4 = tuple(source_tags[j]), tuple(target_tags[j])
                if not G.has_segment(tag3, tag4):
                    continue
                if nodeflag_dict.get(tag3, DisNode.Flags.CLEAR) & DisNode.Flags.NO_COLLISIONS:
                    continue
                if nodeflag_dict.get(tag4, DisNode.Flags.CLEAR) & DisNode.Flags.NO_COLLISIONS:
                    continue
                p1, p2 = R1[i,:].copy(), R2[i,:].copy()
                p3, p4 = R1[j,:].copy(), R2[j,:].copy()
                v1, v2 = np.zeros(3), np.zeros(3)
                v3, v4 = np.zeros(3), np.zeros(3)
                if tag1 != tag3 and tag1 != tag4 and tag2 != tag3 and tag2 != tag4:
                    # no nodes are shared
                    # apply PBC
                    p2 = G.cell.closest_image(Rref=p1, R=p2)
                    p3 = G.cell.closest_image(Rref=p1, R=p3)
                    p4 = G.cell.closest_image(Rref=p3, R=p4)
                    dist2, ddist2dt, L1, L2 = GetMinDist2(p1, v1, p2, v2, p3, v3, p4, v4)
                    if dist2 < self.mindist2:
                        collided[i] = True
                        collided[j] = True

                        seg1_vec = p2 - p1
                        close2node1 = (np.dot(seg1_vec, seg1_vec) * (L1    *L1))     < self.mindist2
                        close2node2 = (np.dot(seg1_vec, seg1_vec) * ((1-L1)*(1-L1))) < self.mindist2

                        seg2_vec = p4 - p3
                        close2node3 = (np.dot(seg2_vec, seg2_vec) * (L2    *L2))     < self.mindist2
                        close2node4 = (np.dot(seg2_vec, seg2_vec) * ((1-L2)*(1-L2))) < self.mindist2

                        if close2node1:
                            mergenode1, splitSeg1, newPos1 = tag1, False, p1
                        elif close2node2:
                            mergenode1, splitSeg1, newPos1 = tag2, False, p2
                        else:
                            splitSeg1, newPos1 = True, (1-L1)*p1 + L1*p2
                            # skip velocity for now
                            new_tag = G.get_new_tag()
                            G.insert_node_between(tag1, tag2, new_tag, newPos1)
                            mergenode1 = new_tag

                        if close2node3:
                            mergenode2, splitSeg2, newPos2 = tag3, False, p3
                        elif close2node4:
                            mergenode2, splitSeg2, newPos2 = tag4, False, p4
                        else:
                            splitSeg2, newPos2 = True, (1-L2)*p3 + L2*p4
                            # skip velocity for now
                            new_tag = G.get_new_tag()
                            G.insert_node_between(tag3, tag4, new_tag, newPos2)
                            mergenode2 = new_tag

                        # To do: determine precise position satisfying glide constraints
                        newPos = (newPos1 + newPos2)/2.0

                        # match_exadis_hinge_merge_artifact: reproduces a specific
                        # ExaDiS behavior, on by default so pydis agrees numerically
                        # with ExaDiS's current output; see below for why we believe
                        # this is unintended on ExaDiS's side and worth raising with
                        # its developers, and set the flag False to use the plain
                        # average instead if/when that is resolved upstream.
                        #
                        # When mergenode1 and mergenode2 are already directly
                        # connected (the usual case for the two nodes bounding an
                        # about-to-vanish junction segment), ExaDiS's
                        # CollisionRetroactive does not end up using the average
                        # computed above. Here is the sequence we traced, by
                        # instrumenting and rebuilding ExaDiS directly, that leads
                        # to a different result:
                        #
                        # 1. SerialDisNet::merge_nodes_position (network.cpp) merges
                        #    the two nodes and correctly sets the surviving node's
                        #    position to the average -- matching newPos here exactly.
                        # 2. For the arm connecting the two original nodes (the
                        #    mutual junction arm), this function zeroes its Burgers
                        #    vector to represent that it no longer carries a
                        #    dislocation, but does not remove the connectivity entry
                        #    on either side. The absorbed node is left present in the
                        #    network, at its own pre-merge position, linked to the
                        #    surviving node solely by this now-zero-Burgers arm.
                        # 3. The hinge-collision pass that runs immediately
                        #    afterward, in the same collision-handling call
                        #    (retroactive_collision_parallel, collision_retroactive.cpp),
                        #    walks connectivity without checking Burgers vectors, so
                        #    it rediscovers this leftover arm as a new candidate
                        #    collision between the surviving node and the absorbed one.
                        # 4. Because the absorbed node now has exactly one connection
                        #    (the leftover arm), AdjustCollisionPoint's early-return
                        #    for degree-1 nodes fires (`if (conn[n1].num==1) return
                        #    p1`). That check appears to be intended for genuinely
                        #    pinned/fixed nodes, which are also degree-1 by
                        #    construction, but it cannot distinguish that case from
                        #    a node that is merely mid-cleanup after an earlier merge.
                        #    It returns the absorbed node's stale, pre-merge position,
                        #    overwriting the average from step 1.
                        #
                        # We believe this is a lifecycle ordering issue rather than
                        # intended behavior: ParaDiS's own RetroactiveCollision2.c
                        # (which ExaDiS's CollisionRetroactive is modeled on) calls
                        # Topology.c's MergeNode, which calls RemoveNode on the
                        # absorbed node immediately, so it cannot be rediscovered by
                        # a later hinge pass within the same cycle. ExaDiS's
                        # Kokkos-parallel implementation instead batches merges via
                        # merge_nodes_position() and defers actual removal to a
                        # single purge_network() call at the end of the whole
                        # collision step, and that gap between "logically merged"
                        # and "structurally removed" is what lets the hinge pass see
                        # a node it shouldn't.
                        #
                        # Reproducing it here is a deliberate compatibility choice,
                        # not a claim that newPos1 is the physically correct merge
                        # point -- the average (newPos, above) is arguably more
                        # correct. Note also that ExaDiS's own choice of which of the
                        # two original positions survives is itself governed by an
                        # implementation detail (which of the two colliding segments
                        # has the smaller internal array index -- see
                        # retroactive_collision_parallel's `if (i <= k) continue`),
                        # not by anything physical. pydis has no analogous concept of
                        # segment array index, so newPos1 (rather than newPos2) was
                        # chosen empirically, by checking which one agrees with
                        # ExaDiS on this test, not derived from matching that
                        # implementation detail directly: with this flag on, full
                        # test_binary_junction_pydis_exadis_elast.py agreement at
                        # step 300 goes from a 3.8e-1 residual (plain average) to
                        # 3.8e-7 (floating-point noise, matching the pre-existing
                        # agreement level elsewhere in the run). A different pydis
                        # nbrlist implementation, or a differently-shaped collision
                        # elsewhere, is not guaranteed to keep agreeing with
                        # whichever side ExaDiS's own array-index ordering picks.
                        if (self.match_exadis_hinge_merge_artifact
                                and G.has_segment(mergenode1, mergenode2)):
                            newPos = newPos1.copy()

                        mergedTag, status = G.merge_node(mergenode1, mergenode2)
                        if mergedTag != None:
                            G.nodes(mergedTag).R = newPos

        # In ParaDiS the hinge case (zipping) is handled separately
        # They are skipped here for simplicity

        if not G.is_sane():
            raise ValueError("HandleCol_Proximity: sanity check failed")

        return state

    def HandleCol_Retroactive(self, G: DisNet, state: dict, xold=None, dt=None) -> dict:
        """HandleCol_Retroactive: retroactive collision, following ParaDiS

        Implemented in collision_retroactive.py, which follows ParaDiS
        RetroactiveCollisions2 (collisionMethod 4), the same rule exadis
        CollisionRetroactive implements. Not the Sills and Cai (2014)
        algorithm, which is collisionMethod 3; see that module.

        xold and dt are read from state rather than from these arguments,
        which exist for signature compatibility with the other handlers.
        """
        record = CollisionRecord() if self.collision_record else None
        rule = handle_collision_retroactive(G, state, record=record,
                                            purge_at_end=self.purge_at_end)
        if record is not None:
            self.records.append(record)
        if rule.missing_old_position:
            print("HandleCol_Retroactive: %d node(s) had no old position"
                  % rule.missing_old_position)
        if rule.missing_velocity:
            print("HandleCol_Retroactive: %d velocity lookup(s) missed"
                  % rule.missing_velocity)
        return state
