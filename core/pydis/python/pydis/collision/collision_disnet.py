"""@package docstring
Collision_DisNet: class for detecting and handing collisions between segments

Provide collision handling functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode
from framework.disnet_manager import DisNetManager

try:
    from .getmindist2_paradis import GetMinDist2_paradis as GetMinDist2
except ImportError:
    # use python version instead
    print("pydis_lib not found, using python version for GetMinDist2")
    from pydis.collision.getmindist2_python  import GetMinDist2_python as GetMinDist2

class Collision:
    """Collision: class for detecting and handling collisions

    """
    def __init__(self, state: dict={}, collision_mode: str='Proximity') -> None:
        self.collision_mode = collision_mode
        self.mindist2 = state.get("rann", np.sqrt(1.0e-3))**2

        self.HandleCol_Functions = {
            'Proximity': self.HandleCol_Proximity }
        
    def HandleCol(self, DM: DisNetManager, state: dict) -> None:
        """HandleCol: handle collision according to collision_mode
        """
        G = DM.get_disnet(DisNet)
        oldpos_dict = state.get('oldpos_dict', None)
        dt = state.get('dt', 0.0)
        xold = np.array(list(oldpos_dict.values())) if oldpos_dict != None else None
        self.HandleCol_Functions[self.collision_mode](G, xold, dt)
        return state

    def HandleCol_Proximity(self, G: DisNet, xold=None, dt=None) -> None:
        """HandleCol_Proximity: handle collision using Proximity criterion
           This is a much simplified version of Proximity collision handling in ParaDiS
        """
        # loop through all segment pairs to check for collision
        segments = G.seg_list()
        midpoints = G.get_segments_midpoint(segments)
        G.cell_list.sort_points_to_list(midpoints)

        collided = np.zeros(len(segments), dtype=bool)
        for i, j in G.cell_list.iterate_nbr_pairs(use_cell_list=all(G.cell.is_periodic)):
                seg1 = segments[i]
                if collided[i]:
                    continue
                tag1, tag2 = seg1["edge"][0], seg1["edge"][1]
                if not G.has_edge(tag1, tag2):
                    continue
                if G.nodes[tag1].get('flag') & DisNode.Flags.NO_COLLISIONS:
                    continue
                if G.nodes[tag2].get('flag') & DisNode.Flags.NO_COLLISIONS:
                    continue

                seg2 = segments[j]
                if collided[i] or collided[j]:
                    continue
                if not G.has_edge(tag1, tag2):
                    continue
                tag3, tag4 = seg2["edge"][0], seg2["edge"][1]
                if not G.has_edge(tag3, tag4):
                    continue
                if G.nodes[tag3].get('flag') & DisNode.Flags.NO_COLLISIONS:
                    continue
                if G.nodes[tag4].get('flag') & DisNode.Flags.NO_COLLISIONS:
                    continue
                p1, p2 = np.array(seg1["R1"]), np.array(seg1["R2"])
                p3, p4 = np.array(seg2["R1"]), np.array(seg2["R2"])
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
                            G.insert_node(tag1, tag2, new_tag, newPos1)
                            mergenode1 = new_tag

                        if close2node3:
                            mergenode2, splitSeg2, newPos2 = tag3, False, p3
                        elif close2node4:
                            mergenode2, splitSeg2, newPos2 = tag4, False, p4
                        else:
                            splitSeg2, newPos2 = True, (1-L2)*p3 + L2*p4
                            # skip velocity for now
                            new_tag = G.get_new_tag()
                            G.insert_node(tag3, tag4, new_tag, newPos2)
                            mergenode2 = new_tag

                        # To do: determine precise position satisfying glide constraints
                        newPos = (newPos1 + newPos2)/2.0
                        mergedTag, status = G.merge_node(mergenode1, mergenode2)
                        if mergedTag != None:
                            G.nodes[mergedTag]['R'] = newPos

        # In ParaDiS the hinge case (zipping) is handled separately
        # They are skipped here for simplicity

        if not G.is_sane():
            raise ValueError("HandleCol_Proximity: sanity check failed")

    def HandleCol_Proximity_ParaDiS(self, G: DisNet, xold=None, dt=None) -> None:
        """HandleCol_Proximity: using ProximityCollision of ParaDiS
        """
        raise NotImplementedError("HandleCol_Proximity_ParaDiS: not implemented yet")

    def HandleCol_Retroactive_ParaDiS(self, G: DisNet, xold=None, dt=None) -> None:
        """HandleCol_Retroactive: using RetroactiveCollision of ParaDiS
        """
        raise NotImplementedError("HandleCol_Retroactive_ParaDiS: not implemented yet")
