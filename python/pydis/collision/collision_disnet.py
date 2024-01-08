"""@package docstring
Collision_DisNet: class for detecting and handing collisions between segments

Provide collision handling functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet

from .getmindist2_paradis import GetMinDist2

class Collision:
    """Collision: class for detecting and handling collisions

    """
    def __init__(self, collision_mode: str='Proximity', **kwargs) -> None:
        self.collision_mode = collision_mode
        self.mindist2 = kwargs.get('mindist2', 1.0e-3)

        self.HandleCol_Functions = {
            'Proximity': self.HandleCol_Proximity }
        
    def HandleCol(self, G: DisNet) -> None:
        """HandleCol: handle collision according to collision_mode
        """
        return self.HandleCol_Functions[self.collision_mode](G)

    def HandleCol_Proximity(self, G: DisNet) -> None:
        """HandleCol_Proximity: handle collision using Proximity criterion
           This is a much simplified version of Proximity collision handling in ParaDiS
        """
        # loop through all segment pairs to check for collision
        # To do: implement cell list to speed up
        segments = G.seg_list()
        collided = np.zeros(len(segments), dtype=bool)
        nseg = len(segments)
        for i in range(nseg):
            seg1 = segments[i]
            if collided[i]:
                continue
            tag1, tag2 = seg1["edge"][0], seg1["edge"][1]
            if not G.has_edge(tag1, tag2):
                continue
            for j in range(i+1, nseg):
                seg2 = segments[j]
                if collided[i] or collided[j]:
                    continue
                if not G.has_edge(tag1, tag2):
                    continue
                tag3, tag4 = seg2["edge"][0], seg2["edge"][1]
                if not G.has_edge(tag3, tag4):
                    continue
                p1, p2 = np.array(seg1["R1"]), np.array(seg1["R2"])
                p3, p4 = np.array(seg2["R1"]), np.array(seg2["R2"])
                v1, v2 = np.zeros(3), np.zeros(3)
                v3, v4 = np.zeros(3), np.zeros(3)
                if tag1 != tag3 and tag1 != tag4 and tag2 != tag3 and tag2 != tag4:
                    # no nodes are shared
                    # To do: apply PBC here
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
                            G.nodes()[mergedTag]['R'] = newPos
                        else:
                            print('mergedTag = None  status = %s', status)

        # In ParaDiS the hinge case (zipping) is handled separately
        # They are skipped here for simplicity

        if not G.is_sane():
            raise ValueError("HandleCol_Proximity: sanity check failed")

    def HandleCol_Proximity_ParaDiS(self, G: DisNet) -> None:
        """HandleCol_Proximity: using ProximityCollision of ParaDiS
        """
        raise NotImplementedError("HandleCol_Proximity_ParaDiS: not implemented yet")

    def HandleCol_Retroactive_ParaDiS(self, G: DisNet) -> None:
        """HandleCol_Retroactive: using RetroactiveCollision of ParaDiS
        """
        raise NotImplementedError("HandleCol_Retroactive_ParaDiS: not implemented yet")
