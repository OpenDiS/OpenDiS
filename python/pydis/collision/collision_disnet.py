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
        self.col_dist = kwargs.get('col_dist', 1.0e-3)

        self.HandleCol_Functions = {
            'Proximity': self.HandleCol_Proximity }
        
    def HandleCol(self, G: DisNet) -> None:
        """HandleCol: handle collision according to collision_mode
        """
        return self.HandleCol_Functions[self.collision_mode](G)

    def HandleCol_Proximity(self, G: DisNet) -> None:
        """HandleCol_Proximity: handle collision using Proximity criterion
        """
        print("HandleCol_Proximity: being implemented right now...")
        # loop through all segment pairs to check for collision
        # To do: implement cell list to speed up
        segments = G.seg_list()
        nseg = len(segments)
        for i in range(nseg):
            for j in range(i+1, nseg):
                seg1 = segments[i]
                seg2 = segments[j]
                p1, p2 = np.array(seg1["R1"]), np.array(seg1["R2"])
                p3, p4 = np.array(seg2["R1"]), np.array(seg2["R2"])
                v1, v2 = np.zeros(3), np.zeros(3)
                v3, v4 = np.zeros(3), np.zeros(3)
                tag1, tag2 = seg1["edge"][0], seg1["edge"][1]
                tag3, tag4 = seg2["edge"][0], seg2["edge"][1]
                if tag1 != tag3 and tag1 != tag4 and tag2 != tag3 and tag2 != tag4:
                    # no nodes are shared
                    # To do: apply PBC here
                    dist2, ddist2dt, L1, L2 = GetMinDist2(p1, v1, p2, v2, p3, v3, p4, v4)
                    if dist2 < self.col_dist:
                        print("collision detected: %s %s %s %s dist2 = %e" % (tag1, tag2, tag3, tag4, dist2))
                        # To do: handle collision
                elif tag1 == tag3 and tag2 != tag4:
                    # hinge case
                    pass
                elif tag1 == tag4 and tag2 != tag3:
                    # hinge case
                    pass
                elif tag2 == tag3 and tag1 != tag4:
                    # hinge case
                    pass
                elif tag2 == tag4 and tag1 != tag3:
                    # hinge case
                    pass

    def HandleCol_Proximity_ParaDiS(self, G: DisNet) -> None:
        """HandleCol_Proximity: using ProximityCollision of ParaDiS
        """
        print("HandleCol_Proximity_ParaDiS: not implemented yet")

    def HandleCol_Retroactive_ParaDiS(self, G: DisNet) -> None:
        """HandleCol_Retroactive: using RetroactiveCollision of ParaDiS
        """
        print("HandleCol_Retroactive_ParaDiS: not implemented yet")
