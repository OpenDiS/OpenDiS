"""@package docstring
Mobility_DisNet: class for defining Mobility Laws

Provide mobility law functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet

class MobilityLaw:
    """MobilityLaw: class for mobility laws

    """
    def __init__(self, mobility_law: str='Relax', **kwargs) -> None:
        self.mobility_law = mobility_law

        self.Mobility_Functions = {
            'Relax': self.Mobility_Relax,
            'SimpleGlide': self.Mobility_SimpleGlide }
        
    def Mobility(self, G: DisNet, nodeforce_dict: dict):
        """Mobility: calculate node velocity according to mobility law function
        """
        return self.Mobility_Functions[self.mobility_law](G, nodeforce_dict)

    def Mobility_Relax(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_Relax: node velocity equal node force
        """
        vel_dict = nodeforce_dict.copy()
        # set velocity of pinned nodes to zero
        for tag in G.nodes:
            if G.nodes[tag]['flag'] == 7:
                vel_dict[tag] = np.zeros(3)
        return vel_dict

    def Mobility_SimpleGlide(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_SimpleGlide: node velocity equal node force divided by sum of arm length / 2
           To do: add glide constraints
        """
        vel_dict = nodeforce_dict.copy()
        for tag in G.nodes:
            node1 = G.nodes[tag]
            # set velocity of pinned nodes to zero
            if node1['flag'] == 7:
                vel_dict[tag] = np.zeros(3)
            else:
                R1 = node1["R"]
                Lsum = 0.0
                for nbr_tag in G.neighbors(tag):
                    node2 = G.nodes[nbr_tag]
                    R2 = node2["R"]
                    # To do: apply PBC here
                    Lsum += np.linalg.norm(R2-R1)
                vel_dict[tag] = vel_dict[tag] / (Lsum/2.0)
        return vel_dict
