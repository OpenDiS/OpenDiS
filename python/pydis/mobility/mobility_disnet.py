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
            'Relax': self.Mobility_Relax }
        
    def Mobility(self, G: DisNet, nodeforce_dict: dict):
        """Mobility: calculate node velocity according to mobility law function
        """
        return self.Mobility_Functions[self.mobility_law](G, nodeforce_dict)

    def Mobility_Relax(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_Relax: node velocity equal node force
        """
        vel_dict = nodeforce_dict.copy()
        return vel_dict
