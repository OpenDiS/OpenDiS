"""@package docstring
Remesh_DisNet: class for defining Remesh functions

Provide remesh functions given a DisNet object
"""

import numpy as np
from disnet import DisNet

class Remesh:
    """Remesh: class for remeshing dislocation network

    """
    def __init__(self, remesh_rule: str='LengthBased', **kwargs) -> None:
        self.remesh_rule = remesh_rule

        self.Remesh_Functions = {
            'LengthBased': self.Remesh_LengthBased }
        
    def Remesh(self, G: DisNet) -> None:
        """Remesh: remesh dislocation network according to remesh_rule
        """
        return self.Remesh_Functions[self.remesh_rule](G)

    def Remesh_LengthBased(self, G: DisNet) -> None:
        """Remesh_LengthBased: remesh dislocation network according to segment length
        """
        print("Remesh_LengthBased: not implemented yet")
