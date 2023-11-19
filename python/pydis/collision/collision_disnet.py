"""@package docstring
Collision_DisNet: class for detecting and handing collisions between segments

Provide collision handling functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet

class Collision:
    """Collision: class for detecting and handling collisions

    """
    def __init__(self, collision_mode: str='Proximity', **kwargs) -> None:
        self.collision_mode = collision_mode

        self.HandleCol_Functions = {
            'Proximity': self.HandleCol_Proximity }
        
    def HandleCol(self, G: DisNet) -> None:
        """HandleCol: handle collision according to collision_mode
        """
        return self.HandleCol_Functions[self.collision_mode](G)

    def HandleCol_Proximity(self, G: DisNet) -> None:
        """HandleCol_Proximity: handle collision using Proximity criterion
        """
        print("HandleCol_Proximity: not implemented yet")
