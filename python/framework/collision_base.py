"""@package docstring
Collision_Base: base class for Collision

Defines interface for the Collision module
"""

from abc import ABC, abstractmethod 
from framework.disnet_manager import DisNetManager

from typing import Tuple
Tag = Tuple[int, int]

class Collision_Base(ABC):
    """Collision_Base: base class for Collision

    Defines interface for the Collision module
    """

    @abstractmethod
    def HandleCol(self, DM: DisNetManager, state: dict) -> dict:
        """HandleCol: handle collision
        """
        pass
