"""@package docstring
CalForce_Base: base class for Calforce

Defines interface for the CalForce module
"""

from abc import ABC, abstractmethod 
from framework.disnet_manager import DisNetManager

from typing import Tuple
Tag = Tuple[int, int]

class CalForce_Base(ABC):
    """CalForce_Base: base class for CalForce

    Defines interface for the CalForce module
    """

    @abstractmethod
    def NodeForce(self, DM: DisNetManager, state: dict) -> dict:
        """NodeForce: return nodal forces in a dictionary
        """
        pass

    def PreCompute(self, DM: DisNetManager, state: dict) -> dict:
        """PreCompute: pre-compute some data for force calculation
        """
        pass

    def OneNodeForce(self, DM: DisNetManager, state: dict, tag: Tag, update_state: bool=True) -> dict:
        """OneNodeForce: compute force calculation on one node
        """
        pass
