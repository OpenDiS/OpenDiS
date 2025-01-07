"""@package docstring
CalForce_Base: base class for Calforce

Defines the interface for the CalForce module
"""

from abc import ABC, abstractmethod 
from framework.disnet_manager import DisNetManager

import numpy as np
from typing import Tuple
Tag = Tuple[int, int]

class CalForce_Base(ABC):
    """CalForce_Base: base class for CalForce

    Defines the interface for the CalForce module
    """

    @abstractmethod
    def NodeForce(self, N: DisNetManager, state: dict, pre_compute: bool=True) -> dict:
        """NodeForce: compute all nodal forces and store them in the state dictionary
        """
        pass

    @abstractmethod
    def PreCompute(self, N: DisNetManager, state: dict) -> dict:
        """PreCompute: pre-compute some data for force calculation
        """
        pass

    @abstractmethod
    def OneNodeForce(self, N: DisNetManager, state: dict, tag: Tag, update_state: bool=True) -> np.array:
        """OneNodeForce: compute and return the force on one node specified by its tag
        """
        pass
