"""@package docstring
MobilityLaw_Base: base class for MobilityLaw

Defines the interface for the MobilityLaw module
"""

from abc import ABC, abstractmethod 
from framework.disnet_manager import DisNetManager

from typing import Tuple
Tag = Tuple[int, int]

class MobilityLaw_Base(ABC):
    """MobilityLaw_Base: base class for MobilityLaw

    Defines the interface for the MobilityLaw module
    """

    @abstractmethod
    def Mobility(self, N: DisNetManager, state: dict) -> dict:
        """Mobility: compute all nodal velocities and store them in the state dictionary
        """
        pass

    @abstractmethod
    def OneNodeMobility(self, N: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
        """OneNodeMobility: compute and return the mobility of one node specified by its tag
        """
        pass
