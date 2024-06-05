"""@package docstring
TimeInt_DisNet: class for Time Integration

Provide time integration functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet
from framework.disnet_manager import DisNetManager

class TimeIntegration:
    """TimeIntegration: class for time integration

    """
    def __init__(self, params: dict={}, integrator: str='EulerForward',
                 dt: float=1e-8) -> None:
        self.integrator = integrator
        self.dt = dt

        self.Update_Functions = {
            'EulerForward': self.Update_EulerForward }
        
    def Update(self, DM: DisNetManager, vel_dict: dict, applied_stress: np.ndarray) -> None:
        """TimeIntegration: update node position given velocity
        """
        G = DM.get_disnet(DisNet)
        return self.Update_Functions[self.integrator](G, vel_dict, applied_stress)

    def Update_EulerForward(self, G: DisNet, vel_dict: dict, applied_stress: np.ndarray) -> None:
        """TimeIntegration_EulerForward: Euler forward time integration
        """
        for tag, vel in vel_dict.items():
            G.nodes[tag]["R"] += vel * self.dt
