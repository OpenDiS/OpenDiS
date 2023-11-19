"""@package docstring
TimeInt_DisNet: class for Time Integration

Provide time integration functions given a DisNet object
"""

import numpy as np
from disnet import DisNet

class TimeIntegration:
    """TimeIntegration: class for time integration

    """
    def __init__(self, integrator: str='EulerForward', 
                 dt: float=1e-8, **kwargs) -> None:
        self.integrator = integrator
        self.dt = dt

        self.Update_Functions = {
            'EulerForward': self.Update_EulerForward }
        
    def Update(self, G: DisNet, vel_dict: dict) -> None:
        """TimeIntegration: update node position given velocity
        """
        return self.Update_Functions[self.integrator](G, vel_dict)

    def Update_EulerForward(self, G: DisNet, vel_dict: dict) -> None:
        """TimeIntegration_EulerForward: Euler forward time integration
        """
        for node, vel in vel_dict.items():
            G.nodes()[node]["R"] += vel * self.dt
