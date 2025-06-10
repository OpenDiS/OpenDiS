"""@package docstring
SimulationDriver: class for simulating dislocation network with boundary value problem (BVP)

Provide simulation functions based on other utlitity classes
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
import pickle
from pydis.disnet import DisNet
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.simulate.sim_disnet import SimulateNetwork
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.visualize.vis_disnet import VisualizeNetwork
from framework.disnet_manager import DisNetManager

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

class SimulationDriver(SimulateNetwork):
    """SimulationDriver: class for simulating dislocation network with boundary value problem (BVP)
           extends SimulateNetwork from PyDiS
    """
    def __init__(self, state: dict,
                 image_stress=None, image_force=None, surface_mobility=None, surface_topology=None,
                 **kwargs) -> None:
        super().__init__(state, **kwargs)

        self.image_stress = image_stress
        self.image_force = image_force
        self.surface_mobility = surface_mobility
        self.surface_topology = surface_topology

    def step_begin(self, DM: DisNetManager, state: dict):
        """step_begin: invoked at the begining of each time step
        """
        if self.image_stress is not None:
            state = self.image_stress.CalImageStress(DM, state)

    def step_integrate(self, DM: DisNetManager, state: dict):
        """step_integrate: invoked for time-integration at each time step
        """
        #self.save_old_nodes(DM, state)

        state = self.calforce.NodeForce(DM, state)
        state = self.mobility.Mobility(DM, state)
        state = self.timeint.Update(DM, state)
        #self.plastic_strain(DM, state)

    def step_topological_operations(self, DM: DisNetManager, state: dict):
        """step_topological_operations: invoked for handling topological events at each time step
        """
        if self.cross_slip is not None:
            self.cross_slip.Handle(DM, state)

        # The order of topology vs collision is opposite to ExaDiS
        if self.topology is not None:
            self.topology.Handle(DM, state)

        #ToDo: may need to write a combined Topology class for both bulk and surface nodes
        if self.surface_topology is not None:
            self.surface_topology.Handle(DM, state)        

        if self.collision is not None:
            self.collision.HandleCol(DM, state)

        if self.remesh is not None:
            self.remesh.Remesh(DM, state)
