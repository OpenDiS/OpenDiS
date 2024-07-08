"""@package docstring
Sim_DisNet: class for simulating dislocation network

Provide simulation functions based on other utlitity classes
"""

import numpy as np
from ..disnet import DisNet
from ..calforce.calforce_disnet import CalForce
from ..mobility.mobility_disnet import MobilityLaw
from ..timeint.timeint_disnet import TimeIntegration
from ..visualize.vis_disnet import VisualizeNetwork
from framework.disnet_manager import DisNetManager

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

class SimulateNetwork:
    """SimulateNetwork: class for simulating dislocation network

    """
    def __init__(self, state: dict, calforce=None,
                 mobility=None, timeint=None, topology=None,
                 collision=None, remesh=None, vis=None,
                 dt0: float=1.0e-8,
                 max_step: int=10,
                 loading_mode: str=None,
                 applied_stress: np.ndarray=None,
                 print_freq: int=None,
                 plot_freq: int=None,
                 plot_pause_seconds: float=None,
                 write_freq: int=None,
                 write_dir: str=".",
                 **kwargs) -> None:
        self.calforce = calforce
        self.mobility = mobility
        self.timeint = timeint
        self.topology = topology
        self.collision = collision
        self.remesh = remesh
        self.vis = vis
        self.dt0 = dt0
        self.max_step = max_step
        self.loading_mode = loading_mode
        self.applied_stress = np.array(applied_stress)
        self.print_freq = print_freq
        self.plot_freq = plot_freq
        self.plot_pause_seconds = plot_pause_seconds
        self.write_freq = write_freq
        self.write_dir = write_dir

        state["applied_stress"] = np.array(applied_stress)
        state["sim"] = self

    def step(self, DM: DisNetManager, state: dict):
        """step: take a time step of DD simulation on DisNet G
        """
        state = self.calforce.NodeForce(DM, state)

        state = self.mobility.Mobility(DM, state)

        # using a constant time step (for now)
        state = self.timeint.Update(DM, state)

        state = self.topology.Handle(DM, state)

        if self.collision is not None:
            state = self.collision.HandleCol(DM, state)

        if self.remesh is not None:
            state = self.remesh.Remesh(DM, state)

        return state

    def run(self, DM: DisNetManager, state: dict):
        G = DM.get_disnet(DisNet)
        if self.plot_freq != None:
            try: 
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)

        for tstep in range(self.max_step):
            self.step(DM, state)

            if self.print_freq != None:
                if tstep % self.print_freq == 0:
                    print("step = %d dt = %e"%(tstep, self.timeint.dt))

            G = DM.get_disnet(DisNet)
            if self.plot_freq != None:
                if tstep % self.plot_freq == 0:
                    self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

        # plot final configuration
        G = DM.get_disnet(DisNet)
        self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)

        return state
