"""@package docstring
Sim_DisNet: class for simulating dislocation network

Provide simulation functions based on other utlitity classes
"""

import numpy as np
from disnet import DisNet
from calforce_disnet import CalForce
from mobility_disnet import MobilityLaw
from timeint_disnet import TimeIntegration
from vis_disnet import VisualizeNetwork

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
    def __init__(self, calforce,
                 mobility=None, timeint=None, collision=None, vis=None,
                 dt0: float=1.0e-8,
                 max_step: int=10,
                 print_freq: int=None,
                 plot_freq: int=None,
                 plot_pause_seconds: float=None,
                 **kwargs) -> None:
        self.calforce = calforce
        if mobility != None:
            self.mobility = mobility
        if timeint != None:
            self.timeint = timeint
        if collision != None:
            self.collision = collision
        if vis != None:
            self.vis = vis
        self.dt0 = dt0
        self.max_step = max_step
        if print_freq != None:
            self.print_freq = print_freq
        if plot_freq != None:
            self.plot_freq = plot_freq
        if plot_pause_seconds != None:
            self.plot_pause_seconds = plot_pause_seconds
        pass

    def step(self, G: DisNet):
        """step: take a time step of DD simulation on DisNet G
        """
        nodeforce_dict = self.calforce.NodeForce(G)

        vel_dict = self.mobility.Mobility(G, nodeforce_dict)

        #for node, vel in vel_dict.items():
        #    G.nodes()[node]["R"] += vel * self.dt

        # using a constant time step (for now)
        self.timeint.dt = self.dt0
        self.timeint.Update(G, vel_dict)

        if hasattr(self, 'collision'):
            self.collision.HandleCol(G)

    def run(self, G: DisNet):
        if self.plot_freq != None:
            try: 
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)

        for tstep in range(self.max_step):
            self.step(G)

            if self.print_freq != None:
                if tstep % self.print_freq == 0:
                    print("step = %d dt = %e"%(tstep, self.timeint.dt))

            if self.plot_freq != None:
                if tstep % self.plot_freq == 0:
                    self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

        # plot final configuration
        self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)
