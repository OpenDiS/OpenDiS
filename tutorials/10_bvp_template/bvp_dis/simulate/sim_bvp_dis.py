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

#class SimulationDriver(SimulateNetwork_Base):
class SimulationDriver:
    """SimulateNetwork: class for simulating dislocation network with boundary value problem (BVP)

    """
    def __init__(self, state: dict, calforce=None,
                 mobility=None, timeint=None, topology=None,
                 collision=None, remesh=None, cross_slip=None, vis=None,
                 image_stress=None, image_force=None, surface_mobility=None, surface_topology=None,
                 dt0: float=1.0e-8,
                 max_step: int=10,
                 loading_mode: str=None,
                 applied_stress: np.ndarray=None,
                 print_freq: int=None,
                 plot_freq: int=None,
                 plot_pause_seconds: float=None,
                 write_freq: int=None,
                 write_dir: str=".",
                 save_state: bool=False,
                 **kwargs) -> None:
        self.calforce = calforce
        self.mobility = mobility
        self.timeint = timeint
        self.topology = topology
        self.collision = collision
        self.remesh = remesh
        self.cross_slip = cross_slip
        self.vis = vis

        self.image_stress = image_stress
        self.image_force = image_force
        self.surface_mobility = surface_mobility
        self.surface_topology = surface_topology

        self.dt0 = dt0
        self.max_step = max_step
        self.loading_mode = loading_mode
        self.applied_stress = np.array(applied_stress)
        self.print_freq = print_freq
        self.plot_freq = plot_freq
        self.plot_pause_seconds = plot_pause_seconds
        self.write_freq = write_freq
        self.write_dir = write_dir
        self.save_state = save_state

        state["applied_stress"] = np.array(applied_stress)

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

        if self.image_force is not None:
            state = self.image_force.AddImageForce(DM, state)

        state = self.mobility.Mobility(DM, state)

        # Adding mobility law implementation for surface nodes
        state = self.surface_mobility.Mobility(DM, state)

        state = self.timeint.Update(DM, state)
        #self.plastic_strain(DM, state)

    def step_post_integrate(self, DM: DisNetManager, state: dict):
        """step_post_integrate: invoked after time-integration of each time step
        """
        pass

    def step_topological_operations(self, DM: DisNetManager, state: dict):
        """step_topological_operations: invoked for handling topological events at each time step
        """
        if self.cross_slip is not None:
            self.cross_slip.Handle(DM, state)

        # The order of topology vs collision is opposite to ExaDiS
        if self.topology is not None:
            self.topology.Handle(DM, state)

        # Adding topology handling for dislocation intersecting surface
        if self.surface_topology is not None:
            self.surface_topology.Handle(DM, state)        

        if self.collision is not None:
            self.collision.HandleCol(DM, state)

        if self.remesh is not None:
            self.remesh.Remesh(DM, state)

    def step_update_response(self, DM: DisNetManager, state: dict):
        """step_update_response: update applied stress and rotation if needed
        """
        if self.loading_mode != 'stress':
            raise ValueError("invalide loading_mode in PyDiS SimulateNetwork")

        return state

    def step_write_files(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                DM.write_json(os.path.join(self.write_dir, f'disnet_{istep}.json'))
                if self.save_state:
                    with open(os.path.join(self.write_dir, f'state_{istep}.pickle'), 'wb') as file:
                        pickle.dump(state, file)

    def step_print_info(self, DM: DisNetManager, state: dict):
        if self.print_freq != None:
            istep = state['istep']
            if istep % self.print_freq == 0:
                print("step = %d dt = %e"%(istep, self.timeint.dt))

    def step_visualize(self, DM: DisNetManager, state: dict):
        G = DM.get_disnet(DisNet)
        if self.plot_freq != None:
            istep = state['istep']
            if istep % self.plot_freq == 0:
                self.vis.plot_disnet(G, fig=self.fig, ax=self.ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

    def step_end(self, DM: DisNetManager, state: dict):
        """step_end: invoked at the end of each time step
        """
        pass

    def step(self, DM: DisNetManager, state: dict):
        """step: take a time step of DD simulation on DisNetManager DM
        """
        # Step begin
        self.step_begin(DM, state)

        # Step time-integrate
        self.step_integrate(DM, state)

        # Step post-integrate
        self.step_post_integrate(DM, state)

        # Step topological operations
        self.step_topological_operations(DM, state)

        # Step update response
        self.step_update_response(DM, state)

        self.step_write_files(DM, state)
        self.step_print_info(DM, state)
        self.step_visualize(DM, state)

        # Step end
        self.step_end(DM, state)

        return state

    def run(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            os.makedirs(self.write_dir, exist_ok=True)

        G = DM.get_disnet(DisNet)
        if self.plot_freq != None:
            try: 
                self.fig = plt.figure(figsize=(8,8))
                self.ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.vis.plot_disnet(G, fig=self.fig, ax=self.ax, trim=True, block=False)

        for istep in range(self.max_step):
            state['istep'] = istep
            self.step(DM, state)


        # plot final configuration
        if self.plot_freq != None:
            G = DM.get_disnet(DisNet)
            self.vis.plot_disnet(G, fig=self.fig, ax=self.ax, trim=True, block=False)

        return state
