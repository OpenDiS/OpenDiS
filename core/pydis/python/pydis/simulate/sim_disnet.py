"""@package docstring
Sim_DisNet: class for simulating dislocation network

Provide simulation functions based on other utlitity classes
"""

import numpy as np
import os, pickle, time
from ..disnet import DisNet
from ..calforce.calforce_disnet import CalForce
from ..mobility.mobility_disnet import MobilityLaw
from ..timeint.timeint_disnet import TimeIntegration
from ..visualize.vis_disnet import VisualizeNetwork
from ..util.state_access import old_positions_by_tag
from ..util.fcc_slip_sys import unit
from ..util.plastic_strain import (PlasticStrainAccum, plastic_strain_from_data,
                                   density_from_data, sym_to_voigt, asym_to_array,
                                   rotate_loading, uniaxial_projectors,
                                   uniaxial_increment, von_mises)
from framework.disnet_manager import DisNetManager

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

#class SimulateNetwork(SimulateNetwork_Base):
class SimulateNetwork:
    """SimulateNetwork: class for simulating dislocation network

    loading_mode, as in exadis' SimulateNetwork: 'stress' keeps
    applied_stress, 'stress_rate' adds srate per unit time, 'strain_rate'
    loads uniaxially along edir at erate. Every mode computes the plastic
    strain of each step, the topological operations' share included, and
    the density, which takes burgmag (default state["burgmag"]). rotation
    turns edir and the applied stress by the plastic spin; max_strain, when
    given, ends the run in place of max_step.
    """
    def __init__(self, state: dict, calforce=None,
                 mobility=None, timeint=None, topology=None,
                 collision=None, remesh=None, cross_slip=None, vis=None,
                 collide_before_split: bool=True,
                 dt0: float=1.0e-8,
                 max_step: int=10,
                 loading_mode: str='stress',
                 applied_stress: np.ndarray=np.zeros(6),
                 print_freq: int=None,
                 plot_freq: int=None,
                 plot_pause_seconds: float=None,
                 write_freq: int=None,
                 write_dir: str=".",
                 save_state: bool=False,
                 srate: np.ndarray=np.zeros(6),
                 erate: float=0.0,
                 edir: np.ndarray=np.array([0., 0., 1.]),
                 burgmag: float=None,
                 rotation: bool=False,
                 max_strain: float=None,
                 **kwargs) -> None:
        self.calforce = calforce
        self.mobility = mobility
        self.timeint = timeint
        self.topology = topology
        self.collision = collision
        self.remesh = remesh
        self.cross_slip = cross_slip
        self.vis = vis
        self.collide_before_split = collide_before_split
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
        self.srate = np.array(srate, dtype=float)
        self.erate = erate
        self.rotation = rotation
        self.max_strain = max_strain
        self.burgmag = state.get("burgmag", 1.0) if burgmag is None else burgmag
        self.results = []
        self.t0 = time.perf_counter()

        state["applied_stress"] = np.array(applied_stress, dtype=float)
        state["edir"] = unit(edir)
        state["Etot"] = np.zeros(6)
        state["strain"], state["stress"], state["density"] = 0.0, 0.0, 0.0
        state["dEp"], state["dWp"] = np.zeros(6), np.zeros(3)
        # this step's dEp and dWp; collision, topology and remesh add to it
        state["plast_acc"] = PlasticStrainAccum()

    def step_begin(self, DM: DisNetManager, state: dict):
        """step_begin: invoked at the begining of each time step
        """
        pass

    def save_old_nodes(self, DM: DisNetManager, state: dict):
        """save_old_nodes: record positions as they are at the start of the step

        Retroactive collision compares the path swept between these and the
        positions after integration, so without them it has nothing to test
        against. Keyed by tag rather than by array offset, because the
        topological operations reindex the nodes.
        """
        G = DM.get_disnet(DisNet)
        state['oldpos_dict'] = {tag: G.nodes(tag).R.copy()
                                for tag in G.all_nodes_tags()}
        return state

    def plastic_strain(self, DM: DisNetManager, state: dict):
        """plastic_strain: dEp, dWp and density swept by this step's time integration

        Runs after time integration and before any topological operation, as
        in exadis' driver: it restarts state["plast_acc"] from this sweep,
        collision, topology and remesh add what their own node moves sweep,
        and step_update_response reads the total. It reads whichever network
        is current, so a mixed pydis/exadis run converts nothing extra.
        """
        oldpos = old_positions_by_tag(state)
        if oldpos is None:
            raise ValueError("plastic_strain: no start-of-step positions in state")
        data, cell = DM.export_data(), DM.cell
        dEp, dWp = plastic_strain_from_data(data, cell, oldpos)
        acc = state.get("plast_acc")
        if acc is not None:
            acc.zero()
            acc.add(dEp, dWp)
        state["dEp"], state["dWp"] = sym_to_voigt(dEp), asym_to_array(dWp)
        state["density"] = density_from_data(data, cell, self.burgmag)
        return state

    def step_integrate(self, DM: DisNetManager, state: dict):
        """step_integrate: invoked for time-integration at each time step
        """
        self.save_old_nodes(DM, state)
        state = self.calforce.NodeForce(DM, state)
        state = self.mobility.Mobility(DM, state)
        state = self.timeint.Update(DM, state)
        self.plastic_strain(DM, state)

    def step_post_integrate(self, DM: DisNetManager, state: dict):
        """step_post_integrate: invoked after time-integration of each time step
        """
        pass

    def step_topological_operations(self, DM: DisNetManager, state: dict):
        """step_topological_operations: invoked for handling topological events at each time step
        """
        if self.cross_slip is not None:
            self.cross_slip.Handle(DM, state)

        # collide_before_split=True matches ExaDiS, which collides before splitting so that a
        # multi-node formed by collision can be split within the same step; collide_before_split=
        # False matches ParaDiS's split-then-collide order instead.
        if self.collide_before_split:
            if self.collision is not None:
                self.collision.HandleCol(DM, state)

            if self.topology is not None:
                self.topology.Handle(DM, state)
        else:
            if self.topology is not None:
                self.topology.Handle(DM, state)

            if self.collision is not None:
                self.collision.HandleCol(DM, state)

        if self.remesh is not None:
            self.remesh.Remesh(DM, state)

    def step_update_response(self, DM: DisNetManager, state: dict):
        """step_update_response: update applied stress and rotation if needed

        As in exadis' driver. The step's dEp and dWp, topological operations
        included, turn edir and the applied stress when rotation is on;
        'strain_rate' then loads along edir (ParaDiS LoadCurve, loadType 1)
        and state["stress"] is the axial stress, while 'stress_rate' adds
        srate dt and state["stress"] is the von Mises stress, as for 'stress'.
        """
        acc = state.get("plast_acc")
        if acc is not None:
            state["dEp"], state["dWp"] = sym_to_voigt(acc.dEp), asym_to_array(acc.dWp)
        if self.rotation:
            state["edir"], state["applied_stress"] = rotate_loading(
                state["edir"], state["applied_stress"], state["dWp"])

        if self.loading_mode == 'strain_rate':
            dstress, dstrain = uniaxial_increment(
                state["dEp"], state["dt"], self.erate, state["mu"], state["nu"],
                state["edir"])
            A, A2 = uniaxial_projectors(state["edir"])
            state["applied_stress"] += dstress * A
            state["Etot"] += dstrain * A
            state["strain"] = np.dot(state["Etot"], A2)
            state["stress"] = np.dot(state["applied_stress"], A2)
        elif self.loading_mode in ('stress', 'stress_rate'):
            if self.loading_mode == 'stress_rate':
                state["applied_stress"] += self.srate * state["dt"]
            state["strain"] = 0.0
            state["stress"] = von_mises(state["applied_stress"])
        else:
            raise ValueError("invalid loading_mode in PyDiS SimulateNetwork")

        return state

    def write_results(self):
        """write_results: stress_strain_dens.dat in write_dir, in the format of exadis
        """
        if len(self.results) > 0:
            if self.write_dir:
                os.makedirs(self.write_dir, exist_ok=True)
            with open(os.path.join(self.write_dir, 'stress_strain_dens.dat'), 'w') as f:
                f.write('# Step Strain Stress Density Walltime\n')
                np.savetxt(f, np.array(self.results), fmt='%d %e %e %e %e')

    def step_write_files(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                DM.write_json(os.path.join(self.write_dir, f'disnet_{istep}.json'))
                if self.save_state:
                    with open(os.path.join(self.write_dir, f'state_{istep}.pickle'), 'wb') as file:
                        pickle.dump(state, file)
                # the results so far, as exadis does
                if self.print_freq != None:
                    self.write_results()

    def step_print_info(self, DM: DisNetManager, state: dict):
        if self.print_freq != None:
            istep = state['istep']
            if istep % self.print_freq == 0:
                dt = self.timeint.dt if self.timeint else 0.0
                elapsed = time.perf_counter() - self.t0
                if self.loading_mode == 'strain_rate':
                    print("step = %d, nodes = %d, dt = %e, strain = %e, elapsed = %.1f sec"
                          % (istep, DM.num_nodes(), dt, state["strain"], elapsed))
                else:
                    print("step = %d, nodes = %d, dt = %e, time = %e, elapsed = %.1f sec"
                          % (istep, DM.num_nodes(), dt, state.get("time", 0.0), elapsed))
                self.results.append([istep, state["strain"], state["stress"],
                                     state["density"], elapsed])

    def step_visualize(self, DM: DisNetManager, state: dict):
        if self.vis != None and self.plot_freq != None:
            istep = state['istep']
            if istep % self.plot_freq == 0:
                self.vis.plot_disnet(DM, fig=self.fig, ax=self.ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

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

    def iterate(self, istep: int, state: dict) -> bool:
        """iterate: whether to take time step istep

        max_strain, when given, decides in place of max_step, as in exadis.
        """
        if self.max_strain is not None:
            return np.abs(state["strain"]) < self.max_strain
        return istep < self.max_step

    def run(self, DM: DisNetManager, state: dict):
        self.t0 = time.perf_counter()
        if self.write_freq != None:
            os.makedirs(self.write_dir, exist_ok=True)

        if self.vis != None and self.plot_freq != None:
            try: 
                self.fig = plt.figure(figsize=(8,8))
                self.ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.vis.plot_disnet(DM, fig=self.fig, ax=self.ax, trim=True, block=False)

        istep = 0
        while self.iterate(istep, state):
            state['istep'] = istep
            self.step(DM, state)
            istep += 1

        if self.print_freq != None:
            self.write_results()

        # plot final configuration
        if self.vis != None and self.plot_freq != None:
            self.vis.plot_disnet(DM, fig=self.fig, ax=self.ax, trim=True, block=False)

        return state
