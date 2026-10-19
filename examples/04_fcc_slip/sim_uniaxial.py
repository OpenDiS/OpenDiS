"""Constant-strain-rate uniaxial loading on top of the stock PyDiS stepper.

SimulateNetworkUniaxial adds what pydis SimulateNetwork does not have yet:
the plastic strain of each step, the loadType 1 stress update that consumes
it, and the stress-strain output. Written to be folded into SimulateNetwork
later, so it overrides only the hooks it needs.

plastic_strain() runs at the end of step_integrate, right after time
integration and before any topological operation, which is where the
commented-out call in SimulateNetwork.step_integrate sits. It zeroes
state["plast_acc"] and adds the integration sweep; collision, topology, and
remesh then add what their own node moves sweep, and the load update reads
the total. topological_strain=False leaves the accumulator out, so only the
integration sweep counts, which is the leak this is there to close.
"""

import os
import time
import numpy as np

from pydis import DisNet, SimulateNetwork
from pydis.util.state_access import old_positions_by_tag
from pydis.util.fcc_slip_sys import unit
from pydis.util.plastic_strain import PlasticStrainAccum
import plast_strain_util as psu


class SimulateNetworkUniaxial(SimulateNetwork):
    """SimulateNetworkUniaxial: loading_mode 'strain_rate' along edir"""

    def __init__(self, state: dict, erate: float=1e3, edir=(0.0, 0.0, 1.0),
                 rotation: bool=False, topological_strain: bool=True,
                 **kwargs) -> None:
        super().__init__(state=state, loading_mode='strain_rate', **kwargs)
        self.erate = erate
        self.rotation = rotation
        self.results = []
        self.history = []
        self.t0 = time.perf_counter()
        state["edir"] = unit(edir)
        state["Etot"] = np.zeros(6)
        state["dEp"], state["dWp"] = np.zeros(6), np.zeros(3)
        state["strain"], state["stress"], state["density"] = 0.0, 0.0, 0.0
        state["dEp_integration"] = np.zeros(6)
        state["plast_acc"] = PlasticStrainAccum() if topological_strain else None

    def step_integrate(self, DM, state: dict):
        super().step_integrate(DM, state)
        self.plastic_strain(DM, state)

    def plastic_strain(self, DM, state: dict):
        """plastic_strain: dEp, dWp, and density swept by this step's integration"""
        G = DM.get_disnet(DisNet)
        oldpos = old_positions_by_tag(state)
        if oldpos is None:
            raise ValueError("plastic_strain: no start-of-step positions in state")
        dEp, dWp = psu.plastic_strain_increment(G, oldpos)
        state["dEp_integration"] = psu.sym_to_voigt(dEp)
        acc = state["plast_acc"]
        if acc is not None:
            acc.zero()
            acc.add(dEp, dWp)
        state["dEp"] = psu.sym_to_voigt(dEp)
        state["dWp"] = psu.asym_to_array(dWp)
        state["density"] = psu.dislocation_density(G, state["burgmag"])
        return state

    def step_update_response(self, DM, state: dict):
        """step_update_response: rotate if asked, then the loadType 1 update"""
        acc = state["plast_acc"]
        if acc is not None:
            state["dEp"] = psu.sym_to_voigt(acc.dEp)
            state["dWp"] = psu.asym_to_array(acc.dWp)
        if self.rotation:
            state["edir"], state["applied_stress"] = psu.rotate_loading(
                state["edir"], state["applied_stress"], state["dWp"])
        dstress, dstrain = psu.uniaxial_increment(
            state["dEp"], state["dt"], self.erate, state["mu"], state["nu"],
            state["edir"])
        A, A2 = psu.uniaxial_projectors(state["edir"])
        state["applied_stress"] += dstress * A
        state["Etot"] += dstrain * A
        state["strain"] = np.dot(state["Etot"], A2)
        state["stress"] = np.dot(state["applied_stress"], A2)
        return state

    def step_print_info(self, DM, state: dict):
        if self.print_freq is not None and state['istep'] % self.print_freq == 0:
            elapsed = time.perf_counter() - self.t0
            print("step = %d, nodes = %d, strain = %e, stress = %e, elapsed = %.1f sec"
                  % (state['istep'], DM.get_disnet(DisNet).num_nodes(),
                     state["strain"], state["stress"], elapsed))
            self.results.append([state['istep'], state["strain"], state["stress"],
                                 state["density"], elapsed])

    def step_end(self, DM, state: dict):
        self.history.append({k: np.copy(state[k]) for k in
                             ("istep", "dEp", "dWp", "dEp_integration",
                              "applied_stress", "Etot",
                              "strain", "stress", "density", "edir")})

    def write_results(self):
        """write_results: stress_strain_dens.dat, in the ExaDiS format"""
        os.makedirs(self.write_dir, exist_ok=True)
        with open(os.path.join(self.write_dir, 'stress_strain_dens.dat'), 'w') as f:
            f.write('# Step Strain Stress Density Walltime\n')
            np.savetxt(f, np.array(self.results).reshape(-1, 5), fmt='%d %e %e %e %e')

    def run(self, DM, state: dict):
        state = super().run(DM, state)
        self.write_results()
        return state
