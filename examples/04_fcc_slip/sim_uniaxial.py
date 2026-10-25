"""Constant-strain-rate uniaxial loading, with the history the tests read.

pydis SimulateNetwork computes each step's plastic strain and, with
loading_mode='strain_rate', the loadType 1 stress update that consumes it.
SimulateNetworkUniaxial only records, at the end of every step, what the
tests check, among it the integration sweep alone (dEp_integration), taken
before collision, topology, and remesh add theirs.

topological_strain=False leaves the accumulator out, so only the
integration sweep counts, which is the leak the accounting is there to
close.
"""

import numpy as np

from pydis import SimulateNetwork


class SimulateNetworkUniaxial(SimulateNetwork):
    """SimulateNetworkUniaxial: loading_mode 'strain_rate', and a per-step history"""

    def __init__(self, state: dict, topological_strain: bool=True,
                 **kwargs) -> None:
        super().__init__(state=state, loading_mode='strain_rate', **kwargs)
        self.history = []
        self.dEp_integration = np.zeros(6)
        if not topological_strain:
            state["plast_acc"] = None

    def plastic_strain(self, DM, state: dict):
        state = super().plastic_strain(DM, state)
        self.dEp_integration = state["dEp"].copy()
        return state

    def step_end(self, DM, state: dict):
        record = {k: np.copy(state[k]) for k in
                  ("istep", "dEp", "dWp", "applied_stress", "Etot",
                   "strain", "stress", "density", "edir")}
        record["dEp_integration"] = self.dEp_integration
        self.history.append(record)
