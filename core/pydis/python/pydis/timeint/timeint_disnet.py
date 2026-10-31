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
    def __init__(self, state: dict={}, integrator: str='EulerForward',
                 dt: float=1e-8, force=None, mobility=None, multi: int=0) -> None:
        self.integrator = integrator
        self.dt = dt

        self.Update_Functions = {
            'EulerForward': self.Update_EulerForward }

        if integrator == 'Trapezoid':
            if force is None or mobility is None:
                raise ValueError("TimeIntegration Trapezoid: give the force and mobility modules")
            self.force, self.mobility, self.multi = force, mobility, multi
            self.nextdt = state.get("nextdt", 1e-12)
            self.maxdt = state.get("maxdt", 1e-7)
            self.rtol = state.get("rtol", -1.0)
            if self.rtol <= 0.0:
                self.rtol = 0.25 * state["a"]

    def Update(self, DM: DisNetManager, state: dict) -> None:
        """TimeIntegration: update node position given velocity
        """
        if self.integrator == 'Trapezoid':
            return self.Update_Trapezoid(DM, state)
        G = DM.get_disnet(DisNet)
        if "nodevels" in state and "nodeveltags" in state:
            DisNet.convert_nodevel_array_to_dict(state)

        vel_dict = state["vel_dict"]
        applied_stress = state["applied_stress"]
        self.Update_Functions[self.integrator](G, vel_dict, applied_stress)
        state["dt"] = self.dt
        t = state.get('time', 0.0)
        state["time"] = t + state["dt"]
        return state

    def Update_EulerForward(self, G: DisNet, vel_dict: dict, applied_stress: np.ndarray) -> None:
        """TimeIntegration_EulerForward: Euler forward time integration
        """
        # Fold into the primary cell to reduce round-off error, as ParaDiS
        # does (FoldBox in Util.c), also consistent with exadis
        for tag, vel in vel_dict.items():
            G.nodes(tag).R = G.cell.fold(G.nodes(tag).R + vel * self.dt)

    def Update_Trapezoid(self, DM: DisNetManager, state: dict) -> dict:
        """Update_Trapezoid: trapezoid steps with an error control on the time step, multi of
           them in a row (as ExaDiS' IntegratorMulti) until maxdt; the start-of-step positions
           stay those of the first, for the plastic strain and collisions
        """
        if "nodevels" in state and "nodeveltags" in state:
            DisNet.convert_nodevel_array_to_dict(state)
        dt, steps = 0.0, 0
        while dt < self.maxdt:
            dt += self.Trapezoid_Step(DM, state)
            steps += 1
            if steps >= max(self.multi, 1):
                break
        self.dt = dt
        state["dt"] = dt
        state["time"] = state.get('time', 0.0) + dt
        return state

    def Trapezoid_Step(self, DM: DisNetManager, state: dict) -> float:
        """Trapezoid_Step: one step of ParaDiS' TrapezoidIntegrator (TrapezoidIntegrator.c):
           predict the positions, then correct them with the velocities there, halving the
           time step until the error is below rtol; returns the time step taken. The prediction
           is ExaDiS' (integrator_trapezoid.h), from the current velocities only.
        """
        G = DM.get_disnet(DisNet)
        newdt = min(self.maxdt, self.nextdt)
        if newdt <= 0.0:
            newdt = self.maxdt
        oldpos = {tag: G.nodes(tag).R.copy() for tag in G.all_nodes_tags()}
        currvel = {tag: state["vel_dict"][tag].copy() for tag in oldpos}

        convergent, incr_delta = False, True
        while not convergent:
            for tag, R in oldpos.items():
                G.nodes(tag).R = G.cell.fold(R + currvel[tag] * newdt)
            for iteration in range(2):
                state = self.force.NodeForce(DM, state)
                state = self.mobility.Mobility(DM, state)
                vel = state["vel_dict"]
                err = {tag: G.nodes(tag).R - G.cell.closest_image(Rref=G.nodes(tag).R, R=R)
                       - (vel[tag] + currvel[tag]) * 0.5 * newdt for tag, R in oldpos.items()}
                if any(np.isnan(e).any() for e in err.values()):
                    raise ValueError("TimeIntegration Trapezoid: NaN in the nodes' positions")
                errmax = max((np.max(np.abs(e)) for e in err.values()), default=0.0)
                if errmax < self.rtol:
                    convergent = True
                    break
                incr_delta = False
                if iteration == 0:
                    for tag, e in err.items():
                        G.nodes(tag).R = G.cell.fold(G.nodes(tag).R - e)
            if not convergent:
                newdt *= 0.5
                if newdt < 1.0e-20:
                    raise ValueError("TimeIntegration Trapezoid: the time step dropped below "
                                     "1e-20 s")

        self.nextdt = min(self.maxdt, newdt * 1.2) if incr_delta else newdt
        return newdt
