"""Plastic strain leaked by topological operations, measured on the A3 source.

The strain-rate run of test_uniaxial_A3_pydis.py, twice: with the
topological accounting on and off. A run leaks nothing exactly when the
plastic strain summed so far, P = dEp + dWp, equals F(start) - F(now) at
every step, F the network's area functional
(plast_strain_util.network_area_dyad). With the accounting on that must hold
to round-off; with it off, the failure is the leak, reported together with
what it does to the stress-strain curve.

F is checked only until the expanding loop could reach its own periodic
image, at about step 594 here: the A3 plane holds the lattice vector
(1000, 1000, 0), 1414 b long, and once the loop reconnects with its image F
is no longer defined, though every merge is still accounted.

    python3 test_leak_measurement.py [nsteps]
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import os
import numpy as np
from framework.testing import report, verdict
from pydis import DisNet, FCC_slip_sys
import fcc_frank_read_util as fr
import plast_strain_util as psu
from sim_uniaxial import SimulateNetworkUniaxial

LBOX, ARM, ERATE = 1000.0, 500.0, 1.0e3
MARGIN = 50.0   # past 2 rann plus a step's growth
NSTEP = int(sys.argv[1]) if len(sys.argv) > 1 else 600


class Tracked(SimulateNetworkUniaxial):
    """Tracked: also records the area functional F after every step"""

    def step_end(self, DM, state):
        super().step_end(DM, state)
        G = DM.get_disnet(DisNet)
        self.history[-1]["extent"] = psu.largest_extent(G)
        self.history[-1]["F"] = (psu.network_area_dyad(G)
                                 if self.history[-1]["extent"] < LBOX - MARGIN
                                 else np.full((3, 3), np.nan))
        self.history[-1]["nodes"] = G.num_nodes()


def run(topological):
    ss = FCC_slip_sys()
    net = fr.init_frank_read_src_on_slip_system(ss, ss.name_to_sysID('A3'),
                                                arm_length=ARM, box_length=LBOX)
    state = fr.default_state(LBOX)
    F0 = psu.network_area_dyad(net.get_disnet())
    sim = Tracked(state=state, erate=ERATE, max_step=NSTEP, print_freq=50,
                  topological_strain=topological,
                  write_dir=os.path.join('output', 'leak_%s' % ('on' if topological else 'off')),
                  **fr.dd_components(net, state))
    sim.run(net, state)
    return F0, sim


def checkable_steps(sim):
    """checkable_steps: steps before any piece could touch its own periodic image"""
    reach = [i for i, h in enumerate(sim.history) if h["extent"] >= LBOX - MARGIN]
    return reach[0] if reach else len(sim.history)


def identity_residual(F0, sim, nsteps):
    """identity_residual: worst and last |sum P - (F0 - F)| over nsteps, relative to |F0 - F|"""
    P = np.cumsum([psu.voigt_to_sym(h["dEp"]) + psu.array_to_asym(h["dWp"])
                   for h in sim.history[:nsteps]], axis=0)
    dF = np.array([F0 - h["F"] for h in sim.history[:nsteps]])
    scale = np.abs(dF).max()
    return np.abs(P - dF).max() / scale, np.abs(P[-1] - dF[-1]).max() / scale


def topological_share(history):
    """topological_share: what the topological operations added, relative to the total"""
    topo = sum(np.abs(h["dEp"] - h["dEp_integration"]).max() for h in history)
    total = sum(np.abs(h["dEp"]).max() for h in history)
    return topo / total


def main():
    print("topological plastic strain leak, A3 source, %d steps" % NSTEP)
    F0_on, on = run(topological=True)
    F0_off, off = run(topological=False)
    n_on, n_off = checkable_steps(on), checkable_steps(off)
    ok = report("loop stays clear of its periodic image for %d of %d steps"
                % (n_on, NSTEP), n_on >= 400)
    worst_on, _ = identity_residual(F0_on, on, n_on)
    worst_off, end_off = identity_residual(F0_off, off, n_off)
    ok &= report("accounting on: sum P = F(start) - F(now) at every one of them, "
                 "worst relative %.1e" % worst_on, worst_on < 1e-9)
    print("    topological operations carried %.2e of |dEp| before step %d, %.2e after"
          % (topological_share(on.history[:n_on]), n_on,
             topological_share(on.history[n_on:])))
    print("    accounting off, same steps: worst residual %.2e, %.2e at step %d"
          % (worst_off, end_off, n_off))
    s_on, s_off = on.history[-1]["stress"], off.history[-1]["stress"]
    print("    final stress %.4f MPa on, %.4f MPa off; nodes %d on, %d off"
          % (s_on/1e6, s_off/1e6, on.history[-1]["nodes"], off.history[-1]["nodes"]))
    return ok


if __name__ == "__main__":
    passed = verdict("test_leak_measurement", main())
    sys.exit(0 if passed else 1)
