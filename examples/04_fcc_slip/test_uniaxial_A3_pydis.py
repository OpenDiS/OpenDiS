"""Frank-Read source on A3 under constant strain rate along [001].

SimulateNetworkUniaxial with erate = 1e3 /s, ParaDiS loadType 1. Checks:

1. A network that cannot move gives dEp = 0 exactly and sigma = E eps.
2. With the source free, sigma = E (eps - eps_p) at every step, eps_p
   summed from e.dEp.e as a tensor contraction, independent of the Voigt
   arithmetic the update uses; and the curve falls measurably below E eps.
3. Every step's dEp is a pure A3 glide increment: traceless, Schmid ratio
   e.dEp.e / (2 n.dEp.b) = 1/sqrt(6), and n.dWp.b = n.dEp.b.
4. Translating the loop across two box faces changes no step's dEp beyond
   round-off.
5. stress_strain_dens.dat is written in the ExaDiS format.

The source bows, and so sweeps area, from the first loaded step; there is
no stage where a free source gives exactly zero plastic strain.

    python3 test_uniaxial_A3_pydis.py
"""

import sys, os
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict
from pydis import FCC_slip_sys
from pydis.util.fcc_slip_sys import unit
import fcc_frank_read_util as fr
import plast_strain_util as psu
from sim_uniaxial import SimulateNetworkUniaxial

LBOX, ARM, ERATE = 1000.0, 500.0, 1.0e3
EDIR = np.array([0.0, 0.0, 1.0])
NSTEP, NSTEP_SHORT = 300, 30
SHIFT = np.array([0.0, -100.0, 450.0])   # the loop then crosses the y and z faces
OUTDIR = 'output'


def run(nstep, shift=np.zeros(3), pin_all=False, write_dir=OUTDIR):
    """run: a strain-rate run on the A3 source; returns (net, state, sim)"""
    ss = FCC_slip_sys()
    net = fr.init_frank_read_src_on_slip_system(
        ss, ss.name_to_sysID('A3'), arm_length=ARM, box_length=LBOX,
        shift=shift, pin_all=pin_all)
    state = fr.default_state(LBOX)
    sim = SimulateNetworkUniaxial(state=state, erate=ERATE, edir=EDIR,
                                  max_step=nstep, print_freq=10,
                                  write_dir=write_dir, **fr.dd_components(net, state))
    sim.run(net, state)
    return net, state, sim


def series(sim, key):
    return np.array([h[key] for h in sim.history])


def dEp_tensors(sim):
    return np.array([psu.voigt_to_sym(v) for v in series(sim, 'dEp')])


def youngs(state):
    return 2.0 * state["mu"] * (1.0 + state["nu"])


def check_immobile():
    _, state, sim = run(50, pin_all=True, write_dir=os.path.join(OUTDIR, 'pinned'))
    ok = report_close("immobile network: dEp", series(sim, 'dEp'),
                      np.zeros((50, 6)), atol=0.0)
    strain, stress = series(sim, 'strain'), series(sim, 'stress')
    return ok & report_close("immobile network: sigma / (E eps)",
                             stress / (youngs(state) * strain), np.ones(50),
                             atol=1e-13)


def check_identity(state, sim):
    eps_p = np.cumsum([EDIR @ T @ EDIR for T in dEp_tensors(sim)])
    strain, stress = series(sim, 'strain'), series(sim, 'stress')
    E = youngs(state)
    ok = report_close("sigma = E (eps - eps_p) at every step, relative",
                      stress / (E * strain), 1.0 - eps_p / strain, atol=1e-12)
    print("    eps_p / eps at the end = %.4f, sigma %.2f MPa vs E eps %.2f MPa"
          % (eps_p[-1] / strain[-1], stress[-1] / 1e6, E * strain[-1] / 1e6))
    ok &= report("plastic strain from the first loaded steps",
                 eps_p[1] > 0.0 and np.all(np.diff(eps_p[1:]) > 0.0))
    return ok & report("curve falls more than 1%% below E eps by the end",
                       stress[-1] < 0.99 * E * strain[-1])


def check_a3_glide(sim):
    n, b = (unit(v) for v in FCC_slip_sys().sysID_to_nb(1))
    T = dEp_tensors(sim)
    W = np.array([psu.array_to_asym(w) for w in series(sim, 'dWp')])
    # steps where integration moved the arm; at step 0 the stress is still zero
    # and dEp holds only the round-off of remesh's midpoint insertions, ~1e-21
    moving = np.array([np.abs(v).max() > 0.0 for v in series(sim, 'dEp_integration')])
    scale = np.array([np.abs(t).max() for t in T])[moving]
    T, W = T[moving], W[moving]
    nEb = np.array([n @ t @ b for t in T])
    ok = report("steps with motion: %d of %d" % (moving.sum(), len(moving)),
                moving.sum() >= len(moving) - 1)
    # round-off here goes as ulp(coordinates) / displacement, and the first
    # loaded steps move the arm by about 0.005 b
    ok &= report_close("dEp traceless, relative",
                       np.array([np.trace(t) for t in T]) / scale,
                       np.zeros(len(T)), atol=1e-10)
    ok &= report_close("Schmid ratio e.dEp.e / (2 n.dEp.b) = 1/sqrt(6)",
                       np.array([EDIR @ t @ EDIR for t in T]) / (2.0 * nEb),
                       np.full(len(T), 1.0 / np.sqrt(6.0)), atol=1e-10)
    return ok & report_close("n.dWp.b = n.dEp.b, relative",
                             np.array([n @ w @ b for w in W]) / nEb,
                             np.ones(len(T)), atol=1e-10)


def check_straddle(sim):
    _, _, sim_s = run(NSTEP_SHORT, shift=SHIFT, write_dir=os.path.join(OUTDIR, 'shifted'))
    a = series(sim, 'dEp')[:NSTEP_SHORT]
    s = series(sim_s, 'dEp')
    return report_close("loop across two box faces: dEp per step, relative",
                        s / np.abs(a).max(), a / np.abs(a).max(), atol=1e-10)


def check_density(state, sim):
    rho0 = 4.0 * ARM / LBOX**3 / state["burgmag"]**2
    return report_close("density at step 0 = 4 L / V / b^2, relative",
                        series(sim, 'density')[0] / rho0, 1.0, atol=1e-12)


def check_output_file():
    path = os.path.join(OUTDIR, 'stress_strain_dens.dat')
    with open(path) as f:
        header = f.readline().strip()
    rows = np.loadtxt(path)
    return report("%s: ExaDiS header and %d rows" % (path, NSTEP // 10),
                  header == '# Step Strain Stress Density Walltime'
                  and rows.shape == (NSTEP // 10, 5))


def main():
    print("Frank-Read source on A3, erate %.0e /s along [001]" % ERATE)
    ok = check_immobile()
    net, state, sim = run(NSTEP)
    ok &= report("final network sane", net.is_sane())
    ok &= report("free arm still on A3",
                 fr.free_arm_sysIDs(net.get_disnet(), FCC_slip_sys()) == {1})
    ok &= check_identity(state, sim)
    ok &= check_a3_glide(sim)
    ok &= check_density(state, sim)
    ok &= check_output_file()
    ok &= check_straddle(sim)
    return ok


if __name__ == "__main__":
    passed = verdict("test_uniaxial_A3_pydis", main())
    sys.exit(0 if passed else 1)
