"""pydis SimulateNetwork's loading modes, on the A3 source with every node pinned.

Nothing can move, so each step's plastic strain is zero bitwise, and what the
load does follows from the mode alone:

1. 'stress': the applied stress stays as given; state["stress"] is its von
   Mises stress, here the axial stress, the load being uniaxial; the strain
   stays zero; the density is the network's line length / V / b^2.
2. 'stress_rate': the applied stress grows by srate dt each step.
3. 'strain_rate' with max_strain: the run ends at the first step whose strain
   reaches max_strain, long before max_step, with sigma = E eps throughout.
4. In every mode, stress_strain_dens.dat holds one row per print_freq steps,
   in the ExaDiS format, matching the state at those steps.

    python3 test_loading_modes_pydis.py
"""

import sys, os
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict
from pydis import SimulateNetwork, FCC_slip_sys
import fcc_frank_read_util as fr

LBOX, ARM, DT = 1000.0, 500.0, 1.0e-8
SIGMA, SRATE, ERATE = 1.0e8, 1.0e14, 1.0e3
NSTEP, PRINT_FREQ = 20, 5
EZZ = np.array([0.0, 0.0, 1.0])
OUTDIR = 'output'


class Recorded(SimulateNetwork):
    """Recorded: SimulateNetwork keeping each step's state"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.history = []

    def step_end(self, DM, state):
        self.history.append({k: np.copy(state[k]) for k in
                             ("istep", "dEp", "dWp", "applied_stress",
                              "strain", "stress", "density")})


def run(mode, max_step=NSTEP, **kwargs):
    """run: the pinned A3 source in loading_mode mode; returns (net, state, sim)"""
    ss = FCC_slip_sys()
    net = fr.init_frank_read_src_on_slip_system(
        ss, ss.name_to_sysID('A3'), arm_length=ARM, box_length=LBOX, pin_all=True)
    state = fr.default_state(LBOX)
    sim = Recorded(state=state, loading_mode=mode, max_step=max_step,
                   print_freq=PRINT_FREQ, write_dir=os.path.join(OUTDIR, mode),
                   **fr.dd_components(net, state, dt=DT), **kwargs)
    sim.run(net, state)
    return net, state, sim


def series(sim, key):
    return np.array([h[key] for h in sim.history])


def line_density(net, burgmag):
    """line_density: segment lengths summed over V burgmag^2, from DisNet.seg_vector"""
    G = net.get_disnet()
    length = sum(np.linalg.norm(G.seg_vector(t1, t2)) for t1, t2 in G.all_segments_dict())
    return length / G.cell.volume() / burgmag**2


def check_no_plastic_strain(mode, sim):
    """check_no_plastic_strain: a network that cannot move sweeps nothing, bitwise"""
    dEp, dWp = series(sim, 'dEp'), series(sim, 'dWp')
    return (report_close("%s: dEp zero at every step" % mode, dEp, np.zeros_like(dEp),
                         atol=0.0)
            & report_close("%s: dWp zero at every step" % mode, dWp, np.zeros_like(dWp),
                           atol=0.0))


def check_results_file(mode, sim):
    """check_results_file: 4. one row per print_freq steps, as the state was then"""
    path = os.path.join(OUTDIR, mode, 'stress_strain_dens.dat')
    with open(path) as f:
        header = f.readline().strip()
    rows = np.loadtxt(path, ndmin=2)
    printed = [h for h in sim.history if int(h['istep']) % PRINT_FREQ == 0]
    want = np.array([[h['istep'], h['strain'], h['stress'], h['density']] for h in printed])
    ok = report("%s: results header in the ExaDiS format" % mode,
                header == '# Step Strain Stress Density Walltime')
    ok &= report("%s: one results row per %d steps" % (mode, PRINT_FREQ),
                 rows.shape == (len(printed), 5), suffix=' %s' % (rows.shape,))
    if rows.shape == (len(printed), 5):
        # relative to each column's largest value, %e keeping 7 digits
        scale = np.where(np.abs(want).max(axis=0) > 0, np.abs(want).max(axis=0), 1.0)
        ok &= report_close("%s: rows match the state, relative to each column" % mode,
                           rows[:, :4] / scale, want / scale, atol=1e-6)
    return ok


def check_stress():
    """check_stress: 1. constant stress"""
    applied = fr.uniaxial_voigt(SIGMA, EZZ)
    net, state, sim = run('stress', applied_stress=applied)
    ok = check_no_plastic_strain('stress', sim)
    ok &= report_close("stress: applied stress unchanged", series(sim, 'applied_stress'),
                       np.tile(applied, (NSTEP, 1)), atol=0.0)
    ok &= report_close("stress: von Mises stress / sigma", series(sim, 'stress') / SIGMA,
                       np.ones(NSTEP), atol=1e-15)
    ok &= report_close("stress: strain", series(sim, 'strain'), np.zeros(NSTEP), atol=0.0)
    ok &= report_close("stress: density / line length over V b^2",
                       series(sim, 'density') / line_density(net, state["burgmag"]),
                       np.ones(NSTEP), atol=1e-15)
    return ok & check_results_file('stress', sim)


def check_stress_rate():
    """check_stress_rate: 2. the applied stress grows by srate dt per step"""
    srate = fr.uniaxial_voigt(SRATE, EZZ)
    _, _, sim = run('stress_rate', srate=srate)
    scale = NSTEP * SRATE * DT
    n = np.arange(1, NSTEP + 1)
    ok = check_no_plastic_strain('stress_rate', sim)
    ok &= report_close("stress_rate: applied stress / (NSTEP srate dt)",
                       series(sim, 'applied_stress') / scale,
                       np.outer(n, srate) * DT / scale, atol=1e-14)
    ok &= report_close("stress_rate: von Mises stress / (NSTEP srate dt)",
                       series(sim, 'stress') / scale, n * SRATE * DT / scale, atol=1e-14)
    return ok & check_results_file('stress_rate', sim)


def check_max_strain():
    """check_max_strain: 3. strain_rate stops at max_strain, before max_step"""
    max_strain = 4.5 * ERATE * DT
    _, state, sim = run('strain_rate', max_step=1000, erate=ERATE, edir=EZZ,
                        max_strain=max_strain)
    eps, sigma = series(sim, 'strain'), series(sim, 'stress')
    E = 2.0 * state["mu"] * (1.0 + state["nu"])
    ok = check_no_plastic_strain('strain_rate', sim)
    ok &= report("strain_rate: run ends after 5 steps, at the first strain past max_strain",
                 len(eps) == 5 and eps[-1] >= max_strain > eps[-2],
                 suffix=' %d steps' % len(eps))
    ok &= report_close("strain_rate: sigma / (E eps)", sigma / (E * eps), np.ones(len(eps)),
                       atol=1e-14)
    return ok & check_results_file('strain_rate', sim)


def main():
    print("pydis SimulateNetwork loading modes, A3 source with every node pinned")
    return all([check_stress(), check_stress_rate(), check_max_strain()])


if __name__ == "__main__":
    passed = verdict("test_loading_modes_pydis", main())
    sys.exit(0 if passed else 1)
