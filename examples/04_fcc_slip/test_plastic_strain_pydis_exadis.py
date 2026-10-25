"""PyDiS' plastic strain against ExaDiS', on the same motion and in the same run.

1. One prescribed Euler step of random loops crossing a box face, in both
   codes: the moved positions agree, bitwise in the strict tier, and dEp,
   dWp and the density agree to round-off four ways: PyDiS; ExaDiS' python
   driver (exadis_plastic_strain=False); ExaDiS' C++, System::plastic_strain;
   and PyDiS' formula read off the ExaDiS network, as PyDiS' driver reads it
   when an ExaDiS integrator runs inside it.
2. The A3 source under constant strain rate along [001], through PyDiS'
   driver and through ExaDiS', each with its own LineTension force,
   SimpleGlide mobility and Euler integrator, and no collision, topology, or
   remesh; the free arm is refined first, so the source can bow. Positions,
   strain, stress, dEp, dWp and density agree at every step, with the loading
   axis fixed and with the axis and the stress turned by the plastic spin.

With topological operations the two codes cannot agree exactly: ExaDiS
accounts no plastic spin there, counts the segment between two merging nodes
twice, and gives none to the segment a split creates (see
test_topological_plastic_strain.py).

    python3 test_plastic_strain_pydis_exadis.py
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict, tier_tol
from framework.disnet_manager import DisNetManager
import pydis
from pydis import DisNet, FCC_slip_sys
from pydis.util.plastic_strain import (plastic_strain_from_data, density_from_data,
                                       sym_to_voigt, asym_to_array)
import pyexadis
import pyexadis_base as ex
import fcc_frank_read_util as fr

LBOX, ARM, DT, ERATE = 1000.0, 500.0, 1.0e-8, 1.0e3
EZZ = np.array([0.0, 0.0, 1.0])
NLOOP, DISP = 10, 40.0      # check 1: random loops, rms nodal displacement
NSTEP = 300                 # check 2: time steps of each run

# check 1: the moved positions, and the four paths, relative to max |dEp|
STEP_POS_TOL = tier_tol(0.0, 1e-12)
STEP_TOL = tier_tol(1e-14, 1e-12)
# check 2: two codes' force, mobility and rotation, whose last bits differ
RUN_TOL = {"positions": tier_tol(2e-11, 2e-8), "strain": tier_tol(1e-14, 1e-11),
           "stress": tier_tol(1e-13, 1e-10), "dEp": tier_tol(1e-12, 1e-9),
           "dWp": tier_tol(1e-12, 1e-9), "density": tier_tol(1e-13, 1e-10),
           "edir": tier_tol(1e-14, 1e-11)}


def positions(DM):
    """positions: tag -> position of the current network"""
    data = DM.export_data()
    return {tuple(t): p.copy() for t, p in zip(np.asarray(data["nodes"]["tags"]),
                                              np.asarray(data["nodes"]["positions"]))}


def pydis_strain(DM, oldpos, burgmag):
    """pydis_strain: (dEp, dWp, density) of the current network, by PyDiS' formula"""
    data, cell = DM.export_data(), DM.cell
    dEp, dWp = plastic_strain_from_data(data, cell, oldpos)
    return sym_to_voigt(dEp), asym_to_array(dWp), density_from_data(data, cell, burgmag)


def one_step(G, vels, state):
    """one_step: one Euler step at vels in both codes

    Returns (PyDiS positions, ExaDiS positions, {path: (dEp, dWp, density)}).
    """
    tags = np.array(list(G.all_nodes_tags()))
    oldpos = {tuple(t): G.nodes(tuple(t)).R.copy() for t in tags}
    burgmag = state["burgmag"]

    st_py = dict(state, nodevels=vels.copy(), nodeveltags=tags.copy(),
                 applied_stress=np.zeros(6))
    py = DisNetManager(G.copy())
    pydis.TimeIntegration(integrator='EulerForward', dt=DT, state=st_py).Update(py, st_py)

    st_ex = dict(st_py, nodevels=vels.copy(), nodeveltags=tags.copy())
    exa = DisNetManager(G.copy())
    driver = ex.SimulateNetwork(state=st_ex, burgmag=burgmag, exadis_plastic_strain=False)
    driver.save_old_nodes(exa, st_ex)
    ex.TimeIntegration(integrator='EulerForward', dt=DT, state=st_ex).Update(exa, st_ex)
    E, W, rho = exa.get_disnet(ex.ExaDisNet).net.get_plastic_strain()
    cpp = (sym_to_voigt(np.array(E).reshape(3, 3)),
           asym_to_array(np.array(W).reshape(3, 3)), rho)
    driver.plastic_strain(exa, st_ex)

    paths = {"PyDiS": pydis_strain(py, oldpos, burgmag),
             "ExaDiS python driver": (st_ex["dEp"], st_ex["dWp"], st_ex["density"]),
             "ExaDiS C++": cpp,
             "PyDiS formula on the ExaDiS network": pydis_strain(exa, oldpos, burgmag)}
    return positions(py), positions(exa), paths


def check_one_step():
    """check_one_step: 1. one prescribed Euler step of random loops"""
    rng = np.random.default_rng(5)
    state = fr.default_state(LBOX)
    pos_py, pos_ex, rows = [], [], {}
    for _ in range(NLOOP):
        G = fr.random_loop(rng, box_length=LBOX)
        vels = rng.normal(scale=DISP/DT, size=(G.num_nodes(), 3))
        p_py, p_ex, paths = one_step(G, vels, state)
        pos_py += [p_py[t] for t in sorted(p_py)]
        pos_ex += [p_ex[t] for t in sorted(p_py)]
        scale = np.abs(paths["PyDiS"][0]).max()
        for name, (dEp, dWp, rho) in paths.items():
            rows.setdefault(name, []).append(
                np.concatenate((dEp / scale, dWp / scale, [rho / paths["PyDiS"][2]])))
    ok = report_close("one step: positions moved by ExaDiS vs PyDiS",
                      np.array(pos_ex), np.array(pos_py), atol=STEP_POS_TOL)
    ref = np.array(rows.pop("PyDiS"))
    for name, values in rows.items():
        ok &= report_close("one step: %s vs PyDiS, dEp, dWp, density, relative" % name,
                           np.array(values), ref, atol=STEP_TOL)
    return ok


def refined_source():
    """refined_source: the A3 source, its free arm bisected until no segment exceeds maxseg"""
    ss = FCC_slip_sys()
    net = fr.init_frank_read_src_on_slip_system(ss, ss.name_to_sysID('A3'),
                                                arm_length=ARM, box_length=LBOX)
    state = fr.default_state(LBOX)
    remesh = pydis.Remesh(remesh_rule='LengthBased', state=state)
    for _ in range(10):
        n = net.num_nodes()
        remesh.Remesh(net, state)
        if net.num_nodes() == n:
            break
    return net.get_disnet(DisNet)


def record(history, DM, state):
    """record: this step's state and positions"""
    history.append({k: np.copy(state[k]) for k in
                    ("strain", "stress", "dEp", "dWp", "density", "edir")})
    history[-1]["pos"] = positions(DM)


class PyDiSRecorded(pydis.SimulateNetwork):
    """PyDiSRecorded: PyDiS' driver, keeping each step's state"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.history = []

    def step_end(self, DM, state):
        record(self.history, DM, state)


class ExaDiSRecorded(ex.SimulateNetwork):
    """ExaDiSRecorded: ExaDiS' driver, keeping each step's state"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.history = []

    def step_end(self, N, state):
        record(self.history, N, state)


def run(code, G0, rotation):
    """run: NSTEP steps at ERATE along [001] in code ('pydis' or 'exadis'); returns the history"""
    lib, sim_class = (pydis, PyDiSRecorded) if code == 'pydis' else (ex, ExaDiSRecorded)
    state = fr.default_state(LBOX)
    net = DisNetManager(G0.copy())
    sim = sim_class(state=state, loading_mode='strain_rate', erate=ERATE, edir=EZZ,
                    rotation=rotation, burgmag=state["burgmag"], max_step=NSTEP, write_dir='',
                    calforce=lib.CalForce(force_mode='LineTension', state=state),
                    mobility=lib.MobilityLaw(mobility_law='SimpleGlide', state=state),
                    timeint=lib.TimeIntegration(integrator='EulerForward', dt=DT, state=state))
    sim.run(net, state)
    return sim.history


def check_run(G0, rotation):
    """check_run: 2. the A3 strain-rate run, step by step, PyDiS against ExaDiS"""
    label = "run, axis turned by the plastic spin" if rotation else "run, axis fixed"
    hp, he = run('pydis', G0, rotation), run('exadis', G0, rotation)
    ok = report("%s: %d steps in both codes, same nodes throughout" % (label, NSTEP),
                len(hp) == len(he) == NSTEP
                and all(a["pos"].keys() == b["pos"].keys() for a, b in zip(hp, he)))
    if not ok:
        return ok

    def series(h, key, norm=None):
        return np.array([x[key] / (1.0 if norm is None else norm(x)) for x in h])

    def dEp_max(x):
        return np.abs(x["dEp"]).max()

    tags = sorted(hp[0]["pos"])
    pos = [np.array([[x["pos"][t] for t in tags] for x in h]) for h in (hp, he)]
    ok &= report_close("%s: positions" % label, pos[1], pos[0], atol=RUN_TOL["positions"])
    for key in ("strain", "stress", "density"):
        ok &= report_close("%s: %s, relative" % (label, key),
                           series(he, key) / series(hp, key), np.ones(NSTEP),
                           atol=RUN_TOL[key])
    for key in ("dEp", "dWp"):
        # PyDiS' max |dEp| of each step as the scale, or 1 where nothing moved
        scale = np.array([dEp_max(x) for x in hp])[:, None]
        scale = np.where(scale > 0.0, scale, 1.0)
        ok &= report_close("%s: %s, relative to each step's max |dEp|" % (label, key),
                           series(he, key) / scale, series(hp, key) / scale,
                           atol=RUN_TOL[key])
    ok &= report_close("%s: loading axis" % label, series(he, "edir"), series(hp, "edir"),
                       atol=RUN_TOL["edir"])
    print("    final stress %.6f MPa in PyDiS, %.6f MPa in ExaDiS; plastic strain %.3e"
          % (hp[-1]["stress"]/1e6, he[-1]["stress"]/1e6,
             sum(np.dot(x["dEp"], [0, 0, 1, 0, 0, 0]) for x in hp)))
    return ok


def main():
    print("PyDiS' plastic strain against ExaDiS'")
    ok = check_one_step()
    G0 = refined_source()
    print("    refined A3 source: %d nodes" % G0.num_nodes())
    return ok & check_run(G0, rotation=False) & check_run(G0, rotation=True)


if __name__ == "__main__":
    pyexadis.initialize()
    passed = verdict("test_plastic_strain_pydis_exadis", main())
    pyexadis.finalize()
    sys.exit(0 if passed else 1)
