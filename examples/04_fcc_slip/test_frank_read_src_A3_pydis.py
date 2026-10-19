"""Frank-Read source on slip system A3 under constant uniaxial [001] stress.

The edge source of examples/02_frank_read_src, rotated onto
A3 = (-1,1,1)[101]/2, whose Schmid factor on [001] is 1/sqrt(6). At 400 MPa,
below the 650 MPa that test_uniaxial_A3_pydis.py reaches, the source must
operate: its free arm grows past the critical semicircle, of length
pi L / 2, and every segment that moves stays on A3.

    python3 test_frank_read_src_A3_pydis.py
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

LBOX, ARM = 1000.0, 500.0
SIGMA, NSTEP = 4.0e8, 200


def check_frame(slip_sys, sysID):
    Q = fr.slip_frame(slip_sys, sysID)
    return (report_close("slip frame orthonormal", Q.T @ Q, np.eye(3), atol=1e-15)
            & report_close("slip frame right-handed, det Q", np.linalg.det(Q), 1.0,
                           atol=1e-15))


def check_initial_labels(net, slip_sys):
    G = net.get_disnet()
    labels = sorted(int(slip_sys.nb_to_sysID(e.plane_normal, e.burg_vec_from(t1))[0])
                    for (t1, t2), e in G.all_segments_dict().items())
    return (report("initial network sane", net.is_sane())
            & report("segments: three on A3, two side arms off {111}",
                     labels == [1, 1, 1, 12, 12], suffix=' %s' % labels)
            & report("free arm on A3", fr.free_arm_sysIDs(G, slip_sys) == {1}))


def run(net, state, slip_sys):
    """run: NSTEP constant-stress steps; returns the free-arm labels seen"""
    sim = SimulateNetwork(state=state, max_step=NSTEP, loading_mode='stress',
                          applied_stress=fr.uniaxial_voigt(SIGMA, [0, 0, 1]),
                          **fr.dd_components(net, state))
    seen = set()
    for istep in range(NSTEP):
        state['istep'] = istep
        sim.step(net, state)
        if istep % 10 == 0:
            seen |= fr.free_arm_sysIDs(net.get_disnet(), slip_sys)
            print("step %3d: free length %7.1f, %d nodes"
                  % (istep, fr.free_line_length(net.get_disnet()),
                     net.get_disnet().num_nodes()))
    return seen


def main():
    print("Frank-Read source on A3, constant %.0f MPa along [001]" % (SIGMA/1e6))
    ss = FCC_slip_sys()
    a3 = ss.name_to_sysID('A3')
    net = fr.init_frank_read_src_on_slip_system(ss, a3, arm_length=ARM, box_length=LBOX)
    ok = check_frame(ss, a3) & check_initial_labels(net, ss)

    seen = run(net, fr.default_state(LBOX), ss)
    length = fr.free_line_length(net.get_disnet())
    ok &= report("free arm stayed on A3 throughout", seen == {1}, suffix=' %s' % seen)
    ok &= report("source operates: free length %.1f > pi L/2 = %.1f"
                 % (length, np.pi*ARM/2), length > np.pi*ARM/2)
    ok &= report("final network sane", net.is_sane())

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_A3_final.json')
    return ok


if __name__ == "__main__":
    passed = verdict("test_frank_read_src_A3_pydis", main())
    sys.exit(0 if passed else 1)
