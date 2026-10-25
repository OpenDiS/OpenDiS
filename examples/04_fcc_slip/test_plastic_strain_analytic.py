"""Plastic strain from swept area: closed form, PBC, sign, rotation, and MLmat.

No force or mobility model is involved; positions are prescribed.

1. A straight infinite dislocation translated rigidly sweeps A = L dx, so
   dEp_xy = 0.5 A b / V and nothing else (Orowan). Inputs are exact binary
   values, so the formula is checked to the last few bits.
2. The same move with the line next to a box face, folded across it as the
   integrator folds, gives the same result.
3. Reversing the move negates dEp and dWp bitwise.
4. spin_matrix equals scipy's extrinsic xyz rotation; rotate_loading keeps a
   uniaxial stress uniaxial along the turned axis.
5. On random loops crossing box faces, the port matches MLmat's own
   calc_plastic_strain_from_disl_network, including its rotation angles.

    python3 test_plastic_strain_analytic.py
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from scipy.spatial.transform import Rotation
from framework.testing import report, report_close, verdict
from pydis import DisNode, DisNet, Cell
import plast_strain_util as psu
import fcc_frank_read_util as fr
import disl_analysis_on_dict_util as mlmat

FREE = DisNode.Constraints.UNCONSTRAINED
LBOX = 1000.0


def network(cell, pos, burg):
    """network: closed loop through pos with Burgers vector burg"""
    N = pos.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        dr = cell.closest_image(Rref=pos[i], R=pos[(i+1) % N]) - pos[i]
        pn = np.cross(burg, dr)
        links[i] = np.concatenate(([i, (i+1) % N], burg, pn/np.linalg.norm(pn)))
    rn = np.column_stack((pos, np.full(N, FREE)))
    return DisNet(cell=cell, rn=rn, links=links)


def snapshot(G):
    """snapshot: tag -> position, as save_old_nodes records it"""
    return {tag: G.nodes(tag).R.copy() for tag in G.all_nodes_tags()}


def move(G, disp):
    """move: displace every node by disp (one row per node), folding as EulerForward"""
    for tag, d in zip(list(G.all_nodes_tags()), disp):
        G.nodes(tag).R = G.cell.fold(G.nodes(tag).R + d)


def straight_line(x0, nseg=8):
    """straight_line: infinite edge line along z at x = x0, b = x, plane normal y"""
    cell = Cell(h=LBOX*np.eye(3), is_periodic=[True, True, True])
    z = np.arange(nseg) * (LBOX / nseg)
    pos = np.column_stack((np.full(nseg, x0), np.full(nseg, 512.0), z))
    return network(cell, pos, np.array([1.0, 0.0, 0.0]))


def translated_increment(x0, dx):
    """translated_increment: (dEp, dWp) for the straight line moved by dx along x"""
    G = straight_line(x0)
    old = snapshot(G)
    move(G, np.tile([dx, 0.0, 0.0], (G.num_nodes(), 1)))
    return psu.plastic_strain_increment(G, old)


def check_closed_form():
    dx = 4.0
    dEp, dWp = translated_increment(512.0, dx)
    gamma = LBOX * dx * 1.0 / LBOX**3
    want = np.zeros((3, 3))
    want[0, 1] = want[1, 0] = 0.5 * gamma
    ok = report_close("rigid translation: dEp is 0.5 A b / V in xy only",
                      dEp / gamma, want / gamma, atol=1e-15)
    ok &= report_close("rigid translation: dWp array [yz, xz, xy] = [0, 0, -0.5 A b/V]",
                       psu.asym_to_array(dWp) / gamma, [0.0, 0.0, -0.5], atol=1e-15)
    return ok


def check_pbc_fold():
    """the line next to the x face, folded across it by the move"""
    dEp_a, dWp_a = translated_increment(512.0, 4.0)
    dEp_b, dWp_b = translated_increment(998.0, 4.0)
    return (report_close("folded across a box face: dEp unchanged",
                         dEp_b, dEp_a, atol=0.0)
            & report_close("folded across a box face: dWp unchanged",
                           dWp_b, dWp_a, atol=0.0))


def check_sign_reversal():
    dEp_p, dWp_p = translated_increment(512.0, 4.0)
    dEp_m, dWp_m = translated_increment(512.0, -4.0)
    return (report_close("reversed move: dEp negated", dEp_m, -dEp_p, atol=0.0)
            & report_close("reversed move: dWp negated", dWp_m, -dWp_p, atol=0.0))


def check_rotation():
    rng = np.random.default_rng(7)
    worst = max(np.abs(psu.spin_matrix(phi)
                       - Rotation.from_euler('xyz', phi).as_matrix()).max()
                for phi in rng.normal(scale=0.5, size=(50, 3)))
    ok = report("spin_matrix == scipy from_euler('xyz'): max error %.1e" % worst,
                worst < 1e-15)
    edir = np.array([0.0, 0.0, 1.0])
    sigma = 3.0e8
    S0 = sigma * psu.sym_to_voigt(np.outer(edir, edir))
    dWp = np.array([2e-3, -1e-3, 5e-4])
    e1, S1 = psu.rotate_loading(edir, S0, dWp)
    ok &= report_close("rotated stress stays sigma e'(x)e'",
                       S1 / sigma, psu.sym_to_voigt(np.outer(e1, e1)), atol=1e-15)
    R = Rotation.from_euler('xyz', [1.0, -1.0, 1.0] * dWp).as_matrix()
    ok &= report_close("rotated axis is R(+yz, -xz, +xy) e", e1, R @ edir, atol=1e-15)
    return ok


def to_mlmat_nodes(G, positions):
    """to_mlmat_nodes: MLmat's tag -> [position + constraint, {nbr: b + plane}]"""
    nodes = {}
    for tag in G.all_nodes_tags():
        arms = {nbr: np.concatenate((e.burg_vec_from(tag), e.plane_normal))
                for nbr, e in G.neighbor_segments_dict(tag).items()}
        nodes[tag] = [np.append(positions[tag], G.nodes(tag).constraint), arms]
    return nodes


def mlmat_increment(G, old, new):
    """mlmat_increment: (d_pl_strain_sum, d_pl_rotation_sum, phispn) from MLmat"""
    p = mlmat.oprec_player(save_data_path=None)
    p.nodes, p.nodes_old = to_mlmat_nodes(G, new), to_mlmat_nodes(G, old)
    p.ctrl_info = {}
    p.calc_plastic_strain_from_disl_network(np.full(3, LBOX))
    W = p.d_pl_rotation_sum
    return p.d_pl_strain_sum, W, np.array([-W[1, 2], W[0, 2], -W[0, 1]])


def check_against_mlmat():
    rng = np.random.default_rng(11)
    err_E = err_W = err_phi = 0.0
    for _ in range(20):
        G = fr.random_loop(rng)
        old = snapshot(G)
        move(G, rng.normal(scale=40.0, size=(G.num_nodes(), 3)))
        dEp, dWp = psu.plastic_strain_increment(G, old)
        E_ml, W_ml, phi_ml = mlmat_increment(G, old, snapshot(G))
        scale = np.abs(dEp).max()
        err_E = max(err_E, np.abs(dEp - E_ml).max() / scale)
        err_W = max(err_W, np.abs(dWp + W_ml).max() / scale)
        err_phi = max(err_phi, np.abs([1, -1, 1] * psu.asym_to_array(dWp) - phi_ml).max() / scale)
    ok = report("vs MLmat on 20 random loops: dEp, max relative error %.1e" % err_E,
                err_E < 1e-12)
    ok &= report("vs MLmat: dWp == -d_pl_rotation_sum, max relative error %.1e" % err_W,
                 err_W < 1e-12)
    ok &= report("vs MLmat: rotation angles agree, max relative error %.1e" % err_phi,
                 err_phi < 1e-12)
    return ok


def main():
    print("plastic strain from swept area, examples/04_fcc_slip")
    checks = [check_closed_form(), check_pbc_fold(), check_sign_reversal(),
              check_rotation(), check_against_mlmat()]
    return all(checks)


if __name__ == "__main__":
    passed = verdict("test_plastic_strain_analytic", main())
    sys.exit(0 if passed else 1)
