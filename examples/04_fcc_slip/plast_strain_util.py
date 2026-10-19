"""Plastic strain from swept area, and the uniaxial strain-rate load update.

Ported from calc_plastic_strain_from_disl_network and update_applied_stress
in disl_analysis_on_dict_util.py (this folder), which follow ParaDiS
DeltaPlasticStrain_FCC and LoadCurve loadType 1.

The spin follows ExaDiS, dWp = 0.5 (P - P^T). MLmat and ParaDiS (after its
2015 sign fix) store -dWp and flip the rotation angles to match, so all
three turn the loading axis the same way.

Voigt order throughout: [xx, yy, zz, yz, xz, xy]; spins as [yz, xz, xy].
"""

import numpy as np

# the one swept-area formula, shared with topological operations
from pydis.util.plastic_strain import swept_area_vector, segment_dyads


def sym_skew(P):
    """sym_skew: (0.5 (P + P^T), 0.5 (P - P^T))"""
    return 0.5 * (P + P.T), 0.5 * (P - P.T)


def old_positions_for(tags, oldpos):
    """old_positions_for: rows of oldpos for tags; a node without one is an error"""
    missing = [tuple(t) for t in tags if tuple(t) not in oldpos]
    if missing:
        raise ValueError("no start-of-step position for nodes %s; was a node "
                         "added after save_old_nodes?" % missing[:5])
    return np.array([oldpos[tuple(t)] for t in tags])


def plastic_strain_increment(G, oldpos):
    """plastic_strain_increment: (dEp, dWp) swept since the positions oldpos

    oldpos maps node tag to its start-of-step position, as
    pydis.util.state_access.old_positions_by_tag returns it. 3x3 tensors.
    """
    segs = G.get_segs_data_with_positions()
    if segs["burgers"].shape[0] == 0:
        return np.zeros((3, 3)), np.zeros((3, 3))
    P = segment_dyads(G.cell, segs["R1"], segs["R2"],
                      old_positions_for(segs["tag1"], oldpos),
                      old_positions_for(segs["tag2"], oldpos), segs["burgers"])
    return sym_skew(P.sum(axis=0))


def dislocation_density(G, burgmag):
    """dislocation_density: total line length / V / burgmag^2, in 1/m^2"""
    segs = G.get_segs_data_with_positions()
    length = np.linalg.norm(segs["R2"] - segs["R1"], axis=1).sum()
    return length / G.cell.volume() / burgmag**2


def sym_to_voigt(T):
    """sym_to_voigt: symmetric 3x3 -> [xx, yy, zz, yz, xz, xy]"""
    return np.array([T[0, 0], T[1, 1], T[2, 2], T[1, 2], T[0, 2], T[0, 1]])


def voigt_to_sym(v):
    """voigt_to_sym: [xx, yy, zz, yz, xz, xy] -> symmetric 3x3"""
    return np.array([[v[0], v[5], v[4]],
                     [v[5], v[1], v[3]],
                     [v[4], v[3], v[2]]])


def asym_to_array(T):
    """asym_to_array: antisymmetric 3x3 -> [yz, xz, xy]"""
    return np.array([T[1, 2], T[0, 2], T[0, 1]])


def array_to_asym(w):
    """array_to_asym: [yz, xz, xy] -> antisymmetric 3x3"""
    return np.array([[0.0, w[2], w[1]],
                     [-w[2], 0.0, w[0]],
                     [-w[1], -w[0], 0.0]])


def spin_matrix(phi):
    """spin_matrix: Rz(phi[2]) Ry(phi[1]) Rx(phi[0]), as ParaDiS SpinMatrix

    With exact sin and cos, where ParaDiS uses sin(a) = a, cos(a) = 1 - a^2/2.
    """
    Ss, Cs = np.sin(phi), np.cos(phi)
    R = np.zeros((3, 3))
    R[0, 0] = Cs[2]*Cs[1]
    R[1, 1] = Cs[2]*Cs[0] + Ss[2]*Ss[1]*Ss[0]
    R[2, 2] = Cs[1]*Cs[0]
    R[0, 1] = -Ss[2]*Cs[0] + Cs[2]*Ss[0]*Ss[1]
    R[1, 2] = -Cs[2]*Ss[0] + Ss[2]*Ss[1]*Cs[0]
    R[2, 0] = -Ss[1]
    R[0, 2] = Ss[2]*Ss[0] + Cs[2]*Cs[0]*Ss[1]
    R[1, 0] = Ss[2]*Cs[1]
    R[2, 1] = Cs[1]*Ss[0]
    return R


def rotate_loading(edir, applied_stress, dWp):
    """rotate_loading: carry edir and the applied stress through the plastic spin

    dWp as [yz, xz, xy]; angles [1, -1, 1] * dWp, as ExaDiS. The stress turns
    with the axis, as in ExaDiS; ParaDiS LoadCurve turns only the axis.
    """
    R = spin_matrix(np.array([1.0, -1.0, 1.0]) * dWp)
    e = R @ edir
    S = voigt_to_sym(applied_stress)
    return e / np.linalg.norm(e), sym_to_voigt(R @ S @ R.T)


def uniaxial_projectors(edir):
    """uniaxial_projectors: (A, A2), Voigt e (x) e without and with doubled shears

    dot(v, A2) is the component along edir of a stress or strain in Voigt form.
    """
    A = sym_to_voigt(np.outer(edir, edir))
    return A, A * np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


def uniaxial_increment(dEp, dt, erate, mu, nu, edir):
    """uniaxial_increment: (dstress, dstrain) along edir for one step, loadType 1

    dEp in Voigt form; dstress = E (erate dt - e.dEp.e), E = 2 mu (1 + nu).
    """
    _, A2 = uniaxial_projectors(edir)
    dstrain = erate * dt
    Eyoung = 2.0 * mu * (1.0 + nu)
    return Eyoung * (dstrain - np.dot(dEp, A2)), dstrain


def unwrapped_pieces(G):
    """unwrapped_pieces: one {tag: position} per connected piece, laid out unwrapped

    Breadth-first from one node of each piece, every neighbour at its image
    nearest the node it is reached from.
    """
    pieces, seen = [], set()
    for root in G.all_nodes_tags():
        if root in seen:
            continue
        pos = {root: G.nodes(root).R.copy()}
        queue = [root]
        while queue:
            t = queue.pop()
            for nbr in G.neighbors_tags(t):
                if nbr not in pos:
                    pos[nbr] = G.cell.closest_image(Rref=pos[t], R=G.nodes(nbr).R)
                    queue.append(nbr)
        seen |= pos.keys()
        pieces.append(pos)
    return pieces


def unwrapped_positions(G):
    """unwrapped_positions: tag -> position over all pieces, or None if one wraps the box

    A piece wraps when some segment disagrees with its own minimum image.
    """
    pos = {}
    for piece in unwrapped_pieces(G):
        pos.update(piece)
    tol = 1e-6 * np.abs(G.cell.h).max()
    for t1, t2 in G.all_segments_dict():
        if np.abs(pos[t2] - pos[t1] - G.seg_vector(t1, t2)).max() > tol:
            return None
    return pos


def largest_extent(G):
    """largest_extent: widest Cartesian span of any one piece, unwrapped

    A piece narrower than the box by more than the collision distance cannot
    touch its own periodic image.
    """
    return max((np.ptp(np.array(list(p.values())), axis=0).max()
                for p in unwrapped_pieces(G)), default=0.0)


def network_area_dyad(G):
    """network_area_dyad: F = sum over segments of outer(0.5 r1 x r2, b) / V

    For a network conserving Burgers vector at every node, the P = dEp + dWp
    of any accounted sequence of moves is F(before) - F(after). Raises if a
    piece wraps the box, where F is undefined; and it stops meaning anything
    once a piece has met its own periodic image, see largest_extent.
    """
    pos = unwrapped_positions(G)
    if pos is None:
        raise ValueError("network_area_dyad: a piece of the network wraps the box")
    F = np.zeros((3, 3))
    for (t1, t2), e in G.all_segments_dict().items():
        F += np.outer(0.5 * np.cross(pos[t1], pos[t2]), e.burg_vec_from(t1))
    return F / G.cell.volume()
