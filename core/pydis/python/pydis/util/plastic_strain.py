"""Plastic strain from swept area, for time integration and topology alike,
and the load update it drives.

A segment 1->2 with Burgers vector b, moving from (r1_old, r2_old) to
(r1, r2), sweeps the area vector

    Sn = 0.5 (r2 - r3) x (r1 - r4),   r3, r4 the old ends imaged near r1

and contributes P = outer(Sn, b) / V, so dEp += 0.5 (P + P^T) and
dWp += 0.5 (P - P^T). This is ParaDiS DeltaPlasticStrain and ExaDiS
System::plastic_strain.

A topological operation moves one node while its neighbours stay put. That
is the same formula with one end held fixed, which is ExaDiS
SerialDisNet::update_node_plastic_strain. Both paths here call
segment_dyads, so they agree bitwise. Unlike ExaDiS, the spin is
accumulated for topological moves as well.

The load update is ParaDiS LoadCurve, loadType 1, and ExaDiS'
SimulateNetwork.step_update_response (core/exadis/python/pyexadis_base.py).
Voigt order throughout: [xx, yy, zz, yz, xz, xy]; spins as [yz, xz, xy].
The spin follows ExaDiS; ParaDiS (after its 2015 sign fix) and MLmat store
-dWp and flip the rotation angles to match, so all three turn the loading
axis the same way.
"""

import numpy as np


def swept_area_vector(r1, r2, r3, r4):
    """swept_area_vector: area vector swept by segment 1->2 moving (r3, r4) -> (r1, r2)

    Half the cross product of the diagonals of the quadrilateral r3, r4, r2,
    r1. Rows may be stacked.
    """
    return 0.5 * np.cross(r2 - r3, r1 - r4)


def segment_dyads(cell, r1, r2, r1_old, r2_old, burgers):
    """segment_dyads: outer(Sn, b) / V for each segment, minimum image throughout

    One (3,) row per segment in every argument, b directed 1 -> 2, ends now
    at (r1, r2) and before at (r1_old, r2_old). Returns (N, 3, 3).
    """
    # asarray: an exadis cell returns lists
    r2 = np.asarray(cell.closest_image(Rref=r1, R=r2))
    r3 = np.asarray(cell.closest_image(Rref=r1, R=r1_old))
    r4 = np.asarray(cell.closest_image(Rref=r3, R=r2_old))
    Sn = swept_area_vector(r1, r2, r3, r4)
    return (1.0 / cell.volume()) * (Sn[:, :, None] * burgers[:, None, :])


def segment_positions(data):
    """segment_positions: (tags1, tags2, r1, r2, burgers), one row per segment

    data as export_data returns it, from a pydis DisNet, an ExaDisNet, or a
    DisNetManager, which exports whichever network is current. b is
    directed 1 -> 2, and r2 is as stored, not yet imaged near r1.
    """
    nodes, segs = data["nodes"], data["segs"]
    ids = np.asarray(segs["nodeids"], dtype=int).reshape(-1, 2)
    tags = np.asarray(nodes["tags"], dtype=int).reshape(-1, 2)
    pos = np.asarray(nodes["positions"], dtype=float).reshape(-1, 3)
    burgers = np.asarray(segs["burgers"], dtype=float).reshape(-1, 3)
    return tags[ids[:, 0]], tags[ids[:, 1]], pos[ids[:, 0]], pos[ids[:, 1]], burgers


def old_positions_for(tags, oldpos):
    """old_positions_for: rows of oldpos for tags; a node without one is an error"""
    missing = [tuple(t) for t in tags if tuple(t) not in oldpos]
    if missing:
        raise ValueError("no start-of-step position for nodes %s; was a node "
                         "added after save_old_nodes?" % missing[:5])
    return np.array([oldpos[tuple(t)] for t in tags])


def sym_skew(P):
    """sym_skew: (0.5 (P + P^T), 0.5 (P - P^T))"""
    return 0.5 * (P + P.T), 0.5 * (P - P.T)


def plastic_strain_from_data(data, cell, oldpos):
    """plastic_strain_from_data: (dEp, dWp), 3x3, swept since the positions oldpos

    data as segment_positions takes it; oldpos maps node tag to its
    start-of-step position, as pydis.util.state_access.old_positions_by_tag
    returns it.
    """
    tags1, tags2, r1, r2, burgers = segment_positions(data)
    if burgers.shape[0] == 0:
        return np.zeros((3, 3)), np.zeros((3, 3))
    P = segment_dyads(cell, r1, r2, old_positions_for(tags1, oldpos),
                      old_positions_for(tags2, oldpos), burgers)
    return sym_skew(P.sum(axis=0))


def density_from_data(data, cell, burgmag):
    """density_from_data: total line length / V / burgmag^2, in 1/m^2 for burgmag in m"""
    _, _, r1, r2, _ = segment_positions(data)
    if r1.shape[0] == 0:
        return 0.0
    r2 = np.asarray(cell.closest_image(Rref=r1, R=r2))
    return np.linalg.norm(r2 - r1, axis=1).sum() / cell.volume() / burgmag**2


def plastic_strain_increment(N, oldpos):
    """plastic_strain_increment: plastic_strain_from_data of a network or DisNetManager N"""
    return plastic_strain_from_data(N.export_data(), N.cell, oldpos)


def dislocation_density(N, burgmag):
    """dislocation_density: density_from_data of a network or DisNetManager N"""
    return density_from_data(N.export_data(), N.cell, burgmag)


def nodes_move_dyad(G, moves):
    """nodes_move_dyad: summed P of every segment touching a moved node

    moves maps tag -> (pos_old, pos_new); every other node stays put. A
    segment joining two moved nodes is counted once, with both ends moving:
    counting one end at a time against the other's old position counts a
    collapsing segment twice, which is what ExaDiS' merge_nodes_position
    does. Segments are taken as stored, b from the first endpoint, exactly as
    DisNet.get_segs_data_with_positions presents them, so a segment gives the
    same bits here as in a time-integration sum.
    """
    P = np.zeros((3, 3))
    seen = set()
    for tag in moves:
        for nbr, attr, tag_first in G._arms_with_direction(tag):
            first, second = (tag, nbr) if tag_first else (nbr, tag)
            if (first, second) in seen:
                continue
            seen.add((first, second))
            old1, new1 = moves.get(first, (G.nodes(first).R,) * 2)
            old2, new2 = moves.get(second, (G.nodes(second).R,) * 2)
            rows = (new1, new2, old1, old2, attr.burg_vec_from(first))
            P += segment_dyads(G.cell, *(np.atleast_2d(r) for r in rows))[0]
    return P


class PlasticStrainAccum:
    """PlasticStrainAccum: running dEp, dWp of one step, as 3x3 tensors"""

    def __init__(self):
        self.dEp = np.zeros((3, 3))
        self.dWp = np.zeros((3, 3))

    def zero(self):
        self.dEp[...] = 0.0
        self.dWp[...] = 0.0

    def add(self, dEp, dWp):
        self.dEp += dEp
        self.dWp += dWp

    def add_dyad(self, P):
        self.dEp += 0.5 * (P + P.T)
        self.dWp += 0.5 * (P - P.T)

    def add_moves(self, G, moves):
        """add_moves: account for nodes moving, {tag: (pos_old, pos_new)}; see nodes_move_dyad"""
        self.add_dyad(nodes_move_dyad(G, moves))


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


def von_mises(stress):
    """von_mises: von Mises equivalent of a stress in Voigt form, as ExaDiS reports it"""
    S = voigt_to_sym(stress)
    Sdev = S - np.trace(S)/3.0*np.eye(3)
    return np.sqrt(3.0/2.0*np.dot(Sdev.ravel(), Sdev.ravel()))
