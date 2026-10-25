"""Test helpers for the plastic strain of this folder.

The formulas were ported here from calc_plastic_strain_from_disl_network
and update_applied_stress in disl_analysis_on_dict_util.py (this folder),
which follow ParaDiS DeltaPlasticStrain_FCC and LoadCurve loadType 1; they
now live in pydis.util.plastic_strain, which pydis SimulateNetwork uses, and
are imported below under the names the tests use. What stays here is the
network's area functional F, which the leak-free identity is checked
against, and the unwrapping it needs.

Voigt order throughout: [xx, yy, zz, yz, xz, xy]; spins as [yz, xz, xy].
"""

import numpy as np

from pydis.util.plastic_strain import (
    swept_area_vector, segment_dyads, sym_skew, old_positions_for,
    plastic_strain_increment, dislocation_density, sym_to_voigt,
    voigt_to_sym, asym_to_array, array_to_asym, spin_matrix,
    rotate_loading, uniaxial_projectors, uniaxial_increment)


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
