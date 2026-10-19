"""Plastic strain from swept area, for time integration and topology alike.

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
    r2 = cell.closest_image(Rref=r1, R=r2)
    r3 = cell.closest_image(Rref=r1, R=r1_old)
    r4 = cell.closest_image(Rref=r3, R=r2_old)
    Sn = swept_area_vector(r1, r2, r3, r4)
    return (1.0 / cell.volume()) * (Sn[:, :, None] * burgers[:, None, :])


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
