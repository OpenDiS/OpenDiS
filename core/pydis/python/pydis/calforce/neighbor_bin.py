"""neighbor_bin: ExaDiS' segment neighbor bins, as its default node_force uses them

ExaDiS' serial node_force, with match_global left False, takes the segments it
pairs with a node's arms from one bin query at the node position
(core/exadis/src/force_types/force_segseglist.h). The query is not a pure
cutoff: it keeps every segment whose folded mid-point lies within cutoff +
maxseg of the node, as reached by NeighborBin's walk over the +-1 neighbour bins
(core/exadis/src/neighbor_types/neighbor_bin.cpp). This module reproduces that
walk, so pydis' OneNodeForce can select the same pairs.

The walk is reproduced as written rather than reduced to a minimum-image
distance test, because the two part where the bins do: with one or two bins in
a direction a segment can be reached through more than one image, and ExaDiS
then counts it more than once.
"""

import numpy as np

MAX_BINS = 50


def bin_layout(cell, radius):
    """bin_layout: (bins per direction, inverse bin matrix) for this cell

    As NeighborBin_t::initialize: floor(cell width / radius) bins per
    direction, clamped to [1, MAX_BINS]. The columns of cell.h are the box
    vectors, as in ExaDiS' H.
    """
    cbox = [cell.h[:, i] for i in range(3)]
    perp = [np.cross(cbox[1], cbox[2]), np.cross(cbox[2], cbox[0]),
            np.cross(cbox[0], cbox[1])]
    dim = np.array([min(max(int(np.floor(abs(np.dot(cbox[i], p/np.linalg.norm(p)))
                                         / radius)), 1), MAX_BINS)
                    for i, p in enumerate(perp)])
    bin_h = np.column_stack([cbox[i]/dim[i] for i in range(3)])
    return dim, np.linalg.inv(bin_h)


def bin_coords(cell, dim, bin_hinv, points):
    """bin_coords: integer bin coordinates, wrapped if periodic, else clamped"""
    c = np.floor(np.dot(bin_hinv, (np.atleast_2d(points) - cell.origin).T).T)
    c = c.astype(int)
    for k in range(3):
        if cell.is_periodic[k]:
            c[:, k] %= dim[k]
        else:
            c[:, k] = np.clip(c[:, k], 0, dim[k] - 1)
    return c


def flat_index(dim, c):
    """flat_index: ExaDiS' bin index, x fastest"""
    return c[..., 2]*dim[0]*dim[1] + c[..., 1]*dim[0] + c[..., 0]


def walk(cell, dim, center_bin):
    """walk: (bin index, image shift) for each bin a query visits, in order

    The +-1 neighbourhood with the z offset fastest, as iterator::next steps
    it. A periodic direction wraps and adds a box vector to the shift; a
    non-periodic one skips the bin.
    """
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                current, shift = center_bin + np.array([dx, dy, dz]), np.zeros(3)
                for k in range(3):
                    if current[k] < 0 or current[k] >= dim[k]:
                        if not cell.is_periodic[k]:
                            break
                        wraps = current[k] if current[k] < 0 else current[k] - dim[k] + 1
                        shift += wraps*cell.h[:, k]
                        current[k] = dim[k] - 1 if current[k] < 0 else 0
                else:
                    yield int(flat_index(dim, current)), shift


def query(cell, mid, center, radius):
    """query: indices of the segments a query at center returns, in order

    mid holds every segment's folded mid-point, in segment order. Within a bin
    the most recently inserted segment comes first, ExaDiS keeping each bin as
    a linked list built by insertion at the head. A segment is kept when its
    shifted mid-point lies within radius of center, and, when any direction has
    a single bin, only through the image nearest center.
    """
    dim, bin_hinv = bin_layout(cell, radius)
    seg_bin = flat_index(dim, bin_coords(cell, dim, bin_hinv, mid))
    center_bin = bin_coords(cell, dim, bin_hinv, center)[0]
    unique_image = bool(np.any(dim == 1))

    found = []
    for b, shift in walk(cell, dim, center_bin):
        members = np.where(seg_bin == b)[0][::-1]
        delta = mid[members] - center + shift
        d2 = delta[:, 0]*delta[:, 0] + delta[:, 1]*delta[:, 1] + delta[:, 2]*delta[:, 2]
        keep = d2 <= radius*radius
        if unique_image:
            keep &= np.all(np.abs(np.dot(cell.hinv, delta.T).T) < 0.5, axis=1)
        found.extend(members[keep].tolist())
    return found
