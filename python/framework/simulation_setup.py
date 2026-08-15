"""@package docstring
simulation_setup: pre-simulation sanity checks and configuration conditioning

These helpers run on plain (rn, links) arrays, before either a PyDiS DisNet or an ExaDiS
ExaDisNet is constructed, so a single implementation serves both codes and nothing under
core/exadis/ has to change.

Two things are provided, both addressing the same underlying requirement: the segment-pair
cutoff used by PyDiS' Elasticity_* modes and by ExaDiS' CUTOFF_MODEL is only exact when the
simulation parameters and the configuration are mutually consistent.

1. check_cutoff_maxseg() -- verifies  cutoff + maxseg <= d/3  for every cell direction.

   ExaDiS bins segments by their mid-point to build the neighbor list, and the number of bins
   per direction is floor(d / (cutoff + maxseg)), clamped up to a minimum of 3
   (core/exadis/src/neighbor_types/neighbor_box.h). Once that clamp engages, the bins are
   smaller than the search radius, and the +-1 bin scan compares mid-points using a periodic
   shift quantized to the bin grid instead of the true minimum image. Pairs whose real
   separation is inside the cutoff are then silently dropped. Keeping cutoff + maxseg <= d/3
   keeps the bin count at 3 or more by construction, where the algorithm is exact.

2. remesh_initial_config() -- subdivides every segment longer than maxseg.

   The neighbor-list bound cutoff + maxseg is only valid when no segment is longer than
   maxseg. LengthBased remesh does NOT guarantee this: it refines a segment only if at least
   one endpoint is free, in both codes (PyDiS remesh_disnet.py, ExaDiS remesh.h), so a
   segment pinned at both ends is never split however long it is. Conditioning the initial
   configuration closes that gap for the whole run, because a pinned-pinned segment never
   moves and so never grows back.

   Segments pinned at both ends are subdivided here, with the inserted nodes marked PINNED so
   that the geometry and the boundary condition are both preserved exactly.

Note that subdividing a segment changes results slightly: the far field is unchanged, but the
non-singular self force depends on segment length, and node/segment counts change.
"""

import numpy as np

UNCONSTRAINED = 0
PINNED_NODE = 7


def cell_widths(h):
    """cell_widths: perpendicular width of the cell along each direction

    Matches the quantity ExaDiS bins against, |cbox[i] . perpVecs[i]|
    (core/exadis/src/neighbor_types/neighbor_box.h). For a cubic cell h = L*I this is
    (L, L, L). h is column-based, i.e. h[:, i] is the i-th cell vector.
    """
    h = np.asarray(h, dtype=float)
    widths = np.zeros(3)
    for i in range(3):
        e1, e2 = h[:, (i+1) % 3], h[:, (i+2) % 3]
        perp = np.cross(e1, e2)
        n = np.linalg.norm(perp)
        if n < 1e-30:
            raise ValueError("cell_widths: degenerate cell")
        widths[i] = abs(np.dot(h[:, i], perp/n))
    return widths


def check_cutoff_maxseg(h, cutoff, maxseg, raise_on_fail=True, verbose=True):
    """check_cutoff_maxseg: verify cutoff + maxseg <= d/3 for every cell direction

    Returns True if the condition holds. If it does not, raises ValueError (default) or
    returns False when raise_on_fail=False. cutoff=None is treated as "no cutoff" and passes.
    """
    if cutoff is None:
        if verbose:
            print("check_cutoff_maxseg: no cutoff, nothing to check")
        return True

    widths = cell_widths(h)
    eff = cutoff + maxseg
    nbin = np.floor(widths/eff).astype(int)
    ok = bool(np.all(nbin >= 3))

    if verbose or not ok:
        print("check_cutoff_maxseg: cutoff = %g, maxseg = %g, cutoff+maxseg = %g" %
              (cutoff, maxseg, eff))
        print("                     cell widths = %s, allowed max = %s" %
              (np.array2string(widths, precision=2),
               np.array2string(widths/3.0, precision=2)))
        print("                     neighbor bins per direction = %s (must be >= 3)" %
              np.array2string(nbin))
    if not ok:
        msg = ("check_cutoff_maxseg: cutoff + maxseg = %g exceeds d/3 = %g for at least one "
               "cell direction. The ExaDiS neighbor list silently drops segment pairs in this "
               "regime; reduce cutoff or maxseg, or enlarge the cell." %
               (eff, float(np.min(widths)/3.0)))
        if raise_on_fail:
            raise ValueError(msg)
        print("WARNING: " + msg)
    elif verbose:
        print("                     OK")
    return ok


def remesh_initial_config(rn, links, maxseg, h, is_periodic=(True, True, True), verbose=True):
    """remesh_initial_config: subdivide every segment longer than maxseg

    rn    : (N,3) or (N,4) array; column 3, if present, is the node constraint
    links : (M,K) array; columns 0,1 are the node indices, columns 2: are per-segment
            attributes (Burgers vector, glide plane, ...) inherited by the sub-segments
    maxseg: no segment longer than this survives
    h     : (3,3) cell matrix, column-based
    Returns (rn, links) with the same column layout.

    A segment is split into ceil(L/maxseg) equal pieces, L being measured under the minimum
    image convention. Inserted nodes are PINNED if and only if both endpoints of the parent
    segment are PINNED, so a fixed segment stays fixed.
    """
    rn = np.asarray(rn, dtype=float).copy()
    links = np.asarray(links, dtype=float).copy()
    h = np.asarray(h, dtype=float)
    hinv = np.linalg.inv(h)
    pbc = np.asarray(is_periodic, dtype=bool)
    has_constraint = rn.shape[1] >= 4

    def min_image(dr):
        s = hinv @ dr
        s = np.where(pbc, s - np.rint(s), s)
        return h @ s

    def fold(r):
        s = hinv @ r
        s = np.where(pbc, s - np.floor(s), s)
        return h @ s

    pts = [rn[i].copy() for i in range(rn.shape[0])]
    segs = []
    n_split = 0
    for k in range(links.shape[0]):
        n1, n2 = int(links[k, 0]), int(links[k, 1])
        attrs = links[k, 2:]
        p1 = pts[n1][:3]
        dr = min_image(rn[n2, :3] - p1)
        length = np.linalg.norm(dr)
        nsub = int(np.ceil(length/maxseg - 1.0e-12))
        if nsub <= 1:
            segs.append(np.concatenate(([n1, n2], attrs)))
            continue

        n_split += 1
        if has_constraint:
            both_pinned = (rn[n1, 3] == PINNED_NODE) and (rn[n2, 3] == PINNED_NODE)
            new_con = PINNED_NODE if both_pinned else UNCONSTRAINED
        prev = n1
        for m in range(1, nsub):
            r_new = fold(p1 + dr*(float(m)/nsub))
            node = np.zeros(rn.shape[1])
            node[:3] = r_new
            if has_constraint:
                node[3] = new_con
            pts.append(node)
            segs.append(np.concatenate(([prev, len(pts)-1], attrs)))
            prev = len(pts)-1
        segs.append(np.concatenate(([prev, n2], attrs)))

    rn_out = np.array(pts)
    links_out = np.array(segs)
    if verbose:
        print("remesh_initial_config: maxseg = %g, split %d of %d segments; "
              "nodes %d -> %d, segments %d -> %d" %
              (maxseg, n_split, links.shape[0], rn.shape[0], rn_out.shape[0],
               links.shape[0], links_out.shape[0]))
    return rn_out, links_out
