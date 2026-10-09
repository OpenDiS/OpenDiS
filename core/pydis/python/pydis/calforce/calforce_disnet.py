"""@package docstring
CalForce_DisNet: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object

Cutoff convention for the Elasticity_* force modes
--------------------------------------------------
When CalForce is constructed with a 'cutoff', the segment-segment elastic interaction is
truncated using the same two-stage criterion as the CUTOFF_MODEL of exadis. A pair of
segments (i, j) contributes only if BOTH hold:

  1. the distance between the two segment MID-POINTS is <= cutoff + maxseg
  2. the minimum distance between the two SEGMENTS is < cutoff

Stage 1 exists because exadis bins segments by their mid-point to build the neighbor list
(src/neighbor_types/neighbor_box.h: p = 0.5*(r1+r2), and cutoff += params.maxseg for
NeiSeg), and stage 2 is the test actually applied when the segment-pair list is assembled
(src/force_types/force_segseglist.h, via get_min_dist2_segseg).

Stage 1 is exact only when no segment is longer than maxseg, which remesh normally
guarantees. Where that does not hold -- pinned edges that remesh never refines, for
instance -- stage 1 discards pairs whose true segment-segment distance is below the cutoff,
and the truncation is then slightly stronger than the nominal cutoff implies. That is
reproduced rather than corrected, so that both codes select the same pairs.

This criterion has no ParaDiS counterpart, and it is worth saying so before anyone goes
looking for one. ParaDiS does not truncate the pair sum at a radius at all: it splits near
from far by CELLS, building pair lists between a cell and its neighbours with CellPriority
deciding ownership (external/paradis/src/LocalSegForces.c), and handles everything beyond
through the fast multipole method (RemoteForceOneSeg, ibid.). Its rc is the core radius,
"a = param->rc", not a cutoff. The mid-point test and the widening by maxseg are a way of
bounding a neighbor search, and they are here for one reason: to make a pair-for-pair
comparison possible.

Everything below that reproduces this convention refers back to these paragraphs rather
than restating the provenance, so the comparison is documented in one place.

Distances are evaluated under the minimum image convention. The largest separation
attainable that way in a cubic cell of side L is the half diagonal sqrt(3)/2*L, so any
cutoff at or above that value reproduces the untruncated result. cutoff=None (the default)
skips both stages and keeps every pair.
"""

import numpy as np
from typing import Tuple
from ..disnet import DisNet, Tag
from ..nbrlist.nbrlist import CellList
from . import neighbor_bin
from framework.disnet_manager import DisNetManager
from framework.calforce_base import CalForce_Base

try:
    from .compute_stress_force_analytic_paradis import compute_selfforce_list
    from .compute_stress_force_analytic_paradis import compute_segseg_force_list, compute_segseg_force
    from .compute_stress_force_analytic_paradis import compute_segseg_force_batch
    from .compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_vec, compute_segseg_force_SBN1
    from .compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_SBA
    from .compute_stress_analytic_paradis       import compute_seg_stress_coord_dep, compute_seg_stress_coord_indep
except ImportError:
    # No compiled library: the self force falls back to its numpy twin, which
    # computes the same expressions in the same grouping. The other kernels
    # have no fallback wired yet.
    from .compute_stress_force_analytic_python import (
        python_selfforce_vec as compute_selfforce_list)
    # use python version instead
    # To do: put import commands here
    print("pydis_lib not found, using python version for force calculation")

from .compute_stress_force_analytic_python  import python_segseg_force_vec
from .bitrepro_math import (ENABLED as BITREPRO_MATH, cross3, dot3,
                            matvec3, norm3)

try:
    from ..collision.getmindist2_paradis import GetMinDist2_paradis as GetMinDist2
except ImportError:
    # No compiled library: the self force falls back to its numpy twin, which
    # computes the same expressions in the same grouping. The other kernels
    # have no fallback wired yet.
    from .compute_stress_force_analytic_python import (
        python_selfforce_vec as compute_selfforce_list)
    # use python version instead
    from ..collision.getmindist2_python  import GetMinDist2_python as GetMinDist2

_ZERO_VEL = np.zeros(3)

def min_dist2_segseg(p1, p2, p3, p4, connected, eps0=1.0e-12):
    """min_dist2_segseg: squared minimum distance between segments (p1,p2) and (p3,p4)

    The same quantity ParaDiS computes in MinSegSegDist
    (external/paradis/src/RetroactiveCollision2.c:219), written to match
    get_min_dist2_segseg() in exadis (src/functions.h) with hinge=0, so that stage 2
    of the cutoff agrees pair for pair. Return convention:
      - returns -1.0 if either segment is degenerate, so that a "dist2 >= 0" test rejects it
      - returns 0.0 for segments sharing a node, so that they are always within any cutoff
    The endpoints are expected to have been mapped to a common periodic image already.
    """
    if np.dot(p2-p1, p2-p1) < eps0 or np.dot(p4-p3, p4-p3) < eps0:
        return -1.0
    if connected:
        return 0.0
    return GetMinDist2(p1, _ZERO_VEL, p2, _ZERO_VEL, p3, _ZERO_VEL, p4, _ZERO_VEL)[0]

def min_dist2_segseg_vec(p1, p2, p3, p4, connected=None, degenerate=None,
                         eps0=1.0e-12, epsM=1e-6):
    """min_dist2_segseg_vec: vectorized squared minimum distance between segment pairs

    Batched form of min_dist2_segseg, following the same branch structure as
    GetMinDist2_python (pydis/collision/getmindist2_python.py) but computing only dist2,
    which is all the cutoff test needs. Every branch is evaluated for the whole batch and
    selected with where(), so there is no data-dependent control flow.

    p1..p4 are (N,3); connected and degenerate are optional (N,) boolean masks applied last,
    reproducing the conventions of min_dist2_segseg: 0.0 for segments sharing a node, -1.0
    for a degenerate segment.
    """
    p1, p2, p3, p4 = (np.asarray(x, dtype=float) for x in (p1, p2, p3, p4))
    r1mr3 = p1 - p3
    r2mr1 = p2 - p1
    r4mr3 = p4 - p3

    dot = lambda u, v: np.einsum('ij,ij->i', u, v)
    A  = dot(r2mr1, r2mr1)          # M[0,0]
    Cm = -dot(r4mr3, r2mr1)         # M[1,0] == M[0,1]
    E  = dot(r4mr3, r4mr3)          # M[1,1]
    rhs0 = -dot(r2mr1, r1mr3)
    rhs1 =  dot(r4mr3, r1mr3)

    def d2(L1, L2):
        v = p1 + r2mr1*L1[:, None] - p3 - r4mr3*L2[:, None]
        return dot(v, v)

    with np.errstate(divide='ignore', invalid='ignore'):
        detM = 1.0 - Cm*Cm/(A*E)
        detM2 = detM * A * E
        s0 = ( E*rhs0 - Cm*rhs1) / detM2
        s1 = (-Cm*rhs0 + A*rhs1) / detM2
        A_safe = np.where(A > eps0, A, 1.0)
        E_safe = np.where(E > eps0, E, 1.0)
        # the four clipped trials of the general out-of-range case
        t = [(np.zeros_like(A),                            np.clip(rhs1/E_safe, 0, 1)),
             (np.ones_like(A),                             np.clip((rhs1 - Cm)/E_safe, 0, 1)),
             (np.clip(rhs0/A_safe, 0, 1),                  np.zeros_like(A)),
             (np.clip((rhs0 - Cm)/A_safe, 0, 1),           np.ones_like(A))]

    pointA = A < eps0                                   # segment 1 is a point
    pointE = (~pointA) & (E < eps0)                     # segment 2 is a point
    par    = (~pointA) & (~pointE) & (detM < epsM)      # parallel
    gen    = ~(pointA | pointE | par)
    inrng  = gen & (s0 >= 0) & (s0 <= 1) & (s1 >= 0) & (s1 <= 1)

    L1 = np.zeros_like(A)
    L2 = np.zeros_like(A)
    L2 = np.where(pointA, np.where(E > eps0, rhs1/E_safe, 0.0), L2)
    L1 = np.where(pointE, rhs0/A_safe, L1)
    L1 = np.where(inrng, s0, L1)
    L2 = np.where(inrng, s1, L2)
    dist2 = d2(np.clip(L1, 0, 1), np.clip(L2, 0, 1))

    # parallel: the minimum is attained at a pair of endpoints
    corners = np.minimum.reduce([d2(np.zeros_like(A), np.zeros_like(A)),
                                 d2(np.zeros_like(A), np.ones_like(A)),
                                 d2(np.ones_like(A),  np.zeros_like(A)),
                                 d2(np.ones_like(A),  np.ones_like(A))])
    dist2 = np.where(par, corners, dist2)

    # general case with the solution outside the unit square: best of the four trials
    trials = np.minimum.reduce([d2(a, b) for a, b in t])
    dist2 = np.where(gen & ~inrng, trials, dist2)

    if connected is not None:
        dist2 = np.where(connected, 0.0, dist2)
    if degenerate is not None:
        dist2 = np.where(degenerate, -1.0, dist2)
    return dist2

def voigt_vector_to_tensor(voigt_vector):
    return np.array([[voigt_vector[0], voigt_vector[5], voigt_vector[4]],
                     [voigt_vector[5], voigt_vector[1], voigt_vector[3]],
                     [voigt_vector[4], voigt_vector[3], voigt_vector[2]]])

def pkforcevec(sigext, segs_data):
    # return Peach-Koehler force vector for each segment
    # half of it should be assigned to each node
    Nseg = segs_data["nodeids"].shape[0]
    burg_vecs = segs_data["burgers"]
    R1 = segs_data["R1"]
    R2 = segs_data["R2"]
    fpk = np.zeros((Nseg, 3))
    for i in range(Nseg):
        sigb = matvec3(sigext, burg_vecs[i]) if BITREPRO_MATH else sigext @ burg_vecs[i]
        dR = R2[i,:] - R1[i,:]
        fpk[i,:] = cross3(sigb, dR) if BITREPRO_MATH else np.cross(sigb, dR)
    return fpk

def selfforcevec_LineTension(MU, NU, Ec, segs_data, eps_L=1e-6):
    # To do: vectorize the calculations
    Nseg = segs_data["nodeids"].shape[0]
    burg_vecs = segs_data["burgers"]
    R1 = segs_data["R1"]
    R2 = segs_data["R2"]
    fs0 = np.zeros((Nseg, 3))
    fs1 = np.zeros((Nseg, 3))
    omninv = 1.0/(1.0-NU)
    for i in range(Nseg):
        dR = R2[i,:] - R1[i,:]
        L = norm3(dR) if BITREPRO_MATH else np.linalg.norm(dR)
        if L < eps_L:
            continue
        t = dR / L
        bs = dot3(burg_vecs[i,:], t) if BITREPRO_MATH else np.dot(burg_vecs[i,:], t)
        bs2 = bs*bs
        bev = burg_vecs[i] - bs*t
        be2 = dot3(bev, bev) if BITREPRO_MATH else np.sum(bev*bev)
        if BITREPRO_MATH:
            # ExaDiS' core_force (force_types/force_core.h), operation for
            # operation: divide by (1-NU) at each use rather than multiplying
            # by a precomputed reciprocal, and assemble as ft*be + fL*t. The
            # same expression either way, rounded differently, and matching it
            # is what takes the two codes from about 16 ulp apart to the last
            # bit; see tests/unit_tests/test1_node_force.
            fL = -Ec*(bs2 + be2/(1.0-NU))
            ft = Ec*2*bs*NU/(1.0-NU)
            fs1[i,:] = ft*bev + fL*t
        else:
            Score = 2.0*NU*omninv*Ec*bs
            LTcore = (bs2+be2*omninv)*Ec
            fs1[i,:] = Score*bev - LTcore*t
    fs0 = -fs1
    return fs0, fs1

def has_length(p1, p2):
    """has_length: whether segments (p1, p2) are long enough to interact

    ExaDiS' segseg_force returns zero force for a pair with either squared
    length below 1e-20 (core/exadis/src/force_types/force_iso.h), which is what
    keeps the zero-length segment between two coincident trial-split nodes out
    of the pair sum. Works on single points or on arrays of them.
    """
    d = np.asarray(p2) - np.asarray(p1)
    return np.sum(d*d, axis=-1) >= 1.e-20


class CalForce(CalForce_Base):
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}, Ec: float=None, cutoff: float=None,
                 force_mode: str='Elasticity_SBA',
                 self_force_mode: str='SelfForce',
                 use_cell_list: bool=True, batch_size: int=100000,
                 force_kernel: str='batch',
                 torch_device=None, torch_dtype=None) -> None:
        self.mu = state.get("mu", 1.0)
        self.nu = state.get("nu", 0.3)
        self.a =  state.get("a", 0.01)
        # Default core energy from the core radius a, which is ParaDiS rc, and not
        # from the segment-pair cutoff set below. Same expression as ParaDiS
        # Initialize.c:625 and exadis force_core.h CoreDefault. It goes negative
        # for a < 0.1, which is a property of the expression rather than an error.
        self.Ec = self.mu/4.0/np.pi*np.log(self.a/0.1) if Ec is None else Ec
        self.force_mode = force_mode

        # How the elasticity modes get the i == j term and the core-energy term.
        #
        # 'SelfForce'   ParaDiS SelfForceIsotropic (NodeForce.c:2628), one
        #               closed form covering both. This is also what ExaDiS
        #               computes, as self_force plus core_force. The default.
        # 'SegSegSelf'  the earlier pydis route: the segment-segment kernel
        #               against the segment itself for the elastic part, with
        #               the line-tension core term added on top. Kept as an
        #               option so results produced with it stay reproducible.
        #               It matches neither ParaDiS nor ExaDiS.
        if self_force_mode not in ('SelfForce', 'SegSegSelf'):
            raise ValueError("CalForce: unknown self_force_mode '%s', expected "
                             "'SelfForce' or 'SegSegSelf'" % self_force_mode)
        self.self_force_mode = self_force_mode

        # Two-stage segment-pair cutoff; see the module
        # docstring for the convention and why the mid-point stage is reproduced rather than
        # corrected. Segment self forces (i == j) are never cut off: ParaDiS applies
        # SelfForce (external/paradis/src/NodeForce.c:2688) to every segment
        # unconditionally.
        self.cutoff = cutoff
        self.cutoff2 = None if cutoff is None or cutoff < 0.0 else cutoff*cutoff
        self.maxseg = state.get("maxseg", None)
        if self.cutoff2 is not None and self.maxseg is None:
            raise KeyError("CalForce: state must define 'maxseg' when a cutoff is used, "
                           "since the mid-point stage of the cutoff is applied at "
                           "cutoff + maxseg (see the module docstring)")

        # Cell list used to avoid testing every segment pair. Only a candidate filter:
        # within_cutoff remains the authoritative test, evaluated on true positions under
        # the minimum image convention, so the bin-derived-distance defect that affects
        # the neighbor-list defect noted in the module docstring is not reproduced here.
        # Built lazily in _nbrlist_for(),
        # because the cell can change between calls.
        self.use_cell_list = use_cell_list
        self._nbrlist = None
        self._nbrlist_key = None

        # Pair forces are evaluated in batches rather than all at once. The pair count grows
        # as Nseg^2 before the cell list prunes it, and the vectorized kernel materialises
        # several hundred temporaries of the batch shape, so peak memory is a large multiple
        # of the batch itself. Tune per machine: a GPU node wants a larger value than a laptop.
        self.batch_size = batch_size
        # 'batch'  compiled ParaDiS kernel via SegSegForceList: one ctypes call for the
        #          whole batch, looping inside C. Same numbers as 'scalar', without paying
        #          ~8 us of marshalling per pair. The default.
        # 'scalar' the same kernel called once per pair; kept as the reference.
        # 'vec'    numpy, batched (python_segseg_force_vec).
        # 'torch'  the same kernel on a GPU (torch_segseg_force_vec). Needs torch;
        #          steer the device with torch_device or PYDIS_TORCH_DEVICE.
        #
        # 'vec' and 'torch' share a lineage and give the same answers as each other. Both
        # differ from 'batch'/'scalar' by more than rounding once a pair comes within
        # 1-c^2 < 1e-6 of parallel, up to 2e-4 relative at 1e-12; that is a property of the
        # numpy kernel rather than of torch, and it is why 'batch' remains the default.
        if force_kernel not in ('batch', 'scalar', 'vec', 'torch'):
            raise ValueError("CalForce: force_kernel must be one of 'batch', "
                             "'scalar', 'vec', 'torch'; got %r" % (force_kernel,))
        self.force_kernel = force_kernel
        self.torch_device = torch_device
        self.torch_dtype = torch_dtype

        self.NodeForce_Functions = {
            'LineTension': self.NodeForce_LineTension,
            'Elasticity_SBA': self.NodeForce_Elasticity_SBA,
            'Elasticity_SBN1_SBA': self.NodeForce_Elasticity_SBN1_SBA }
        self.OneNodeForce_Functions = {
            'LineTension': self.OneNodeForce_LineTension,
            'Elasticity_SBA': self.OneNodeForce_Elasticity_SBA,
            'Elasticity_SBN1_SBA': self.OneNodeForce_Elasticity_SBN1_SBA }

    def _nbrlist_for(self, cell, Nseg):
        """_nbrlist_for: cell list sized for this cell and cutoff, or None for all pairs

        The division count is n_div = floor(width / (cutoff + maxseg)), so that
        anything within the search radius lies in an adjacent cell. Returns None when a cell
        list cannot be used or would not help, in which case the caller falls back to
        iterating every pair:
          - no cutoff, so every pair is wanted anyway
          - the cell is not periodic in all three directions (sort_points_to_list only
            builds the list under full PBC)
          - fewer than 4 divisions. Below 3 the +-1 neighbour scan would reach the same
            cell by more than one offset and double count; at exactly 3 it reaches every
            cell, so nothing is filtered and the bookkeeping is pure overhead.
        """
        if self.cutoff2 is None or not self.use_cell_list:
            return None
        if not all(cell.is_periodic):
            return None

        h = np.asarray(cell.h, dtype=float)
        eff = self.cutoff + self.maxseg
        n_div = []
        for i in range(3):
            e1, e2 = h[:, (i+1) % 3], h[:, (i+2) % 3]
            perp = np.cross(e1, e2)
            width = abs(np.dot(h[:, i], perp/np.linalg.norm(perp)))
            n_div.append(int(np.floor(width/eff)))
        if min(n_div) < 4:
            return None

        key = (tuple(n_div), h.tobytes())
        if self._nbrlist_key != key:
            self._nbrlist = CellList(cell=cell, n_div=n_div)
            self._nbrlist_key = key
        return self._nbrlist

    def candidate_pairs(self, G, R1, R2, Nseg):
        """candidate_pairs: iterate (i, j), i < j, of segment pairs worth testing

        Uses the cell list when it applies, otherwise every pair. Segments are binned by
        mid-point, as the cutoff convention and the collision handler both require.
        """
        nbrlist = self._nbrlist_for(G.cell, Nseg)
        if nbrlist is None:
            for i in range(Nseg):
                for j in range(i+1, Nseg):
                    yield i, j
        else:
            nbrlist.sort_points_to_list(0.5*(R1 + R2))
            for i, j in nbrlist.iterate_nbr_pairs(use_cell_list=True):
                yield i, j

    def _pair_forces(self, P1, P2, P3, P4, B12, B34):
        """_pair_forces: forces on the four endpoints of each pair, in bounded batches

        Yields (start, stop, f1, f2, f3, f4) so the caller can scatter each batch without
        ever holding the whole pair set of results at once.
        """
        n = P1.shape[0]
        step = max(1, int(self.batch_size))
        for start in range(0, n, step):
            stop = min(start + step, n)
            sl = slice(start, stop)
            if self.force_kernel == 'batch':
                f1, f2, f3, f4 = compute_segseg_force_batch(
                    P1[sl], P2[sl], P3[sl], P4[sl], B12[sl], B34[sl],
                    self.mu, self.nu, self.a)
            elif self.force_kernel == 'torch':
                from .compute_stress_force_analytic_torch import (
                    torch_segseg_force_vec)
                f1, f2, f3, f4 = torch_segseg_force_vec(
                    P1[sl], P2[sl], P3[sl], P4[sl], B12[sl], B34[sl],
                    self.mu, self.nu, self.a,
                    device=self.torch_device, dtype=self.torch_dtype)
            elif self.force_kernel == 'scalar':
                f1 = np.empty((stop-start, 3)); f2 = np.empty_like(f1)
                f3 = np.empty_like(f1); f4 = np.empty_like(f1)
                for k in range(stop-start):
                    f1[k], f2[k], f3[k], f4[k] = compute_segseg_force(
                        P1[sl][k], P2[sl][k], P3[sl][k], P4[sl][k],
                        B12[sl][k], B34[sl][k], self.mu, self.nu, self.a)
            else:
                f1, f2, f3, f4 = python_segseg_force_vec(
                    P1[sl], P2[sl], P3[sl], P4[sl], B12[sl], B34[sl],
                    self.mu, self.nu, self.a)
            yield start, stop, f1, f2, f3, f4

    def select_pairs(self, G, R1, R2, source_tags, target_tags, Nseg):
        """select_pairs: (i, j) index arrays of the segment pairs that interact

        Same two-stage criterion as within_cutoff, applied to the whole candidate set at
        once instead of one pair at a time. Candidates come from the cell list when it
        applies, so this is O(N) rather than O(N^2) for large networks. Returns (idx_i,
        idx_j) with idx_i < idx_j, plus the PBC-mapped endpoints of each surviving pair so
        the caller does not have to image them again.
        """
        pairs = np.fromiter((k for ij in self.candidate_pairs(G, R1, R2, Nseg) for k in ij),
                            dtype=np.int64)
        if pairs.size == 0:
            e = np.zeros((0, 3))
            return np.zeros(0, np.int64), np.zeros(0, np.int64), e, e, e, e
        i, j = pairs[0::2], pairs[1::2]

        # PBC: image segment j next to segment i, matching the scalar path exactly
        p1 = R1[i]
        p2 = G.cell.closest_image(Rref=p1, R=R2[i])
        p3 = G.cell.closest_image(Rref=p1, R=R1[j])
        p4 = G.cell.closest_image(Rref=p3, R=R2[j])

        if self.cutoff2 is None:
            return i, j, p1, p2, p3, p4

        # stage 1: mid-point separation <= cutoff + maxseg
        c12 = 0.5*(p1 + p2)
        c34 = G.cell.closest_image(Rref=c12, R=0.5*(p3 + p4))
        keep = np.einsum('ij,ij->i', c34-c12, c34-c12) <= (self.cutoff + self.maxseg)**2
        i, j, p1, p2, p3, p4 = (x[keep] for x in (i, j, p1, p2, p3, p4))
        if i.size == 0:
            return i, j, p1, p2, p3, p4

        # stage 2: minimum segment-segment distance < cutoff
        eps0 = 1.0e-12
        d12 = np.einsum('ij,ij->i', p2-p1, p2-p1)
        d34 = np.einsum('ij,ij->i', p4-p3, p4-p3)
        degenerate = (d12 < eps0) | (d34 < eps0)
        connected = ((source_tags[i] == source_tags[j]).all(axis=1) |
                     (source_tags[i] == target_tags[j]).all(axis=1) |
                     (target_tags[i] == source_tags[j]).all(axis=1) |
                     (target_tags[i] == target_tags[j]).all(axis=1))
        dist2 = min_dist2_segseg_vec(p1, p2, p3, p4,
                                     connected=connected, degenerate=degenerate)
        keep = (dist2 >= 0.0) & (dist2 < self.cutoff2)
        return (i[keep], j[keep], p1[keep], p2[keep], p3[keep], p4[keep])

    def within_cutoff(self, cell, p1, p2, p3, p4, tags12, tags34) -> bool:
        """within_cutoff: whether segments (p1,p2) and (p3,p4) interact under self.cutoff

        Applies the two-stage criterion documented at the top of this module: mid-point separation <= cutoff + maxseg, then minimum segment-segment
        distance < cutoff. Endpoints must already be mapped to a common periodic image by
        the caller.
        """
        # stage 1: mid-point test, reproducing the segment binning described in the module
        # docstring. Unlike a bound based on the true segment lengths this can discard pairs
        # when a segment is longer than maxseg.
        c12, c34 = 0.5*(p1+p2), 0.5*(p3+p4)
        c34 = cell.closest_image(Rref=c12, R=c34)
        dc = np.linalg.norm(c34 - c12)
        if dc > self.cutoff + self.maxseg:
            return False

        # stage 2: exact minimum distance between the two segments
        connected = bool(set(tags12) & set(tags34))
        dist2 = min_dist2_segseg(p1, p2, p3, p4, connected)
        return dist2 >= 0.0 and dist2 < self.cutoff2

    def NodeForce(self, DM: DisNetManager, state: dict, pre_compute: bool=True) -> dict:
        """NodeForce: return nodal forces in a dictionary

        Using different force calculation functions depending on force_mode
        """
        applied_stress = state["applied_stress"]
        G = DM.get_disnet(DisNet)
        nodeforce_dict, segforce_dict = self.NodeForce_Functions[self.force_mode](G, applied_stress)
        state["nodeforce_dict"] = nodeforce_dict
        state["segforce_dict"] = segforce_dict

        # prepare nodeforces and nodeforce_tags arrays for compatibility with exadis
        state = DisNet.convert_nodeforce_dict_to_array(state)

        return state

    def PreCompute(self, DM: DisNetManager, state: dict) -> dict:
        """PreCompute: pre-compute some data for force calculation
        """
        #G = DM.get_disnet(DisNet)
        #segs_data_with_positions = G.get_segs_data_with_positions()
        #state["segs_data_with_positions"] = segs_data_with_positions
        return state

    def self_and_core_force(self, G, R1, R2, burg_vecs, segs_data, rows):
        """self_and_core_force: the i == j term and the core-energy term

        Returns (f0, f1), the contribution to each segment's two nodes, zero on
        the segments not named in rows. Which route is taken is set by
        self_force_mode; see __init__.
        """
        nseg = burg_vecs.shape[0]
        f0 = np.zeros((nseg, 3))
        f1 = np.zeros((nseg, 3))
        idx = np.asarray(list(rows), dtype=int)
        if idx.size == 0:
            return f0, f1

        if self.self_force_mode == 'SelfForce':
            p1 = np.array([R1[i, :] for i in idx])
            p2 = np.array([G.cell.closest_image(Rref=R1[i, :], R=R2[i, :])
                           for i in idx])
            a0, a1 = compute_selfforce_list(burg_vecs[idx, :], p1, p2,
                                            self.mu, self.nu, self.a, self.Ec)
            f0[idx] = a0
            f1[idx] = a1
            return f0, f1

        # 'SegSegSelf': the segment against itself, then the core term on top
        for i in idx:
            p1 = R1[i, :].copy()
            p2 = G.cell.closest_image(Rref=p1, R=R2[i, :].copy())
            b12 = burg_vecs[i, :].copy()
            g1, g2, _, _ = compute_segseg_force(p1, p2, p1, p2, b12, b12,
                                                self.mu, self.nu, self.a)
            f0[i] += g1
            f1[i] += g2
        fs0, fs1 = selfforcevec_LineTension(self.mu, self.nu, self.Ec, segs_data)
        f0[idx] += fs0[idx]
        f1[idx] += fs1[idx]
        return f0, f1

    def OneNodeForce(self, DM: DisNetManager, state: dict, tag: Tag,
                     update_state: bool=True, match_global: bool=False) -> dict:
        """OneNodeForce: compute force calculation on one node

        match_global=True selects the segment pairs exactly as NodeForce does,
        so the result equals that node's NodeForce entry to rounding. The
        default reproduces ExaDiS' default node_force instead, which Topology
        is meant to use: one bin query at the node position, kept for every
        arm, with no segment-to-segment distance test (see neighbor_bin). The
        two coincide when there is no cutoff.
        """
        applied_stress = state["applied_stress"]
        G = DM.get_disnet(DisNet)
        f = self.OneNodeForce_Functions[self.force_mode](G, applied_stress, tag,
                                                         match_global=match_global)
        # update force dictionary if needed
        if update_state:
            if "nodeforces" in state and "nodeforcetags" in state:
                nodeforcetags = state["nodeforcetags"]
                ind = np.where((nodeforcetags[:,0]==tag[0])&(nodeforcetags[:,1]==tag[1]))[0]
                if ind.size == 1:
                    state["nodeforces"][ind[0]] = f
                else:
                    state["nodeforces"] = np.vstack((state["nodeforces"], f))
                    state["nodeforcetags"] = np.vstack((state["nodeforcetags"], tag))
            else:
                state["nodeforces"] = np.array([f])
                state["nodeforcetags"] = np.array([tag])

        return f

    def OneNodeForce_LineTension(self, G: DisNet, applied_stress: np.ndarray, tag,
                                 match_global: bool=False) -> float:
        """OneNodeForce_LineTension: return force on one node from line tension

        match_global is accepted for the common signature; with no pair sum
        there is nothing for it to select.
        """
        # To do: refactor this into a function
        Nseg = G.out_degree(tag) # This line is different
        nodeids = np.zeros((Nseg, 2), dtype=int)
        tag1 = np.zeros((Nseg, 2), dtype=int)
        tag2 = np.zeros((Nseg, 2), dtype=int)
        burgers = np.zeros((Nseg, 3))
        planes = np.zeros((Nseg, 3))
        R1 = np.zeros((Nseg, 3))
        R2 = np.zeros((Nseg, 3))
        i = 0
        for nbr_tag, edge_attr in G.neighbor_segments_dict(tag).items(): # This line is different
            nodeids[i,:] = -1, -1 # This line is different
            source = tag          # This line is different
            target = nbr_tag      # This line is different
            tag1[i,:] = source
            tag2[i,:] = target
            burgers[i,:] = edge_attr.burg_vec_from(source).copy()
            planes[i,:] = getattr(edge_attr, "plane_normal", np.zeros(3)).copy()
            r1_local = G.nodes(source).R  # This line is different
            r2_local = G.nodes(target).R  # This line is different
            # apply PBC
            r2_local = G.cell.closest_image(Rref=r1_local, R=r2_local) # This line is different
            R1[i,:] = r1_local
            R2[i,:] = r2_local
            i += 1
        segs_data_with_positions = {
            "nodeids": nodeids,
            "tag1": tag1,
            "tag2": tag2,
            "burgers": burgers,
            "planes":  planes,
            "R1": R1,
            "R2": R2
        }
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]

        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fs0, fs1 = selfforcevec_LineTension(self.mu, self.nu, self.Ec, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5 + fs0, fpk*0.5 + fs1))

        f = np.zeros(3)

        for i in range(Nseg):
            tag1 = tuple(source_tags[i])
            tag2 = tuple(target_tags[i])
            if tag == tag1:
                f += fseg[i, 0:3]
            elif tag == tag2:
                f += fseg[i, 3:6]

        return f

    def OneNodeForce_Elasticity_SBA(self, G: DisNet, applied_stress: np.ndarray, tag,
                                    match_global: bool=False) -> float:
        """OneNodeForce_Elasticity_SBA: force on one node, elastic interactions

        Translated from SetOneNodeForce in ParaDiS
        (external/paradis/src/NodeForce.c, the FULL_N2_FORCES branch, which is
        the one that matches this force mode: a full pair sum with no remote or
        FMM contribution). For each arm of the node ParaDiS computes the
        segment self force, the Peach-Koehler force from the external stress,
        and then the interaction of that arm with every other segment, adding
        only to that arm; the nodal force is the sum over the node's arms.
        The osmotic and FEM terms in the C are not reproduced, having no
        counterpart here.

        The result must equal the NodeForce_Elasticity_SBA entry for this node
        to rounding, since both are the same sum of the same terms. That is
        what makes it usable by Topology, and
        tests/unit_tests/test1_node_force/test_onenodeforce_pydis.py pins it.
        Getting there requires selecting the same pairs and imaging them the
        same way, so the pair set comes from self.select_pairs rather than from
        a second implementation of the cutoff test: two ways of choosing pairs
        would be two things to keep in agreement.

        Cost. ParaDiS loops over every segment for each arm, and so does this:
        the saving over NodeForce is that forces are accumulated for the node's
        arms only, not for all segments, not that fewer pairs are examined.
        With the cell list the candidate scan is O(Nseg) rather than O(Nseg^2),
        which is what makes it usable inside the trial splits Topology performs
        per multi-arm node per step.
        """
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        burg_vecs = segs_data_with_positions["burgers"]

        # the arms of this node, i.e. the segments ParaDiS loops over as
        # armIndex. A segment can appear with the node at either end, and a
        # self-loop appears at both.
        tag_arr = np.array(tag)
        at_source = (source_tags == tag_arr).all(axis=1)
        at_target = (target_tags == tag_arr).all(axis=1)
        arms = np.where(at_source | at_target)[0]
        if arms.size == 0:
            return np.zeros(3)

        # ExtPKForce: the applied stress contributes to each arm, split evenly
        # between its two endpoints. Computed for every segment because
        # pkforcevec is vectorized over the whole list; only the arm rows are
        # ever read.
        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        # Self force and core force, on the arms only. Never cut off, matching
        # NodeForce_Elasticity_SBA.
        sf0, sf1 = self.self_and_core_force(G, R1, R2, burg_vecs,
                                            segs_data_with_positions, arms)
        fseg[:, 0:3] += sf0
        fseg[:, 3:6] += sf1

        # ComputeForces against the other segments
        if match_global or self.cutoff2 is None:
            self._add_global_pairs(G, R1, R2, burg_vecs, source_tags,
                                   target_tags, Nseg, arms, fseg)
        else:
            self._add_node_query_pairs(G, R1, R2, burg_vecs, tag, arms, fseg)

        # the nodal force is the sum of this node's arm forces. Both ends are
        # added when a segment starts and finishes at this node, matching the
        # unconditional accumulation in NodeForce_Elasticity_SBA.
        f = np.zeros(3)
        for i in arms:
            if at_source[i]:
                f += fseg[i, 0:3]
            if at_target[i]:
                f += fseg[i, 3:6]
        return f

    def _add_global_pairs(self, G, R1, R2, burg_vecs, source_tags, target_tags,
                          Nseg, arms, fseg):
        """_add_global_pairs: add the arms' share of the pairs NodeForce computes

        select_pairs returns i < j, so an arm can appear on either side of a
        pair and takes (f1, f2) or (f3, f4) accordingly. Pairs touching no arm
        are dropped before the kernel runs, which is the whole saving over
        NodeForce.
        """
        idx_i, idx_j, P1, P2, P3, P4 = self.select_pairs(G, R1, R2, source_tags,
                                                         target_tags, Nseg)
        if idx_i.size:
            is_arm = np.zeros(Nseg, dtype=bool)
            is_arm[arms] = True
            touches = is_arm[idx_i] | is_arm[idx_j]
            idx_i, idx_j = idx_i[touches], idx_j[touches]
            P1, P2, P3, P4 = P1[touches], P2[touches], P3[touches], P4[touches]
        if not idx_i.size:
            return
        f_a, f_b = fseg[:, 0:3], fseg[:, 3:6]
        for start, stop, f1, f2, f3, f4 in self._pair_forces(
                P1, P2, P3, P4, burg_vecs[idx_i], burg_vecs[idx_j]):
            ii, jj = idx_i[start:stop], idx_j[start:stop]
            np.add.at(f_a, ii, f1)
            np.add.at(f_b, ii, f2)
            np.add.at(f_a, jj, f3)
            np.add.at(f_b, jj, f4)

    def node_query_partners(self, G, R1, R2, tag) -> np.ndarray:
        """node_query_partners: segments ExaDiS' default node_force pairs with tag

        One query at the node's own position, unfolded as ExaDiS leaves it,
        over the segments' folded mid-points at cutoff + maxseg; see
        neighbor_bin. R2 is already imaged next to R1.
        """
        mid = G.cell.fold(0.5*(R1 + R2))
        return np.array(neighbor_bin.query(G.cell, mid, G.nodes(tag).R,
                                           self.cutoff + self.maxseg),
                        dtype=np.int64)

    def _add_node_query_pairs(self, G, R1, R2, burg_vecs, tag, arms, fseg):
        """_add_node_query_pairs: add each arm's pairs with the queried segments

        Every arm k against every queried segment but itself, k taken first
        and the partner imaged next to it, as ExaDiS' segseg_force(k, n) does.
        Only k's end forces are kept.
        """
        near = self.node_query_partners(G, R1, R2, tag)
        K, N = np.repeat(arms, near.size), np.tile(near, arms.size)
        K, N = K[K != N], N[K != N]
        if not K.size:
            return
        P1 = R1[K]
        P2 = G.cell.closest_image(Rref=P1, R=R2[K])
        P3 = G.cell.closest_image(Rref=P1, R=R1[N])
        P4 = G.cell.closest_image(Rref=P3, R=R2[N])
        live = has_length(P1, P2) & has_length(P3, P4)
        if not live.any():
            return
        K, N, P1, P2, P3, P4 = (x[live] for x in (K, N, P1, P2, P3, P4))
        f_a, f_b = fseg[:, 0:3], fseg[:, 3:6]
        for start, stop, f1, f2, _, _ in self._pair_forces(
                P1, P2, P3, P4, burg_vecs[K], burg_vecs[N]):
            np.add.at(f_a, K[start:stop], f1)
            np.add.at(f_b, K[start:stop], f2)

    @staticmethod
    def _sbn1_partners(i, Nseg, near):
        """_sbn1_partners: (a, b) segment pairs for arm i, as the kernel takes them

        near None: every segment, in the order NodeForce evaluated each pair, so
        the imaging matches. Otherwise arm i's self term, then i first against
        each queried segment, as ExaDiS orders a pair in node_force.
        """
        if near is None:
            return [(i, j) if j >= i else (j, i) for j in range(Nseg)]
        return [(i, i)] + [(i, int(n)) for n in near if n != i]

    def OneNodeForce_Elasticity_SBN1_SBA(self, G: DisNet, applied_stress: np.ndarray, tag,
                                         match_global: bool=False) -> float:
        """OneNodeForce_Elasticity_SBN1_SBA: force on one node, SBN1 quadrature

        The SBN1 counterpart of OneNodeForce_Elasticity_SBA, and translated
        from the same ParaDiS function, SetOneNodeForce in
        external/paradis/src/NodeForce.c: per arm, the segment self force, the
        Peach-Koehler force from the external stress, and the interaction of
        that arm with every other segment, accumulated onto that arm alone.

        It must reproduce the NodeForce_Elasticity_SBN1_SBA entry for this node
        to rounding, so it follows that function's scalar path rather than the
        batched one used by the SBA mode: same imaging, same within_cutoff
        test, same kernel. Two details are load-bearing for that:

          - the self term is the j == i case, and is never cut off
          - a pair is evaluated in the order NodeForce evaluated it. Its loop
            runs j from i upward, so a partner with j < i was computed as the
            pair (j, i) with the imaging anchored on segment j, and this arm
            took (f3, f4) from it. Recomputing it here as (i, j) would image it
            differently and give a different answer in the last digits.
        """
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        burg_vecs = segs_data_with_positions["burgers"]

        tag_arr = np.array(tag)
        at_source = (source_tags == tag_arr).all(axis=1)
        at_target = (target_tags == tag_arr).all(axis=1)
        arms = np.where(at_source | at_target)[0]
        if arms.size == 0:
            return np.zeros(3)

        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        # must match NodeForce_Elasticity_SBN1_SBA, which hardcodes
        # force_nint = 3
        quad_points = np.array([-0.774596669241483, 0.0, 0.774596669241483])
        weights = np.array([0.555555555555556, 0.888888888888889,
                            0.555555555555556])

        near = (None if match_global or self.cutoff2 is None
                else self.node_query_partners(G, R1, R2, tag))
        for i in arms:
            for a_idx, b_idx in self._sbn1_partners(i, Nseg, near):
                p1 = R1[a_idx,:].copy()
                p2 = R2[a_idx,:].copy()
                p3 = R1[b_idx,:].copy()
                p4 = R2[b_idx,:].copy()
                b12 = burg_vecs[a_idx,:].copy()
                b34 = burg_vecs[b_idx,:].copy()

                # apply PBC
                p2 = G.cell.closest_image(Rref=p1, R=p2)
                p3 = G.cell.closest_image(Rref=p1, R=p3)
                p4 = G.cell.closest_image(Rref=p3, R=p4)

                tag1, tag2 = tuple(source_tags[a_idx]), tuple(target_tags[a_idx])
                tag3, tag4 = tuple(source_tags[b_idx]), tuple(target_tags[b_idx])
                if (near is None and self.cutoff2 is not None
                        and a_idx != b_idx):
                    if not self.within_cutoff(G.cell, p1, p2, p3, p4,
                                              (tag1, tag2), (tag3, tag4)):
                        continue
                if (near is not None and a_idx != b_idx
                        and not (has_length(p1, p2) and has_length(p3, p4))):
                    continue

                f1, f2, f3, f4 = compute_segseg_force_SBN1_SBA(
                    p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a,
                    quad_points, weights)

                if a_idx == b_idx:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                elif a_idx == i:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                else:
                    fseg[i, 0:3] += f3
                    fseg[i, 3:6] += f4

        f = np.zeros(3)
        for i in arms:
            if at_source[i]:
                f += fseg[i, 0:3]
            if at_target[i]:
                f += fseg[i, 3:6]
        return f

    def NodeForce_LineTension(self, G: DisNet, applied_stress: np.ndarray) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from line tension in a dictionary

        Only Peach-Koehler force from external stress and line tension forces
        (from DDLab/src/segforcevec.m)
        Note: assuming G.seg_list already accounts for PBC
        """
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]

        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fs0, fs1 = selfforcevec_LineTension(self.mu, self.nu, self.Ec, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5 + fs0, fpk*0.5 + fs1))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.all_nodes_tags():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})

        for i in range(Nseg):
            tag1 = tuple(source_tags[i])
            tag2 = tuple(target_tags[i])
            nodeforce_dict[tag1] += fseg[i, 0:3]
            nodeforce_dict[tag2] += fseg[i, 3:6]
            segforce_dict[(tag1, tag2)] = fseg[i, :]

        return nodeforce_dict, segforce_dict

    def NodeForce_Elasticity_SBA(self, G: DisNet, applied_stress: np.ndarray) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        Note: assuming G.get_segs_data_with_positions already accounts for PBC
        """
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        burg_vecs = segs_data_with_positions["burgers"]

        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.all_nodes_tags():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})

        # Self force and core force on every segment. Never cut off: ParaDiS
        # applies SelfForce unconditionally (NodeForce.c:2688), and the cell
        # list never yields i == j anyway, since it yields only i < j.
        sf0, sf1 = self.self_and_core_force(G, R1, R2, burg_vecs,
                                           segs_data_with_positions,
                                           range(Nseg))
        fseg[:, 0:3] += sf0
        fseg[:, 3:6] += sf1

        # pair forces (i < j). Selection runs on the whole candidate set at once, then the
        # kernel is evaluated in bounded batches and scattered back onto the segments.
        idx_i, idx_j, P1, P2, P3, P4 = self.select_pairs(G, R1, R2,
                                                         source_tags, target_tags, Nseg)
        if idx_i.size:
            B12, B34 = burg_vecs[idx_i], burg_vecs[idx_j]
            f_a, f_b = fseg[:, 0:3], fseg[:, 3:6]
            for start, stop, f1, f2, f3, f4 in self._pair_forces(P1, P2, P3, P4, B12, B34):
                ii, jj = idx_i[start:stop], idx_j[start:stop]
                np.add.at(f_a, ii, f1)
                np.add.at(f_b, ii, f2)
                np.add.at(f_a, jj, f3)
                np.add.at(f_b, jj, f4)

        # nodal forces are the accumulated segment forces; equivalent to updating the dict
        # inside the loops above, but done once now that the loops are batched
        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            nodeforce_dict[tag1] += fseg[i, 0:3]
            nodeforce_dict[tag2] += fseg[i, 3:6]
            segforce_dict[(tag1, tag2)] = fseg[i, :]

        return nodeforce_dict, segforce_dict

    def NodeForce_Elasticity_SBN1_SBA(self, G: DisNet, applied_stress: np.ndarray) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        Note: assuming G.seg_list already accounts for PBC
        """
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        burg_vecs = segs_data_with_positions["burgers"]

        # To do: need to run test case for this function
        sigext = voigt_vector_to_tensor(applied_stress)
        fpk = pkforcevec(sigext, segs_data_with_positions)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.all_nodes_tags():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            nodeforce_dict[tag1] += fseg[i, 0:3]
            nodeforce_dict[tag2] += fseg[i, 3:6]

        """ hardcode force_nint = 3
        """
        quad_points = np.array([-0.774596669241483, 0.0, 0.774596669241483])
        weights = np.array([0.555555555555556, 0.888888888888889, 0.555555555555556])

        for i in range(Nseg):
            for j in range(i, Nseg):
                p1 = R1[i,:].copy()
                p2 = R2[i,:].copy()
                p3 = R1[j,:].copy()
                p4 = R2[j,:].copy()
                b12 = burg_vecs[i,:].copy()
                b34 = burg_vecs[j,:].copy()

                # apply PBC
                p2 = G.cell.closest_image(Rref=p1, R=p2)
                p3 = G.cell.closest_image(Rref=p1, R=p3)
                p4 = G.cell.closest_image(Rref=p3, R=p4)
                tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
                tag3, tag4 = tuple(source_tags[j]), tuple(target_tags[j])
                # apply the spherical cutoff (self forces, i == j, are never cut off)
                if self.cutoff2 is not None and i != j:
                    if not self.within_cutoff(G.cell, p1, p2, p3, p4, (tag1, tag2), (tag3, tag4)):
                        continue
                f1, f2, f3, f4 = compute_segseg_force_SBN1_SBA(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a, quad_points, weights)
                if i == j:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    nodeforce_dict[tag1] += f1
                    nodeforce_dict[tag2] += f2
                else:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    fseg[j, 0:3] += f3
                    fseg[j, 3:6] += f4
                    nodeforce_dict[tag1] += f1
                    nodeforce_dict[tag2] += f2
                    nodeforce_dict[tag3] += f3
                    nodeforce_dict[tag4] += f4

        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            segforce_dict[(tag1, tag2)] = fseg[i, :]

        return nodeforce_dict, segforce_dict
