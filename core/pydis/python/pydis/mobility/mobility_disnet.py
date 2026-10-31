"""@package docstring
Mobility_DisNet: class for defining Mobility Laws

Provide mobility law functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager
from framework.mobility_base import MobilityLaw_Base

from typing import Tuple


def matrix33_inverse(A: np.ndarray) -> np.ndarray:
    """matrix33_inverse: the inverse of a 3x3 matrix, by its adjugate, as ParaDiS' Matrix33Invert
       (Matrix.c)
    """
    C = np.array([[A[1,1]*A[2,2] - A[1,2]*A[2,1], A[1,2]*A[2,0] - A[1,0]*A[2,2], A[1,0]*A[2,1] - A[2,0]*A[1,1]],
                  [A[2,1]*A[0,2] - A[0,1]*A[2,2], A[2,2]*A[0,0] - A[2,0]*A[0,2], A[2,0]*A[0,1] - A[2,1]*A[0,0]],
                  [A[0,1]*A[1,2] - A[0,2]*A[1,1], A[0,2]*A[1,0] - A[1,2]*A[0,0], A[0,0]*A[1,1] - A[0,1]*A[1,0]]])
    det = A[0,0]*C[0,0] + A[0,1]*C[0,1] + A[0,2]*C[0,2]
    if abs(det) < 1e-20:
        raise ValueError("matrix33_inverse: the matrix cannot be inverted")
    return C.T / det


def orthogonalize(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """orthogonalize: a less its component along b, as ParaDiS' Orthogonalize (Util.c); a if b is zero
    """
    b2 = np.dot(b, b)
    return a if b2 <= 0.0 else a - (np.dot(a, b) / b2) * b


def independent(vectors: list, eps: float=0.05) -> list:
    """independent: each vector orthogonalized against the ones before it, and zero when less than
       eps of its square is left (ParaDiS' FFACTOR_ORTH)
    """
    found = []
    for v in vectors:
        for u in found:
            v = orthogonalize(v, u)
        found.append(v if np.dot(v, v) >= eps else np.zeros(3))
    return found


class MobilityLaw(MobilityLaw_Base):
    """MobilityLaw: class for mobility laws
    """
    def __init__(self, state: dict={}, mobility_law: str='Relax', vmax: float=None,
                 Medge: float=None, Mscrew: float=None, Mclimb: float=None) -> None:
        self.mobility_law = mobility_law
        self.mob = state.get("mob", 1.0)
        # vmax caps the speed: in b/s for SimpleGlide, in m/s, smoothly, for BCC_0B and FCC_0 (as ExaDiS)
        self.vmax = (-1.0 if mobility_law in ('BCC_0B', 'FCC_0') else 1e9) if vmax is None else vmax
        self.burgmag = state.get("burgmag", 1.0)

        self.NodeMobility_Functions = {
            'Relax': self.NodeMobility_Relax,
            'SimpleGlide': self.NodeMobility_SimpleGlide,
            'BCC_0B': self.NodeMobility_BCC_0B,
            'FCC_0': self.NodeMobility_FCC_0
        }
        if mobility_law == 'FCC_0':
            if state.get("crystal", "fcc").lower() != "fcc":
                raise ValueError("MobilityLaw FCC_0: the crystal must be fcc")
            if not all(M is not None and M > 0.0 for M in (Medge, Mscrew)):
                raise ValueError("MobilityLaw FCC_0: Medge and Mscrew must be positive")
            self.Medge, self.Mscrew = Medge, Mscrew
        if mobility_law == 'BCC_0B':
            if state.get("crystal", "bcc").lower() != "bcc":
                raise ValueError("MobilityLaw BCC_0B: the crystal must be bcc")
            if not all(M is not None and M > 0.0 for M in (Medge, Mscrew, Mclimb)):
                raise ValueError("MobilityLaw BCC_0B: Medge, Mscrew and Mclimb must be positive")
            self.Bedge, self.Bscrew, self.Beclimb = 1.0/Medge, 1.0/Mscrew, 1.0/Mclimb
            self.Bscrew2, self.Beclimb2 = self.Bscrew*self.Bscrew, self.Beclimb*self.Beclimb
            self.invBscrew2 = 1.0/(self.Bscrew*self.Bscrew)
            self.invBedge2 = 1.0/(self.Bedge*self.Bedge)
            self.Bline = 1.0e-2 * min(self.Bscrew, self.Bedge)
            self.BlmBsc, self.BlmBecl = self.Bline - self.Bscrew, self.Bline - self.Beclimb
            self.split3node = state.get("split3node", 1)
        
    def Mobility(self, DM: DisNetManager, state: dict) -> dict:
        """Mobility: calculate node velocity according to mobility law function
        """
        G = DM.get_disnet(DisNet)
        if "nodeforces" in state and "nodeforcetags" in state:
            DisNet.convert_nodeforce_array_to_dict(state)

        nodeforce_dict = state["nodeforce_dict"]
        vel_dict = nodeforce_dict.copy()
        for tag in G.all_nodes_tags():
            f = vel_dict[tag].copy()
            vel_dict[tag] = self.NodeMobility_Functions[self.mobility_law](G, tag, f)
        state["vel_dict"] = vel_dict

        # prepare nodeforces and nodeforce_tags arrays for compatibility with exadis
        state = DisNet.convert_nodevel_dict_to_array(state)
        return state
        
    def OneNodeMobility(self, DM: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
        """OneNodeMobility: compute and return the mobility of one node specified by its tag
        """
        G = DM.get_disnet(DisNet)
        v = self.NodeMobility_Functions[self.mobility_law](G, tag, f)
        # update velocity dictionary if needed
        if update_state:
            if "nodevels" in state and "nodeveltags" in state:
                nodeveltags = state["nodeveltags"]
                ind = np.where((nodeveltags[:,0]==tag[0])&(nodeveltags[:,1]==tag[1]))[0]
                if ind.size == 1:
                    state["nodevels"][ind[0]] = v
                else:
                    state["nodevels"] = np.vstack((state["nodevels"], v))
                    state["nodeveltags"] = np.vstack((state["nodeveltags"], tag))
            else:
                state["nodevels"] = np.array([v])
                state["nodeveltags"] = np.array([tag])
        return v
    
    def NodeMobility_Relax(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_Relax: node velocity equal node force
        """
        vel = 1.0*f
        # set velocity of pinned nodes to zero
        for tag in G.all_nodes_tags():
            if G.nodes(tag).constraint == DisNode.Constraints.PINNED_NODE:
                vel = np.zeros(3)
        return vel
        
    @staticmethod
    def ortho_vel_glide_planes(vel: np.ndarray, normals: np.ndarray, eps_normal=0.05) -> np.ndarray:
        """ortho_vel_glide_planes: project velocity onto glide planes

        RULE CHANGE: eps_normal used to be 1e-10, compared against the raw
        norm -- effectively "discard only if exactly zero". ParaDiS instead
        uses FFACTOR_ORTH=0.05 (external/paradis/include/Constants.h:65),
        compared against the *squared* norm after Gram-Schmidt, in every
        mobility law that does this projection (MobilityLaw_FCC_0.c,
        MobilityLaw_BCC_glide_0.c, MobilityLaw_Relax.c, ...); exadis
        (core/exadis/src/mobility_types/mobility_glide.h, glide_constraints)
        carries the same 0.05 cutoff forward. 0.05 is a physical judgment
        call -- a glide plane whose normal is mostly cancelled by earlier
        ones is treated as redundant (not a real independent constraint)
        rather than kept and renormalized to unit length. pydis's tight
        1e-10 kept such near-redundant normals as full constraints, which
        over-restricts velocity at multi-arm junction nodes; found via a
        cross-code trajectory divergence in tests/full_runs/03_binary_junction,
        where a junction node's third glide-plane normal had squared norm
        0.046875 -- below 0.05 (discarded by paradis/exadis) but far above
        1e-10 (kept by pydis).
        eps_normal is named for a norm but compared against normals[i]'s
        squared value, matching FFACTOR_ORTH's own convention, not renamed
        to avoid disturbing existing callers that pass it positionally.
        """
        # first orthogonalize glide plane normals among themselves
        for i in range(normals.shape[0]):
            for j in range(i):
                normals[i] -= np.dot(normals[i], normals[j]) * normals[j]
            if np.dot(normals[i], normals[i]) < eps_normal:
                normals[i] = np.array([0.0, 0.0, 0.0])
            else:
                normals[i] /= np.linalg.norm(normals[i])

        # then orthogonalize velocity with glide plane normals
        vel -= np.dot( np.dot(vel, normals.T), normals )
        return vel
    
    @staticmethod
    def bcc_binary_junction(G: DisNet, tag: Tag) -> tuple:
        """bcc_binary_junction: the junction arm's direction at a binary junction node, and
           whether the junction lies in its parent planes, as ExaDiS' BCC_binary_junction_node;
           (None, False) for any other node
        """
        segs = G.neighbor_segments_dict(tag)
        if len(segs) != 3:
            return None, False
        R1 = G.nodes(tag).R
        tjunc, num_junc, num_glide, normals = None, 0, 0, []
        for nbr_tag, edge in segs.items():
            t = G.cell.closest_image(Rref=R1, R=G.nodes(nbr_tag).R) - R1
            if np.dot(t, t) < 1e-12:
                continue
            t = t / np.linalg.norm(t)
            b = edge.burg_vec_from(tag)
            if np.linalg.norm(b) > 1.0 + 1e-12:
                tjunc, num_junc = t, num_junc + 1
            else:
                n = np.cross(b, t)
                if np.dot(n, n) > 1.0e-2:
                    normals.append(n)
                num_glide += 1
        if num_junc != 1 or num_glide != 2:
            return None, False
        angletol = 5.0*np.pi/180.0
        if len(normals) == 2:
            tdir = np.cross(normals[0], normals[1])
            cos_angle = np.dot(tdir, tjunc) / np.linalg.norm(tdir) / np.linalg.norm(tjunc)
            planar = abs(abs(cos_angle) - 1.0) < 1.0 - np.cos(angletol)
        elif len(normals) == 1:
            planar = abs(np.dot(normals[0], tjunc)) < np.sin(angletol)
        else:
            planar = True  # both arms screw
        return tjunc, planar

    def NodeMobility_BCC_0B(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_BCC_0B: velocity from the drag matrix of the arms' edge, screw and climb
           mobilities, as ParaDiS' Mobility_BCC_0b (MobilityLaw_BCC_0b.c), with ExaDiS' handling
           of binary junction nodes and its velocity cap (mobility_bcc0b.h)
        """
        node = G.nodes(tag)
        segs = G.neighbor_segments_dict(tag)
        if node.constraint == DisNode.Constraints.PINNED_NODE or len(segs) < 2:
            return np.zeros(3)

        linejunc, tjunc = False, None
        if self.split3node:
            tjunc, planar = self.bcc_binary_junction(G, tag)
            # the node moves freely along a planar junction, or one that it unzips
            linejunc = tjunc is not None and (planar or np.dot(f, tjunc) > 0.0)

        eps = 1.0e-12
        I = np.eye(3)
        Btotal = np.zeros((3, 3))
        num_nonzero = 0
        R1 = node.R
        for nbr_tag, edge in segs.items():
            b = edge.burg_vec_from(tag)
            dr = G.cell.closest_image(Rref=R1, R=G.nodes(nbr_tag).R) - R1
            mag = np.linalg.norm(dr)
            if mag < eps:
                continue
            num_nonzero += 1
            dr = dr / mag
            half = 0.5 * mag
            if abs(b[0]*b[1]*b[2]) < eps:
                # a [0 0 1] arm moves less readily
                if linejunc:
                    Btotal += half * self.Bline * I
                elif len(segs) == 2:
                    Btotal += half * self.Beclimb * I
                else:
                    Btotal += half * (self.BlmBecl * np.outer(dr, dr) + self.Beclimb * I)
                continue
            # the drag of a screw arm, then corrected for its edge part
            Btotal += half * (self.BlmBsc * np.outer(dr, dr) + self.Bscrew * I)
            bMag2 = np.dot(b, b)
            costheta = np.dot(dr, b)
            costheta2 = (costheta*costheta) * (1.0/bMag2)
            if 1.0 - costheta2 > eps:
                n = np.cross(b, dr) * (1.0/np.sqrt((1.0 - costheta2) * bMag2))
                m = np.cross(n, dr)
                if linejunc:
                    n = n - np.dot(n, tjunc) * tjunc
                    n = n / np.linalg.norm(n)
                Bglide = 1.0 / np.sqrt(self.invBedge2 + (self.invBscrew2 - self.invBedge2) * costheta2)
                Bclimb = np.sqrt(self.Beclimb2 + (self.Bscrew2 - self.Beclimb2) * costheta2)
                Btotal += half * ((Bclimb - self.Bscrew) * np.outer(n, n) + (Bglide - self.Bscrew) * np.outer(m, m))
        if num_nonzero == 0:
            return np.zeros(3)

        return self.speed_cap(np.dot(matrix33_inverse(Btotal), f))

    def speed_cap(self, vel: np.ndarray) -> np.ndarray:
        """speed_cap: the velocity with its speed (m/s) capped smoothly at vmax, as ExaDiS'
           apply_velocity_cap (mobility.h); unchanged when vmax is not positive
        """
        speed = np.linalg.norm(vel) * self.burgmag
        if self.vmax > 0.0 and speed >= 1e-5:
            vel = vel * (1.0 / (1.0 + (speed/self.vmax)**10)**0.1)
        return vel

    def NodeMobility_FCC_0(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_FCC_0: the force over the arms' drag, from Medge to Mscrew by their
           character, kept in their glide planes, and along an arm off the {111} planes, as
           ParaDiS' Mobility_FCC_0 (MobilityLaw_FCC_0.c), with ExaDiS' velocity cap (mobility_fcc0.h)
        """
        node = G.nodes(tag)
        if node.constraint == DisNode.Constraints.PINNED_NODE:
            return np.zeros(3)
        R1 = node.R
        arms = [(G.cell.closest_image(Rref=R1, R=G.nodes(nbr_tag).R) - R1, edge.burg_vec_from(tag),
                 np.array(edge.plane_normal, dtype=float))
                for nbr_tag, edge in G.neighbor_segments_dict(tag).items()]
        if all(np.dot(dr, dr) <= 1.0e-12 for dr, _, _ in arms):
            return np.zeros(3)

        normals, lines = [], []
        for dr, _, n in arms:
            off111 = abs(abs(n[0]) - abs(n[1])) > 1.0e-4 or abs(abs(n[1]) - abs(n[2])) > 1.0e-4
            line = dr if off111 else np.zeros(3)
            nmag, lmag = np.linalg.norm(n), np.linalg.norm(line)
            normals.append(n / nmag if nmag > 0.0 else n)
            lines.append(line / lmag if lmag > 0.0 else line)
        planes = independent(normals)
        lines = [line for line in independent(lines) if line.any()]

        LtimesB = 0.0
        for dr, b, _ in arms:
            lr = np.linalg.norm(dr)
            if lr == 0.0:
                continue
            dangle = abs(np.dot(b / np.linalg.norm(b), dr / lr))
            LtimesB += lr / (self.Medge + (self.Mscrew - self.Medge) * dangle)
        LtimesB /= 2.0

        vel = f / LtimesB
        if sum(n.any() for n in planes) >= 3:
            return np.zeros(3)
        for n in planes:
            vel = orthogonalize(vel, n)
        if len(lines) >= 2:
            return np.zeros(3)
        if len(lines) == 1:
            # the arm off the {111} planes moves along itself, and only if it lies in every plane
            lc = lines[0]
            if any(abs(np.dot(n, lc)) > 0.05 for n in normals):
                return np.zeros(3)
            vel = np.dot(vel, lc) * lc
        return self.speed_cap(vel)

    def NodeMobility_SimpleGlide(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_SimpleGlide: node velocity equal node force divided by sum of arm length / 2
           To do: add glide constraints
        """
        node1 = G.nodes(tag)
        # set velocity of pinned nodes to zero
        if node1.constraint == DisNode.Constraints.PINNED_NODE:
            vel = np.zeros(3)
        else:
            R1 = node1.R.copy()
            Lsum = 0.0
            for nbr_tag, node2 in G.neighbors_dict(tag).items():
                R2 = node2.R.copy()
                # apply PBC
                R2 = G.cell.closest_image(Rref=R1, R=R2)
                Lsum += np.linalg.norm(R2-R1)
            # reciprocal-then-multiply, not divide-then-multiply: ExaDiS'
            # MobilityGlide::node_velocity (mobility_types/mobility_glide.h)
            # computes vi = P*(1.0/LtimesB * fi), and a/b is not always
            # bit-identical to (1.0/b)*a in IEEE754. Measured directly
            # against ExaDiS' own OneNodeMobility for
            # test5_topology_mode's step 194: this form, not f/(Lsum/2)*mob,
            # is what closes the gap (plan_topology_mode.md section 10).
            vel = (1.0/(Lsum/2.0)) * f * self.mob
            normals = np.array([edge.plane_normal for edge in G.neighbor_segments_dict(tag).values()])
            #print("Mobility_SimpleGlide: tag = %s, vel = %s, normals = %s"%(tag, str(vel), str(normals)))
            vel = self.ortho_vel_glide_planes(vel, normals)
            vel_norm = np.linalg.norm(vel)
            if vel_norm > self.vmax:
                vel *= self.vmax / vel_norm
        return vel
