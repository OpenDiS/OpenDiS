"""@package docstring
Mobility_DisNet: class for defining Mobility Laws

Provide mobility law functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager
from framework.mobility_base import MobilityLaw_Base

from typing import Tuple


class MobilityLaw(MobilityLaw_Base):
    """MobilityLaw: class for mobility laws
    """
    def __init__(self, state: dict={}, mobility_law: str='Relax', vmax: float=1e9) -> None:
        self.mobility_law = mobility_law
        self.mob = state.get("mob", 1.0)
        self.vmax = vmax

        self.NodeMobility_Functions = {
            'Relax': self.NodeMobility_Relax,
            'SimpleGlide': self.NodeMobility_SimpleGlide
        }
        
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
