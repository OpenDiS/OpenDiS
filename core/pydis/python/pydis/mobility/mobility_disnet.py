"""@package docstring
Mobility_DisNet: class for defining Mobility Laws

Provide mobility law functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode
from framework.disnet_manager import DisNetManager

# To do:
# 1. Implement OneNodeMobility in pydis, to use exadis_topology with pydis_mobility

class MobilityLaw:
    """MobilityLaw: class for mobility laws

    """
    def __init__(self, state: dict={}, mobility_law: str='Relax', vmax: float=1e9) -> None:
        self.mobility_law = mobility_law
        self.mob = state.get("mob", 1.0)
        self.vmax = vmax

        self.Mobility_Functions = {
            'Relax': self.Mobility_Relax,
            'SimpleGlide': self.Mobility_SimpleGlide }
        
    def Mobility(self, DM: DisNetManager, state: dict) -> dict:
        """Mobility: calculate node velocity according to mobility law function
        """
        G = DM.get_disnet(DisNet)
        if "nodeforces" in state and "nodeforcetags" in state:
            DisNet.convert_nodeforce_array_to_dict(state)

        nodeforce_dict = state["nodeforce_dict"]
        vel_dict = self.Mobility_Functions[self.mobility_law](G, nodeforce_dict)
        state["vel_dict"] = vel_dict

        # prepare nodeforces and nodeforce_tags arrays for compatibility with exadis
        state = DisNet.convert_nodevel_dict_to_array(state)
        return state

    def Mobility_Relax(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_Relax: node velocity equal node force
        """
        vel_dict = nodeforce_dict.copy()
        # set velocity of pinned nodes to zero
        for tag in G.all_nodes_tags():
            if G.nodes[tag].get('constraint') == DisNode.Constraints.PINNED_NODE:
                vel_dict[tag] = np.zeros(3)
        return vel_dict

    @staticmethod
    def ortho_vel_glide_planes(vel: np.ndarray, normals: np.ndarray, eps_normal = 1.0e-10) -> np.ndarray:
        """ortho_vel_glide_planes: project velocity onto glide planes
        """
        # first orthogonalize glide plane normals among themselves
        for i in range(normals.shape[0]):
            for j in range(i):
                normals[i] -= np.dot(normals[i], normals[j]) * normals[j]
            if np.linalg.norm(normals[i]) < eps_normal:
                normals[i] = np.array([0.0, 0.0, 0.0])
            else:
                normals[i] /= np.linalg.norm(normals[i])

        # then orthogonalize velocity with glide plane normals
        vel -= np.dot( np.dot(vel, normals.T), normals )
        return vel

    def Mobility_SimpleGlide(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_SimpleGlide: node velocity equal node force divided by sum of arm length / 2
           To do: add glide constraints
        """
        vel_dict = nodeforce_dict.copy()
        for tag in G.all_nodes_tags():
            node1 = G.nodes(tag)
            # set velocity of pinned nodes to zero
            if node1.constraint == DisNode.Constraints.PINNED_NODE:
                vel_dict[tag] = np.zeros(3)
            else:
                R1 = node1.R.copy()
                Lsum = 0.0
                for nbr_tag, node2 in G.neighbors_dict(tag).items():
                    R2 = node2.R.copy()
                    # apply PBC
                    R2 = G.cell.closest_image(Rref=R1, R=R2)
                    Lsum += np.linalg.norm(R2-R1)
                vel = vel_dict[tag] / (Lsum/2.0) * self.mob
                normals = np.array([edge.plane_normal for edge in G.neighbor_segments_dict(tag).values()])
                #print("Mobility_SimpleGlide: tag = %s, vel = %s, normals = %s"%(tag, str(vel), str(normals)))
                vel = self.ortho_vel_glide_planes(vel, normals)
                vel_norm = np.linalg.norm(vel)
                if vel_norm > self.vmax:
                    vel *= self.vmax / vel_norm
                vel_dict[tag] = vel
        return vel_dict
