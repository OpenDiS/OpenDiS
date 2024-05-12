"""@package docstring
Mobility_DisNet: class for defining Mobility Laws

Provide mobility law functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode
from framework.disnet_manager import DisNetManager

class MobilityLaw:
    """MobilityLaw: class for mobility laws

    """
    def __init__(self, params: dict={}, mobility_law: str='Relax') -> None:
        self.mobility_law = mobility_law
        self.mob = params.get("mob", 1.0)

        self.Mobility_Functions = {
            'Relax': self.Mobility_Relax,
            'SimpleGlide': self.Mobility_SimpleGlide }
        
    def Mobility(self, DM: DisNetManager, nodeforce_dict: dict):
        """Mobility: calculate node velocity according to mobility law function
        """
        G = DM.get_disnet(DisNet)
        return self.Mobility_Functions[self.mobility_law](G, nodeforce_dict)

    def Mobility_Relax(self, G: DisNet, nodeforce_dict: dict) -> dict:
        """Mobility_Relax: node velocity equal node force
        """
        vel_dict = nodeforce_dict.copy()
        # set velocity of pinned nodes to zero
        for tag in G.nodes:
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
        for tag in G.nodes:
            node1 = G.nodes[tag]
            # set velocity of pinned nodes to zero
            if node1.get('constraint') == DisNode.Constraints.PINNED_NODE:
                vel_dict[tag] = np.zeros(3)
            else:
                R1 = node1["R"]
                Lsum = 0.0
                for nbr_tag in G.neighbors(tag):
                    node2 = G.nodes[nbr_tag]
                    R2 = node2["R"]
                    # apply PBC
                    R2 = G.cell.map_to(R2, R1)
                    Lsum += np.linalg.norm(R2-R1)
                vel = vel_dict[tag] / (Lsum/2.0) * self.mob
                normals = np.array([G.edges[id]["plane_normal"] for id in G.edges(tag)])
                vel_dict[tag] = self.ortho_vel_glide_planes(vel, normals)
        return vel_dict
