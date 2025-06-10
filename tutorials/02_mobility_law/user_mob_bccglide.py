import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell
from pydis import Topology
from pydis.disnet import Tag
from framework.mobility_base import MobilityLaw_Base


class My_MobilityLaw(MobilityLaw_Base):
    """MobilityLaw: class for mobility laws
    """
    def __init__(self, state: dict={}, mobility_law: str='Relax', vmax: float=1e9) -> None:
        self.mobility_law = mobility_law
        self.mob = state.get("mob", 1.0)
        self.vmax = vmax
        self.mob_edge = state.get("mob_edge", 1.0)
        self.mob_screw = state.get("mob_screw", 0.1)

        self.NodeMobility_Functions = {
            'BCCGlide': self.NodeMobility_BCCGlide
        }

    # Implement method required by MobilityLaw_Base
    def Mobility(self, DM: DisNetManager, state: dict) -> dict:
        """Mobility: compute all nodal velocities and store them in the state dictionary
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

    # Implement method required by MobilityLaw_Base
    def OneNodeMobility(self, N: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
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

    @staticmethod
    def ortho_vel_glide_planes(vel: np.ndarray, normals: np.ndarray, eps_normal=1.0e-10) -> np.ndarray:
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

    def NodeMobility_BCCGlide(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_BCCGlide: node velocity equal local mobility times node force divided by sum of arm length / 2
               local mobility depends on local dislocation character angle (for 1-arm and 2-arm nodes)
               multi-arm nodes still use self.mob as SimpleGlide
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
            if G.out_degree(tag) > 2:
                # for multi-arm nodes, use self.mob (same as in SimpleGlide)
                local_mob = self.mob
            else:
                # compute local mobility value based on dislocation character angle (only for nodes of degree <= 2)
                cos_theta_sq_avg = 0
                for nbr_tag, edge_attr in G.neighbor_segments_dict(tag).items(): # This line is different
                    local_burg = edge_attr.burg_vec_from(tag).copy()
                    R2 = G.nodes(nbr_tag).R
                    R2 = G.cell.closest_image(Rref=R1, R=R2) # handle PBC
                    local_line_vec = R2 - R1
                    dot_product = np.dot(local_line_vec, local_burg)
                    burg_sq = np.sum(np.dot(local_burg, local_burg))
                    L_sq = np.sum(np.dot(local_line_vec, local_line_vec))
                    local_L = np.sqrt(L_sq)
                    cos_theta_sq = dot_product**2 / (burg_sq * L_sq)
                    cos_theta_sq_avg += cos_theta_sq * local_L / Lsum
                local_mob = (1 - cos_theta_sq_avg) * self.mob_edge + cos_theta_sq_avg * self.mob_screw

            vel = f / (Lsum/2.0) * local_mob
            normals = np.array([edge.plane_normal for edge in G.neighbor_segments_dict(tag).values()])
            #print("Mobility_BCCGlide: tag = %s, vel = %s, normals = %s"%(tag, str(vel), str(normals)))
            vel = self.ortho_vel_glide_planes(vel, normals)
            vel_norm = np.linalg.norm(vel)
            if vel_norm > self.vmax:
                vel *= self.vmax / vel_norm
        return vel
    

class My_Topology(Topology):
    """Topology: class for selecting and handling multi node splitting
    """
    def __init__(self, state: dict={}, split_mode: str='MaxDiss', **kwargs) -> None:
        self.split_mode = split_mode
        self.force = kwargs.get('force')
        if not self.force.__module__.split('.')[0] in ['pydis', 'pyexadis_base']:
            print("Topology: force.__module__ = ", self.force.__module__)
            raise ValueError("Topology: force must come compatible modules")
        self.mobility = kwargs.get('mobility')
        if not self.mobility.__module__.split('.')[0] in ['pydis', 'user_mob_bccglide']:
            print("Topology: mobility.__module__ = ", self.mobility.__module__)
            raise ValueError("Topology: mobility must come compatible modules")
        self.Handle_Functions = {
            'MaxDiss': self.Handle_MaxDiss }
