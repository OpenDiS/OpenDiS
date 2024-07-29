"""@package docstring
CalForce_DisNet: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object
"""

import numpy as np
from typing import Tuple
from ..disnet import DisNet
from framework.disnet_manager import DisNetManager

try:
    from .compute_stress_force_analytic_paradis import compute_segseg_force_vec, compute_segseg_force
    from .compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_vec, compute_segseg_force_SBN1
    from .compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_SBA
    from .compute_stress_analytic_paradis       import compute_seg_stress_coord_dep, compute_seg_stress_coord_indep
except ImportError:
    # use python version instead
    # To do: put import commands here
    print("pydis_lib not found, using python version for force calculation")

from .compute_stress_force_analytic_python  import python_segseg_force_vec

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
        sigb = sigext @ burg_vecs[i]
        dR = R2[i,:] - R1[i,:]
        fpk[i,:] = np.cross(sigb, dR)
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
        L = np.linalg.norm(dR)
        if L < eps_L:
            continue
        t = dR / L
        bs = np.dot(burg_vecs[i,:], t)
        bs2 = bs*bs
        bev = burg_vecs[i] - bs*t
        be2 = np.sum(bev*bev)
        Score = 2.0*NU*omninv*Ec*bs
        LTcore = (bs2+be2*omninv)*Ec
        fs1[i,:] = Score*bev - LTcore*t
    fs0 = -fs1
    return fs0, fs1

class CalForce:
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}, Ec: float=None,
                 force_mode: str='Elasticity_SBA') -> None:
        self.mu = state.get("mu", 1.0)
        self.nu = state.get("nu", 0.3)
        self.a =  state.get("a", 0.01)
        self.Ec = self.mu/4.0/np.pi*np.log(self.a/0.1) if Ec is None else Ec
        self.force_mode = force_mode

        self.NodeForce_Functions = {
            'LineTension': self.NodeForce_LineTension,
            'Elasticity_SBA': self.NodeForce_Elasticity_SBA,
            'Elasticity_SBN1_SBA': self.NodeForce_Elasticity_SBN1_SBA }

    def NodeForce(self, DM: DisNetManager, state: dict) -> dict:
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

    def NodeForce_from_SegForce(self, G: DisNet, state: dict) -> dict:
        """NodeForce_from_SegForce: return nodal forces by assembling segment forces
        """
        segforce_dict = state["segforce_dict"]
        nodeforce_dict = {}
        for tag in G.all_nodes():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for segment in segforce_dict:
            tag1, tag2 = segment
            nodeforce_dict[tag1] += segforce_dict[segment][0:3]
            nodeforce_dict[tag2] += segforce_dict[segment][3:6]
        state["nodeforce_dict"] = nodeforce_dict

        # prepare nodeforces and nodeforce_tags arrays for compatibility with exadis
        state = DisNet.convert_nodeforce_dict_to_array(state)

        return state

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
        for tag in G.all_nodes_dict():
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
        for tag in G.all_nodes():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            nodeforce_dict[tag1] += fseg[i, 0:3]
            nodeforce_dict[tag2] += fseg[i, 3:6]

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
                f1, f2, f3, f4 = compute_segseg_force(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a)
                tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
                tag3, tag4 = tuple(source_tags[j]), tuple(target_tags[j])
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
        for tag in G.all_nodes():
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            nodeforce_dict[tag1] += fseg[i, 0:3]
            nodeforce_dict[tag2] += fseg[i, 3:6]

        """ hardcode force_nint = 3
        """
        quad_points = np.array([-0.774596669241483, 0.0, 0.774596669241483])
        weights = np.array([0.555555555555556, 0.888888888888889, 0.555555555555556])
        print('quad_points = ', quad_points)
        print('weights = ', weights)

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
                f1, f2, f3, f4 = compute_segseg_force_SBN1_SBA(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a, quad_points, weights)
                tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
                tag3, tag4 = tuple(source_tags[j]), tuple(target_tags[j])
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
