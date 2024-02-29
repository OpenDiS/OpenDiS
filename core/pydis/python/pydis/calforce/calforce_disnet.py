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

def pkforcevec(sigext, segments):
    # return Peach-Koehler force vector for each segment
    # half of it should be assigned to each node
    nseg = len(segments)
    fpk = np.zeros((nseg, 3))
    for idx, segment in enumerate(segments):
        sigb = sigext @ segment["burg_vec"]
        dR = segment["R2"] - segment["R1"]
        fpk[idx] = np.cross(sigb, dR)
    return fpk

def selfforcevec_LineTension(MU, NU, Ec, segments, eps_L=1e-6):
    # to do: vectorize the calculations
    nseg = len(segments)
    fs0 = np.zeros((nseg, 3))
    fs1 = np.zeros((nseg, 3))
    omninv = 1.0/(1.0-NU)
    for idx, segment in enumerate(segments):
        dR = segment["R2"] - segment["R1"]
        L = np.linalg.norm(dR)
        if L < eps_L:
            continue
        t = dR / L
        bs = np.dot(segment["burg_vec"], t)
        bs2 = bs*bs
        bev = segment["burg_vec"] - bs*t
        be2 = np.sum(bev*bev)
        Score = 2.0*NU*omninv*Ec*bs
        LTcore = (bs2+be2*omninv)*Ec
        fs1[idx] = Score*bev - LTcore*t
    fs0 = -fs1
    return fs0, fs1


class CalForce:
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, mu: float, nu: float, a: float, Ec: float=0,
                 applied_stress: np.ndarray=np.zeros(6),
                 force_mode: str='Elasticity_SBA',
                 **kwargs) -> None:
        self.mu = mu
        self.nu = nu
        self.a = a
        self.Ec = Ec
        self.applied_stress = applied_stress
        self.force_mode = force_mode

        self.NodeForce_Functions = {
            'LineTension': self.NodeForce_LineTension,
            'Elasticity_SBA': self.NodeForce_Elasticity_SBA,
            'Elasticity_SBN1_SBA': self.NodeForce_Elasticity_SBN1_SBA }

    def NodeForce(self, DM: DisNetManager) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces in a dictionary

        Using different force calculation functions depending on force_mode
        """
        G = DM.get_disnet(DisNet)
        return self.NodeForce_Functions[self.force_mode](G)

    def NodeForce_from_SegForce(self, G: DisNet, segforce_dict: dict) -> dict:
        """NodeForce_from_SegForce: return nodal forces by assembling segment forces
        """
        nodeforce_dict = {}
        for tag in G.nodes:
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for segment in segforce_dict:
            tag1, tag2 = segment
            nodeforce_dict[tag1] += segforce_dict[segment][0:3]
            nodeforce_dict[tag2] += segforce_dict[segment][3:6]

        return nodeforce_dict

    def NodeForce_LineTension(self, G: DisNet) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from line tension in a dictionary

        Only Peach-Koehler force from external stress and line tension forces
        (from DDLab/src/segforcevec.m)
        Note: assuming G.seg_list already accounts for PBC
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fs0, fs1 = selfforcevec_LineTension(self.mu, self.nu, self.Ec, segments)
        fseg = np.hstack((fpk*0.5 + fs0, fpk*0.5 + fs1))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.nodes:
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            tag1 = segment["edge"][0]
            tag2 = segment["edge"][1]
            nodeforce_dict[tag1] += fseg[idx, 0:3]
            nodeforce_dict[tag2] += fseg[idx, 3:6]
            segforce_dict[segment["edge"]] = fseg[idx, :]

        return nodeforce_dict, segforce_dict

    def NodeForce_Elasticity_SBA(self, G: DisNet) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        Note: assuming G.seg_list already accounts for PBC
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.nodes:
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            tag1 = segment["edge"][0]
            tag2 = segment["edge"][1]
            nodeforce_dict[tag1] += fseg[idx, 0:3]
            nodeforce_dict[tag2] += fseg[idx, 3:6]

        nseg = len(segments)
        for i in range(nseg):
            for j in range(i, nseg):
                seg1 = segments[i]
                seg2 = segments[j]
                p1 = np.array(seg1["R1"])
                p2 = np.array(seg1["R2"])
                p3 = np.array(seg2["R1"])
                p4 = np.array(seg2["R2"])
                b12 = np.array(seg1["burg_vec"])
                b34 = np.array(seg2["burg_vec"])
                # apply PBC
                p2 = G.cell.map_to(p2, p1)
                p3 = G.cell.map_to(p3, p1)
                p4 = G.cell.map_to(p4, p3)
                f1, f2, f3, f4 = compute_segseg_force(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a)
                tag1 = seg1["edge"][0]
                tag2 = seg1["edge"][1]
                tag3 = seg2["edge"][0]
                tag4 = seg2["edge"][1]
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

        for idx, segment in enumerate(segments):
            segforce_dict[segment["edge"]] = fseg[idx, :]

        return nodeforce_dict, segforce_dict

    def NodeForce_Elasticity_SBN1_SBA(self, G: DisNet) -> Tuple[dict, dict]:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        Note: assuming G.seg_list already accounts for PBC
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fseg = np.hstack((fpk*0.5, fpk*0.5))

        nodeforce_dict, segforce_dict = {}, {}
        for tag in G.nodes:
            nodeforce_dict.update({tag: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            tag1 = segment["edge"][0]
            tag2 = segment["edge"][1]
            nodeforce_dict[tag1] += fseg[idx, 0:3]
            nodeforce_dict[tag2] += fseg[idx, 3:6]

        """ hardcode force_nint = 3
        """
        quad_points = np.array([-0.774596669241483, 0.0, 0.774596669241483])
        weights = np.array([0.555555555555556, 0.888888888888889, 0.555555555555556])
        print('quad_points = ', quad_points)
        print('weights = ', weights)

        nseg = len(segments)
        for i in range(nseg):
            for j in range(i, nseg):
                seg1 = segments[i]
                seg2 = segments[j]
                p1 = np.array(seg1["R1"])
                p2 = np.array(seg1["R2"])
                p3 = np.array(seg2["R1"])
                p4 = np.array(seg2["R2"])
                b12 = np.array(seg1["burg_vec"])
                b34 = np.array(seg2["burg_vec"])
                # apply PBC
                p2 = G.cell.map_to(p2, p1)
                p3 = G.cell.map_to(p3, p1)
                p4 = G.cell.map_to(p4, p3)
                f1, f2, f3, f4 = compute_segseg_force_SBN1_SBA(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a, quad_points, weights)
                tag1 = seg1["edge"][0]
                tag2 = seg1["edge"][1]
                tag3 = seg2["edge"][0]
                tag4 = seg2["edge"][1]
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

        for idx, segment in enumerate(segments):
            segforce_dict[segment["edge"]] = fseg[idx, :]

        return nodeforce_dict, segforce_dict
