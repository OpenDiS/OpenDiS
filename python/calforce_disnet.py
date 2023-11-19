"""@package docstring
CalForce_DisNet: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object
"""

import numpy as np
from disnet import DisNet

from compute_stress_force_analytic_paradis import compute_segseg_force_vec, compute_segseg_force
from compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_vec, compute_segseg_force_SBN1
from compute_stress_force_analytic_paradis import compute_segseg_force_SBN1_SBA
from compute_stress_force_analytic_python  import python_segseg_force_vec
from compute_stress_analytic_paradis       import compute_seg_stress_coord_dep, compute_seg_stress_coord_indep

def voigt_vector_to_tensor(voigt_vector):
    return np.array([[voigt_vector[0], voigt_vector[5], voigt_vector[4]],
                     [voigt_vector[5], voigt_vector[1], voigt_vector[3]],
                     [voigt_vector[4], voigt_vector[3], voigt_vector[2]]])

def pkforcevec(sigext, segments):
    nseg = len(segments)
    fpk = np.zeros((nseg, 3))
    for idx, segment in enumerate(segments):
        sigb = sigext @ segment["burg_vec"]
        dR = segment["R2"] - segment["R1"]
        fpk[idx] = np.cross(sigb, dR)
    return fpk

def selfforcevec_LineTension(MU, NU, Ec, segments):
    # to do: vectorize the calculations
    nseg = len(segments)
    fs0 = np.zeros((nseg, 3))
    fs1 = np.zeros((nseg, 3))
    omninv = 1.0/(1.0-NU)
    for idx, segment in enumerate(segments):
        dR = segment["R2"] - segment["R1"]
        L = np.linalg.norm(dR)
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

    def NodeForce(self, G: DisNet) -> dict:
        """NodeForce: return nodal forces in a dictionary

        Using different force calculation functions depending on force_mode
        """
        return self.NodeForce_Functions[self.force_mode](G)

    def NodeForce_LineTension(self, G: DisNet) -> dict:
        """NodeForce: return nodal forces from line tension in a dictionary

        Only Peach-Koehler force from external stress and line tension forces
        (from DDLab/src/segforcevec.m)
        Note: PBC not handled
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fs0, fs1 = selfforcevec_LineTension(self.mu, self.nu, self.Ec, segments)
        fseg = np.hstack((fpk + fs0, fpk + fs1))

        nodeforce_dict = {}
        for node in G.nodes():
            nodeforce_dict.update({node: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            node1 = segment["edge"][0]
            node2 = segment["edge"][1]
            nodeforce_dict[node1] += fseg[idx, 0:3]
            nodeforce_dict[node2] += fseg[idx, 3:6]

        return nodeforce_dict

    def NodeForce_Elasticity_SBA(self, G: DisNet) -> dict:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        To do: fseg can be return if needed
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fseg = np.hstack((fpk, fpk))

        nodeforce_dict = {}
        for node in G.nodes():
            nodeforce_dict.update({node: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            node1 = segment["edge"][0]
            node2 = segment["edge"][1]
            nodeforce_dict[node1] += fseg[idx, 0:3]
            nodeforce_dict[node2] += fseg[idx, 3:6]

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
                f1, f2, f3, f4 = compute_segseg_force(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a)
                node1 = seg1["edge"][0]
                node2 = seg1["edge"][1]
                node3 = seg2["edge"][0]
                node4 = seg2["edge"][1]
                if i == j:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    nodeforce_dict[node1] += f1
                    nodeforce_dict[node2] += f2
                else:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    fseg[j, 0:3] += f3
                    fseg[j, 3:6] += f4
                    nodeforce_dict[node1] += f1
                    nodeforce_dict[node2] += f2
                    nodeforce_dict[node3] += f3
                    nodeforce_dict[node4] += f4

        return nodeforce_dict

    def NodeForce_Elasticity_SBN1_SBA(self, G: DisNet) -> dict:
        """NodeForce: return nodal forces from external stress and elastic interactions

        (from ParaDiS)
        To do: fseg can be return if needed
        """
        segments = G.seg_list()
        sigext = voigt_vector_to_tensor(self.applied_stress)
        fpk = pkforcevec(sigext, segments)
        fseg = np.hstack((fpk, fpk))

        nodeforce_dict = {}
        for node in G.nodes():
            nodeforce_dict.update({node: np.array([0.0,0.0,0.0])})
        for idx, segment in enumerate(segments):
            node1 = segment["edge"][0]
            node2 = segment["edge"][1]
            nodeforce_dict[node1] += fseg[idx, 0:3]
            nodeforce_dict[node2] += fseg[idx, 3:6]

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
                f1, f2, f3, f4 = compute_segseg_force_SBN1_SBA(p1, p2, p3, p4, b12, b34, self.mu, self.nu, self.a, quad_points, weights)
                node1 = seg1["edge"][0]
                node2 = seg1["edge"][1]
                node3 = seg2["edge"][0]
                node4 = seg2["edge"][1]
                if i == j:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    nodeforce_dict[node1] += f1
                    nodeforce_dict[node2] += f2
                else:
                    fseg[i, 0:3] += f1
                    fseg[i, 3:6] += f2
                    fseg[j, 0:3] += f3
                    fseg[j, 3:6] += f4
                    nodeforce_dict[node1] += f1
                    nodeforce_dict[node2] += f2
                    nodeforce_dict[node3] += f3
                    nodeforce_dict[node4] += f4

        return nodeforce_dict

