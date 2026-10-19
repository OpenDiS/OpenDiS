'''''''''''''''''''''''''''''''''''''''''''''''''''''
    Slip system analysis for FCC metals 
    (Hanfeng Zhai | hzhai@stanford.edu)
    (Wei Cai | caiwei@stanford.edu)
'''''''''''''''''''''''''''''''''''''''''''''''''''''
import numpy as np
import sys

class FCC_slip_sys():
    def print_info(self):
        print(self)

    def __str__(self):
        return ('    ##########################################################\n'
                '    #                                                        #\n'
                '    #                 FCC Slip System Analysis               #\n'
                '    #                                                        #\n'
                '    ##########################################################\n')

    def __init__(self):

        self.SB_index_names = {
            0: 'A2', 1: 'A3', 2: 'A6', 3: 'B2', 4:  'B4',  5: 'B5',
            6: 'C1', 7: 'C3', 8: 'C5', 9: 'D1', 10: 'D4', 11: 'D6'}

        # https://link.springer.com/article/10.1186/s41313-020-00024-y/tables/1
        slip_system_data = [([-1,  1, 1], [ 0, -1, 1], 'A2'),
                            ([-1,  1, 1], [ 1,  0, 1], 'A3'),
                            ([-1,  1, 1], [ 1,  1, 0], 'A6'),
                            ([ 1,  1, 1], [ 0, -1, 1], 'B2'),
                            ([ 1,  1, 1], [-1,  0, 1], 'B4'),
                            ([ 1,  1, 1], [-1,  1, 0], 'B5'),
                            ([-1, -1, 1], [ 0,  1, 1], 'C1'),
                            ([-1, -1, 1], [ 1,  0, 1], 'C3'),
                            ([-1, -1, 1], [-1,  1, 0], 'C5'),
                            ([ 1, -1, 1], [ 0,  1, 1], 'D1'),
                            ([ 1, -1, 1], [-1,  0, 1], 'D4'),
                            ([ 1, -1, 1], [ 1,  1, 0], 'D6')] # TODO: put 0.5 inside | add comment for n, b, SB index

        self.slipSystem = {'slipPlane': {}, 'burgVect': {}, 'SB_index': {}}

        for i, (slip_plane, burg_vect, sb_index) in enumerate(slip_system_data):
            self.slipSystem['slipPlane'][i] = slip_plane
            self.slipSystem['burgVect'][i] = 0.5 * np.asarray(burg_vect) # TODO: no need for 0.5
            self.slipSystem['SB_index'][i] = sb_index

        self.plane_normals = {i: np.array(vec).reshape(-1, 1) for i, vec in enumerate([
            [-1.0,  1.0, 1.0], [1.0,  1.0, 1.0], 
            [-1.0, -1.0, 1.0], [1.0, -1.0, 1.0]])}

        self.burgers_vectors = {i: np.array(vec).reshape(-1, 1) * 0.5 for i, vec in enumerate([
            [0.0,  1.0, 1.0], [0.0, -1.0, 1.0], [1.0,  0.0, 1.0], 
            [-1.0, 0.0, 1.0], [-1.0, 1.0, 0.0], [1.0,  1.0, 0.0]])}

        plane_burgers_indices = np.array([[0, 1], [0, 2], [0, 5], [1, 1], [1, 3], [1, 4],
                                          [2, 0], [2, 2], [2, 4], [3, 0], [3, 3], [3, 5]])
        self.sysID_to_plane_burgers = {
            i: (self.plane_normals[plane], self.burgers_vectors[burg])
            for i, (plane, burg) in enumerate(plane_burgers_indices)}
        self.sysID_to_plane_burgers[len(plane_burgers_indices)] = (None, None)

        self.nb_pairs = [self.sysID_to_nb(i) for i in range(12)]

        self.ref_perm_tab = np.array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11],
                                      [ 3,  4,  5,  0,  1,  2,  9, 10, 11,  6,  7,  8],
                                      [ 6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5],
                                      [ 9, 10, 11,  6,  7,  8,  3,  4,  5,  0,  1,  2],
                                      [ 8,  6,  7,  5,  3,  4, 11,  9, 10,  2,  0,  1],
                                      [ 2,  0,  1, 11,  9, 10,  5,  3,  4,  8,  6,  7],
                                      [11,  9, 10,  2,  0,  1,  8,  6,  7,  5,  3,  4],
                                      [ 5,  3,  4,  8,  6,  7,  2,  0,  1, 11,  9, 10],
                                      [10, 11,  9,  4,  5,  3,  1,  2,  0,  7,  8,  6],
                                      [ 1,  2,  0,  7,  8,  6, 10, 11,  9,  4,  5,  3],
                                      [ 4,  5,  3, 10, 11,  9,  7,  8,  6,  1,  2,  0],
                                      [ 7,  8,  6,  1,  2,  0,  4,  5,  3, 10, 11,  9],
                                      [ 3,  5,  4,  0,  2,  1,  6,  8,  7,  9, 11, 10],
                                      [ 1,  0,  2,  7,  6,  8,  4,  3,  5, 10,  9, 11],
                                      [ 2,  1,  0, 11, 10,  9,  8,  7,  6,  5,  4,  3],
                                      [ 8,  7,  6,  5,  4,  3,  2,  1,  0, 11, 10,  9],
                                      [10,  9, 11,  4,  3,  5,  7,  6,  8,  1,  0,  2],
                                      [ 0,  2,  1,  3,  5,  4,  9, 11, 10,  6,  8,  7],
                                      [ 9, 11, 10,  6,  8,  7,  0,  2,  1,  3,  5,  4],
                                      [ 6,  8,  7,  9, 11, 10,  3,  5,  4,  0,  2,  1],
                                      [ 5,  4,  3,  8,  7,  6, 11, 10,  9,  2,  1,  0],
                                      [11, 10,  9,  2,  1,  0,  5,  4,  3,  8,  7,  6],
                                      [ 4,  3,  5, 10,  9, 11,  1,  0,  2,  7,  6,  8],
                                      [ 7,  6,  8,  1,  0,  2, 10,  9, 11,  4,  3,  5]])

        self.interaction_pairs = {'S': [(i, i) for i in range(12)],
                                  'CP': [(0, 1), (0, 2), (1, 2), (3, 4 ), (3, 5 ), (4, 5),
                                         (6, 7), (6, 8), (7, 8), (9, 10), (9, 11), (10, 11)],
                                  'CL': [(0, 3), (2, 11), (4, 10), (5, 8), (6, 9), (1, 7)],
                                  'G': [(0, 4), (0, 5), (0, 7), (0, 11), (1, 3), (1, 6),
                                        (1, 8), (1, 11), (2, 3), (2, 7), (2, 9), (2, 10),
                                        (3, 8), (3, 10), (4, 8), (4, 9), (4, 11), (5, 6),
                                        (5, 7), (5, 10), (6, 10), (6, 11), (7, 9), (8, 9)],
                                  'H': [(0, 6), (0, 9), (1, 4), (1, 10), (2, 5), (2, 8),
                                        (3, 6), (3, 9), (4, 7), (5, 11), (7, 10), (8, 11)],
                                  'L': [(0, 8), (0, 10), (1, 5), (1, 9), (2, 4), (2, 6),
                                        (3, 7), (3, 11), (4, 6), (5, 9), (7, 11), (8, 10)]
                                  }

    def get_slip_Snb(self, edotDir):
        # creates 'SchmidFactor' in the slipSystem dict
        # add 12th entries in the burgVect & slipPlane for junctions
        # returns the domSch
        self.schmid_factor = 0; domSch = -1
        self.slipSystem['SchmidFactor'] = {}

        for i in range(12):
            UnitedotDir  = np.asarray(edotDir).reshape(-1,3) / np.linalg.norm(edotDir)
            UnitBurgVect = np.asarray(self.slipSystem['burgVect'][i]).reshape(-1,3) / np.linalg.norm(self.slipSystem['burgVect'][i])
            UnitNormVect = np.asarray(self.slipSystem['slipPlane'][i]).reshape(-1,3) / np.linalg.norm(self.slipSystem['slipPlane'][i])

            if ( i != self.nb_to_sysID(UnitNormVect, UnitBurgVect) ):
                sys.exit("mismatch in slip system")

            self.slipSystem['SchmidFactor'][i] = float( abs( np.inner(UnitedotDir, UnitBurgVect ) * \
                                                             np.inner(UnitedotDir, UnitNormVect ) ) )
            if self.slipSystem['SchmidFactor'][i]>domSch:
                self.schmid_factor = self.slipSystem['SchmidFactor'][i]

        self.slipSystem['SchmidFactor'][12] = 1      # Junctions # DANGER: why is it 1 ? Check if you use this value anywhere. Max( Sch ) returns 1, which will cause a bug
        self.slipSystem['slipPlane'][12] = [0, 0, 0] # Junctions
        self.slipSystem['burgVect'][12]  = [0, 0, 0] # Junctions
        return self.slipSystem, domSch

    def sysID_to_nb(self, slip_sys):
        if slip_sys in self.sysID_to_plane_burgers:
            self.slip_plane, self.burg_vec = self.sysID_to_plane_burgers[slip_sys]
        else:
            raise ValueError("Invalid slip system")
        return self.slip_plane, self.burg_vec

    def nb_to_sysID(self, nPlane, burg, tol=1e-4):
        # Finds the slip system of a FCC dislocation
        # Based on GetSlipSystem in WriteParaview.c (Hamedakh)
        nPlane = nPlane.reshape( -1, 3 )
        burg = burg.reshape( -1, 3 )
        sysID = np.zeros(( burg.shape[0] ))

        for i in range( burg.shape[0] ):
            nid = self.find_collinear(nPlane[i, :], self.plane_normals, tol)
            bid = self.find_collinear(burg[i, :], self.burgers_vectors, tol)

            sysID_index = [(0, 1, 0), (0, 2, 1),  (0, 5, 2),
                           (1, 1, 3), (1, 3, 4),  (1, 4, 5),
                           (2, 0, 6), (2, 2, 7),  (2, 4, 8),
                           (3, 0, 9), (3, 3, 10), (3, 5, 11)]

            self.nb2sysID_tab = self.get_nb2sysID_tab(12, sysID_index)

            if 0 <= int(nid) < 4 and 0 <= int(bid) < 6:
                sys_id = self.nb2sysID_tab[int(nid), int(bid)]
            else:
                sys_id = 12  # Default value for junctions

            sysID[i] = sys_id
        return sysID

    def find_collinear(self, vector, candidates, tol):
        for idx, _ in enumerate(candidates):
            if self.isCollinear(vector, candidates[idx].T[0], tol=tol):
                return idx
        return -1

    @staticmethod
    def get_nb2sysID_tab(default_val, assignments):
        nb2sysID_tab = np.full((4, 6), default_val)
        for nid, bid, sysID in assignments:
            nb2sysID_tab[nid, bid] = sysID
        return nb2sysID_tab

    @staticmethod
    def isCollinear(v1, v2, tol=1e-4):
        dot_v12_norm = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return abs(abs(dot_v12_norm) - 1.0) < tol