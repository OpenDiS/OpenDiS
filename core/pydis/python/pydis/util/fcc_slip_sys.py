"""FCC slip systems: Schmid-Boas naming, plane/Burgers lookup, Schmid factors.

Slip system analysis for FCC metals, by Hanfeng Zhai (hzhai@stanford.edu) and
Wei Cai (caiwei@stanford.edu). Ported into pydis from the MLmat analysis code
with its tables unchanged.

The 12 systems are indexed 0..11 in Schmid-Boas order, A2 through D6, after
https://link.springer.com/article/10.1186/s41313-020-00024-y/tables/1. Index
12 stands for any (plane, Burgers) pair that is not one of them, such as a
junction.

Plain arrays only; nothing here knows about dislocation networks.
"""

import numpy as np

N_SLIP_SYSTEMS = 12
JUNCTION_ID = 12

SB_INDEX_NAMES = {
    0: 'A2', 1: 'A3', 2: 'A6', 3: 'B2', 4: 'B4', 5: 'B5',
    6: 'C1', 7: 'C3', 8: 'C5', 9: 'D1', 10: 'D4', 11: 'D6'}

# (plane normal, Burgers direction, name); the Burgers vector is half of it
SLIP_SYSTEM_DATA = [([-1,  1, 1], [ 0, -1, 1], 'A2'),
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
                    ([ 1, -1, 1], [ 1,  1, 0], 'D6')]

PLANE_NORMALS = [[-1.0,  1.0, 1.0], [1.0,  1.0, 1.0],
                 [-1.0, -1.0, 1.0], [1.0, -1.0, 1.0]]

BURGERS_DIRECTIONS = [[ 0.0, 1.0, 1.0], [0.0, -1.0, 1.0], [1.0, 0.0, 1.0],
                      [-1.0, 0.0, 1.0], [-1.0, 1.0, 0.0], [1.0, 1.0, 0.0]]

# (plane index, Burgers index) of slip system i, into the two lists above
PLANE_BURGERS_INDICES = [[0, 1], [0, 2], [0, 5], [1, 1], [1, 3], [1, 4],
                         [2, 0], [2, 2], [2, 4], [3, 0], [3, 3], [3, 5]]


def unit(v):
    """unit: v as a flat unit vector"""
    v = np.ravel(np.asarray(v, dtype=float))
    return v / np.linalg.norm(v)


def schmid_factor(edir, n, b):
    """schmid_factor: |cos(e,b) cos(e,n)| for uniaxial loading along edir"""
    e = unit(edir)
    return float(abs(np.dot(e, unit(b)) * np.dot(e, unit(n))))


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
        self.SB_index_names = dict(SB_INDEX_NAMES)

        self.slipSystem = {'slipPlane': {}, 'burgVect': {}, 'SB_index': {}}
        for i, (slip_plane, burg_vect, sb_index) in enumerate(SLIP_SYSTEM_DATA):
            self.slipSystem['slipPlane'][i] = slip_plane
            self.slipSystem['burgVect'][i] = 0.5 * np.asarray(burg_vect)
            self.slipSystem['SB_index'][i] = sb_index

        self.plane_normals = {i: np.array(vec).reshape(-1, 1)
                              for i, vec in enumerate(PLANE_NORMALS)}
        self.burgers_vectors = {i: np.array(vec).reshape(-1, 1) * 0.5
                                for i, vec in enumerate(BURGERS_DIRECTIONS)}

        self.sysID_to_plane_burgers = {
            i: (self.plane_normals[plane], self.burgers_vectors[burg])
            for i, (plane, burg) in enumerate(PLANE_BURGERS_INDICES)}
        self.sysID_to_plane_burgers[JUNCTION_ID] = (None, None)

        self.nb2sysID_tab = self.get_nb2sysID_tab(
            JUNCTION_ID, [(nid, bid, i) for i, (nid, bid)
                          in enumerate(PLANE_BURGERS_INDICES)])

        self.nb_pairs = [self.sysID_to_nb(i) for i in range(N_SLIP_SYSTEMS)]

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

    def name_to_sysID(self, name):
        """name_to_sysID: index of the slip system with Schmid-Boas name, e.g. 'A3'"""
        for i, sb_name in self.SB_index_names.items():
            if sb_name == name:
                return i
        raise ValueError("name_to_sysID: no slip system named %r" % (name,))

    def get_slip_Snb(self, edotDir):
        """get_slip_Snb: Schmid factors of the 12 systems for loading along edotDir

        Fills slipSystem['SchmidFactor'], adds the junction entries at index
        12, and sets self.schmid_factor to the largest Schmid factor.
        Returns (slipSystem, domSch), domSch being the index of that system;
        on a tie the lowest index wins.
        """
        factors = {}
        for i in range(N_SLIP_SYSTEMS):
            n = self.slipSystem['slipPlane'][i]
            b = self.slipSystem['burgVect'][i]
            if self.nb_to_sysID(n, b)[0] != i:
                raise ValueError("get_slip_Snb: mismatch in slip system %d" % i)
            factors[i] = schmid_factor(edotDir, n, b)

        domSch = max(range(N_SLIP_SYSTEMS), key=factors.get)
        self.schmid_factor = factors[domSch]

        # junction placeholder, deliberately kept out of the search above
        factors[JUNCTION_ID] = 1
        self.slipSystem['SchmidFactor'] = factors
        self.slipSystem['slipPlane'][JUNCTION_ID] = [0, 0, 0]
        self.slipSystem['burgVect'][JUNCTION_ID] = [0, 0, 0]
        return self.slipSystem, domSch

    def sysID_to_nb(self, sysID):
        """sysID_to_nb: (plane normal, Burgers vector) of slip system sysID

        Both as (3,1) columns; the junction index 12 gives (None, None).
        """
        if sysID not in self.sysID_to_plane_burgers:
            raise ValueError("sysID_to_nb: invalid slip system %r" % (sysID,))
        return self.sysID_to_plane_burgers[sysID]

    def nb_to_sysID(self, nPlane, burg, tol=1e-4):
        """nb_to_sysID: slip system index of each (plane, Burgers) row pair

        Matches up to sign. A pair that is not one of the 12 systems, such as
        a junction, maps to 12. Based on GetSlipSystem in WriteParaview.c
        (Hamedakh).
        """
        nPlane = np.asarray(nPlane, dtype=float).reshape(-1, 3)
        burg = np.asarray(burg, dtype=float).reshape(-1, 3)
        sysID = np.zeros(burg.shape[0])
        for i in range(burg.shape[0]):
            nid = self.find_collinear(nPlane[i, :], self.plane_normals, tol)
            bid = self.find_collinear(burg[i, :], self.burgers_vectors, tol)
            sysID[i] = (self.nb2sysID_tab[nid, bid] if nid >= 0 and bid >= 0
                        else JUNCTION_ID)
        return sysID

    def find_collinear(self, vector, candidates, tol):
        """find_collinear: key of the first candidate parallel to vector, else -1"""
        for idx, cand in candidates.items():
            if self.isCollinear(vector, cand.T[0], tol=tol):
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
