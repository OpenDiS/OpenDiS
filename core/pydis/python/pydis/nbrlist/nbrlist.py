"""@package docstring
NbrList: class for constructing neighbor list

Provide neighbor list iterator given a set of points (and cell)
"""

import numpy as np
from ..disnet import DisNet, Cell

class CellList:
    """CellList: class for cell list

    Implements cell list for efficient neighbor search
    Note: only works for PBC in all 3 directions
    """
    def __init__(self, cell: Cell=None, n_div: np.array=[3,3,3]) -> None:
        # To do: comput ndiv from cell and cutoff
        self.cell = Cell() if cell is None else cell
        self.n_div = n_div
        self._cell_list = [[[ [] for n2 in range(self.n_div[2])] for n1 in range(self.n_div[1])] for n0 in range(self.n_div[0])]
        self._cell_indices = None

    def get_cell_index(self, R: np.ndarray) -> np.ndarray:
        """get_cell_index: get cell index of a position
        """
        s = np.dot(self.cell.hinv, R.T - np.broadcast_to(self.cell.center(), shape=R.shape).T).T
        s -= np.round(s)
        ind = np.floor((s+0.5)*self.n_div).astype(int)
        ind = np.mod(ind, np.array(self.n_div))
        return ind

    def sort_points_to_list(self, R: np.ndarray) -> None:
        """build: build cell list

        Note: only do this when PBC is applied in all 3 directions
        """
        self._cell_list = [[[ [] for n2 in range(self.n_div[2])] for n1 in range(self.n_div[1])] for n0 in range(self.n_div[0])]
        if all(self.cell.is_periodic):
            self._cell_indices = self.get_cell_index(R)
            for i, ind in enumerate(self._cell_indices):
                self._cell_list[ind[0]][ind[1]][ind[2]].append(i)
        else:
            self._cell_indices = [None]*len(R)
        return

    def get_objs_in_cell(self, ind: np.ndarray) -> list:
        """get_indices_in_cell: get indices in a cell
        """
        return self._cell_list[ind[0]][ind[1]][ind[2]]

    def get_objs_in_nbr_cells(self, obj_id: int) -> list:
        """get_indices_in_cell: get indices in the same cell of a point and all neighbor cells
        """
        ind = self._cell_indices[obj_id]
        nbr_ids = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    ind_nbr = np.array([ind[0]+i, ind[1]+j, ind[2]+k])
                    ind_nbr = np.mod(ind_nbr, self.n_div)
                    nbr_ids.extend(self.get_objs_in_cell(ind_nbr))
        return nbr_ids

    def iterate_nbr_pairs(self, use_cell_list: bool=True):
        """iterate_nbr_pairs: iterate over all pairs of segments in the same cell and all neighbor cells

        Note: self pairs (i, i) not included
        """
        if not use_cell_list:
            n = len(self._cell_indices)
            for i in range(n):
                for j in range(i+1, n):
                    yield i, j
        else:
            n = len(self._cell_indices)
            for i in range(n):
                nbr_ids = self.get_objs_in_nbr_cells(i)
                for j in nbr_ids:
                    if i < j: yield i, j
