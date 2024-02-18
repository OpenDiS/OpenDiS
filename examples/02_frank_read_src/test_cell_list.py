import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])

from pydis.disnet import DisNet, DisNode, Cell, CellList

L = 10
R = (np.random.rand(1000,3)-0.15) * 1.1 * L

cell = Cell(h=L*np.eye(3), is_periodic=[True,True,True])
n_div = [5, 9, 3]
cell_list = CellList(cell=cell, n_div=n_div)

cell_list.sort_points_to_list(R)

R_mapped = cell.map_to(R, np.zeros(3))

points_within_bounds = True

for i in range(n_div[0]):
    for j in range(n_div[1]):
        for k in range(n_div[2]):
            idx = cell_list.get_objs_in_cell([i,j,k])
            xmin, xmax = (np.array([1.0*i/n_div[0], 1.0*(i+1)/n_div[0]]) - 0.5)*L
            ymin, ymax = (np.array([1.0*j/n_div[1], 1.0*(j+1)/n_div[1]]) - 0.5)*L
            zmin, zmax = (np.array([1.0*k/n_div[2], 1.0*(k+1)/n_div[2]]) - 0.5)*L
            R_mapped_in_cell = R_mapped[idx]
            points_within_bounds = points_within_bounds and all(np.min(R_mapped_in_cell, axis=0) >= [xmin, ymin, zmin])
            points_within_bounds = points_within_bounds and all(np.max(R_mapped_in_cell, axis=0) <  [xmax, ymax, zmax])

print('list of points around R[0] =', cell_list.get_objs_in_nbr_cells(0))

print('iterate over all neighbor pairs:')
for i, j in cell_list.iterate_nbr_pairs(use_cell_list=True):
    print(i, j)

print("points_within_bounds = ", points_within_bounds)
if points_within_bounds:
    print("test" + '\033[32m' + " PASSED" + '\033[0m')
else:
    print("test" + '\033[31m' + " FAILED" + '\033[0m')