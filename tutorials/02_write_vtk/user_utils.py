import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

def save_DisNet_to_vtp(DM: DisNetManager, filename: str):
    # Adapted from save_graph_for_paraview() written by Mychul Kim (@mckim2023) and Hanfeng Zhai (@hanfengzhai2)
    # ToDo: refactor this to be a member function of DisNet class
    import numpy as np
    from vtk import (vtkPoints, vtkCellArray, vtkPolyData, vtkLine, 
                    vtkFloatArray, vtkXMLPolyDataWriter, vtkIdList)
    import time

    start_time = time.time()
    if not filename.endswith('.vtp'): filename += '.vtp'

    G = DM.get_disnet(DisNet)
    # ToDo: refactor code to not use networkx as an intermediate data structure
    nx_digraph = G.to_networkx()
    node_list = list(nx_digraph.nodes())
    node_to_index = {node: i for i, node in enumerate(node_list)}
    positions = {node: nx_digraph.nodes[node]['R'] for node in node_list}
    
    # Create points for nodes
    points = vtkPoints()
    for node in node_list:
        pos = positions[node]
        points.InsertNextPoint(pos[0], pos[1], pos[2])
    
    # Process edges correctly with PBC
    point_count = len(node_list)  # Starting index for periodic image points
    
    # Store all points in a single array
    all_points = vtkPoints()
    for node in node_list:
        pos = positions[node]
        all_points.InsertNextPoint(pos[0], pos[1], pos[2])
    
    # Store all lines in a single array
    all_lines = vtkCellArray()
    
    for u, v in nx_digraph.edges():
        if u > v:
            R1 = positions[u]
            R2 = positions[v]
            image_pos = G.cell.closest_image(Rref=R1, R=R2)
            all_points.InsertNextPoint(image_pos[0], image_pos[1], image_pos[2])
            
            line = vtkLine()
            line.GetPointIds().SetId(0, node_to_index[u])
            line.GetPointIds().SetId(1, point_count)
            all_lines.InsertNextCell(line)
            
            point_count += 1

    # Create polydata with all points and lines
    polydata = vtkPolyData()
    polydata.SetPoints(all_points)
    polydata.SetLines(all_lines)
    
    # Write the file
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()

    elapsed_time = time.time() - start_time
    file_size_mb = os.path.getsize(filename) / (1024 * 1024)
    print(f"  Saved to {filename} ({file_size_mb:.2f} MB) in {elapsed_time:.2f} seconds")

class My_SimulateNetwork(SimulateNetwork):
    def step_write_files(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                #DM.write_json(os.path.join(self.write_dir, f'disnet_{istep}.json'))
                save_DisNet_to_vtp(DM, os.path.join(self.write_dir, f'disnet_{istep}.vtp'))
                if self.save_state:
                    with open(os.path.join(self.write_dir, f'state_{istep}.pickle'), 'wb') as file:
                        pickle.dump(state, file)
