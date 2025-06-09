import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

from pydis.disnet import Tag

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

def save_DisNet_to_vtp(G, filename):
    # Adapted from save_graph_for_paraview() written by Mychul Kim (@mckim2023) and Hanfeng Zhai (@hanfengzhai2)
    # ToDo: refactor this to be a member function of DisNet class
    import numpy as np
    from vtk import (vtkPoints, vtkCellArray, vtkPolyData, vtkLine, 
                    vtkFloatArray, vtkXMLPolyDataWriter, vtkIdList)
    import time

    start_time = time.time()
    if not filename.endswith('.vtp'): filename += '.vtp'

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
    def run(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            os.makedirs(self.write_dir, exist_ok=True)

        G = DM.get_disnet(DisNet)
        if self.plot_freq != None:
            try: 
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)

        for tstep in range(self.max_step):
            self.step(DM, state)

            if self.write_freq != None:
                if tstep % self.write_freq == 0:
                    #DM.write_json(os.path.join(self.write_dir, f'disnet_{tstep}.json'))
                    print("I would like to write vtk file instead of json")
                    save_DisNet_to_vtp(G, os.path.join(self.write_dir, f'disnet_{tstep}.vtp'))
                    if self.save_state:
                        with open(os.path.join(self.write_dir, f'state_{tstep}.pickle'), 'wb') as file:
                            pickle.dump(state, file)

            if self.print_freq != None:
                if tstep % self.print_freq == 0:
                    print("step = %d dt = %e"%(tstep, self.timeint.dt))

            G = DM.get_disnet(DisNet)
            if self.plot_freq != None:
                if tstep % self.plot_freq == 0:
                    self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

        # plot final configuration
        if self.plot_freq != None:
            G = DM.get_disnet(DisNet)
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)
