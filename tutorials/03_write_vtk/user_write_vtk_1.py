import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from framework.disnet_manager import DisNetManager
from pydis import DisNet, SimulateNetwork


def save_DisNet_to_vtp(DM: DisNetManager, filename: str, pbc_tol=1e-10):
    import xml.etree.ElementTree as ET
    import time

    start_time = time.time()
    if not filename.endswith('.vtp'): filename += '.vtp'

    G = DM.get_disnet(DisNet)
    node_list = list(G.all_nodes_tags())
    node_to_index = {node: i for i, node in enumerate(node_list)}
    positions = {node: G.nodes(node).view()['R'] for node in node_list}
    
    # Process edges correctly with PBC
    point_count = len(node_list)  # Starting index for periodic image points
    
    # Store all points in a single array
    all_points = [ positions[node] for node in node_list ]
    
    # Store all lines in a single array
    connectivity = []

    for u, v in G.all_segments_tags():
        R1 = positions[u]
        R2 = positions[v]
        image_pos = G.cell.closest_image(Rref=R1, R=R2)

        if np.max(np.abs(image_pos - R2)) > pbc_tol:
            all_points.append(image_pos)
            line = [node_to_index[u], point_count]
            point_count += 1
        else:
            line = [node_to_index[u], node_to_index[v]]

        connectivity.extend(line)

    num_points = len(all_points)
    num_lines = G.num_segments()
    offsets = list(range(2, 2*num_lines+2, 2))

    # Format points for DataArray text content
    points_str = " ".join(f"{x} {y} {z}" for x, y, z in all_points)
    connectivity_str = " ".join(map(str, connectivity))
    offsets_str = " ".join(map(str, offsets))

    # Root element
    vtk_file = ET.Element("VTKFile")
    vtk_file.set("type", "PolyData")
    vtk_file.set("version", "1.0")
    vtk_file.set("byte_order", "LittleEndian")
    vtk_file.set("header_type", "UInt64")

    poly_data = ET.SubElement(vtk_file, "PolyData")
    piece = ET.SubElement(poly_data, "Piece")
    piece.set("NumberOfPoints", str(num_points))
    piece.set("NumberOfLines", str(num_lines))

    # Points section
    points_element = ET.SubElement(piece, "Points")
    points_data_array = ET.SubElement(points_element, "DataArray")
    points_data_array.set("type", "Float64")
    points_data_array.set("Name", "Points")
    points_data_array.set("NumberOfComponents", "3")
    points_data_array.set("format", "ascii")
    points_data_array.text = points_str

    # Lines section
    lines_element = ET.SubElement(piece, "Lines")
    connectivity_data_array = ET.SubElement(lines_element, "DataArray")
    connectivity_data_array.set("type", "Int64")
    connectivity_data_array.set("Name", "connectivity")
    connectivity_data_array.set("format", "ascii")
    connectivity_data_array.text = connectivity_str

    offsets_data_array = ET.SubElement(lines_element, "DataArray")
    offsets_data_array.set("type", "Int64")
    offsets_data_array.set("Name", "offsets")
    offsets_data_array.set("format", "ascii")
    offsets_data_array.text = offsets_str


    # Write to file
    tree = ET.ElementTree(vtk_file)
    try:
        ET.indent(tree, '  ')
        tree.write(filename, encoding="utf-8", xml_declaration=True)
        print(f"Successfully wrote {num_lines} line segments to {filename}")
    except IOError as e:
        print(f"Error writing file {filename}: {e}")

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
