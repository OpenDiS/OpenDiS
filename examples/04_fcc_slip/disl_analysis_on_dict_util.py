# Python utility functions for manipulating dislocation network structure
'''''''''''''''''''''''''''''''''''''''''''''''''''''
    Analyzing DDD simulation data on slip systems from oprec files
        (Nicolas Bertin | bertin1@llnl.gov)
        (Hamed Akhondzadeh | hamedakh@alumni.stanford.edu)
        (Wei Cai | caiwei@stanford.edu)
        (Hanfeng Zhai | hzhai@stanford.edu)
'''''''''''''''''''''''''''''''''''''''''''''''''''''
import matplotlib.pyplot as plt
import numpy as np
import tqdm
import copy, tabulate, re, os
from fcc_slip_sys import FCC_slip_sys

class TerminalColors:
    """A class to add colors to terminal output."""

    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    ITALIC = '\033[3m'
    UNDERLINE = '\033[4m'
    BLINK = '\033[5m'
    REVERSE = '\033[7m'
    HIDDEN = '\033[8m'

    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'

    BG_BLACK = '\033[40m'
    BG_RED = '\033[41m'
    BG_GREEN = '\033[42m'
    BG_YELLOW = '\033[43m'
    BG_BLUE = '\033[44m'
    BG_MAGENTA = '\033[45m'
    BG_CYAN = '\033[46m'
    BG_WHITE = '\033[47m'
    BG_BRIGHT_BLACK = '\033[100m'
    BG_BRIGHT_RED = '\033[101m'
    BG_BRIGHT_GREEN = '\033[102m'
    BG_BRIGHT_YELLOW = '\033[103m'
    BG_BRIGHT_BLUE = '\033[104m'
    BG_BRIGHT_MAGENTA = '\033[105m'
    BG_BRIGHT_CYAN = '\033[106m'
    BG_BRIGHT_WHITE = '\033[107m'

    @staticmethod
    def colored(text, *colors):
        """Applies the given colors to the text."""
        color_codes = "".join(colors)
        return f"{color_codes}{text}{TerminalColors.RESET}"


def get_edotdir_init(loadDirs):
    loadDirInterp = loadDirs[2]
    unitDir1 = loadDirs[0] / np.linalg.norm( loadDirs[0] )
    unitDir2 = loadDirs[1] / np.linalg.norm( loadDirs[1] )
    edotdir_init = (loadDirInterp/10)*unitDir1+(1-loadDirInterp/10)*unitDir2
    edotdir_init = edotdir_init / np.linalg.norm(edotdir_init)
    return edotdir_init


def get_ddd_database_folder_name(loadDirs):
    loadDirInterp = loadDirs[2]
    return 'changeOrientation_%d%d%d_%d%d%d_v1/%d/'%(loadDirs[0][0], loadDirs[0][1], loadDirs[0][2],
                                                     loadDirs[1][0], loadDirs[1][1], loadDirs[1][2],
                                                     loadDirInterp)


def pbc_position(r1, r2, box_size):
    if len(box_size) != 3:
        raise ValueError("box_size need to be an array of length 3")
    ds = (r2-r1)/box_size
    ds = ds - ds.round()
    return r1 + ds*box_size


def nodes_dicts_are_equiv(dict_1, dict_2, pos_tol=0.0, arm_tol=0.0, ignore_planes=False):
    are_equiv = True
    if len(dict_1.keys()) != len(dict_2.keys()):
        print(f'  {TerminalColors.RED}dictionaries have different number of keys{TerminalColors.RESET}')
        return False
    for tag in dict_1.keys():
        pos_1, arms_dict_1 = dict_1[tag][0], dict_1[tag][1]
        if tag not in dict_2:
            print(f'  {TerminalColors.RED}node {tag} in dict_1 not found in dict_2{TerminalColors.RESET}')
            return False
        pos_2, arms_dict_2 = dict_2[tag][0], dict_2[tag][1]
        if np.max(np.abs(pos_1 - pos_2)) > pos_tol:
            print(f'  {TerminalColors.BRIGHT_MAGENTA}node {tag} position mismatch:{TerminalColors.RESET}')
            print(f'    pos_1 {pos_1}')
            print(f'    pos_2 {pos_2}')
            print(f'    difference = {np.max(np.abs(pos_1 - pos_2))}')
            are_equiv = False
        if len(arms_dict_1.keys()) != len(arms_dict_2.keys()):
            print(f'  {TerminalColors.RED}node {tag} have different number arms{TerminalColors.RESET}')
            return False
        for nbr_tag in arms_dict_1:
            if not ignore_planes:
                arm_data_1 = arms_dict_1[nbr_tag]
                arm_data_2 = arms_dict_2[nbr_tag]
            else:
                arm_data_1 = arms_dict_1[nbr_tag][0:3]
                arm_data_2 = arms_dict_2[nbr_tag][0:3]
            if np.max(np.abs((arm_data_1 - arm_data_2))) > arm_tol:
                print(f'  {TerminalColors.BRIGHT_MAGENTA}arm {tag} - {nbr_tag} data mismatch:{TerminalColors.RESET}')
                print(f'    arm_data_1 {arm_data_1}')
                print(f'    arm_data_2 {arm_data_2}')
                print(f'    difference = {np.max(np.abs((arm_data_1 - arm_data_2)))}')
                are_equiv = False
    return are_equiv


def plot_nodes_dict(nodes, bounds, plot_links=True, trim=False, fig=None, ax=None, block=False, pause_seconds=0.01):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    if fig==None:
        fig = plt.figure(figsize=(8,8))
    if ax==None:
        ax = plt.axes(projection='3d')

    box_size = np.array(bounds[1] - bounds[0])

    rn = np.array([nodes[i][0].tolist() for i in list(nodes.keys())])
    p_link = np.empty((0,6))

    plt.cla()
    if plot_links:
        for my_tag in list(nodes.keys()):
            my_coords = nodes[my_tag][0][0:3]
            arms_dict = nodes[my_tag][1]
            for nbr_tag in list(arms_dict.keys()):
                if my_tag < nbr_tag:
                    r_link = np.zeros((2,3))
                    if not nbr_tag in nodes.keys():
                        print("plot_nodes_dict(): cannot find neighbor")
                        print([my_tag, nbr_tag])
                        continue
                    nbr_coords = nodes[nbr_tag][0][0:3]
                    r_link[0,:] = my_coords
                    r_link[1,:] = pbc_position(my_coords, nbr_coords, box_size)
                    if (not trim) or np.max(np.absolute(r_link)) <= box_size.max()/2:
                        p_link = np.append(p_link, [r_link[0,:], r_link[1,:]])

    ls = p_link.reshape((-1,2,3))
    lc = Line3DCollection(ls, linewidths=0.5, colors='b')
    ax.add_collection(lc)

    ax.scatter(rn[:,0], rn[:,1], rn[:,2], c='r', s=4)
    ax.set_xlim(bounds[0][0], bounds[1][0])
    ax.set_ylim(bounds[0][1], bounds[1][1])
    ax.set_zlim(bounds[0][2], bounds[1][2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_box_aspect([1,1,1])

    plt.draw()
    plt.show(block=block)
    plt.pause(pause_seconds)

    return fig, ax


def spin_matrix(phi):
    # Ss, Cs = phi, 1 - 0.5*phi**2 # small angle approximation
    Ss, Cs = np.sin(phi), np.cos(phi)

    Rspin = np.zeros([3,3])

    Rspin[0,0] =  Cs[2]*Cs[1]
    Rspin[1,1] =  Cs[2]*Cs[0] + Ss[2]*Ss[1]*Ss[0];
    Rspin[2,2] =  Cs[1]*Cs[0];
    Rspin[0,1] = -Ss[2]*Cs[0] + Cs[2]*Ss[0]*Ss[1];
    Rspin[1,2] = -Cs[2]*Ss[0] + Ss[2]*Ss[1]*Cs[0];
    Rspin[2,0] = -Ss[1];
    Rspin[0,2] =  Ss[2]*Ss[0] + Cs[2]*Cs[0]*Ss[1];
    Rspin[1,0] =  Ss[2]*Cs[1];
    Rspin[2,1] =  Cs[1]*Ss[0];

    return Rspin


def read_nodes_to_dict(datafile, print_func=print):
# Read dislocation nodes from ParaDiS data file into a python dictionary
# Works for restart files generated from parallel runs
#
# Syntax
#  nodes, bounds = read_nodes_to_dict(datafile)
#
# Input
#  datafile: path to the ParaDiS data file
#
# Outputs
#  nodes: dictionary indexed by tag
#         each value of the dictionary is a list [ my_coord, arms_dict ]
#         where my_coord is numpy array, arms_dict is a dictionary indexed by tag whose values are numpy arrays
#  bounds: [[xMin, yMin, zMin],[xMax, yMax, zMax]]
    print_func("read_nodes_to_dict = %s"%(datafile))
    try:
        with open(datafile,'r') as fid:
            line = fid.readline(); dataFileVersion = int(line.split()[2])
            print_func("  dataFileVersion = %d"%(dataFileVersion))
            line = fid.readline(); numFileSegments = int(line.split()[2])
            print_func("  numFileSegments = %d"%(numFileSegments))

            line = fid.readline()
            line = fid.readline(); xMin = float(line.split()[0])
            line = fid.readline(); yMin = float(line.split()[0])
            line = fid.readline(); zMin = float(line.split()[0])
            line = fid.readline()
            line = fid.readline()
            line = fid.readline(); xMax = float(line.split()[0])
            line = fid.readline(); yMax = float(line.split()[0])
            line = fid.readline(); zMax = float(line.split()[0])
            line = fid.readline()

            bounds = np.array([[xMin,yMin,zMin], [xMax,yMax,zMax]])

            line = fid.readline(); nodeCount = int(line.split()[2])
            print_func("  nodeCount = %d"%(nodeCount))
            line = fid.readline(); dataDecompType = int(line.split()[2])
            print_func("  dataDecompType = %d"%(dataDecompType))

            # initialize an empty dictionary
            nodes = {}
            while True:
                line = fid.readline()
                if line.find('nodalData =') != -1:
                    line = fid.readline(); skiptwo = line.find('#') != -1
                    break
                if not line:
                    print_func("Cannot find nodalData in" + datafile)
                    return
            print_func("  skiptwo = %r"%(skiptwo))

            if skiptwo:
                line = fid.readline()
            while True:
                if not skiptwo and len(nodes) == 0:
                    data_list = line.split()
                else:
                    line = fid.readline(); data_list = line.split()
                if len(data_list) == 0:
                    break
                if len(data_list) != 6:
                    print_func('Unknown data file format')
                    return nodes, bounds

                my_tag = tuple([ int(i) for i in data_list[0].split(',') ])
                my_coord = np.array([float(data_list[1]), float(data_list[2]), float(data_list[3]), float(data_list[5])])

                num_arms = int(data_list[4])
                arms_dict = {}
                for i in range(num_arms):
                    line = fid.readline(); arm_data_list = line.split()
                    if len(arm_data_list) != 4:
                        print_func('Unknown data file format')
                        return nodes, bounds
                    nbr_tag = tuple([ int(i) for i in arm_data_list[0].split(',')])
                    line = fid.readline(); plane_data_list = line.split()
                    if len(plane_data_list) != 3:
                        print_func('Unknown data file format')
                        return nodes, bounds
                    arm_data = np.array([ float(arm_data_list[1]),float(arm_data_list[2]),float(arm_data_list[3]),
                                          float(plane_data_list[0]),float(plane_data_list[1]),float(plane_data_list[2])])
                    arms_dict[nbr_tag] = arm_data
                nodes[my_tag] = [my_coord, arms_dict]
            fid.close()
    except:
        print_func('Cannot open data file ' + datafile)
        raise

    return nodes, bounds


def clear_arms_extra_dict(nodes):
    for my_tag in list(nodes.keys()):
        node_attr = nodes[my_tag]
        arms_dict = node_attr[1]
        arms_extra_dict = {}
        for nbr_tag in list(arms_dict.keys()):
            arms_extra_dict[nbr_tag] = {}
        if len(node_attr) > 2:
            node_attr[2] = arms_extra_dict
        else:
            node_attr.append(arms_extra_dict)
    return nodes


def make_arms_extra_dict_key_consistent(nodes):
    # make nodes have same keys for arms_dict and arms_extra_dict
    for my_tag in list(nodes.keys()):
        node_attr = nodes[my_tag]
        arms_dict = node_attr[1]
        if len(node_attr) > 2:
            arms_extra_dict = node_attr[2]
            for nbr_tag in list(arms_extra_dict.keys()):
                if not nbr_tag in arms_dict:
                    arms_extra_dict.pop(nbr_tag)
            for nbr_tag in list(arms_dict.keys()):
                if not nbr_tag in arms_extra_dict:
                    arms_extra_dict[nbr_tag] = {}
    return nodes


def check_arms_extra_dict_keys(nodes):
    # verify that nodes has same keys for arms_dict and arms_extra_dict
    for my_tag in list(nodes.keys()):
        node_attr = nodes[my_tag]
        if set(node_attr[1].keys()) != set(node_attr[2].keys()):
            print("find_arm_sysID_for_nodes_dict:")
            print("my_tag: ", my_tag)
            print("node_attr[my_tag][1].keys: ", set(node_attr[1].keys()))
            print("node_attr[my_tag][2].keys: ", set(node_attr[2].keys()))
            raise ValueError("nodes have different keys for arms_dict and arms_extra_dict !")


def find_armsysID_for_nodes_dict(nodes, slip_sys):
    # create list for keys since we are changing the dictionary
    for my_tag in list(nodes.keys()):
        node_attr = nodes[my_tag]
        arms_dict = node_attr[1]
        if len(node_attr) > 2:
            arms_extra_dict = node_attr[2]
        else:
            arms_extra_dict = {nbr_tag: {} for nbr_tag in arms_dict}
        # create list for keys since we are changing the dictionary
        for nbr_tag in list(arms_dict.keys()):
            arm_attr = arms_dict[nbr_tag]
            burg = arm_attr[0:3]
            nPlane = arm_attr[3:6]
            armsysID = int(slip_sys.nb_to_sysID(nPlane, burg)[0])
            if nbr_tag in arms_extra_dict:
                arms_extra_dict[nbr_tag]['sysID'] = armsysID
            else:
                arms_extra_dict[nbr_tag] = {'sysID': armsysID}
        # in case the memory location of arms_extra_dict has changed
        if len(node_attr) > 2:
            nodes[my_tag][2] = arms_extra_dict
        else:
            nodes[my_tag].append(arms_extra_dict)

    check_arms_extra_dict_keys(nodes)
    return nodes


def find_segment_length_for_nodes_dict(nodes, bounds, dt=None): #***
    box_size = np.array(bounds[1] - bounds[0])
    for my_tag, node_attr in nodes.items():
        if len(node_attr) > 2:
            arms_extra_dict = node_attr[2]
        else:
            arms_extra_dict = {nbr_tag: {} for nbr_tag in arms_dict}
        
        my_coords = node_attr[0][0:3]
        arms_dict = node_attr[1]
        # arms_extra_dict = node_attr[2]
        
        for nbr_tag, arm_attr in arms_dict.items():
            nbr_coords = nodes[nbr_tag][0][0:3]
            nbr_coords = pbc_position(my_coords, nbr_coords, box_size)
            length = np.linalg.norm(nbr_coords - my_coords)
            if nbr_tag in arms_extra_dict:
                arms_extra_dict[nbr_tag]['length'] = length
                if dt!=None: arms_extra_dict[nbr_tag]['l_dot'] = length / dt
            else:
                arms_extra_dict[nbr_tag] = {'length': length}
                if dt!=None: arms_extra_dict[nbr_tag] = {'l_dot': length / dt}

        if len(node_attr) > 2:
            nodes[my_tag][2] = arms_extra_dict
        else:
            nodes[my_tag].append(arms_extra_dict)

    check_arms_extra_dict_keys(nodes)
    return nodes


def find_swept_area_for_nodes_dict(stats, nodes_old, nodes, bounds):
    box_size = np.array(bounds[1] - bounds[0])
    for my_tag, node_attr in nodes.items():
        if len(node_attr) > 2:
            arms_extra_dict = node_attr[2]
        else:
            arms_extra_dict = {nbr_tag: {} for nbr_tag in arms_dict}
        my_coords = node_attr[0][0:3]
        arms_dict = node_attr[1]

        for nbr_tag, arm_attr in arms_dict.items():
            # nbr_coords = nodes[nbr_tag][0][0:3]
            # nbr_coords = pbc_position(my_coords, nbr_coords, box_size)
            # V = L**3
            link_info = arms_dict[nbr_tag]
            r1 = nodes[my_tag][0][0:3]
            r2 = pbc_position(r1, nodes[nbr_tag][0][0:3], box_size)
            r3 = pbc_position(r1, nodes_old[my_tag][0][0:3], box_size)
            r4 = pbc_position(r3, nodes_old[nbr_tag][0][0:3], box_size)
            Sn = 0.5*np.cross(r2-r3, r1-r4) # Sn_i
            swept_area = np.linalg.norm(Sn)

            burg_vec = arms_dict[nbr_tag][0:3]
            # sym_nb_prod = (1/2) * (np.outer(plane_normal, burg_vec) - np.outer(plane_normal, burg_vec).T)
            Snb = np.outer(Sn, burg_vec)
            d_pl_strain = (1/2) * (Snb + Snb.T) / (box_size[0]*box_size[1]*box_size[2])

            # plane_normal = link_info[:3]; burg_vec = link_info[3:7];
            # P = np.multiply(n[ [0,0,0,1,1,1,2,2,2] ], burg_vec[ [0,1,2,0,1,2,0,1,2] ] )
            # S_nb = P[ [1,2,5]] - P[ [3,6,7]]; #print('snb',S_nb,P)
            # swept_area = 0.5 * ( P[ [0,4,8,1,2,5]]+P[ [0,4,8,3,6,7]] )#/V
            if nbr_tag in arms_extra_dict:
                arms_extra_dict[nbr_tag]['d_pl_strain'] = d_pl_strain
                arms_extra_dict[nbr_tag]['swept_area'] = swept_area
            else:
                arms_extra_dict[nbr_tag] = {'d_pl_strain': d_pl_strain}
                arms_extra_dict[nbr_tag] = {'swept_area': swept_area}
        if len(node_attr) > 2:
            nodes[my_tag][2] = arms_extra_dict
        else:
            nodes[my_tag].append(arms_extra_dict)
    check_arms_extra_dict_keys(nodes)
    return nodes


def nodes_dict_to_graph(nodes, DiGraph=False):
    # convert nodes dictionary to networkx DiGraph
    import networkx as ntx
    if DiGraph: G = ntx.DiGraph()
    else: G = ntx.Graph()
    G.add_nodes_from(nodes)
    for node, node_attr in nodes.items():
        G.nodes[node]['R'] = node_attr[0][0:3]
        G.nodes[node]['flag']     = node_attr[0][3]
        arms_data = node_attr[1]
        for nbr_node, arm_attr in arms_data.items():
            G.add_edge( node, nbr_node )
            G.edges[node, nbr_node]['burg_vec'] = arm_attr[0:3]
            G.edges[node, nbr_node]['plane_normal'] = arm_attr[3:6]
    return G


def phys_nodes_dict_to_graph(nodes, DiGraph=False):
    # convert physical network in nodes dictionary to networkx DiGraph
    import networkx as ntx
    G = ntx.DiGraph() if DiGraph else ntx.Graph()

    for node, node_attr in nodes.items():
        if isinstance(node_attr[0], np.ndarray) and node_attr[0].size > 0:  
            G.add_node(node)  # Only add the node if it meets the condition
            G.nodes[node]['R'] = node_attr[0]  # Assign the position data to the node

            arms_data = node_attr[1]  # Get arms data
            arms_extra_data = node_attr[2]  # Get additional arms data
            
            for nbr_node, arm_attr in arms_data.items():
                G.add_edge(node, nbr_node)  # Add edge to the graph
                G.edges[node, nbr_node]['burg'] = arm_attr[0:3]  # Assign burg attribute
                G.edges[node, nbr_node]['plane_normal'] = arm_attr[3:6]  # Assign plane normal attribute
            
            for nbr_node, arm_extra_attr in arms_extra_data.items():
                if G.has_edge(node, nbr_node):  # Ensure the edge exists before trying to add attributes
                    G.edges[node, nbr_node]['slip_sys'] = arm_extra_attr['sysID']  # Assign slip_sys attribute
                    G.edges[node, nbr_node]['phys_length'] = arm_extra_attr['physical_length']  # Assign phys_length attribute
    return G


def sanity_check_nodes_dict(nodes, btol=1e-6):
# Perform sanity check on nodes dictionary
# Burgers vector on every arm stored on two end nodes should sum up to zero
# Burgers vector of all arms on the same arm should sum up to zero
#
# Syntax
#  is_sane = sanity_check_node_dict(nodes)
#
# Inputs
#  nodes: dictionary of dislocation network structure organized by nodes
#
# Output
#  is_sane: True or False

    is_sane = True
    for my_tag, my_node in nodes.items():
        arms_dict = my_node[1]

        bsum = np.zeros((1,3))
        for nbr_tag, arms_attr in arms_dict.items():
            b1 = arms_attr[0:3]
            b2 = nodes[nbr_tag][1][my_tag][0:3]
            if np.max(np.abs(b1+b2)) > btol:
                print("b1 + b2 != 0"); print(b1); print(b2); print(b1+b2)
                is_sane = False
                return is_sane
            bsum += b1
        if np.max(np.abs(bsum)) > btol:
            print("bsum = 0"); print(bsum);
            is_sane = False
            return is_sane

    return is_sane


def sanity_check_phys_nodes_dict(phys_nodes):
    # Sanity check: all arms in phys_nodes points to nodes in phys_nodes

    is_sane = True
    for my_tag, node_attr in phys_nodes.items():
        arms_dict = node_attr[2]
        arms_extra_dict = node_attr[2]
        #for nbr_tag, arm_extra_attr in arms_extra_dict.items():
        for nbr_tag in arms_dict.keys():
            arm_extra_attr = arms_extra_dict[nbr_tag]
            if "dest_tag" not in arm_extra_attr:
                is_sane = False
                print("*** sanity_check_phys_nodes_dict ***")
                print("dest_tag not a key in arm_extra_attr")
                print("my_tag:", my_tag)
                print("nbr_tag:", nbr_tag)
                print("arm_extra_attr:", arm_extra_attr)
                print(f"arms_extra_dict[{nbr_tag}]:", arms_extra_dict[nbr_tag])
                print("arms_extra_dict:", arms_extra_dict)
                print("arms_dict:", node_attr[1])
                print("************************************")
                return is_sane
            dest_tag = arm_extra_attr["dest_tag"]
            if not dest_tag in phys_nodes.keys():
                print([my_tag, nbr_tag, dest_tag])
                print(node_attr)
                print("dest_tag not a key in phys_nodes")
                is_sane = False
                return is_sane
    return is_sane


def is_physical_node(nodes, ptol=1e-5):
# Determine if each node in the nodes dictionary are physical node
# A node is a physical node if it has 1 or more than 2 arms,
# or (for a 2-arm node) if its arms are on different slip planes.
# Otherwise, it is a discretization node
#
# Syntax
#  is_phys = is_physical_node(nodes, ptol)
#
# Inputs
#  nodes: dictionary of dislocation network structure organized by nodes
#  ptol:  tolerance of whether two planes are identical
#
# Output
#  is_phys: dictionary of tag -> True (if a node is a physical node) or False (otherwise)

    is_phys = dict( (my_tag, False) for my_tag in list(nodes.keys()) )
    for my_tag, node_attr in nodes.items():
        arms_dict = node_attr[1]
        num_nbrs = len(arms_dict.keys())
        if num_nbrs != 2:
            is_phys[my_tag] = True
        else:
            planes = np.zeros((2,3))
            id = 0
            for nbr_tag in list(arms_dict.keys()):
                planes[id,:] = arms_dict[nbr_tag][3:6]
                if np.linalg.norm(planes[id,:]) > 0:
                    planes[id,:] = planes[id,:] / np.linalg.norm(planes[id,:])
                else:
                    # TODO: check whether it is ok to have zero plane normal vector
                    # print("is_physical_node: warning plane normal is zero")
                    pass
                id += 1
            planes[1,:] = planes[1,:] - np.dot(planes[0,:], planes[1,:])*planes[0,:]
            if np.linalg.norm(planes[1,:]) > ptol:
                is_phys[my_tag] = True
    return is_phys


def extract_link_length(phys_nodes):
# Extract the direct and physical lengths stored in phys_nodes dictionary
# Syntax
#  link_lengths = physical_structure(nodes, bounds, glide)
#
# Inputs
#  phys_nodes: dictionary of dislocation network structure organized by physical nodes
#
# Outputs
#  link_lengths: nx2 numpy array
#    link_lengths[:,0] -- direct lengths
#    link_lengths[:,1] -- physical lengths (contour lengths of curved paths)

    link_lengths = np.empty((0,2))
    for my_tag, node_attr in phys_nodes.items():
        arms_dict = node_attr[1]
        arms_extra_dict = node_attr[2]
        for nbr_tag, arm_extra_attr in arms_extra_dict.items():
            nbr_link_length = np.array([[arm_extra_attr['direct_length'], arm_extra_attr['physical_length']]])
            dest_tag = arm_extra_attr['dest_tag']
            if my_tag < dest_tag:
                link_lengths = np.append(link_lengths, nbr_link_length, axis=0)
    return link_lengths


def extract_link_info(phys_nodes):
# Extract the link information stored in phys_nodes dictionary
# Syntax
#  phys_nodes = physical_structure(nodes, bounds, glide)
#  link_info = extract_link_info(phys_nodes)
#
# Inputs
#  phys_nodes: dictionary of dislocation network structure organized by physical nodes
#
# Outputs
#  link_info: nx3 numpy array
#    link_info[:,0] -- slip system IDs of the links
#    link_info[:,1] -- direct lengths of the links
#    link_info[:,2] -- physical lengths of the links (contour lengths of curved paths)
    link_info = np.empty((0,3))
    for my_tag, node_attr in phys_nodes.items():
        arms_dict = node_attr[1]
        arms_extra_dict = node_attr[2]
        for nbr_tag, arm_attr in arms_dict.items():
            if len(arm_attr) != 6:
                raise ValueError("extract_link_info: arm_attr does not have length 6")
        for nbr_tag, arm_extra_attr in arms_extra_dict.items():
            dest_tag = arm_extra_attr["dest_tag"]
            if my_tag < dest_tag:
                sysID = arm_extra_attr["sysID"]
                direct_length = arm_extra_attr["direct_length"]
                physical_length = arm_extra_attr["physical_length"]
                link_info = np.append(link_info, np.array([[sysID, direct_length, physical_length]]), axis=0)
    return link_info


def physical_structure(nodes, bounds, glide=True, calc_pl_strain=False):
# Reduce the DD structure to physical nodes only by removing all discretization nodes
# Note: physical_structure allows two links to go from the same (source) node to the same destination node
#       the key in the arms_dict is the original nbr_tag of the (source) node
#       the destination node tag is stored as a value 'dest_tag' in the arms_extra_dict
#
# Syntax
#  phys_nodes = physical_structure(nodes, bounds, glide)
#
# Inputs
#  nodes: dictionary of dislocation network structure organized by nodes
#  bounds: bounds of the primary network
#  glide: (optional) if glide = True (default), the glide plane constraints
#         are used and preserved
#
# Outputs
#  phys_nodes: new dictionary of dislocation network containing only physical nodes
#              the keys in the arms_dict for each node of the original neighbors
#              the detination node ('dest_tag') of each link is stored in arms_extra_dict
#              arms_extra_dict also contains the 'direct_length' and 'physical_length' of each link
#

    # determine if each node is a physical node
    is_phys = is_physical_node(nodes)
    # initialize phys_nodes from the sub-set of nodes dictionary
    phys_nodes = {}

    '''old code but assign physical length based on the segment length'''
    for my_tag, node_attr in nodes.items():
        if is_phys[my_tag]:
            node_attr = nodes[my_tag]
            my_coords = node_attr[0][0:3]
            arms_dict = { nbr_tag: node_attr[1][nbr_tag] for nbr_tag in node_attr[1] }
            arms_extra_dict = { nbr_tag: node_attr[2][nbr_tag] for nbr_tag in node_attr[2] } if len(node_attr) > 2 else {}
            # initialize memory for arms_extra_dict[nbr_tag]
            for nbr_tag in list(arms_dict.keys()):
                arms_extra_dict[nbr_tag].update({'dest_tag': None, 'direct_length': None, 'physical_length': None,
                                                 'link_d_pl_strain': None, 'num_segments': None, 'plane_intersect': None})

            phys_nodes[my_tag] = [my_coords.copy(), arms_dict, arms_extra_dict]

    # verify that nodes has same keys for arms_dict and arms_extra_dict
    check_arms_extra_dict_keys(nodes)
    box_size = np.array(bounds[1] - bounds[0])

    # make sure node_attr is a list of length 3
    for my_tag in list(phys_nodes.keys()):
        node_attr = phys_nodes[my_tag]
        if len(node_attr)==2:
            node_attr.append({})

    # create list for keys since we are changing the dictionary
    for my_tag in list(phys_nodes.keys()):
        node_attr = phys_nodes[my_tag]
        my_coords = node_attr[0][0:3]
        arms_dict = node_attr[1]
        arms_extra_dict = node_attr[2]
        for nbr_tag in list(arms_dict.keys()):
            if is_phys[nbr_tag]:
                # the link consists of one segment from my_tag to nbr_tag
                nbr_coords = phys_nodes[nbr_tag][0][0:3]
                nbr_coords = pbc_position(my_coords, nbr_coords, box_size)
                direct_length = np.linalg.norm(nbr_coords - my_coords)
                new_nbr_tag = nbr_tag
                physical_length = arms_extra_dict[nbr_tag]['length']
                if calc_pl_strain: link_d_pl_strain = arms_extra_dict[nbr_tag]['d_pl_strain']
            else:
                # travel along the dislocation line until finding a physical node
                zeroth_nbr_tag = my_tag
                first_nbr_tag = nbr_tag
                physical_length = arms_extra_dict[nbr_tag]['length']
                ref_coords = nodes[zeroth_nbr_tag][0][0:3]
                first_nbr_coords = nodes[first_nbr_tag][0][0:3]
                first_nbr_coords = pbc_position(ref_coords, first_nbr_coords, box_size)
                if calc_pl_strain: link_d_pl_strain = arms_extra_dict[nbr_tag]['d_pl_strain']
                while not is_phys[first_nbr_tag]:
                    second_nbr_tag_list = list( nodes[first_nbr_tag][1].keys() )
                    second_nbr_tag_list.remove(zeroth_nbr_tag)
                    if len(second_nbr_tag_list) != 1:
                        raise ValueError("There should be only one nbr left here!")
                    second_nbr_tag = second_nbr_tag_list[0]
                    zeroth_nbr_tag = first_nbr_tag
                    first_nbr_tag = second_nbr_tag
                    ref_coords = first_nbr_coords.copy()
                    first_nbr_coords = nodes[first_nbr_tag][0][0:3]
                    first_nbr_coords = pbc_position(ref_coords, first_nbr_coords, box_size)
                    physical_length += nodes[zeroth_nbr_tag][2][first_nbr_tag]['length']
                    if calc_pl_strain: link_d_pl_strain += nodes[zeroth_nbr_tag][2][first_nbr_tag]['d_pl_strain']

                # the physical node at the other end is called new_nbr_tag
                new_nbr_tag = first_nbr_tag
                if nbr_tag == new_nbr_tag:
                    print([my_tag, nbr_tag, new_nbr_tag])
                    raise ValueError("physical_structure(): This is not supposed to happen!")

                if not is_phys[new_nbr_tag]:
                    raise ValueError(f"physical_structure(): new_nbr_tag {new_nbr_tag} is not physical!")

                # compute the direct length between node and new_nbr_tag (handles PBC correctly when going along the link)
                direct_length = np.linalg.norm(first_nbr_coords - my_coords)

            if 'dest_tag' in arms_extra_dict[nbr_tag]: arms_extra_dict[nbr_tag]['dest_tag'] = new_nbr_tag
            if 'direct_length' in arms_extra_dict[nbr_tag]: arms_extra_dict[nbr_tag]['direct_length'] = direct_length
            if 'physical_length' in arms_extra_dict[nbr_tag]: arms_extra_dict[nbr_tag]['physical_length'] = physical_length
            if calc_pl_strain:
                if 'link_d_pl_strain' in arms_extra_dict[nbr_tag]: arms_extra_dict[nbr_tag]['link_d_pl_strain'] = link_d_pl_strain

        # in case the memory location of arms_extra_dict has changed
        node_attr[2] = arms_extra_dict

    if not sanity_check_phys_nodes_dict(phys_nodes):
        raise ValueError("physical_structure(): phys_nodes does not pass sanity check")

    return phys_nodes


class oprec_player():
    # class for playing oprec file on nodes dictionary
    def __init__(self, target_step_list=None, target_time_list=None,
                 edotdir_init=None, edotdir=None,
                 separate_step_at_optype=[10,11],
                 voigt_test=False, voigt_test_threshold=1e-16,
                 update_stress_strain=False, mu=None, nu=None, eRate=None,
                 tqdm_target_time_multiplier=1.1,
                 print_freq=1, print_func=print,
                 fig=None, ax=None, plot_freq=100, plot_physical_structure=False,
                 plot_links=True, plot_trim=False,
                 print_analysis_details=False,
                 slip_sys_obj=None, slip_analysis=False,
                 save_nodes_to_file=False, save_data_path="."):
        self.target_step_list = target_step_list
        self.target_step_list_counter = 0
        self.target_time_list = target_time_list
        self.target_time_list_counter = 0
        self.edotdir_init = edotdir_init
        self.edotdir = edotdir
        self.separate_step_at_optype = separate_step_at_optype
        self.voigt_test = voigt_test
        self.voigt_test_threshold = voigt_test_threshold
        self.update_stress_strain = update_stress_strain
        # stress and strain (testing)
        self.d_stress_sum, self.d_strain_sum = np.zeros((3,3)), np.zeros((3,3))

        self.tqdm_target_time_multiplier = tqdm_target_time_multiplier
        self.print_freq = print_freq
        self.print_func = print_func
        self.fig = fig
        self.ax = ax
        self.plot_freq = plot_freq
        self.plot_physical_structure = plot_physical_structure
        self.plot_links = plot_links
        self.plot_trim = plot_trim
        self.print_analysis_details = print_analysis_details
        self.slip_sys_obj = slip_sys_obj
        self.slip_analysis = slip_analysis
        self.save_nodes_to_file = save_nodes_to_file
        self.save_data_path = save_data_path
        if self.save_data_path is not None:
            os.makedirs(self.save_data_path, exist_ok=True)
        self.played_time_list = []

        self.nodes = None
        self.nodes_old = None
        self.bounds = None
        self.arms_dict = None

        self.end_of_file = False
        self.prev_optype = None
        self.prev_step = None
        self.dt = None
        self.time_now = 0
        self.step = 0
        self.io = 0
        self.d_pl_strain_sum, self.d_pl_rotation_sum = np.zeros((3,3)), np.zeros((3,3))
        self.stats = None
        self.colors = TerminalColors()

    def clear_step_io(self):
        self.step = 0
        self.io = 0

    def parse_ctrl_file(self, filename):
        with open(filename) as f:
            ctrl_content = f.read()
        # This regex finds lines with:
        #   key = value
        # where value is either:
        #   - a list enclosed in [ ... ] (non-greedy, DOTALL to allow newlines)
        #   - or a quoted string, or any non-empty value up to the newline.
        #
        # It ignores commented lines (those starting with '#' after optional whitespace).
        pattern = re.compile(
            r'^\s*([^#\s][\w_]*)\s*=\s*'
            r'(\[(?:.|\n)*?\]|".*?"|\S.*)$', 
            re.MULTILINE)

        def convert(s):
            s = s.strip()
            # Remove surrounding quotes if present.
            if s.startswith('"') and s.endswith('"'):
                return s[1:-1]
            try:
                # Convert to int if possible; if the number contains a dot or an exponent, use float.
                if '.' in s or 'e' in s.lower():
                    return float(s)
                else:
                    return int(s)
            except Exception:
                return s
        ctrl_info = {}
        # Iterate over each key-value match.
        for m in pattern.finditer(ctrl_content):
            key = m.group(1)
            value_str = m.group(2).strip()
            # If the value is a bracketed list, remove the brackets and split by any whitespace.
            if value_str.startswith('[') and value_str.endswith(']'):
                inner = value_str[1:-1].strip()
                # Split on any whitespace, including newlines.
                items = re.split(r'\s+', inner)
                ctrl_info[key] = [convert(x) for x in items if x]
            else:
                ctrl_info[key] = convert(value_str)
        self.ctrl_info = ctrl_info
        self.assign_params_from_ctrl_info()
        return ctrl_info

    def assign_params_from_ctrl_info(self):
        self.mu = self.ctrl_info['shearModulus']
        self.nu = self.ctrl_info['pois']
        self.eRate = self.ctrl_info['eRate']
        self.init_step = self.ctrl_info.get('cycleStart', 0)
        self.step = self.init_step if self.step == 0 else self.step
        self.init_time = self.ctrl_info.get('timeNow', 0.0)
        self.time_now = self.init_time

    def set_nodes_data(self, nodes, bounds):
        self.nodes = nodes
        self.bounds = bounds

    def play(self, oprecfile):
    # Input: oprecfile - path of the oprecfile
        print("open oprecfile (" + oprecfile + ")")
        with open(oprecfile) as oprecfile_handle:
            self.play_open_file(oprecfile_handle)
        return

    def play_open_file(self, f):
    # Inputs: f - opened file handler of oprecfile
        self.init_vars_before_play()

        if self.target_step_list is None and self.target_time_list is None:
            raise ValueError("target_step_list and target_time_list are both None")
        if self.target_step_list is not None and self.target_time_list is not None:
            raise ValueError("target_step_list and target_time_list are both not None")

        if self.target_step_list is not None: tqdm_target = self.target_step_list[-1] - self.init_step
        elif self.target_time_list is not None: tqdm_target = (self.target_time_list[-1] - self.init_time) * self.tqdm_target_time_multiplier
        else: tqdm_target = None

        with tqdm.tqdm(total=tqdm_target, desc="Oprec") as pbar:
            self.print_func = tqdm.tqdm.write
            while True:
                self.print_func("status = %d total = %d"%(self.step, tqdm_target))
                self.play_one_block(f)
                if self.need_data_process_after_play_one_block(): self.data_process()
                if self.target_step_list is not None: pbar.n=self.step - self.init_step
                elif self.target_time_list is not None: pbar.n=self.time_now - self.init_time
                if self.update_stress_strain:
                    if self.optype == 11: self.update_applied_stress()
                pbar.refresh()
                if self.need_break_after_play_one_block(): break
        self.print_func = print
        return

    def init_vars_before_play(self):
        if self.fig==None:
            self.fig = plt.figure(figsize=(8,8))
        if self.ax==None:
            self.ax = plt.axes(projection='3d')

        if self.edotdir is None:
            self.edotdir = self.edotdir_init.copy() # self.edotdir_init.copy()
            if self.edotdir is None:
                raise ValueError("edotdir is not set")

        if self.slip_analysis:
            if self.slip_sys_obj==None:
                self.slip_sys_obj = FCC_slip_sys() # set it to be FCC for analysis as default
                self.print_func("set slip_sys_obj = FCC_slip_sys()")
            self.nodes = clear_arms_extra_dict(self.nodes)
            self.dstress_array, self.dstrain_array = [], []

        if self.prev_optype == 10:
            self.nodes_old = copy.deepcopy(self.nodes)

    def play_one_block(self, f):
        # play for one block (multiple operations) of the oprec file
        # note that a block may be a fraction of a step, i.e. a step can have multiple blocks
        while True:
            self.play_one_op(f)
            if self.need_break_after_play_one_op(): break
        return

    def need_data_process_after_play_one_block(self):
        if self.prev_optype not in self.separate_step_at_optype:
            raise ValueError("self.prev_optype = %d not in %s"%(self.prev_optype, str(self.separate_step_at_optype)))

        need_data_process = False
        if self.slip_analysis:
            need_data_process = True

        if self.target_time_list is not None:
            if self.time_now >= self.target_time_list[self.target_time_list_counter]:
                self.target_time_list_counter += 1
                need_data_process = True
        elif self.target_step_list is not None:
            if self.step >= self.target_step_list[self.target_step_list_counter]:
                self.target_step_list_counter += 1
                need_data_process = True
        return need_data_process

    def data_process(self):
        self.played_time_list.append(self.time_now)

        if self.slip_analysis:
            self.data_process_slip_analysis()

        if self.save_nodes_to_file:
            self.data_process_save_nodes_to_file()

    def data_process_slip_analysis(self):
        if hasattr(self, 'dstress'):
            if self.d_normal_stress is not None: self.dstress_array.append(self.d_normal_stress)
        if hasattr(self, 'dstrain'):
            if self.d_normal_strain is not None: self.dstrain_array.append(self.d_normal_strain)

        # should not be needed, already in play_one_op
        #if self.prev_optype == 10:
        #    self.nodes_old = copy.deepcopy(self.nodes)

        if self.prev_optype == 11:
            self.nodes_old = clear_arms_extra_dict(self.nodes_old)
            self.nodes_old = find_armsysID_for_nodes_dict(self.nodes_old, self.slip_sys_obj)
            self.nodes_old = find_segment_length_for_nodes_dict(self.nodes_old, self.bounds, self.dt)

            self.nodes = clear_arms_extra_dict(self.nodes)
            self.nodes = find_armsysID_for_nodes_dict(self.nodes, self.slip_sys_obj)
            self.nodes = find_segment_length_for_nodes_dict(self.nodes, self.bounds, self.dt)
            self.nodes = find_swept_area_for_nodes_dict(self.stats, self.nodes_old, self.nodes, self.bounds)

            # if not sanity_check_nodes_dict(nodes):
            #     raise ValueError("nodes not satisfying sanity check!")

            self.phys_nodes = physical_structure(self.nodes, self.bounds, calc_pl_strain=True)

            self.stats.inc_step = self.prev_step - self.init_step
            self.stats.add_link_length_per_slipsys_hist(self.phys_nodes, print_analysis_details=self.print_analysis_details)
            self.stats.add_dGamma_per_slipsys_hist(self.phys_nodes, time_step=self.dt, slip_sys=self.slip_sys_obj, 
                                                   print_analysis_details=self.print_analysis_details, print_func=self.print_func,
                                                   save_intermediate_graph=self.save_nodes_to_file, save_data_path=self.save_data_path,
                                                   ddd_step=self.step)

    def data_process_save_nodes_to_file(self, ddd_data_dir=None):
        import pickle
        if self.target_step_list is not None:   counter = self.target_step_list_counter
        elif self.target_time_list is not None: counter = self.target_time_list_counter
        else: raise ValueError("both target_step_list and target_time_list are None")
        self.rsNum = int(self.step/100)
        if ddd_data_dir == None:
            load_base = os.path.basename(self.save_data_path)  # 'graph_001_123_9'
            _, v1_str, v2_str, v3_str = load_base.split('_')  # ['graph', '001', '123', '9']
            vec_1 = [int(c) for c in v1_str]
            vec_2 = [int(c) for c in v2_str]
            num = int(v3_str)
            load_info = [vec_1, vec_2, num]
            ddd_data_dir = os.path.join(os.environ['OAK_DIR'], 'shared_data/DislocationDynamics/hamed_ddd_data/15um_1e3_changeOrientation/',
                                        get_ddd_database_folder_name(load_info))
            self.ctrlfile = f'rs{self.rsNum:04d}'
        self.ctrl_info = self.parse_ctrl_file(os.path.join(ddd_data_dir, 'restart/', self.ctrlfile))
        with open(f'{self.save_data_path}/oprec_player_saved_{self.step}.pkl', 'wb') as f:
            phys_nodes = getattr(self, 'phys_nodes', None)
            pickle.dump({'ctrl_info': self.ctrl_info,
                         'nodes': self.nodes, 'dt': self.dt,
                         'step': self.step, 'phys_nodes': phys_nodes}, f)
            if getattr(self, 'phys_nodes', None) is not None:
                self.print_func(f"---> saved nodes (w/ phys_nodes) & info to {self.save_data_path}/oprec_player_saved_{self.step}.pkl")
            else: self.print_func(f"---> saved nodes & info to {self.save_data_path}/oprec_player_saved_{self.step}.pkl")

    def need_break_after_play_one_block(self):
        if self.end_of_file:
            return True

        if self.target_step_list is not None:
            if self.step >= self.target_step_list[-1]:
                self.print_func('step = ' + str(self.step) + ' target_step (' + str(self.target_step_list[-1]) + ') has reached')
                return True

        if self.target_time_list is not None:
            if self.time_now >= self.target_time_list[-1]:
                self.print_func(f"reach target time:\n-> time_now = {self.time_now}, target_time = {self.target_time_list[-1]}")
                return True
        return False

    def need_break_after_play_one_op(self):
        if self.end_of_file:
            return True

        if self.optype in self.separate_step_at_optype:
            return True
        return False

    def plot_nodes_during_play(self):
        if self.plot_physical_structure:
            phys_nodes = physical_structure(nodes, bounds, calc_pl_strain=True)
            self.fig, self.ax = plot_nodes_dict(phys_nodes, bounds, fig=self.fig, ax=self.ax, plot_links=self.plot_links, trim=self.plot_trim)
        else:
            self.fig, self.ax = plot_nodes_dict(nodes, bounds, fig=self.fig, ax=self.ax, plot_links=self.plot_links, trim=self.plot_trim)

    # Note: ParaDiS algorithm in ParadisStep.c 
    # The following functions are called in this sequence in each step
    #    StartStepOpRecList(home, START_STEP);        /* mark with optype 10 */
    #    NodeMoveOpRecList(home);                     /* after moving all nodes, mark with optype 11 */
    #    FixGlideViolations(home,NULL,NULL);          /* not logged in oprec (can cause mismatch with restart file) */
    #    ResetGlidePlanes(home);                      /* not logged in oprec (can cause mismatch with restart file) */
    #    DeltaPlasticStrain(home);
    #    GenerateOutput(home, STAGE_CYCLE);           /* restart file is written here */
    #    StartStepOpRecList(home, START_SPLIT_MULTI); /* mark with optype 12 */
    #    SplitMultiNodes(home);
    #    StartStepOpRecList(home, START_CROSS_SLIP);  /* mark with optype 13 */
    #    CrossSlip(home);
    #    StartStepOpRecList(home, START_COLLISION);   /* mark with optype 14 */
    #    HandleCollisions(home);
    #    StartStepOpRecList(home, START_REMESH);      /* mark with optype 15 */
    #    Remesh(home);
    #    LoadCurve(home, deltaStress);

    def play_one_op(self, f):
    # Replay the oprecfile for one step and update the DD structure stored in the nodes dictionary
    #  (from oprec_player_on_dict)
    #
        box_size = np.array(self.bounds[1] - self.bounds[0])

        self.oprec_line = f.readline()
        if self.oprec_line == '':
            self.print_func('End of oprec file')
            self.end_of_file = True

        self.optype = int(self.oprec_line)
        if self.optype != self.prev_optype:
            if self.optype >= 6:
                self.print_func("io = %d  optype = %d"%(self.io,self.optype));
                self.prev_optype = self.optype

        if self.step is not None:
            if self.step != self.prev_step:
                self.print_func("step: %d"%(self.step))
                self.prev_step = self.step

        if self.optype == 10:
            # START_STEP
            if self.update_stress_strain:
                self.nodes_old = copy.deepcopy(self.nodes)

            self.io_inc = 1

        elif self.optype == 5:
            # RESET_COORD
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]
            my_tag = tuple(np.array(opdata)[0:2].astype(int))
            self.nodes[my_tag][0][0:3] = np.array(opdata[2:5])

            self.io_inc = 2

        elif self.optype == 11:
            # TIME_INTEGRATION
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]
            self.dt     = opdata[0]

            self.time_now += self.dt
            self.step = self.step + 1

            if self.print_freq > 0:
                if self.step % self.print_freq == 0:
                    self.print_func('Step %d, dt = %e, timeNow = %e' % (self.step, self.dt, self.time_now))

            if self.update_stress_strain: # Following lines to calcualte dpstn were copied from oprec_player_dict_gammaidot but was not tested
                self.calc_plastic_strain_from_disl_network(box_size)

            if self.plot_freq is not None and self.plot_freq > 0:
                if self.step % self.plot_freq == 0:
                    self.plot_nodes_during_play()

            # if write_movie # to be written later
            self.io_inc = 2

        elif self.optype == 4:
            # SPLIT_NODE
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]
            my_tag = tuple(np.array(opdata)[0:2].astype(int))

            # create new node with empty arms dictionary
            # print('splitting node', my_tag)
            self.nodes[my_tag] = [np.concatenate([np.array(opdata)[2:5],[0]]), {}]

            self.io_inc = 2

        elif self.optype == 0:
            # CHANGE_CONNECTION
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]

            tag_1 = tuple(np.array(opdata)[0:2].astype(int))
            tag_2 = tuple(np.array(opdata)[2:4].astype(int))
            tag_3 = tuple(np.array(opdata)[4:6].astype(int))

            self.arms_dict = self.nodes[tag_1][1]
            if not tag_2 in list(self.arms_dict.keys()):
                print("oprec_player_on_dict(): cannot find tag_2 in nodes(tag_1)")
                return None
            if tag_3 in list(self.arms_dict.keys()):
                # tag_3 already exist as a neighbor for tag_1
                # merge double links
                self.arms_dict[tag_3][0:3] += self.arms_dict[tag_2][0:3]
                self.arms_dict.pop(tag_2)
            else:
                self.arms_dict[tag_3] = self.arms_dict.pop(tag_2)
            self.io_inc = 2

        elif self.optype == 1:
            # INSERT_ARM
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]

            tag_1 = tuple(np.array(opdata)[0:2].astype(int))
            tag_2 = tuple(np.array(opdata)[2:4].astype(int))

            self.arms_dict = self.nodes[tag_1][1]
            burg_vec = np.array(opdata)[4:7]
            if tag_2 in list(self.arms_dict.keys()):
                # tag_2 is already a neighbor of tag_1, merge double links
                self.arms_dict[tag_2][0:3] += np.array(opdata)[4:7]
            else:
                self.arms_dict[tag_2] = np.array(opdata)[4:10]
            self.io_inc = 2

        elif self.optype == 3:
            # CHANGE_ARM_BURG
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]

            tag_1 = tuple(np.array(opdata)[0:2].astype(int))
            tag_2 = tuple(np.array(opdata)[2:4].astype(int))

            self.arms_dict = self.nodes[tag_1][1]
            burg_vec = np.array(opdata)[4:7]
            if not tag_2 in list(self.arms_dict.keys()):
                if np.linalg.norm(burg_vec) > 0:
                    print("oprec_player_on_dict(): [3] tag_2 not a neighbor of nodes(tag_1)")
                    print(tag_1); print(tag_2); print(nodes[tag_1]); print(nodes[tag_2])
                    raise ValueError("error in optype 3")
            else:
                self.arms_dict[tag_2][0:6] = np.array(opdata)[4:10]
                burg_vec = self.arms_dict[tag_2][0:3]
                if np.linalg.norm(burg_vec) > 0:
                    pass
                else:
                    self.arms_dict.pop(tag_2)
            self.io_inc = 2

        elif self.optype == 9:
            # RESET_GLIDE_PLANE
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]

            tag_1 = tuple(np.array(opdata)[0:2].astype(int))
            tag_2 = tuple(np.array(opdata)[2:4].astype(int))

            self.arms_dict = self.nodes[tag_1][1]
            if not tag_2 in list(self.arms_dict.keys()):
                print("oprec_player_on_dict(): [9] tag_2 not a neighbor of nodes(tag_1)")
                raise ValueError("error encountered in oprec.play_one_op()")

            self.arms_dict[tag_2][3:6] = np.array(opdata)[4:7]
            self.io_inc = 2

        elif self.optype == 2:
            # REMOVE_NODE
            self.oprec_line = f.readline()
            opdata = [float(e) for e in self.oprec_line.split()]

            tag_1 = tuple(np.array(opdata)[0:2].astype(int))

            self.nodes.pop(tag_1)
            self.io_inc = 2

        elif self.optype >= 12 and self.optype <= 15:
            # START_*: start of stages log, just skip it
            # 10: START_STEP,
            # 11: TIME_INTEGRATION,
            # 12: START_SPLIT_MULTI,
            # 13: START_CROSS_SLIP,
            # 14: START_COLLISION,
            # 15: START_REMESH
            self.io_inc = 1

        else: print('WARNING: Unknown operation type ' + str(self.optype) + '. Skipping.')
        self.io = self.io + self.io_inc
        return

    @staticmethod
    def sym_tensor_to_voigt(tensor):
        '''convert 3x3 symmetric tensor to 6x1 voigt notation'''
        return np.array([
                tensor[0,0], tensor[1,1], tensor[2,2],
                tensor[1,2], tensor[0,2], tensor[0,1]])

    @staticmethod
    def asym_tensor_to_array(tensor):
        '''convert 3x3 anti-symmetric tensor to 6x1 vector'''
        return np.array([0.0, 0.0, 0.0, tensor[1,2], tensor[0,2], tensor[0,1]])

    def calc_plastic_strain_from_disl_network(self, box_size):
        # should be consistent with ParaDiS code:
        # https://gitlab.com/micronano/ParaDiS/-/blob/master/src/LoadCurve.c
        # https://gitlab.com/micronano/ParaDiS/-/blob/master/src/DeltaPlasticStrain_FCC.c
        list_keys = list( self.nodes.keys() )
        self.d_pl_strain_sum, self.d_pl_rotation_sum = np.zeros((3,3)), np.zeros((3,3))

        dstn, dspn, seg_num = np.zeros((6,)), np.zeros((3,)), 0
        for counter_my_node, my_node_tag in enumerate(list_keys):
            self.arms_dict = self.nodes[my_node_tag][1]
            for counter, nbr_tag in enumerate(self.arms_dict.keys()):
                link = self.arms_dict[nbr_tag]

                if nbr_tag > my_node_tag: continue

                r1 = self.nodes[my_node_tag][0][0:3]
                r2 = pbc_position(r1, self.nodes[nbr_tag][0][0:3], box_size)
                r3 = pbc_position(r1, self.nodes_old[my_node_tag][0][0:3], box_size)
                r4 = pbc_position(r3, self.nodes_old[nbr_tag][0][0:3], box_size)

                Sn = 0.5*np.cross(r2-r3, r1-r4) # Sn_i
                swept_area = np.linalg.norm(Sn)

                # V = L**3 # generalized to non-cubic box
                V = box_size[0] * box_size[1] * box_size[2]
                burg_vec = self.arms_dict[nbr_tag][0:3]
                Snb = np.outer(Sn, burg_vec) # TODO: figure out why the factor of 2 in Snb (don't think factor of 2 is correct)
                d_pl_strain = (1/2) * (Snb + Snb.T) / V
                d_pl_rotation = (1/2) * (Snb.T - Snb) / V # sign reverse corresponding to ParaDiS bug fix in DeltaPlasticStrain.c (line 41) by Amin (2015)

                self.d_pl_strain_sum += d_pl_strain
                self.d_pl_rotation_sum += d_pl_rotation
                if self.voigt_test:
                    n = 0.5*np.cross(r2-r3, r1-r4)
                    P = np.multiply( n[[0,0,0,1,1,1,2,2,2]], link[[0,1,2,0,1,2,0,1,2]] )
                    dstn += (0.5 / V) * ( P[[0,4,8,1,2,5]] + P[[0,4,8,3,6,7]] )
                    dspn += (0.5 / V) * ( P[[1,2,5]] - P[[3,6,7]] )

        phispn = np.array([-self.d_pl_rotation_sum[1,2], self.d_pl_rotation_sum[0,2], -self.d_pl_rotation_sum[0,1]])
        if self.voigt_test:
            # TODO: figure out why phispn_ref has the reversed order of dspn
            phispn_ref = -1.0 * np.array([-dspn[2], dspn[1], -dspn[0]]) # extra -1 needed due to ParaDiS bug fix in DeltaPlasticStrain.c (line 41) by Amin (2015)
            strain_voigt = self.sym_tensor_to_voigt(self.d_pl_strain_sum)
            dstn_voigt = np.array([dstn[0], dstn[1], dstn[2], dstn[5], dstn[4], dstn[3]])
            strain_err = np.abs(strain_voigt - dstn_voigt)
            rotat_err = np.abs(phispn - phispn_ref)
            if np.max(strain_err)<self.voigt_test_threshold and np.max(rotat_err)<self.voigt_test_threshold:
                self.print_func(f'Strain & rotation tensor -- compare with Hamed implementation {TerminalColors.GREEN}PASS{TerminalColors.RESET}')
            else:
                self.print_func(f'Strain & rotation tensor -- compare with Hamed implementation {TerminalColors.RED}FAIL{TerminalColors.RESET}')
                self.print_func(f'\n strain (Hamed)-> {dstn_voigt} \n -> strain (here) {strain_voigt} \n -> error {strain_err}' )
                self.print_func(f'\n rotation (Hamed)-> {phispn_ref} \n -> rotation (here) {phispn} \n -> error {rotat_err}' )

        self.ctrl_info['cycleStart']        = self.step
        self.ctrl_info['timeNow']           = self.time_now

    def update_applied_stress(self):
        plastic_strain_array = np.array(self.ctrl_info['totpStn']) + self.sym_tensor_to_voigt(self.d_pl_strain_sum)
        plastic_spin_array   = np.array(self.ctrl_info['totpSpn']) + self.asym_tensor_to_array(self.d_pl_rotation_sum)
        self.ctrl_info['totpStn'] = plastic_strain_array.tolist()
        self.ctrl_info['totpSpn'] = plastic_spin_array.tolist()

        if self.ctrl_info['loadType'] == 1:
            phispn = np.array([-self.d_pl_rotation_sum[1,2], self.d_pl_rotation_sum[0,2], -self.d_pl_rotation_sum[0,1]])
            Rspin = spin_matrix(phispn)
            if type(self.edotdir) is not np.ndarray:
                raise TypeError(f"self.edotdir should be a numpy array, but got {type(self.edotdir)}")
            self.edotdir = np.matmul(Rspin, self.edotdir)
            self.edotdir = self.edotdir / np.linalg.norm(self.edotdir)
            # dpstn = edotdir.T @ (self.d_pl_strain_sum @ edotdir) # dpstn = np.matmul(self.d_pl_strain_sum, edotdir).T @ edotdir
            self.dpstn = np.dot(self.edotdir.T, np.dot(self.d_pl_strain_sum, self.edotdir))

            self.print_func(f'loading direction: {self.edotdir}')
            self.d_normal_strain = self.dt * self.eRate
            self.d_normal_stress = (2*self.mu*(1+self.nu)) * (self.d_normal_strain - self.dpstn)
            self.print_func(f' -> {TerminalColors.BLUE}Stress strain analysis: Δε, Δε_pl, Δσ, Δt{TerminalColors.RESET}: {self.d_normal_strain}, {self.dpstn}, {self.d_normal_stress}, {self.dt}')
            applied_stress_array = np.array(self.ctrl_info['appliedStress']) + self.sym_tensor_to_voigt(self.d_normal_stress * np.outer(self.edotdir, self.edotdir))
            total_strain_array   = np.array(self.ctrl_info['totstraintensor']) + self.sym_tensor_to_voigt(self.d_normal_strain * np.outer(self.edotdir, self.edotdir))
        else:
            raise NotImplementedError(f'loadType {self.ctrl_info["loadType"]} not implemented yet')

        self.ctrl_info['edotdir']           = self.edotdir.tolist()
        self.ctrl_info['appliedStress']     = applied_stress_array.tolist()
        self.ctrl_info['totstraintensor']   = total_strain_array.tolist()


class link_statistics():
    '''class for analyzing and storing link statistics (length, dGamma)'''
    def __init__(self, nbins_links=126, max_normalized_linklength=25,
                 nslipsys=13, edotdir_init=None, print_func=print):
        self.nbins_links = nbins_links
        self.max_normalized_linklength = max_normalized_linklength
        self.nslipsys = nslipsys
        self.edotdir_init = edotdir_init
        self.print_func = print_func

    def init_link_length_arrays(self):
        self.lBarPerSystem = np.zeros((self.data_length, self.nslipsys))
        self.numberOfLinks = np.zeros((self.data_length, self.nslipsys))
        self.nLinks_bin = np.zeros((self.data_length, self.nslipsys, self.nbins_links))
        self.lBar_overall = np.zeros((self.data_length))
        self.lBar_dir_overall = np.zeros((self.data_length))
        self.numberOfLinks_overall = np.zeros(self.data_length)
        self.nLinks_bin_overall = np.zeros((self.data_length, self.nbins_links))
        self.nLinks_dir_bin_overall = np.zeros((self.data_length, self.nbins_links))
        self.bins_links = np.linspace(0, self.max_normalized_linklength, self.nbins_links)
        self.dotGammai_all = np.zeros((self.data_length,self.nslipsys))
        self.dotGammai_bin = np.zeros((self.data_length,self.nslipsys,self.nbins_links))
        self.dotEpsilon_bin_overall = np.zeros((self.data_length,self.nbins_links))

    def add_link_length_per_slipsys_hist(self, phys_nodes, print_analysis_details=False):
        step = self.inc_step
        PhysLength_all, DirLength_all, SlipSystem_links = self.linkLength_system(phys_nodes, IncludeJunc=True)
        for sysID in range(self.nslipsys):
            self.lBarPerSystem[step, sysID] = np.mean( PhysLength_all[sysID] )
            self.numberOfLinks[step, sysID] = len( PhysLength_all[sysID] )
            self.nLinks_bin[step, sysID, 0:(self.nbins_links-1) ] = np.histogram( PhysLength_all[sysID] / self.lBarPerSystem[step,sysID], self.bins_links )[0]
            if print_analysis_details:
                # self.print_func('\033[95mSystem ID\033[0m:', sysID+1, '| \033[95mMean Physical Length\033[0m: ', np.mean( PhysLength_all[sysID]))
                self.print_func(f'\033[95mSystem ID\033[0m: {sysID+1} | \033[95mMean Physical Length\033[0m: {self.lBarPerSystem[step, sysID]}')

        # link_lengths_overall = extract_link_length(phys_nodes) # old method
        phy_link_lengths_overall = np.array([j for i in PhysLength_all for j in i])
        dir_link_lengths_overall = np.array([j for i in DirLength_all for j in i])
        link_lengths_overall = np.vstack((dir_link_lengths_overall, phy_link_lengths_overall)).T

        self.lBar_overall[step] = np.mean( link_lengths_overall[:,1] )
        self.lBar_dir_overall[step] = np.mean( link_lengths_overall[:,0] )
        self.numberOfLinks_overall[step] = len( link_lengths_overall[:,1] )
        self.nLinks_bin_overall[step, 0:(self.nbins_links-1) ] = np.histogram( link_lengths_overall[:,1] / self.lBar_overall[step], self.bins_links )[0]
        self.nLinks_dir_bin_overall[step, 0:(self.nbins_links-1) ] = np.histogram( link_lengths_overall[:,0] / self.lBar_dir_overall[step], self.bins_links )[0]

    def add_dGamma_per_slipsys_hist(self, phys_nodes, time_step, slip_sys, IncludeJunc=True, 
                                    test_plotting=False, print_analysis_details=False, print_func=print,
                                    save_intermediate_graph=False, save_data_path=None, ddd_step=0):
        step = self.inc_step
        nRows = 13 if IncludeJunc else 12

        PhysLength_all, DirLength_all, SlipSystem_links = self.linkLength_system(phys_nodes, IncludeJunc=True)
        d_pl_strain_all, dGammai_all = self.linkDplstrain_system(phys_nodes, slip_sys, IncludeJunc=True)

        edotdir = self.edotdir_init
        Schmid_dict = slip_sys.get_slip_Snb(edotdir)[0]['SchmidFactor']
        Schmid = [ Schmid_dict[sysID] for sysID in range(13) ]

        dotGamma_all = 0
        for sysID in range(self.nslipsys):
            if sysID==12:
                break

            self.lBarPerSystem[step, sysID] = np.mean( PhysLength_all[sysID] )
            dotGammai_all_per_link = np.concatenate(dGammai_all[sysID]) / time_step
            self.dotGammai_all[step, sysID] = np.sum(dotGammai_all_per_link)
            if print_analysis_details:
                self.print_func('\033[95mSystem ID\033[0m:',sysID+1, '| \033[95mShear Strain Rate\033[0m: ', self.dotGammai_all[step, sysID])
            for linkID, link_len in enumerate(PhysLength_all[sysID]/self.lBarPerSystem[step,sysID]):
                binID = np.minimum(int(link_len / self.max_normalized_linklength * (self.nbins_links-1)), self.nbins_links-1)
                self.dotGammai_bin[step, sysID, binID] += dotGammai_all_per_link[linkID]

            dotGamma_all += self.dotGammai_all[step, sysID] * Schmid[sysID]
        self.print_func(f'-> {TerminalColors.BLUE}Overall Gamma rate for step{TerminalColors.RESET} {step+1}: {dotGamma_all}')
        if save_intermediate_graph:
            import pickle
            if save_data_path is None: save_data_path = 'tmp'
            with open(f'{save_data_path}/LinkLen_dGamma_IncStep{ddd_step}.pkl', 'wb') as f:
                pickle.dump({'dGammai_all': dGammai_all,
                             'PhysLength_all': PhysLength_all,
                             'DirLength_all': DirLength_all,}, f)
            self.print_func(f'---> saved link length and dGamma data to LinkLen_dGamma_IncStep{step+1}.pkl')

    def linkLength_system(self, phys_nodes, IncludeJunc=True, ignoreLoops=True):
        is_phys = is_physical_node(phys_nodes)
        if not all(is_phys):
            raise ValueError('not all nodes are physical')

        nRows = 13 if IncludeJunc else 12

        PhysLength_all = [ [] for i in range(nRows) ]
        DirLength_all = [ [] for i in range(nRows) ]
        SlipSystem_links = [ [] for i in range(nRows) ]

        for my_tag, node_attr in phys_nodes.items():
            arms_extra_dict = node_attr[2]
            for nbr_tag, arm_extra_attr in arms_extra_dict.items():
                dest_tag = arm_extra_attr['dest_tag']
                if dest_tag > my_tag: continue
                if dest_tag == my_tag and ignoreLoops: continue
                double_counting_factor = 0.5 if (dest_tag == my_tag) else 1.0
                link_sysID = arms_extra_dict[nbr_tag]['sysID']
                if (link_sysID==12 and not IncludeJunc):    # if junction, move on!
                    continue
                # TODO: fix double-counting issue for dislocation loops:
                # instead of two links with (both) half length, we should have one link with full length
                if double_counting_factor != 1.0:
                    raise ValueError('double counting factor is not 1, which means that the link is a loop')
                PhysLength_all[link_sysID].append( arms_extra_dict[nbr_tag]['physical_length'] * double_counting_factor )
                DirLength_all[link_sysID].append( arms_extra_dict[nbr_tag]['direct_length'] * double_counting_factor )
                link_info = {'node_tag': my_tag, 'nbr_tag': nbr_tag, 'dest_tag': dest_tag,
                             'physical_length': PhysLength_all, 'direct_length': DirLength_all, 'sysID': link_sysID }
                SlipSystem_links[link_sysID].append(link_info)
        return PhysLength_all, DirLength_all, SlipSystem_links

    def linkDplstrain_system(self, phys_nodes, slip_sys, IncludeJunc=True, ignoreLoops=True):
        is_phys = is_physical_node(phys_nodes)
        if not all(is_phys):
            raise ValueError('not all nodes are physical')

        nRows = 13 if IncludeJunc else 12
        d_pl_strain_all, dGammai_all = [ [] for i in range(nRows) ], [ [] for i in range(nRows) ]

        for my_tag in list(phys_nodes.keys()):
            node_attr = phys_nodes[my_tag]
            my_coords = node_attr[0][0:3]
            arms_dict, arms_extra_dict = node_attr[1], node_attr[2]

            for nbr_tag, arm_extra_attr in arms_extra_dict.items():
                dest_tag = arm_extra_attr['dest_tag']
                if dest_tag > my_tag: continue
                if dest_tag == my_tag and ignoreLoops: continue
                double_counting_factor = 0.5 if (dest_tag == my_tag) else 1.0
                link_sysID = arms_extra_dict[nbr_tag]['sysID']

                if (link_sysID==12 and not IncludeJunc):    # if junction, move on!
                    continue

                burg_vec, plane_normal = slip_sys.sysID_to_nb(link_sysID);

                if burg_vec is not None and plane_normal is not None:
                    plane_normal_unit, burg_vec_unit = plane_normal / np.linalg.norm(plane_normal), burg_vec / np.linalg.norm(burg_vec)
                    n_dot_d_pl_strain_dot_b = np.dot(np.dot(plane_normal_unit.T, arms_extra_dict[nbr_tag]['link_d_pl_strain']), burg_vec_unit) #arms_extra_dict[nbr_tag]['d_pl_strain']
                    '''definition of gamma: γij = 2εij'''
                    dGammai_all[link_sysID].append( 2 * n_dot_d_pl_strain_dot_b * double_counting_factor )

                    # Implementation for d_pl_strain_all (to be tested)
                    d_pl_strain = arm_extra_attr['link_d_pl_strain'] * double_counting_factor
                    d_pl_strain_all[link_sysID].append(d_pl_strain)

                    '''for testing purpose'''
                    b_test, n_test = burg_vec/np.linalg.norm(burg_vec), plane_normal/np.linalg.norm(plane_normal)
                    b_test, n_test = np.array(b_test).T, np.array(n_test).T
                    sysID_check = int(slip_sys.nb_to_sysID(n_test, b_test)[0])
                    if sysID_check!=12:
                        if link_sysID != sysID_check:
                            self.print_func(f'{TerminalColors.RED}Slip system ID mismatch{TerminalColors.RESET}: {link_sysID}, {sysID_check}')
        return d_pl_strain_all, dGammai_all

    def error_analysis_per_oprec_step(self, data_gen_op3, data_ref_op3, maxstep=5):
        nLinks_bin_ref = data_ref_op3['nLinks_bin'] 
        dotGammai_bin_ref = data_ref_op3['dotGammai_bin']
        lBarPerSystem_ref = data_ref_op3['lBarPerSystem']
        numberOfLinks_ref = data_ref_op3['numberOfLinks']
        bins_links_ref = data_ref_op3['bins_links']
        dotEpsilon_bin_overall_ref = data_ref_op3['dotEpsilon_bin_overall']
        nLinks_bin_overall_ref = data_ref_op3['nLinks_bin_overall']
        lBar_overall_ref = data_ref_op3['lBar_overall']
        nLinks_dir_bin_overall_ref = data_ref_op3['nLinks_dir_bin_overall']
        lBar_dir_overall_ref = data_ref_op3['lBar_dir_overall']

        nLinks_bin = data_gen_op3['nLinks_bin'] 
        dotGammai_bin = data_gen_op3['dotGammai_bin']
        lBarPerSystem = data_gen_op3['lBarPerSystem']
        numberOfLinks = data_gen_op3['numberOfLinks']
        bins_links = data_gen_op3['bins_links']
        dotEpsilon_bin_overall = data_gen_op3['dotEpsilon_bin_overall']
        nLinks_bin_overall = data_gen_op3['nLinks_bin_overall']
        lBar_overall = data_gen_op3['lBar_overall']
        nLinks_dir_bin_overall = data_gen_op3['nLinks_dir_bin_overall']
        lBar_dir_overall = data_gen_op3['lBar_dir_overall']

        color = TerminalColors
        nLinks_bin_error = nLinks_bin[0:maxstep,:,:] - nLinks_bin_ref[0:maxstep,:,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for nLinks_bin:{color.RESET}', np.max(nLinks_bin_error))

        dotGammai_bin_error = dotGammai_bin[0:maxstep,:,:] - dotGammai_bin_ref[0:maxstep,:,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for dotGammai_bin:{color.RESET}', np.max(dotGammai_bin_error))

        lBarPerSystem_error = lBarPerSystem[0:maxstep,:] - lBarPerSystem_ref[0:maxstep,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for lBarPerSystem:{color.RESET}', np.max(lBarPerSystem_error))

        numberOfLinks_error = numberOfLinks[0:maxstep,:] - numberOfLinks_ref[0:maxstep,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for numberOfLinks:{color.RESET}', np.max(numberOfLinks_error))

        dotEpsilon_bin_overall_error = dotEpsilon_bin_overall[0:maxstep,:] - dotEpsilon_bin_overall_ref[0:maxstep,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for dotEpsilon_bin_overall:{color.RESET}', np.max(dotEpsilon_bin_overall_error))

        nLinks_bin_overall_error = nLinks_bin_overall[0:maxstep,:] - nLinks_bin_overall_ref[0:maxstep,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for nLinks_bin_overall:{color.RESET}', np.max(nLinks_bin_overall_error))

        lBar_overall_error = lBar_overall[:maxstep] - lBar_overall_ref[:maxstep]
        print(f'{color.BRIGHT_MAGENTA}Max error for lBar_overall:{color.RESET}', np.max(lBar_overall_error))

        lBar_dir_overall_error = lBar_dir_overall[:maxstep] - lBar_dir_overall_ref[:maxstep]
        print(f'{color.BRIGHT_MAGENTA}Max error for lBar_dir_overall:{color.RESET}', np.max(lBar_dir_overall_error))

        nLinks_dir_bin_overall_error = nLinks_dir_bin_overall[0:maxstep,:] - nLinks_dir_bin_overall_ref[0:maxstep,:]
        print(f'{color.BRIGHT_MAGENTA}Max error for nLinks_dir_bin_overall:{color.RESET}', np.max(nLinks_dir_bin_overall_error))

        ALL_ERRORS = np.array([nLinks_bin_error.max(), dotGammai_bin_error.max(), lBarPerSystem_error.max(), numberOfLinks_error.max(), \
            dotEpsilon_bin_overall_error.max(), nLinks_bin_overall_error.max(), lBar_overall_error.max(),\
                lBar_dir_overall_error.max(), nLinks_dir_bin_overall_error.max()])
        return ALL_ERRORS

    def save_to_mat_old(self, matfile_name):
        from scipy.io import savemat
        # for oprec player for all loading directions
        folder = os.path.dirname(matfile_name)
        if folder and not os.path.isdir(folder):
            os.makedirs(folder, exist_ok=True)
        # -------------------------------------------
        # savemat(matfile_name, vars(self))
        data_to_save = {k: v for k, v in vars(self).items() if v is not None}
        savemat(matfile_name, data_to_save)
        print(f"---> saved to {matfile_name}")
    def save_to_mat(self, matfile_name):
        from scipy.io import savemat
        import os, types
        folder = os.path.dirname(matfile_name)
        if folder and not os.path.isdir(folder): os.makedirs(folder, exist_ok=True)
        # keep only MATLAB-serializable-ish things; drop callables/modules/classes/None
        data_to_save = {k: v for k, v in vars(self).items()
                        if (v is not None)
                        and not callable(v)
                        and not isinstance(v, (types.ModuleType, type))}
        savemat(matfile_name, data_to_save, do_compression=True)
        print(f"---> saved to {matfile_name}")

    def save_mat_link_avg(self, matFileName):
        from scipy.io import savemat
        savemat(matFileName, {"dotGammai_bin_avg":self.dotGammai_bin_avg[:,0,:,:].reshape(9,13,self.nbins_links),#self.dotGammai_bin_avg,
                "nLinks_bin_avg":self.nLinks_bin_avg[:,0,:,:].reshape(9,13,self.nbins_links), 
                "Vi": self.Vi[:,0,:,:].reshape(9,13,self.nbins_links),
                "numberOfLinks_avg":self.numberOfLinks_avg[:,0,:].reshape(9,13), 
                "lBarPerSystem_avg":self.lBarPerSystem_avg[:,0,:].reshape(9,13),
                "Rhoi_avg":self.Rhoi_avg[:,0,:].reshape(9,13), 
                "tau_avg":self.tau_avg[:,0,:].reshape(9,12), 
                "stress_avg":self.stress_avg[:,0].reshape(9), 
                "strain_avg":self.strain_avg[:,0].reshape(9),
                "niL_fitting_coefs":self.niL_fitting_coefs, 
                "niL_fitting_coefs_lasso":self.niL_fitting_coefs_lasso,
                "Schmids":self.Schmids})