import networkx as nx
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt

opendis_lib = __import__('opendis')

def pbc_position_L(r1, r2, L):
    ds = (r2-r1)/L
    ds = ds - ds.round()
    return r1 + ds*L

# DisNetwork is an extension of DiGraph for which each node has attributes such as R, F, V
#  and each edge has attributes such as burg_vec and plane_normal
class DisNetwork(nx.DiGraph):
    def __init__(self, data=None, **attr):
         super(DisNetwork, self).__init__()

    def pos_array(self):
        return np.array([self.nodes[node]['R'] for node in self.nodes])
    
    def insert_node(self, node1, node2, tag, pos=None):
        # insert a new node on the link connecting node1 and node2
        if not self.has_edge(node1, node2) or not self.has_edge(node2, node1):
            raise ValueError('insert_node: no link exists between '+str(node1)+' and '+str(node2)) 

        self.add_node(tag, R = pos)
        link12_attr = self.edges[(node1, node2)]
        self.add_edge(node1, tag, **link12_attr)
        self.add_edge(tag, node2, **link12_attr)
        link21_attr = self.edges[(node2, node1)]
        self.add_edge(node2, tag, **link21_attr)
        self.add_edge(tag, node1, **link21_attr)
        self.remove_edge(node1, node2)
        self.remove_edge(node2, node1)

    def remove_two_arm_node(self, node):
        # remove two-arm node
        if not self.has_node(node):
            raise ValueError('remove_two_arm_node: node '+str(node)+' does not exist')
        if self.out_degree(node) != 2:
            raise ValueError('remove_two_arm_node: node '+str(node)+' does not have 2 arms')
        
        node1 = list(self.neighbors(node))[0]
        node2 = list(self.neighbors(node))[1]
        link12_attr = self.edges[(node, node2)]
        self.add_edge(node1, node2, **link12_attr)
        link21_attr = self.edges[(node, node1)]
        self.add_edge(node2, node1, **link21_attr)
        self.remove_edge(node , node1)
        self.remove_edge(node1, node )
        self.remove_edge(node , node2)
        self.remove_edge(node2, node )
        self.remove_node(node)

    def merge_node(self, node1, node2):
        # merge two nodes into one node
        pass

    def split_node(self, node, partition):
        # split node into two nodes, with neighbors split according to partition
        pass



class DDParam():
    def __init__(self, mu=1.0, nu=0.3, bounds=None):
        self.mu = mu
        self.nu = nu
        self.bounds = bounds

class DDSim():
    def __init__(self, param=None):
        self.disnet = DisNetwork()
        #self.disnet = nx.DiGraph()
        self.param = DDParam()

    def LinkForce(self):
        # compute forces on every link
        pass

    def NodeForce(self):
        # compute forces on every node by calling
        # opendis_lib.SegSegForce(...)
        pass

    def MobilityLaw(self):
        # compute velocity on all nodes
        pass

    def TimeIntegration(self):
        # advance nodes by one time step
        pass

    def plot_disnet(self, plot_links=True, trim=False, fig=None, ax=None, block=False, pause_seconds=0.01):
        if fig==None:
            fig = plt.figure(figsize=(8,8))
        if ax==None:
            ax = plt.axes(projection='3d')

        # to do: extend to non-cubic box
        L = self.param.bounds[1][0] - self.param.bounds[0][0]

        rn = self.disnet.pos_array()
        p_link = np.empty((0,6))

        plt.cla()
        if plot_links:
            for my_tag in list(self.disnet.nodes()):
                my_coords = self.disnet.nodes[my_tag]['R']
                for arm in self.disnet.out_edges(my_tag):
                    nbr_tag = arm[1]
                    if my_tag < nbr_tag:
                        r_link = np.zeros((2,3))
                        nbr_coords = self.disnet.nodes[nbr_tag]['R']
                        r_link[0,:] = my_coords
                        # to do: extend to non-cubic box
                        r_link[1,:] = pbc_position_L(my_coords, nbr_coords, L)
                        if (not trim) or np.max(np.absolute(r_link)) <= L/2:
                            p_link = np.append(p_link, [r_link[0,:], r_link[1,:]])

        ls = p_link.reshape((-1,2,3))
        lc = Line3DCollection(ls, linewidths=0.5, colors='b')
        ax.add_collection(lc)

        ax.scatter(rn[:,0], rn[:,1], rn[:,2], c='r', s=4)
        ax.set_xlim(self.param.bounds[0][0], self.param.bounds[1][0])
        ax.set_ylim(self.param.bounds[0][1], self.param.bounds[1][1])
        ax.set_zlim(self.param.bounds[0][2], self.param.bounds[1][2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_box_aspect([1,1,1])

        plt.draw()
        plt.show(block=block)
        plt.pause(pause_seconds)

        return fig, ax
