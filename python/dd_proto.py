import networkx as nx
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt

opendis_lib = __import__('opendis')

# DisNetwork is an extension of DiGraph for which each node has attributes such as R, F, V
#  and each edge has attributes such as burg_vec and plane_normal
class DisNetwork(nx.DiGraph):
    def __init__(self, data=None, **attr):
         super(DisNetwork, self).__init__(data, **attr)

    def pos_array(self):
        return np.array([self.nodes[node]['R'] for node in self.nodes])
    
    def seg_array(self):
        # construct segment list (each link appear once: node1 < node2)
        return None

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
        # merge two nodes, node1 and node2, into node1
        # redirect all links to/from node2 into node1

        # remove any links between node1 and node2

        # if node1 has double-links to any neighbor, combine them into one (or zero) link
        pass

    def split_node(self, node, partition):
        # split node into two nodes, with neighbors split according to partition
        pass



class DDParam():
    def __init__(self, mu=1.0, nu=0.3, bounds=None, boundary_cond=None, force_mode=None,
                       collision_mode=None, mobility_law=None, time_integrator=None):
        self.mu = mu
        self.nu = nu
        self.bounds = bounds
        self.boundary_cond = boundary_cond
        self.force_mode = force_mode
        self.collision_mode = collision_mode
        self.mobility_law = mobility_law
        self.time_integrator = time_integrator

        self.load_type       = None
        self.applied_stress  = np.zeros(6)

class DDSim():
    def __init__(self, param=None):
        self.disnet = DisNetwork()
        self.param = DDParam()
        self.NodeForce_Functions = {'LineTension': self.NodeForce_LineTension}
        self.Collision_Functions = {'Proximity': self.Collision_Proximity}
        self.MobilityLaw_Functions = {'LineTension': self.NodeForce_LineTension}
        self.TimeIntegration_Functions = {'EulerForward': self.TimeIntegration_EulerForward}

    def pbc_position_L(self, r1, r2, L):
        ds = (r2-r1)/L
        ds = ds - ds.round()
        return r1 + ds*L

    def NodeForce(self):
        self.NodeForce_Functions[self.param.force_mode]()

    def Collision(self):
        self.Collision_Functions[self.param.collision_mode]()

    def MobilityLaw(self):
        self.MobilityLaw_Functions[self.param.mobility_law]()

    def TimeIntegration(self):
        self.TimeIntegration_Functions[self.param.time_integrator]()

    def Step(self):
        self.NodeForce()
        self.MobilityLaw()
        self.TimeIntegration()
        self.Collision()

    def Run(self):
        # for tstep in range(self.param.maxstep):
        #     self.Step()
        pass

    def NodeForce_LineTension(self):
        # line tension forces only
        print('NodeForce_LineTension')
        # to be implemented ...
        pass

    def Collision_Proximity(self):
        # use current node position to detect collision
        print('Collision_Proximity')
        # to be implemented ...
        pass

    def MobilityLaw_FCC0(self):
        # compute velocity on all nodes
        print('MobilityLaw_FCC0')
        # to be implemented ...
        pass

    def TimeIntegration_EulerForward(self):
        # advance nodes by one time step using Eurler-Forward method
        print('TimeIntegration_EulerForward')
        # to be implemented ...
        pass

    def LinkForce(self):
        # compute forces on every link by calling
        # opendis_lib.SegSegForce(...)
        # shall we move it to an extension ?
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
                        r_link[1,:] = self.pbc_position_L(my_coords, nbr_coords, L)
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
