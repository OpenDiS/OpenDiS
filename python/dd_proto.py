import networkx as nx
import numpy as np

class Node:
    def __init__(self, tag, pos):
        self.tag = tag
        self.pos = pos
        self.force = np.array([0.0, 0.0, 0.0])
        self.vel   = np.array([0.0, 0.0, 0.0])

class Link:
    def __init__(self, burg, plane_normal):
        self.burg = burg
        self.plane_normal = plane_normal
        #self.node1 = 0
        #self.node2 = 0

# DisNetwork is a Graph for which each node has attributes specified by class Node
#  and each edge has attributes specified by class Link
class DisNetwork(nx.Graph):
    def __init__(self, data=None, **attr):
        nx.Graph.__init__(self, data, **attr)

class DDParam():
    def __init__(self, mu=1.0, nu=0.3, box=None):
        self.mu = mu
        self.nu = nu
        self.box = box

class DDSim():
    def __init__(self, param=None):
        self.disnet = DisNetwork()
        self.param = DDParam()

    def insert_node(self, node1, node2, pos=None):
        # insert a new node on the link connecting node1 and node2
        pass

    def remove_two_arm_node(self, node):
        # remove two-arm node
        pass

    def merge_node(self, node1, node2):
        # merge two nodes into one node
        pass

    def split_node(self, node, partition):
        # split node into two nodes, with neighbors split according to partition
        pass

    def LinkForce(self):
        # compute forces on every link
        pass

    def NodeForce(self):
        # compute forces on every node
        pass

    def MobilityLaw(self):
        # compute velocity on all nodes
        pass

    def TimeIntegration(self):
        # advance nodes by one time step
        pass