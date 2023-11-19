"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
import networkx as nx
from typing import Tuple

Tag = Tuple[int, int]

class DisNode:
    """DisNode: class for dislocation node

    Defines the basic attributes on a node
    """
    def __init__(self, R: np.ndarray) -> None:
        self.R = R

class DisEdge:
    """DisEdge: class for dislocation edge

    Defines the basic features on a edge
    """
    def __init__(self, burg_vec: np.ndarray, plane_normal: np.ndarray=None) -> None:
        self.burg_vec = burg_vec
        if plane_normal is not None:
            self.plane_normal = plane_normal

class DisNet:
    """DisNet: class for dislocation network

    Implements basic topological operations on dislocation networks
    """
    def __init__(self, data=None, **attr) -> None:
        self._G = nx.DiGraph(data, **attr)

    def neighbors(self, tag: Tag):
        """neighbors: return neighbors (as iterator) of a node
        """
        return self._G.neighbors(tag)

    def nodes(self):
        """nodes: return node view of the network
        """
        return self._G.nodes
    
    def edges(self):
        """edges: return out edge view of the network
        """
        return self._G.edges
    
    def pos_array(self) -> np.ndarray:
        """pos_array: return a numpy array of node positions
        """
        return np.array([node.R for node in self._G.nodes()])
    
    def seg_list(self) -> list:
        """seg_list: return a list of segments

        To do: also return the Burgers vector and plane normal
        """
        return [(node1.R, node2.R) for node1, node2 in self._G.edges() if node1 < node2]
    
    def has_node(self, tag: Tag) -> bool:
        """has_node: check if a node exists in the network
        """
        return self._G.has_node(tag)
    
    def has_edge(self, node1: Tag, node2: Tag) -> bool:
        """has_edge: check if an edge exists in the network
        """
        return self._G.has_edge(node1, node2)
    
    def _add_node(self, tag: Tag, node_attr: DisNode) -> None:
        """add_node: add a node to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self._G.add_node(tag, **vars(node_attr))
    
    def _add_edge(self, node1: Tag, node2: Tag, edge_attr: DisEdge) -> None:
        """add_edge: add an edge to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self._G.add_edge(node1, node2, **vars(edge_attr))
    
    def _remove_edge(self, node1: Tag, node2: Tag) -> None:
        """remove_edge: remove an edge from the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self._G.remove_edge(node1, node2)

    def add_nodes_links_from_list(self, rn, links) -> None:
        """add_nodes_links_from_list: add nodes and edges stored in lists to network
           sanity after this operation depends on the input
        """
        N = rn.shape[0]
        num_links = links.shape[0]
        for i in range(N):
            self._add_node((0,i), DisNode(R=rn[i]))
        for j in range(num_links):
            seg = links[j, :2].astype(int)
            bv  = links[j, 2:]
            node = (0,seg[0])
            nbr_node = (0, seg[1])
            self._add_edge(node, nbr_node, DisEdge(burg_vec= bv))
            self._add_edge(nbr_node, node, DisEdge(burg_vec=-bv))

        if not self.is_sane():
            raise ValueError("add_nodes_links_from_list: sanity check failed")
    
    def insert_node(self, node1: Tag, node2: Tag, tag: Tag, R: np.ndarray) -> None:
        """insert_node: insert a node between two existing nodes
           guarantees sanity after operation
        """
        self._add_node(tag, R=R)
        link12_attr = self.edges[(node1, node2)]
        self._add_edge(node1, tag, **link12_attr)
        self._add_edge(tag, node2, **link12_attr)
        link21_attr = self.edges[(node2, node1)]
        self._add_edge(node2, tag, **link21_attr)
        self._add_edge(tag, node1, **link21_attr)
        self._remove_edge(node1, node2)
        self._remove_edge(node2, node1)
    
    def remove_two_arm_node(self, tag: Tag) -> None:
        """remove_two_arm_node: remove a node with two arms from the network
           guarantees sanity after operation
        """
        if not self.has_node(tag):
            raise ValueError("remove_two_arm_node: Node %s does not exist" % str(tag))
        if self.out_degree(tag) != 2:
            raise ValueError("remove_two_arm_node: Node %s does not have two arms" % str(tag))

        node1, node2 = self.neighbors(node)
        link12_attr = self.edges[(tag, node2)]
        self._add_edge(node1, node2, **link12_attr)
        link21_attr = self.edges[(tag, node1)]
        self._add_edge(node2, node1, **link21_attr)
        self._remove_edge(tag, node1)
        self._remove_edge(node1, tag)
        self._remove_edge(tag, node2)
        self._remove_edge(node2, tag)
        self._remove_node(tag)
    
    def merge_node(self, node1: Tag, node2: Tag) -> None:
        """merge_node: merge two nodes into one
           guarantees sanity after operation

        Remove any links between node1 and node2
        If node1 has double-links to any neighbor, combine them into one (or zero) link
        """
        pass
    
    def split_node(self, node: Tag, partition: list) -> None:
        """split_node: split a node into two nodes with neighbors split according to partition
           guarantees sanity after operation

        """
        pass

    def is_sane(self, tol: float=1e-8) -> bool:
        """is_sane: check if the network is sane
           guarantees sanity after operation

        The two arms connecting two nodes should have opposite Burgers vectors and parallel plane_normal vectors
        The sum of all Burgers vectors of outgoing arms from a node should be zero
        """
        for node, _ in self.nodes().items():
            for nbr_node in self.neighbors(node):
                b12 = self.edges()[(node, nbr_node)]['burg_vec']
                b21 = self.edges()[(nbr_node, node)]['burg_vec']
                if np.max(np.abs(b12 + b21)) > tol:
                    print("Burgers vector of edge (%s, %s) is not opposite to that of edge (%s, %s)" % (str(node), str(nbr_node), str(nbr_node), str(node)))
                    return False

        for node, _ in self.nodes().items():
            b_tot = np.zeros(3)
            for nbr_node in self.neighbors(node):
                b_tot += self.edges()[(node, nbr_node)]['burg_vec']
            if np.max(np.abs(b_tot)) > tol:
                print("Total Burgers vector from node %s is not zero" % (str(node), str(nbr_node), str(nbr_node), str(node)))
                return False

        return True
    
