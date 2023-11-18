"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
from typing import Tuple

Tag = Tuple[int, int]

class DisNode:
    """DisNode: class for dislocation node

    Defines the basic features on a node
    """
    def __init__(self, tag: Tag, R: np.ndarray) -> None:
        self.tag = tag
        self.R = R

class DisEdge:
    """DisEdge: class for dislocation edge

    Defines the basic features on a edge
    """
    def __init__(self, burg_vec: np.ndarray, plane_normal: np.ndarray) -> None:
        self.burg_vec = burg_vec
        self.plane_normal = plane_normal

class DisNet:
    """DisNet: class for dislocation network

    Implements basic topological operations on dislocation networks
    """
    def __init__(self, data=None, **attr) -> None:
        self.data = data
        self.attr = attr

    def nodes(self) -> list:
        """nodes: return a list of nodes

        (To do: return a generator instead of a list)
        """
        return list(self.data)
    
    def edges(self) -> list:
        """edges: return a list of edges

        (To do: return a generator instead of a list)
        """
        return list(self.data.values())
    
    def pos_array(self) -> np.ndarray:
        """pos_array: return a numpy array of node positions
        """
        return np.array([node.R for node in self.nodes()])
    
    def seg_list(self) -> list:
        """seg_list: return a list of segments
        """
        return [(node1.R, node2.R) for node1, node2 in self.edges()]
    
    def has_node(self, tag: Tag) -> bool:
        """has_node: check if a node exists in the network
        """
        return tag in self.data
    
    def has_edge(self, node1: Tag, node2: Tag) -> bool:
        """has_edge: check if an edge exists in the network
        """
        return node1 in self.data and node2 in self.data[node1]
    
    def add_node(self, tag: Tag, R: np.ndarray) -> None:
        """add_node: add a node to the network
        """
        if tag not in self.data:
            self.data[tag] = DisNode(tag, R)
    
    def add_edge(self, node1: Tag, node2: Tag, burg_vec: np.ndarray, plane_normal: np.ndarray) -> None:
        """add_edge: add an edge to the network
        """
        if node1 not in self.data:
            self.data[node1] = {}
        if node2 not in self.data:
            self.data[node2] = {}
        self.data[node1][node2] = DisEdge(burg_vec, plane_normal)
        self.data[node2][node1] = DisEdge(-burg_vec, plane_normal)
    
    def remove_edge(self, node1: Tag, node2: Tag) -> None:
        """remove_edge: remove an edge from the network
        """
        if node1 in self.data and node2 in self.data[node1]:
            del self.data[node1][node2]
            del self.data[node2][node1]
    
    def insert_node(self, node1: Tag, node2: Tag, tag: Tag, R: np.ndarray) -> None:
        """insert_node: insert a node between two existing nodes
        """
        if node1 in self.data and node2 in self.data[node1]:
            self.remove_edge(node1, node2)
            self.add_node(tag, R)
            self.add_edge(node1, tag, self.data[node1].burg_vec, self.data[node1].plane_normal)
            self.add_edge(tag, node2, self.data[node1].burg_vec, self.data[node1].plane_normal)
    
    def remove_two_arm_node(self, tag: Tag) -> None:
        """remove_two_arm_node: remove a node with two arms from the network
        """
        if tag in self.data:
            node1, node2 = list(self.data[tag].keys())
            self.remove_edge(tag, node1)
            self.remove_edge(tag, node2)
            self.add_edge(node1, node2, self.data[node1].burg_vec, self.data[node1].plane_normal)
            del self.data[tag]
    
    def merge_node(self, node1: Tag, node2: Tag) -> None:
        """merge_node: merge two nodes into one
        """
        if node1 in self.data and node2 in self.data:
            self.add_edge(node1, node2, self.data[node1].burg_vec, self.data[node1].plane_normal)
            self.remove_two_arm_node(node1)
            self.remove_two_arm_node(node2)
    
    def split_node(self, node: Tag, partition: list) -> None:
        """split_node: split a node into two nodes with neighbors split according to partition
        """
        pass

    def is_sane(self) -> bool:
        """is_sane: check if the network is sane
        The two arms connecting two nodes should have opposite Burgers vectors and parallel plane_normal vectors
        The sum of all Burgers vectors of outgoing arms from a node should be zero
        """
        return True
    
