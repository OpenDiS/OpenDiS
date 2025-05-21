"""@package docstring
DisNet_Base: base class for DisNet (dislocation network)

Defines interface for the DisNet class
"""

from abc import ABC, abstractmethod 

class DisNet_Base(ABC):
    """DisNet_Base: base class for DisNet (dislocation network)

    Defines interface for the DisNet class
    """

    #@abstractmethod
    # To do: find out if property can be abstract method
    #def cell(self):
    #    """cell: return simulation cell
    #    """
    #    pass

    @abstractmethod
    def export_data(self):
        """export_data: export network to data
        """
        pass
    
    @abstractmethod
    def import_data(self, data):
        """import_data: import network from data
        """
        pass
    
    @abstractmethod
    def num_nodes(self):
        """num_nodes: return total number of nodes
        """
        pass
    
    @abstractmethod
    def num_segments(self):
        """num_segments: return total number of segments
        """
        pass
    
    @abstractmethod
    def is_sane(self):
        """is_sane: check if the network is sane
        """
        pass


class DisNet_Python(DisNet_Base):
    """DisNet_Python: base class for DisNet inherited by PyDiS

    Defines interface for the DisNet class for Python implementation
    """
    @abstractmethod
    def all_nodes_tags(self):
        """nodes: return list of all node tags (keys)
        """
        pass
    
    @abstractmethod
    def nodes(self, tag):
        """nodes: return property of node specified by tag (key)
        """
        pass

    @abstractmethod
    def all_segments_tags(self):
        """segments: return list of all segments (tag pairs)
        """
        pass

    @abstractmethod
    def segments(self, tag_pair):
        """segments: return property of segments specified by tag_pair
        """
        pass
    
    @abstractmethod
    def out_degree(self, tag):
        """out_degree: return out degree of a node specified by tag (key)
        """
        pass

    @abstractmethod
    def has_node(self, tag) -> bool:
        """has_node: check if a node exists in the network
        """
        pass
    
    @abstractmethod
    def has_segment(self, tag1, tag2) -> bool:
        """has_edge: check if an edge exists in the network
        """
        pass
    
    @abstractmethod
    def add_nodes_segments_from_list(self, rn, links) -> None:
        """add_nodes_segments_from_list: add nodes and edges stored in lists to network
        """
        pass
    
    # To do: rename this function to insert_node_between
    @abstractmethod
    def insert_node(self, tag1, tag2, tag):
        """insert_node: insert a node (tag) between two existing nodes (tag1 and tag2)
        """
        pass
    
    @abstractmethod
    def remove_two_arm_node(self, tag):
        """remove_two_arm_node: remove a node with two arms from the network
        """
        pass

    @abstractmethod 
    def merge_node(self, tag1, tag2):
        """merge_node: merge two nodes into one
           return mergedTag (tag1 or tag2) if merge is successful, None otherwise

        Remove any links between node1 and node2
        If merged node has double-links to any neighbor, combine them into one (or zero) link
        """
        pass
    
    @abstractmethod
    def split_node(self, tag, nbrs_to_split):
        """split_node: split a node into two nodes with neighbors split according to partition

        inputs:
           tag: tag of the node to be split
           nbrs_to_split: list of tags of neighbors to be split off to the new node

        outputs:
           split_node1: tag of the node to which all unselected arms are connected
           split_node2: tag of the node to which all selected arms are connected
        """
        pass

    @abstractmethod
    def is_sane(self) -> bool:
        """is_sane: check if the network is sane

        The two arms connecting two nodes should have opposite Burgers vectors and parallel plane_normal vectors
        The sum of all Burgers vectors of outgoing arms from a node should be zero
        """
        pass
    
    # perhaps split_multi_nodes should belong to another (higher level) module
