"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
import networkx as nx
from typing import Tuple
from copy import deepcopy

Tag = Tuple[int, int]

class DisNode:
    """DisNode: class for dislocation node

    Defines the basic attributes on a node
    """
    def __init__(self, R: np.ndarray, F: np.ndarray=None,
                v: np.ndarray=None, flag: np.ndarray=None) -> None:
        self.R = R
        if F is not None:
            self.F = F
        if v is not None:
            self.v = v
        if flag is not None:
            self.flag = flag

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
        self._recycled_tags = []

    def neighbors(self, tag: Tag):
        """neighbors: return neighbors (as iterator) of a node
        """
        return self._G.neighbors(tag)

    @property
    def nodes(self):
        """nodes: return node view of the network
        """
        return self._G.nodes
    
    @property
    def edges(self):
        """edges: return out edge view of the network
        """
        return self._G.edges
    
    def pos_array(self) -> np.ndarray:
        """pos_array: return a numpy array of node positions
        """
        return np.array([node_attr["R"] for node, node_attr in self._G.nodes.items()])
    
    def seg_list(self) -> list:
        """seg_list: return a list of segments

        Each link appear once: tag1 < tag2
        """
        segments = []
        for edge in self._G.edges:
            tag1 = edge[0]
            tag2 = edge[1]
            if tag1 < tag2:
                segments.append({"edge":edge,
                                 "burg_vec":self._G.edges[edge]["burg_vec"],
                                 "R1":self._G.nodes[tag1]["R"],
                                 "R2":self._G.nodes[tag2]["R"]})
        return segments
    
    def has_node(self, tag: Tag) -> bool:
        """has_node: check if a node exists in the network
        """
        return self._G.has_node(tag)
    
    def has_edge(self, tag1: Tag, tag2: Tag) -> bool:
        """has_edge: check if an edge exists in the network
        """
        return self._G.has_edge(tag1, tag2)
    
    def _add_node(self, tag: Tag, node_attr: DisNode) -> None:
        """add_node: add a node to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if self.has_node(tag):
            raise ValueError("_add_node: node %s already exists" % (str(tag)))
        self._G.add_node(tag, **vars(node_attr))
    
    def _add_edge(self, tag1: Tag, tag2: Tag, edge_attr: DisEdge) -> None:
        """add_edge: add an edge to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if self.has_edge(tag1, tag2):
            raise ValueError("_add_edge: Edge (%s, %s) already exists" % (str(tag1), str(tag2)))
        self._G.add_edge(tag1, tag2, **vars(edge_attr))
    
    def _remove_node(self, tag: Tag) -> None:
        """remove_edge: remove a node from the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self._G.remove_node(tag)
        self._recycled_tags.append(tag)

    def _remove_edge(self, tag1: Tag, tag2: Tag) -> None:
        """remove_edge: remove an edge from the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self._G.remove_edge(tag1, tag2)

    def _combine_edge(self, tag1: Tag, tag2: Tag, edge_attr: DisEdge) -> None:
        """combine_edge: combine an edge with an existing edge
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if not self.has_edge(tag1, tag2):
            self._add_edge(tag1, tag2, edge_attr)
            return

        self._G.edges[(tag1, tag2)]["burg_vec"] += edge_attr.burg_vec

        # To do: update glide plane normal

    def add_nodes_links_from_list(self, rn, links) -> None:
        """add_nodes_links_from_list: add nodes and edges stored in lists to network
           sanity after this operation depends on the input
        """
        N = rn.shape[0]
        num_links = links.shape[0]
        for i in range(N):
            if rn.shape[1] > 3:
                node = DisNode(R=rn[i,:3], flag=int(rn[i,3]))
            else:
                node = DisNode(R=rn[i,:3])
            self._add_node((0,i), deepcopy(node))
        for j in range(num_links):
            seg = links[j, :2].astype(int)
            bv  = links[j, 2:5]
            tag = (0,seg[0])
            nbr_tag = (0, seg[1])
            if links.shape[1] > 5:
                pn  = links[j, 5:8]
                self._add_edge(tag, nbr_tag, deepcopy(DisEdge(burg_vec= bv, plane_normal=pn)))
                self._add_edge(nbr_tag, tag, deepcopy(DisEdge(burg_vec=-bv, plane_normal=pn)))
            else:
                self._add_edge(tag, nbr_tag, deepcopy(DisEdge(burg_vec= bv)))
                self._add_edge(nbr_tag, tag, deepcopy(DisEdge(burg_vec=-bv)))

        if not self.is_sane():
            raise ValueError("add_nodes_links_from_list: sanity check failed")
    
    def get_new_tag(self, recycle = True) -> Tag:
        """get_new_tag: return a new tag for a new node
           if recycle == True, then take from list of recycled node tags
           recycle == False makes it easier to debug as node tags are never reused
        """
        if recycle and len(self._recycled_tags) > 0:
            return self._recycled_tags.pop(0)
        else:
            max_tag = max(self.nodes)
            return (max_tag[0], max_tag[1]+1)

    def insert_node(self, tag1: Tag, tag2: Tag, tag: Tag, R: np.ndarray) -> None:
        """insert_node: insert a node between two existing nodes
           guarantees sanity after operation
        """
        # To do: update plastic strain if new node position is not on segment

        self._add_node(tag, deepcopy(DisNode(R=R, flag=0)))
        link12_attr = self.edges[(tag1, tag2)]
        self._add_edge(tag1, tag, deepcopy(DisEdge(**link12_attr)))
        self._add_edge(tag, tag2, deepcopy(DisEdge(**link12_attr)))
        link21_attr = self.edges[(tag2, tag1)]
        self._add_edge(tag2, tag, deepcopy(DisEdge(**link21_attr)))
        self._add_edge(tag, tag1, deepcopy(DisEdge(**link21_attr)))
        self._remove_edge(tag1, tag2)
        self._remove_edge(tag2, tag1)
    
    def remove_two_arm_node(self, tag: Tag) -> None:
        """remove_two_arm_node: remove a node with two arms from the network
           guarantees sanity after operation
        """
        if not self.has_node(tag):
            raise ValueError("remove_two_arm_node: Node %s does not exist" % str(tag))
        if self._G.out_degree(tag) != 2:
            raise ValueError("remove_two_arm_node: Node %s does not have two arms" % str(tag))

        # To do: update plastic strain due to removed node operation

        tag1, tag2 = self.neighbors(tag)
        if not self.has_edge(tag1, tag2):
            link12_attr = self.edges[(tag, tag2)]
            self._add_edge(tag1, tag2, deepcopy(DisEdge(**link12_attr)))
            link21_attr = self.edges[(tag, tag1)]
            self._add_edge(tag2, tag1, deepcopy(DisEdge(**link21_attr)))
            self._remove_edge(tag, tag1)
            self._remove_edge(tag1, tag)
            self._remove_edge(tag, tag2)
            self._remove_edge(tag2, tag)
            self._remove_node(tag)
        else:
            # To do: refactor this section
            link12_attr = self.edges[(tag, tag2)]
            self._combine_edge(tag1, tag2, deepcopy(DisEdge(**link12_attr)))
            link21_attr = self.edges[(tag, tag1)]
            self._combine_edge(tag2, tag1, deepcopy(DisEdge(**link21_attr)))
            self._remove_edge(tag, tag1)
            self._remove_edge(tag1, tag)
            self._remove_edge(tag, tag2)
            self._remove_edge(tag2, tag)
            self._remove_node(tag)

            self.remove_empty_arms(tag1)
            self.remove_empty_arms(tag2)
            if len(self.edges(tag1)) == 0:
                self._remove_node(tag1)
            if len(self.edges(tag2)) == 0:
                self._remove_node(tag2)

    def remove_empty_arms(self, tag: Tag) -> None:
        """remove_empty_arms: remove any zero-Burgers vector arms between a node and any of its neighbors
           guarantees sanity after operation

        Remove neighbor nodes if they become orphaned, but do not remove the node itself
        """
        nbr_list = list(self.neighbors(tag))
        for nbr in nbr_list:
            bv = self.edges[(tag, nbr)]["burg_vec"]
            if np.max(np.abs(bv)) < 1e-8:
                self._remove_edge(tag, nbr)
                self._remove_edge(nbr, tag)

        for nbr in nbr_list:
            if len(self.edges(nbr)) == 0:
                self._remove_node(nbr)
    
    def merge_node(self, tag1: Tag, tag2: Tag):
        """merge_node: merge two nodes into one
           guarantees sanity after operation
           return mergedTag (tag1 or tag2) if merge is successful, None otherwise
                  and status (MERGE_NOT_PERMITTED, MERGE_NODE_ORPHANED or MERGE_NODE_SUCCESS)

        Remove any links between node1 and node2
        If merged node has double-links to any neighbor, combine them into one (or zero) link
        """

        node1Deletable = self.nodes[tag1]["flag"] != 7
        node2Deletable = self.nodes[tag2]["flag"] != 7

        if node1Deletable:
            targetNode, deadNode = tag2, tag1
        elif node2Deletable:
            targetNode, deadNode = tag1, tag2
        else:
            mergedTag = None
            status = 'MERGE_NOT_PERMITTED'
            return mergedTag, status

        # To do: update plastic strain due to merged node operation

        # Remove any links between targetNode and deadNode
        if self.has_edge(targetNode, deadNode):
            self._remove_edge(targetNode, deadNode)
        if self.has_edge(deadNode, targetNode):
            self._remove_edge(deadNode, targetNode)

        # Move all connections from the dead node to the target node
        # and add a new connection from the target node to each of the
        # dead node's neighbors.
        for edge in self.edges(deadNode):
            nbr = edge[1]
            if nbr != targetNode:
                link_attr = self.edges[(deadNode, nbr)]
                self._combine_edge(targetNode, nbr, deepcopy(DisEdge(**link_attr)))

                link_attr = self.edges[(nbr, deadNode)]
                self._combine_edge(nbr, targetNode, deepcopy(DisEdge(**link_attr)))

            # To do: reset seg forces

        self._remove_node(deadNode)

        #if not self.is_sane():
        #    raise ValueError("sanity check failed 2")

        self.remove_empty_arms(targetNode)

        if not self.is_sane():
            raise ValueError("sanity check failed 3")

        # Remove target node if orphaned (i.e. no longer has any arms)
        if len(self.edges(targetNode)) == 0:
            self._remove_node(targetNode)
            targetNode = None
            status = 'MERGE_NODE_ORPHANED'
        else:
            mergedTag = targetNode
            status = 'MERGE_NODE_SUCCESS'

        return mergedTag, status
    
    def split_node(self, tag: Tag, partition: list) -> None:
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
        for tag in self.nodes:
            for nbr_tag in self.neighbors(tag):
                b12 = self.edges[(tag, nbr_tag)]['burg_vec']
                b21 = self.edges[(nbr_tag, tag)]['burg_vec']
                if np.max(np.abs(b12 + b21)) > tol:
                    print("Burgers vector of edge (%s, %s) is not opposite to that of edge (%s, %s)" % (str(tag), str(nbr_tag), str(nbr_tag), str(tag)))
                    return False

        for tag in self.nodes:
            b_tot = np.zeros(3)
            for nbr_tag in self.neighbors(tag):
                b_tot += self.edges[(tag, nbr_tag)]['burg_vec']
            if np.max(np.abs(b_tot)) > tol:
                print("Total Burgers vector from node %s is not zero" % (str(tag)))
                return False

        return True
    
