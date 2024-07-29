"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
from .graph.graph import Graph, Node as Node_bare, Edge as Edge_bare
from typing import Tuple
from enum import IntEnum

from base_classes.disnet_base import DisNet_BASE

Tag = Tuple[int, int]

class DisNode():
    """DisNode: class for dislocation node (properties)

    Defines the basic attributes on a node
    """

    class Constraints(IntEnum):
        """Constraints: enum class for node constraints
        """
        UNCONSTRAINED = 0
        SURFACE_NODE  = 1
        THINFILM_SURFACE_NODE = 5 # To do: resolve this difference from ParaDiS
        CYLINDER_SURFACE_NODE = 6
        PINNED_NODE    = 7

    class Flags(IntEnum):
        """Flags: enum class for node flags
        """
        CLEAR              = 0x00
        NODE_RESET_FORCES  = 0x01
        NODE_OUTSIDE_SURF  = 0x02
        NO_COLLISIONS      = 0x04
        NO_MESH_COARSEN    = 0x08
        NODE_CHK_DBL_LINK  = 0x10

    def __init__(self, R: np.ndarray,
                 constraint: int=Constraints.UNCONSTRAINED) -> None:
        self.R = R
        self.constraint = int(constraint)

    def is_equivalent(self, other) -> bool:
        # check if all attributes in self have the same value in other
        for key, value in vars(self).items():
            if type(value) == np.ndarray:
                if not np.all(value == vars(other)[key]):
                    return False
            else:
                if value != vars(other)[key]:
                    return False
        return True

    def copy(self):
        """copy: return a deep copy of the node attr
        """
        return DisNode(R=self.R.copy(), constraint=self.constraint)

    def view(self):
        """view: return a dictionary view of the node attr
        """
        return vars(self)

class DisEdge():
    """DisEdge: class for dislocation edge (properties)

    Defines the basic features on a edge
    """
    def __init__(self, source_tag: Tag, burg_vec: np.ndarray, plane_normal: np.ndarray=None) -> None:
        self.source_tag = source_tag
        self.burg_vec = burg_vec
        if plane_normal is not None:
            self.plane_normal = plane_normal

    def burg_vec_from(self, from_tag: Tag) -> np.ndarray:
        return self.burg_vec if self.source_tag == from_tag else -self.burg_vec

    def add(self, other) -> None:
        self.burg_vec += other.burg_vec_from(self.source_tag)
        # To do: update glide plane normal

    def is_equivalent(self, other) -> bool:
        # check if burg_vec and plane_normal in self have the equivalent values in other
        burg_vec_eq = np.all(self.burg_vec == other.burg_vec_from(self.source_tag))
        plane_normal_eq = np.all(self.plane_normal == other.plane_normal) or np.all(self.plane_normal == -other.plane_normal)
        return burg_vec_eq and plane_normal_eq

    def copy(self):
        """copy: return a deep copy of the edge attr
        """
        return DisEdge(source_tag=self.source_tag, burg_vec=self.burg_vec.copy(),
                       plane_normal=self.plane_normal.copy() if hasattr(self, "plane_normal") else None)

    def view(self):
        """view: return a dictionary view of the edge attr
        """
        return vars(self)

class Cell:
    """Cell: class for simulation cell in which dislocation network is embedded

    Defines periodic or non-periodic boundary conditions in each direction
    """
    def __init__(self, h: np.ndarray=np.eye(3), origin: np.ndarray=np.zeros(3), is_periodic: list=[False,False,False]) -> None:
        self.h = h
        self.hinv = np.linalg.inv(h)
        self.origin = origin # "left-most" corner of the cell
        self.is_periodic = is_periodic

    def map(self, dr: np.ndarray) -> np.ndarray:
        """map: map a vector to the simulation cell
        """
        if not any(self.is_periodic):
            return dr
        ds = np.dot(self.hinv, dr.T).T
        ds -= self.is_periodic * np.round(ds)
        return np.dot(self.h, ds.T).T

    def closest_image(self, Rref: np.ndarray, R: np.ndarray) -> np.ndarray:
        """map: map R to the nearest image of Rref (if PBC is applied)
        """
        if not any(self.is_periodic):
            return R
        ds = np.dot(self.hinv, (R-Rref).T).T
        ds -= self.is_periodic * np.round(ds)
        return np.dot(self.h, ds.T).T + Rref

    def center(self) -> np.ndarray:
        """center: return the center of the cell
        """
        return self.origin + 0.5*np.dot(self.h, np.ones(3))

    def copy(self):
        """copy: return a deep copy of the cell
        """
        return Cell(h=self.h.copy(), origin=self.origin.copy(), is_periodic=self.is_periodic.copy())

    def view(self):
        """view: return a dictionary view of the cell
        """
        return vars(self)

class DisNet(DisNet_BASE):
    """DisNet: class for dislocation network

    Implements basic topological operations on dislocation networks
    """
    class Node_with_attr(Node_bare):
        def __init__(self, tag: Tag, node_attr: DisNode):
            super().__init__()
            self.tag = tag
            self.attr = node_attr

    class Edge_with_attr(Edge_bare):
        def __init__(self, source: Node_bare, target: Node_bare, edge_attr: DisEdge):
            super().__init__(source, target)
            self.attr = edge_attr

    def __init__(self, cell=None, rn=None, links=None) -> None:
        self._G = Graph()
        # provide a reference from tags back to nodes (with attr)
        self.tags_to_nodes = {}
        self.cell = Cell() if cell is None else cell
        self._recycled_tags = []
        if rn is not None or links is not None:
            self.add_nodes_segments_from_list(rn, links)

    def clear_graph(self):
        """clear: clear the network
        """
        # To do: uncomment _G.clear after implemented in Graph
        self._G.clear()
        self.tags_to_nodes.clear()
        self._recycled_tags.clear()

    def neighbors_tags(self, tag: Tag):
        """neighbors: return neighbor tags (as iterator) of a node
                      (not required in base class)
        """
        connected_edges = self.tags_to_nodes[tag].edges()

        return (edge.source.tag if edge.source.tag != tag else edge.target.tag
                for edge in connected_edges)

    def neighbors_dict(self, tag: Tag):
        """neighbors_dict: return neighbor tags and attributes of a node
        """
        connected_edges = self.tags_to_nodes[tag].edges()

        result = {}
        for edge in connected_edges:
            other_node = edge.source if edge.source.tag != tag else edge.target
            result[other_node.tag] = other_node.attr

        return result

    def neighbor_segments_dict(self, tag: Tag):
        """neighbor_segments_dict: return nbr_tags and attributes of edges
        """
        connected_edges = self.tags_to_nodes[tag].edges()
        return {
           edge.source.tag if edge.source.tag != tag else edge.target.tag : edge.attr
           for edge in connected_edges
        }

    def all_nodes_tags(self):
        """all_nodes_tags: return iterator of all node tags
        """
        return self.tags_to_nodes.keys()

    def all_nodes_mapping(self):
        """all_nodes_mapping: return iterator of all node tags and attributes
        """
        return ( (tag, node.attr) for tag, node in self.tags_to_nodes.items() )

    def all_nodes_dict(self):
        """nodes: return dictionary of all node tags and attributes
        """
        return { tag: node.attr for tag, node in self.tags_to_nodes.items() }
    
    def nodes(self, tag: Tag):
        """nodes: return property of node specified by tag (key)
        """
        return self.tags_to_nodes[tag].attr

    def num_nodes(self):
        """num_nodes: return total number of nodes in the network
        """
        return len(self._G._nodes)

    def num_segments(self):
        """num_segments: return total number of segments in the network
        """
        return len(self._G._edges)

    def all_segments_tags(self):
        """segments: return iterator of all segments (tag pairs)
        """
        return ( (edge.source.tag, edge.target.tag) for edge in self._G.edges() )
    
    def all_segments_mapping(self):
        """all_segments_mapping: return iterator of all segments (tag pairs) and attributes
        """
        return ( ((edge.source.tag, edge.target.tag), edge.attr) for edge in self._G.edges() )

    def all_segments_dict(self):
        """all_segments_dict: return dictionary of all segments (tag pairs) -> attributes
        """
        return { (edge.source.tag, edge.target.tag): edge.attr for edge in self._G.edges() }

    def segments(self, tag_pair: Tuple[Tag, Tag]):
        """segments: return property of segments specified by tag_pair
        """
        node_source = self.tags_to_nodes[tag_pair[0]]
        node_target = self.tags_to_nodes[tag_pair[1]]
        return self._G.edge_between(node_source, node_target).attr

    def out_degree(self, tag):
        """out_degree: return out degree of a node
        """
        return self.tags_to_nodes[tag].num_neighbors()

    def pos_array(self) -> np.ndarray:
        """pos_array: return a numpy array of node positions
        """
        return np.array([node.attr.R for node in self.tags_to_nodes.values()])

    # To do: remove function node_prop_list (after removed from base class)
    def node_prop_list(self) -> list:
        """node_prop_list: return a list of node properties
        """
        raise NotImplementedError("node_prop_list: not implemented")

    # To do: remove function seg_prop_list (after removed from base class)
    def seg_prop_list(self) -> list:
        """seg_prop_list: return a list of segment properties
        """
        raise NotImplementedError("seg_prop_list: not implemented")

    def get_nodes_data(self):
        """get_nodes_data: collect nodes data into a dictionary format
           Returns:
                nodes_data: dictionary of nodes data
                ntags: position of each tag in the nodes array
        """
        Nnode = self.num_nodes()
        tags = np.zeros((Nnode, 2), dtype=int)
        R = np.zeros((Nnode, 3))
        constraints = np.zeros((Nnode, 1), dtype=int)
        ntags = {}
        i = 0
        for tag, node_attr in self.all_nodes_mapping():
            tags[i,:] = tag
            R[i,:] = node_attr.R
            constraints[i,:] = node_attr.constraint
            ntags[tag] = i
            i += 1
        return {
            "tags": tags,
            "positions": R,
            "constraints": constraints
        }, ntags

    def get_segs_data(self, ntags: dict):
        """get_segs_data: collect segments data into a dictionary format
        """
        Nseg = self.num_segments()
        nodeids = np.zeros((Nseg, 2), dtype=int)
        burgers = np.zeros((Nseg, 3))
        planes = np.zeros((Nseg, 3))
        i = 0
        for (source, target), edge_attr in self.all_segments_dict().items():
            nodeids[i,:] = ntags[source], ntags[target]
            burgers[i,:] = edge_attr.burg_vec_from(source)
            planes[i,:] = edge_attr.plane_normal
            i += 1
        return {
            "nodeids": nodeids,
            "burgers": burgers,
            "planes": planes
        }

    def get_segs_data_with_positions(self):
        """get_segs_data_with_positions: collect segments data into a dictionary format
        """
        _, ntags = self.get_nodes_data()
        Nseg = self.num_segments()
        nodeids = np.zeros((Nseg, 2), dtype=int)
        tag1 = np.zeros((Nseg, 2), dtype=int)
        tag2 = np.zeros((Nseg, 2), dtype=int)
        burgers = np.zeros((Nseg, 3))
        planes = np.zeros((Nseg, 3))
        R1 = np.zeros((Nseg, 3))
        R2 = np.zeros((Nseg, 3))
        i = 0
        for (source, target), edge_attr in self.all_segments_dict().items():
            nodeids[i,:] = ntags[source], ntags[target]
            tag1[i,:] = source
            tag2[i,:] = target
            burgers[i,:] = edge_attr.burg_vec_from(source).copy()
            planes[i,:] = getattr(edge_attr, "plane_normal", np.zeros(3)).copy()
            r1_local = self.nodes(source).R
            r2_local = self.nodes(target).R
            # apply PBC
            r2_local = self.cell.closest_image(Rref=r1_local, R=r2_local)
            R1[i,:] = r1_local
            R2[i,:] = r2_local
            i += 1
        return {
            "nodeids": nodeids,
            "tag1": tag1,
            "tag2": tag2,
            "burgers": burgers,
            "planes":  planes,
            "R1": R1,
            "R2": R2
        }

    def export_data(self):
        """export_data: export network to data
        """
        cell = {"h": self.cell.h, "origin": self.cell.origin, "is_periodic": self.cell.is_periodic}
        nodes_data, ntags = self.get_nodes_data()
        segs_data = self.get_segs_data(ntags)
        data = {"cell": cell, "nodes": nodes_data, "segs": segs_data}
        return data
    
    def import_data(self, data):
        """import_data: import network from data
        """
        cell = data.get("cell")
        self.cell = Cell(h=cell.get("h"), origin=cell.get("origin"), is_periodic=cell.get("is_periodic"))
        self.clear_graph()
        nodes_data = data.get("nodes")
        nodes_array = np.hstack((nodes_data["tags"], nodes_data["positions"], nodes_data["constraints"]))
        segs_data = data.get("segs")
        segs_array = np.hstack((segs_data["nodeids"], segs_data["burgers"], segs_data["planes"]))
        self.add_nodes_segments_from_list(nodes_array, segs_array)

    def copy(self):
        """copy: return a deep copy of the network
        """
        result = DisNet(cell=self.cell.copy())
        for tag, node_attr in self.all_nodes_mapping():
            result._add_node(tag, node_attr.copy())

        for (source, target), edge_attr in self.all_segments_dict().items():
            result._add_edge(source, target, edge_attr.copy())

        return result

    def to_networkx(self):
        """to_networkx: convert to networkx graph
        """
        # To do: implement function to_networkx
        import networkx as nx
        nx_graph = nx.DiGraph()
        for tag, node in self.all_nodes_mapping():
            nx_graph.add_node(tag, **vars(node))

        for (source, target), edge_attr in self.all_segments_dict().items():
            edge_dict = vars(edge_attr.copy())
            del edge_dict["source_tag"]
            edge_dict["burg_vec"] = edge_attr.burg_vec_from(source)
            nx_graph.add_edge(source, target, **edge_dict)

            edge_dict = vars(edge_attr.copy())
            del edge_dict["source_tag"]
            edge_dict["burg_vec"] = edge_attr.burg_vec_from(target)
            nx_graph.add_edge(target, source, **edge_dict)

        return nx_graph

    def from_networkx(self, nx_graph):
        """from_networkx: convert from networkx graph
                nx_graph: networkx.DiGraph
        """
        import networkx as nx
        if type(nx_graph) != nx.DiGraph:
            raise ValueError("from_networkx: input must be networkx.DiGraph")
        self.clear_graph()
        for tag, node_attr in nx_graph.nodes(data=True):
            self._add_node(tag, DisNode(**node_attr))

        for source, target, edge_attr in nx_graph.edges(data=True):
            if source < target:
                self._add_edge(source, target, DisEdge(source, edge_attr.get("burg_vec").copy(), edge_attr.get("plane_normal", None).copy()))

    def is_equivalent(self, G_compare):
        for tag, node in self.all_nodes_mapping():
            if not node.is_equivalent(G_compare.nodes(tag)):
                return False
        for (source, target), edge_attr in self.all_segments_dict().items():
            if not edge_attr.is_equivalent(G_compare.segments((source, target))):
                return False
        return True

    def has_node(self, tag: Tag) -> bool:
        """has_node: check if a node exists in the network
        """
        return tag in self.tags_to_nodes

    def has_segment(self, tag1: Tag, tag2: Tag) -> bool:
        """has_segment: check if an edge exists in the network
        """
        node1 = self.tags_to_nodes.get(tag1, None)
        node2 = self.tags_to_nodes.get(tag2, None)

        if node1 is None or node2 is None:
            return False

        return self._G.has_edge_between(node1, node2)

    def _add_node(self, tag: Tag, node_attr: DisNode) -> None:
        """add_node: add a node to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if self.has_node(tag):
            raise ValueError("_add_node: node %s already exists" % (str(tag)))
        node = self.Node_with_attr(tag, node_attr)
        self.tags_to_nodes[tag] = node
        self._G.add_node(node)
    
    def _add_edge(self, tag1: Tag, tag2: Tag, edge_attr: DisEdge) -> None:
        """add_edge: add an edge to the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if self.has_segment(tag1, tag2):
            raise ValueError("_add_edge: Edge (%s, %s) already exists" % (str(tag1), str(tag2)))

        node1 = self.tags_to_nodes[tag1]
        node2 = self.tags_to_nodes[tag2]

        edge = self.Edge_with_attr(node1, node2, edge_attr)
        self._G.add_edge(edge)

    def _remove_node(self, tag: Tag) -> None:
        """remove_edge: remove a node from the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        node = self.tags_to_nodes.pop(tag)
        self._G.remove_node(node)
        self._recycled_tags.append(tag)

    def _remove_edge(self, tag1: Tag, tag2: Tag) -> None:
        """remove_edge: remove an edge from the network
           user is not supposed to call this low level function, does not guarantee sanity
        """
        node1 = self.tags_to_nodes[tag1]
        node2 = self.tags_to_nodes[tag2]
        edge = self._G.edge_between(node1, node2)
        self._G.remove_edge(edge)

    def _combine_edge(self, tag1: Tag, tag2: Tag, edge_attr: DisEdge) -> None:
        """combine_edge: combine an edge with an existing edge
           user is not supposed to call this low level function, does not guarantee sanity
        """
        if not self.has_segment(tag1, tag2):
            self._add_edge(tag1, tag2, edge_attr)
            return

        node1 = self.tags_to_nodes[tag1]
        node2 = self.tags_to_nodes[tag2]
        self._G.edge_between(node1, node2).attr.add(edge_attr)

    def add_nodes_segments_from_list(self, rn, links) -> None:
        """add_nodes_segments_from_list: add nodes and edges stored in lists to network
           sanity after this operation depends on the input
        """
        N = rn.shape[0]
        num_links = links.shape[0]
        for i in range(N):
            if rn.shape[1] == 6:
                tag = (int(rn[i,0]), int(rn[i,1]))
                node = DisNode(R=rn[i,2:5].copy(), constraint=int(rn[i,5]))
            elif rn.shape[1] == 4:
                tag = (0, i)
                node = DisNode(R=rn[i,:3].copy(), constraint=int(rn[i,3]))
            elif rn.shape[1] == 3:
                tag = (0, i)
                node = DisNode(R=rn[i,:3].copy())
            else:
                raise ValueError("add_nodes_segments_from_list: invalid node format")
            self._add_node(tag, node)
        for j in range(num_links):
            seg = links[j, :2].astype(int)
            bv  = links[j, 2:5]
            if rn.shape[1] == 6:
                tag = (int(rn[int(seg[0]),0]), int(rn[int(seg[0]),1]))
                nbr_tag = (int(rn[int(seg[1]),0]), int(rn[int(seg[1]),1]))
            elif rn.shape[1] == 4 or rn.shape[1] == 3:
                tag = (0,seg[0])
                nbr_tag = (0, seg[1])
            if links.shape[1] > 5:
                pn  = links[j, 5:8]
                # Note: now we add edges only once
                self._add_edge(tag, nbr_tag, DisEdge(tag, burg_vec=bv.copy(), plane_normal=pn.copy()))
            else:
                # Note: now we add edges only once
                self._add_edge(tag, nbr_tag, DisEdge(tag, burg_vec=bv.copy()))

        if not self.is_sane():
            raise ValueError("add_nodes_segments_from_list: sanity check failed")
    
    def get_new_tag(self, recycle = True) -> Tag:
        """get_new_tag: return a new tag for a new node
           if recycle == True, then take from list of recycled node tags
           recycle == False makes it easier to debug as node tags are never reused
        """
        if recycle and len(self._recycled_tags) > 0:
            return self._recycled_tags.pop(0)
        else:
            max_tag = max(self.all_nodes_tags())
            return (max_tag[0], max_tag[1]+1)

    def insert_node(self, tag1: Tag, tag2: Tag, new_tag: Tag, R: np.ndarray) -> None:
        insert_node_between(tag1, tag2, new_tag, R)

    def insert_node_between(self, tag1: Tag, tag2: Tag, new_tag: Tag, R: np.ndarray) -> None:
        """insert_node_between: insert a node between two existing nodes
           guarantees sanity after operation
        """
        # To do: update plastic strain if new node position is not on segment
        self._add_node(new_tag, DisNode(R=R.copy()))
        prev_edge_attr = self.segments((tag1, tag2))
        new_edge_attr = DisEdge(tag1, prev_edge_attr.burg_vec_from(tag1).copy(), prev_edge_attr.plane_normal.copy())
        self._add_edge(tag1, new_tag, new_edge_attr)
        new_edge_attr = DisEdge(tag2, prev_edge_attr.burg_vec_from(tag2).copy(), prev_edge_attr.plane_normal.copy())
        self._add_edge(new_tag, tag2, new_edge_attr)
        self._remove_edge(tag1, tag2)
    
    def remove_two_arm_node(self, old_tag: Tag) -> None:
        """remove_two_arm_node: remove a node with two arms from the network
           guarantees sanity after operation
        """
        if not self.has_node(old_tag):
            raise ValueError("remove_two_arm_node: Node %s does not exist" % str(old_tag))
        if self.out_degree(old_tag) != 2:
            raise ValueError("remove_two_arm_node: Node %s does not have two arms" % str(old_tag))

        # To do: update plastic strain due to removed node operation

        tag1, tag2 = self.neighbors_tags(old_tag)
        end_nodes_connected = self.has_segment(tag1, tag2)

        prev_link_attr = self.segments((old_tag, tag2))
        new_link_attr = DisEdge(tag2, prev_link_attr.burg_vec_from(tag2).copy(), prev_link_attr.plane_normal.copy())
        self._combine_edge(tag1, tag2, new_link_attr)
        self._remove_node(old_tag)

        # remove neighbor nodes if they become orphaned
        self.remove_empty_arms(tag1)
        self.remove_empty_arms(tag2)

    def remove_empty_arms(self, tag: Tag, tol = 1e-8) -> None:
        """remove_empty_arms: remove any zero-Burgers vector arms between a node and any of its neighbors
           guarantees sanity after operation

        Remove neighbor nodes if they become orphaned, but do not remove the node itself
        """
        if not self.has_node(tag):
            return

        nbr_list = list(self.neighbors_tags(tag))
        node = self.tags_to_nodes[tag]
        edges_to_remove = []
        for edge in node.edges():
            bv = edge.attr.burg_vec_from(tag)
            if np.max(np.abs(bv)) < tol:
                edges_to_remove.append(edge)
        for edge in edges_to_remove:
            self._G.remove_edge(edge)

        if node.num_neighbors() == 0:
            self._remove_node(node.tag)

        for nbr_tag in nbr_list:
            if self.out_degree(nbr_tag) == 0:
                self._remove_node(nbr_tag)

    def merge_node(self, tag1: Tag, tag2: Tag):
        """merge_node: merge two nodes into one
           guarantees sanity after operation
           return mergedTag (tag1 or tag2) if merge is successful, None otherwise
                  and status (MERGE_NOT_PERMITTED, MERGE_NODE_ORPHANED or MERGE_NODE_SUCCESS)

        Remove any links between node1 and node2
        If merged node has double-links to any neighbor, combine them into one (or zero) link
        """
        node1Deletable = self.nodes(tag1).constraint != DisNode.Constraints.PINNED_NODE
        node2Deletable = self.nodes(tag2).constraint != DisNode.Constraints.PINNED_NODE

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
        if self.has_segment(targetNode, deadNode):
            self._remove_edge(targetNode, deadNode)

        # Move all connections from the dead node to the target node
        # and add a new connection from the target node to each of the
        # dead node's neighbors.
        for nbr_tag, link_attr in self.neighbor_segments_dict(deadNode).items():
            new_link_attr = DisEdge(nbr_tag, link_attr.burg_vec_from(nbr_tag).copy(), link_attr.plane_normal.copy())
            self._combine_edge(targetNode, nbr_tag, new_link_attr)

            # To do: reset seg forces

        self._remove_node(deadNode)

        self.remove_empty_arms(targetNode)

        if self.has_node(targetNode) and self.out_degree(targetNode) > 0:
            mergedTag = targetNode
            status = 'MERGE_NODE_SUCCESS'
        else:
            mergedTag = None
            status = 'MERGE_NODE_ORPHANED'

        return mergedTag, status
    
    def find_precise_glide_plane(self, bv: np.ndarray, dirv: np.ndarray, dot_cutoff=0.9995) -> np.ndarray:
        """find_precise_glide_plane: find glide plane normal given burgers vector and line direction
        """
        # following ParaDiS FindPreciseGlidePlane.c
        if np.dot(dirv, dirv) < 1e-8:
            return np.zeros(3)

        if np.square(np.dot(bv, dirv)) / (np.dot(bv, bv)*np.dot(dirv, dirv)) > dot_cutoff:
            return np.zeros(3)

        pn = np.cross(bv, dirv)
        pn /= np.linalg.norm(pn)

        # To do: implement geometries based on FCC or BCC slip systems

        return pn

    def split_node(self, tag: Tag, pos1: np.ndarray, pos2: np.ndarray, nbrs_to_split: list, eps_b = 1e-3) -> (Tag, Tag):
        """split_node: split a node into two nodes with neighbors split according to partition
           guarantees sanity after operation

        inputs:
           tag: tag of the node to be split
           nbrs_to_split: list of tags of neighbors to be split off to the new node

        outputs:
           split_node1: tag of the node to which all unselected arms are connected
           split_node2: tag of the node to which all selected arms are connected
        """
        split_node1 = tag
        split_node2 = self.get_new_tag()
        self._add_node(split_node2, DisNode(R=pos2.copy()))

        self.nodes(split_node1).R = pos1.copy()

        bv = np.zeros(3)
        for nbr in nbrs_to_split:
            if not self.has_segment(tag, nbr):
                raise ValueError("split_node: Node %s and %s are not connected" % (str(tag), str(nbr)))

            link_attr = self.segments((tag, nbr))
            new_link_attr = DisEdge(nbr, link_attr.burg_vec_from(nbr).copy(), link_attr.plane_normal.copy())
            self._add_edge(split_node2, nbr, new_link_attr)
            bv += link_attr.burg_vec_from(tag)

            self._remove_edge(tag, nbr)

        # ParadiS calls AssignNodeToCell here

        # possibly adding a link between split_node1 and split_node2
        if np.inner(bv, bv) > eps_b:
            dirv = pos2 - pos1
            # To do: apply PBC

            # To do: move two nodes apart along their velocity

            if np.inner(dirv, dirv) < eps_b:
                dirv = np.array([0.0, 0.0, 0.0])
            else:
                dirv = dirv / np.linalg.norm(dirv)
            pn = self.find_precise_glide_plane(bv, dirv)

            self._add_edge(split_node1, split_node2, DisEdge(split_node1, burg_vec=bv.copy(), plane_normal=pn.copy()))

        #print("split_node: original tag = %s -> new tags = %s, %s, bv = %s" % (str(tag), str(split_node1), str(split_node2)))
        return split_node1, split_node2

    def is_sane(self, tol: float=1e-8) -> bool:
        """is_sane: check if the network is sane
           guarantees sanity after operation

        The two arms connecting two nodes should have opposite Burgers vectors and parallel plane_normal vectors
        The sum of all Burgers vectors of outgoing arms from a node should be zero
        """
        for (source, target), edge_attr in self.all_segments_dict().items():
            if edge_attr.source_tag != source and edge_attr.source_tag != target:
                print("source_tag ", edge_attr.source_tag, " is not in link", (source, target))
                return False

        for tag, node in self.tags_to_nodes.items():
            if node.num_neighbors() == 0:
                print("Node %s has no neighbors" % (str(tag)))
                return False
            if node.attr.constraint == DisNode.Constraints.PINNED_NODE: continue
            b_tot = np.zeros(3)
            for edge in node.edges():
                b_tot += edge.attr.burg_vec_from(tag)
            if np.max(np.abs(b_tot)) > tol:
                print(f"Total Burgers vector from node {tag} is not zero {b_tot}")
                return False

        return True

    @staticmethod
    def convert_nodeforce_dict_to_array(state: dict) -> dict:
        nodeforce_dict = state["nodeforce_dict"]
        state["nodeforces"] = np.array([f for f in nodeforce_dict.values()])
        state["nodeforcetags"] = np.array([ [domainID, index] for domainID, index in nodeforce_dict.keys()])
        return state

    @staticmethod
    def convert_nodeforce_array_to_dict(state: dict) -> dict:
        nodeforces = state["nodeforces"]
        nodeforcetags = state["nodeforcetags"]
        nodeforce_dict = {(nodeforcetags[i][0],nodeforcetags[i][1]): nodeforces[i] for i in range(len(nodeforcetags))}
        state["nodeforce_dict"] = nodeforce_dict
        return state
