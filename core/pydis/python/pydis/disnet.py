"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
import networkx as nx
from typing import Tuple
from copy import deepcopy
from enum import IntEnum

from base_classes.disnet_base import DisNet_BASE

Tag = Tuple[int, int]

class DisNode:
    """DisNode: class for dislocation node

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

    def __init__(self, R: np.ndarray, F: np.ndarray=None,
                v: np.ndarray=None, flag: int=Flags.CLEAR,
                constraint: int=Constraints.UNCONSTRAINED) -> None:
        self.R = R
        if F is not None:
            self.F = F
        if v is not None:
            self.v = v
        if flag is not None:
            self.flag = int(flag)
        if constraint is not None:
            self.constraint = int(constraint)

class DisEdge:
    """DisEdge: class for dislocation edge

    Defines the basic features on a edge
    """
    def __init__(self, burg_vec: np.ndarray, plane_normal: np.ndarray=None) -> None:
        self.burg_vec = burg_vec
        if plane_normal is not None:
            self.plane_normal = plane_normal

class Cell:
    """Cell: class for simulation cell in which dislocation network is embedded

    Defines periodic or non-periodic boundary conditions in each direction
    """
    def __init__(self, h: np.ndarray=np.eye(3), is_periodic: list=[False,False,False]) -> None:
        self.h = h
        self.hinv = np.linalg.inv(h)
        self.is_periodic = is_periodic

    def map(self, dr: np.ndarray) -> np.ndarray:
        """map: map a vector to the simulation cell
        """
        if not any(self.is_periodic):
            return dr
        ds = np.dot(self.hinv, dr.T).T
        ds -= self.is_periodic * np.round(ds)
        return np.dot(self.h, ds.T).T

    def map_to(self, r2: np.ndarray, r1: np.ndarray) -> np.ndarray:
        """map: map r2 to the nearest image of r1 (if PBC is applied)
        """
        if not any(self.is_periodic):
            return r2
        ds = np.dot(self.hinv, (r2-r1).T).T
        ds -= self.is_periodic * np.round(ds)
        return np.dot(self.h, ds.T).T + r1

class CellList:
    """CellList: class for cell list

    Implements cell list for efficient neighbor search
    Note: only works for PBC in all 3 directions
    """
    def __init__(self, cell: Cell=None, n_div: np.array=[3,3,3]) -> None:
        # To do: comput ndiv from cell and cutoff
        self.cell = Cell() if cell is None else cell
        self.n_div = n_div
        self._cell_list = [[[ [] for n2 in range(self.n_div[2])] for n1 in range(self.n_div[1])] for n0 in range(self.n_div[0])]
        self._cell_indices = None

    def get_cell_index(self, R: np.ndarray) -> np.ndarray:
        """get_cell_index: get cell index of a position
        """
        s = np.dot(self.cell.hinv, R.T).T
        s -= np.round(s)
        ind = np.floor((s+0.5)*self.n_div).astype(int)
        return ind

    def sort_points_to_list(self, R: np.ndarray) -> None:
        """build: build cell list

        Note: only do this when PBC is applied in all 3 directions
        """
        self._cell_list = [[[ [] for n2 in range(self.n_div[2])] for n1 in range(self.n_div[1])] for n0 in range(self.n_div[0])]
        if all(self.cell.is_periodic):
            self._cell_indices = self.get_cell_index(R)
            for i, ind in enumerate(self._cell_indices):
                self._cell_list[ind[0]][ind[1]][ind[2]].append(i)
        else:
            self._cell_indices = [None]*len(R)
        return

    def get_objs_in_cell(self, ind: np.ndarray) -> list:
        """get_indices_in_cell: get indices in a cell
        """
        return self._cell_list[ind[0]][ind[1]][ind[2]]

    def get_objs_in_nbr_cells(self, obj_id: int) -> list:
        """get_indices_in_cell: get indices in the same cell of a point and all neighbor cells
        """
        ind = self._cell_indices[obj_id]
        nbr_ids = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    ind_nbr = np.array([ind[0]+i, ind[1]+j, ind[2]+k])
                    ind_nbr = np.mod(ind_nbr, self.n_div)
                    nbr_ids.extend(self.get_objs_in_cell(ind_nbr))
        return nbr_ids

    def iterate_nbr_pairs(self, use_cell_list: bool=True):
        """iterate_nbr_pairs: iterate over all pairs of segments in the same cell and all neighbor cells

        Note: self pairs (i, i) not included
        """
        if not use_cell_list:
            n = len(self._cell_indices)
            for i in range(n):
                for j in range(i+1, n):
                    yield i, j
        else:
            n = len(self._cell_indices)
            for i in range(n):
                nbr_ids = self.get_objs_in_nbr_cells(i)
                for j in nbr_ids:
                    if i < j: yield i, j

class DisNet(DisNet_BASE):
    """DisNet: class for dislocation network

    Implements basic topological operations on dislocation networks
    """
    def __init__(self, data=None, cell=None, cell_list=None, rn=None, links=None, **attr) -> None:
        self._G = nx.DiGraph(data, **attr)
        self.cell = Cell() if cell is None else cell
        self.cell_list = CellList(cell=self.cell) if cell_list is None else cell_list
        self._recycled_tags = []
        if rn is not None or links is not None:
            self.add_nodes_links_from_list(rn, links)

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
    
    # To do: rename function edges
    def segments(self):
        """segments: return segment view of the network
        """
        return self.edges()

    def segments(self, tag_pair: Tuple[Tag, Tag]):
        """segments: return segment view of the network
        """
        return self.edges(tag_pair)

    @property
    def out_degree(self):
        """out_degree: return out degree of a node
        """
        return self._G.out_degree

    def pos_array(self) -> np.ndarray:
        """pos_array: return a numpy array of node positions
        """
        return np.array([node_attr["R"] for node, node_attr in self._G.nodes.items()])
    
    # To do: implement function node_prop_list
    def node_prop_list(self) -> list:
        """node_prop_list: return a list of node properties
        """
        return [node_attr for node, node_attr in self._G.nodes.items()]

    # To do: implement function seg_prop_list
    def seg_prop_list(self) -> list:
        """seg_prop_list: return a list of segment properties
        """
        return [edge_attr for edge, edge_attr in self._G.edges.items()]

    def seg_list(self) -> list:
        """seg_list: return a list of segments

        Each link appear once: tag1 < tag2
        """
        segments = []
        for edge in self._G.edges:
            tag1 = edge[0]
            tag2 = edge[1]
            if tag1 < tag2:
                r1 = self._G.nodes[tag1]["R"]
                r2 = self._G.nodes[tag2]["R"]
                # apply PBC
                r2 = self.cell.map_to(r2, r1)
                segments.append({"edge":edge,
                                 "burg_vec":self._G.edges[edge]["burg_vec"],
                                 "R1":r1,
                                 "R2":r2})
        return segments
    
    def get_segments_midpoint(self, segments: list) -> np.ndarray:
        """get_segments_midpoint: return the midpoints of a segment list
        """
        mp = np.zeros((len(segments), 3))
        for i, seg in enumerate(segments):
            r1 = seg["R1"]
            r2 = seg["R2"]
            # apply PBC
            r2 = self.cell.map_to(r2, r1)
            mp[i,:] = 0.5*(r1+r2)
        return mp
        
    def nodes_array(self) -> list:
        """nodes_array: pack nodes into an array for export
        Node format: x,y,z,constraint
        """
        nodes = []
        ntags = {}
        for i, tag in enumerate(self._G.nodes()):
            r = self._G.nodes[tag]["R"]
            constraint = self._G.nodes[tag]["constraint"]
            nodes.append([r[0], r[1], r[2], constraint])
            ntags[tag] = i
        return nodes, ntags
    
    def segs_array(self, ntags: dict) -> list:
        """segs_array: pack segments into an array for export
        Seg format: node1,node2,burg,plane
        """
        segments = []
        for edge in self._G.edges():
            n1, n2 = ntags[edge[0]], ntags[edge[1]]
            if n1 < n2:
                b = self._G.edges[edge]["burg_vec"]
                p = self._G.edges[edge]["plane_normal"]
                segments.append([n1, n2, b[0], b[1], b[2], p[0], p[1], p[2]])
        return segments
        
    def export_data(self):
        """export_data: export network to data
        """
        cell = {"h": self.cell.h, "is_periodic": self.cell.is_periodic}
        nodes, ntags = self.nodes_array()
        segs = self.segs_array(ntags)
        data = {"cell" : cell,
                "nodes": nodes,
                "segs" : segs }
        return data
    
    def import_data(self, data):
        """import_data: import network from data
        """
        cell = data.get("cell")
        self.cell = Cell(h=cell.get("h"), is_periodic=cell.get("is_periodic"))
        self._G.clear()
        self._recycled_tags = []
        rn = data.get("nodes")
        segs = data.get("segs")
        self.add_nodes_links_from_list(rn, segs)

    def sort_segments_to_cell_list(self, segments):
        """sort_segments_to_cell: sort segments to cell list
        """
        cell_list = []
        for edge in self._G.edges:
            tag1 = edge[0]
            tag2 = edge[1]
            r1 = self._G.nodes[tag1]["R"]
            r2 = self._G.nodes[tag2]["R"]
            # apply PBC
            r2 = self.cell.map_to(r2, r1)
            self._G.edges[edge]["R1"] = r1
            self._G.edges[edge]["R2"] = r2
        return cell_list

    def has_node(self, tag: Tag) -> bool:
        """has_node: check if a node exists in the network
        """
        return self._G.has_node(tag)
    
    # To do: rename function has_edge
    def has_segment(self, tag1: Tag, tag2: Tag) -> bool:
        return self.has_edge(tag1, tag2)

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

    # To do: rename function add_nodes_links_from_list
    def add_nodes_segments_from_list(self, rn, links) -> None:
        return self.add_nodes_links_from_list(rn, links)

    def add_nodes_links_from_list(self, rn, links) -> None:
        """add_nodes_links_from_list: add nodes and edges stored in lists to network
           sanity after this operation depends on the input
        """
        N = rn.shape[0]
        num_links = links.shape[0]
        for i in range(N):
            if rn.shape[1] > 3:
                node = DisNode(R=rn[i,:3], constraint=int(rn[i,3]))
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

        self._add_node(tag, deepcopy(DisNode(R=R)))
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
        end_nodes_connected = self.has_edge(tag1, tag2)

        link12_attr = self.edges[(tag, tag2)]
        self._combine_edge(tag1, tag2, deepcopy(DisEdge(**link12_attr)))
        link21_attr = self.edges[(tag, tag1)]
        self._combine_edge(tag2, tag1, deepcopy(DisEdge(**link21_attr)))
        self._remove_edge(tag, tag1)
        self._remove_edge(tag1, tag)
        self._remove_edge(tag, tag2)
        self._remove_edge(tag2, tag)
        self._remove_node(tag)

        if end_nodes_connected:
            # cleaning up is needed if the two end nodes were connected
            self.remove_empty_arms(tag1)
            if len(self.edges(tag1)) == 0:
                self._remove_node(tag1)
            if self.has_node(tag2):
                self.remove_empty_arms(tag2)
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

        node1Deletable = self.nodes[tag1]["constraint"] != DisNode.Constraints.PINNED_NODE
        node2Deletable = self.nodes[tag2]["constraint"] != DisNode.Constraints.PINNED_NODE

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

        self.remove_empty_arms(targetNode)

        # Remove target node if orphaned (i.e. no longer has any arms)
        if len(self.edges(targetNode)) == 0:
            self._remove_node(targetNode)
            targetNode = None
            status = 'MERGE_NODE_ORPHANED'
        else:
            mergedTag = targetNode
            status = 'MERGE_NODE_SUCCESS'

        #print("merge_node: %s + %s -> %s" % (str(tag1), str(tag2), str(mergedTag)))
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
        self._add_node(split_node2, deepcopy(DisNode(R=pos2.copy())))

        self.nodes[split_node1]['R'] = pos1.copy()

        bv = np.zeros(3)
        for nbr in nbrs_to_split:
            if not self.has_edge(tag, nbr):
                raise ValueError("split_node: Node %s and %s are not connected" % (str(tag), str(nbr)))

            link_attr = self.edges[(tag, nbr)]
            self._add_edge(split_node2, nbr, deepcopy(DisEdge(**link_attr)))
            bv += np.array(link_attr['burg_vec'])

            link_attr = self.edges[(nbr, tag)]
            self._add_edge(nbr, split_node2, deepcopy(DisEdge(**link_attr)))

            self._remove_edge(tag, nbr)
            self._remove_edge(nbr, tag)

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

            self._add_edge(split_node1, split_node2, deepcopy(DisEdge(burg_vec= bv, plane_normal=pn)))
            self._add_edge(split_node2, split_node1, deepcopy(DisEdge(burg_vec=-bv, plane_normal=pn)))

        #print("split_node: original tag = %s -> new tags = %s, %s, bv = %s" % (str(tag), str(split_node1), str(split_node2)))
        return split_node1, split_node2

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
            if self.nodes[tag]["constraint"] == DisNode.Constraints.PINNED_NODE: continue
            b_tot = np.zeros(3)
            for nbr_tag in self.neighbors(tag):
                b_tot += self.edges[(tag, nbr_tag)]['burg_vec']
            if np.max(np.abs(b_tot)) > tol:
                print("Total Burgers vector from node %s is not zero" % (str(tag)))
                return False

        return True
