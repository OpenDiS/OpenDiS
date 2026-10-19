"""@package docstring
DisNet: class for dislocation network

Implements basic topological operations on dislocation networks
"""

import numpy as np
from .graph.graph import Graph, Node as Node_bare, Edge as Edge_bare
from typing import Tuple
from enum import IntEnum

from framework.disnet_base import DisNet_Python

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
    def __init__(self, source_tag: Tag, target_tag: Tag, burg_vec: np.ndarray, plane_normal: np.ndarray=None) -> None:
        self.source_tag = source_tag
        self.target_tag = target_tag
        self.burg_vec = burg_vec
        if plane_normal is not None:
            self.plane_normal = plane_normal

    def burg_vec_from(self, from_tag: Tag) -> np.ndarray:
        if from_tag != self.source_tag and from_tag != self.target_tag:
            raise ValueError("burg_vec_from: from_tag not in edge")
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
        return DisEdge(source_tag=self.source_tag, target_tag=self.target_tag,
                       burg_vec=self.burg_vec.copy(),
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

    def is_orthogonal(self) -> bool:
        """is_orthogonal: are the cell vectors axis-aligned

        True when every off-diagonal entry of h is zero, which is the
        negation of exadis Cell::is_triclinic in core/exadis/src/network.h.
        Only an orthogonal cell has a single box edge per component, which
        is what the 'edge' branch of closest_image needs.
        """
        off = self.h - np.diag(np.diag(self.h))
        return not np.any(off)

    def _closest_image_scaled(self, Rref: np.ndarray, R: np.ndarray) -> np.ndarray:
        """_closest_image_scaled: via scaled coordinates, for any cell shape

        Transforms to scaled coordinates and back unconditionally, so R is
        rebuilt from hinv @ (R - Rref) even when the image index is zero and
        does not round-trip exactly.
        """
        ds = np.dot(self.hinv, (R-Rref).T).T
        ds -= self.is_periodic * np.round(ds)
        return np.dot(self.h, ds.T).T + Rref

    def _closest_image_edge(self, Rref: np.ndarray, R: np.ndarray) -> np.ndarray:
        """_closest_image_edge: subtract whole box edges, orthogonal cells only

        Keeps R and removes an exact multiple of the box edge per component,
        so a point already in the nearest image comes back untouched. Caller
        checks is_orthogonal first; there is no per-component edge otherwise.
        """
        out = np.array(R, dtype=float)
        ref = np.asarray(Rref, dtype=float)
        for k in range(3):
            if self.is_periodic[k]:
                edge = self.h[k, k]
                out[..., k] -= np.rint((out[..., k] - ref[..., k])/edge)*edge
        return out

    def closest_image(self, Rref: np.ndarray, R: np.ndarray,
                      method: str='edge') -> np.ndarray:
        """closest_image: map R to the nearest periodic image of Rref

        Same convention as exadis' Cell::pbc_position
        (core/exadis/src/network.h) either way; method selects the
        arithmetic, which changes the last bits but never which image is
        chosen.

        'scaled'
            Always goes out to scaled coordinates and back, rebuilding R
            from hinv @ (R - Rref) even when the image index is zero. The
            long-standing pydis form, kept because every stored reference
            up to 2026-08-24 was produced with it.

        'edge'
            For an orthogonal cell, keeps R and subtracts an exact multiple
            of the box edge one component at a time,

                rpbc.x -= rint((rpbc.x - r0.x)/Lbox.x) * Lbox.x

            so a point already in the nearest image is returned untouched.
            This is exadis' orthogonal branch. For a triclinic cell exadis
            has no such edge either and does the scaled round trip, so
            'edge' falls back to 'scaled' there and the two codes already
            agree.

        The two differ on about 4% of calls in a 1000-wide box, by up to one
        ulp at the coordinate magnitude, 5.7e-14. Where a folded point sits
        near a box face the loss is an ulp of the box edge instead, 1.1e-13
        at 1000. Measured over tests/unit_tests/test4_collision_mode, moving
        this rule to 'edge' takes the worst per-step disagreement with exadis
        from 4.5e-13 to 1.2e-13.

        'edge' is the more accurate of the two and is the default. It was
        made the default on 2026-08-24 after measuring that it moves no
        stored reference far enough to fail: all eight test folders pass
        with it, including the zero-tolerance comparisons that
        tests/unit_tests/test1_node_force asserts on a bitrepro build.

        The np.dot calls in the scaled form are a further wrinkle for
        reproducibility rather than for agreement: a 3x3 matvec goes to BLAS,
        so the result depends on which BLAS numpy was built against. The
        'edge' form has no matrix product and so no such dependence. See
        pydis.calforce.bitrepro_math.
        """
        if not any(self.is_periodic):
            return R
        if method == 'edge' and self.is_orthogonal():
            return self._closest_image_edge(Rref, R)
        if method not in ('scaled', 'edge'):
            raise ValueError("closest_image: method must be 'scaled' or "
                             "'edge', got %r" % (method,))
        return self._closest_image_scaled(Rref, R)

    def fold(self, R: np.ndarray) -> np.ndarray:
        """fold: map a position into the primary cell (if PBC is applied)

        The image of R nearest the cell centre. Corresponds to FoldBox in
        external/ParaDiS.git/src/Util.c.
        """
        return self.closest_image(Rref=self.center(), R=R)

    def center(self) -> np.ndarray:
        """center: return the center of the cell
        """
        return self.origin + 0.5*np.dot(self.h, np.ones(3))

    def volume(self) -> float:
        """volume: return the volume of the cell
        """
        return np.abs(np.linalg.det(self.h))

    def copy(self):
        """copy: return a deep copy of the cell
        """
        return Cell(h=self.h.copy(), origin=self.origin.copy(), is_periodic=self.is_periodic.copy())

    def view(self):
        """view: return a dictionary view of the cell
        """
        return vars(self)

class DisNet(DisNet_Python):
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

    def __init__(self, cell=None, rn=None, links=None,
                 recycle_tags=True, recycle_LIFO=True) -> None:
        self._G = Graph()
        # provide a reference from tags back to nodes (with attr)
        self.tags_to_nodes = {}
        self.cell = Cell() if cell is None else cell
        self._recycled_tags = []
        # highest tag ever handed out, so an index is not reissued just
        # because the node holding it was removed
        self._max_tag = (0, -1)
        # how get_new_tag reuses the tags of removed nodes
        self.recycle_tags = recycle_tags
        self.recycle_LIFO = recycle_LIFO
        if rn is not None or links is not None:
            self.add_nodes_segments_from_list(rn, links)

    def clear_graph(self):
        """clear: clear the network
        """
        # To do: uncomment _G.clear after implemented in Graph
        self._G.clear()
        self.tags_to_nodes.clear()
        self._recycled_tags.clear()
        self._max_tag = (0, -1)

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
                self._add_edge(source, target, DisEdge(source, target, edge_attr.get("burg_vec").copy(), edge_attr.get("plane_normal", None).copy()))

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
        self._max_tag = max(self._max_tag, tag)
    
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

    def segment_as_stored(self, tag1: Tag, tag2: Tag) -> Tuple[Tag, Tag]:
        """segment_as_stored: this segment's two endpoints, in the stored order

        Which way round a segment is stored is no physics, every reader going
        through DisEdge.burg_vec_from, but rules break ties on it: exadis'
        retroactive collision merges into a segment's FIRST endpoint whenever
        the collision point is within rann of it, so on a segment shorter than
        2*rann the stored direction alone picks the survivor. The operations
        below therefore put a surviving or new node into the endpoint slot the
        node it replaces held, and callers naming a segment from a node's point
        of view come through here first.
        """
        node1 = self.tags_to_nodes[tag1]
        node2 = self.tags_to_nodes[tag2]
        edge = self._G.edge_between(node1, node2)
        if edge is None:
            raise ValueError("segment_as_stored: no segment %s-%s"
                             % (str(tag1), str(tag2)))
        return edge.source.tag, edge.target.tag

    def _zero_edge(self, tag1: Tag, tag2: Tag) -> None:
        """zero_edge: set an edge's Burgers vector to zero, keeping the edge
           user is not supposed to call this low level function, does not guarantee sanity
        """
        self.segments((tag1, tag2)).burg_vec[:] = 0.0

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
                self._add_edge(tag, nbr_tag, DisEdge(tag, nbr_tag, burg_vec=bv.copy(), plane_normal=pn.copy()))
            else:
                # Note: now we add edges only once
                self._add_edge(tag, nbr_tag, DisEdge(tag, nbr_tag, burg_vec=bv.copy()))

        if not self.is_sane():
            raise ValueError("add_nodes_segments_from_list: sanity check failed")
    
    def get_new_tag(self, recycle = None, LIFO = None) -> Tag:
        """get_new_tag: return a new tag for a new node
           recycle and LIFO default to the network settings recycle_tags
           and recycle_LIFO, and may be overridden per call
           if recycle == True, then take from list of recycled node tags
           recycle == False makes it easier to debug as node tags are never reused
           LIFO == True takes the tag freed most recently, matching exadis
        """
        if recycle is None: recycle = self.recycle_tags
        if LIFO is None: LIFO = self.recycle_LIFO
        if recycle and len(self._recycled_tags) > 0:
            return self._recycled_tags.pop() if LIFO \
                   else self._recycled_tags.pop(0)
        else:
            # RULE CHANGE (LengthBased): counts up from the highest tag ever
            # issued, not from the highest still in use, so removing the
            # top node does not put its index back in circulation. Matches
            # exadis SerialDisNet::get_new_tag, which counts from maxindex.
            return (self._max_tag[0], self._max_tag[1]+1)

    def insert_node(self, tag1: Tag, tag2: Tag, new_tag: Tag, R: np.ndarray) -> None:
        insert_node_between(tag1, tag2, new_tag, R)

    def move_node(self, tag: Tag, pos: np.ndarray, plast_acc=None) -> None:
        """move_node: place a node at pos, accounting for the area its arms sweep

        Corresponds to exadis SerialDisNet::move_node
        (core/exadis/src/network.cpp), except that pos is stored as given
        rather than folded, as the assignments it replaces did. plast_acc,
        a pydis.util.plastic_strain.PlasticStrainAccum, collects the swept
        area; None skips the accounting.
        """
        pos = np.array(pos, dtype=float)
        if plast_acc is not None:
            plast_acc.add_moves(self, {tag: (self.nodes(tag).R, pos)})
        self.nodes(tag).R = pos

    def insert_node_between(self, tag1: Tag, tag2: Tag, new_tag: Tag, R: np.ndarray,
                            plast_acc=None) -> None:
        """insert_node_between: insert a node between two existing nodes
           guarantees sanity after operation

        With plast_acc, the area swept by bending the segment out to R is
        accounted, as a node created at tag1 and moved to R. It is zero when
        R lies on the segment.
        """
        self._add_node(new_tag, DisNode(R=R.copy()))
        prev_edge_attr = self.segments((tag1, tag2))
        new_edge_attr = DisEdge(tag1, new_tag, prev_edge_attr.burg_vec_from(tag1).copy(), prev_edge_attr.plane_normal.copy())
        self._add_edge(tag1, new_tag, new_edge_attr)
        new_edge_attr = DisEdge(tag2, new_tag, prev_edge_attr.burg_vec_from(tag2).copy(), prev_edge_attr.plane_normal.copy())
        self._add_edge(new_tag, tag2, new_edge_attr)
        self._remove_edge(tag1, tag2)
        if plast_acc is not None:
            plast_acc.add_moves(self, {new_tag: (self.nodes(tag1).R, R)})
    
    def remove_two_arm_node(self, old_tag: Tag, survivor: Tag=None,
                            position: np.ndarray=None, plast_acc=None) -> None:
        """remove_two_arm_node: remove a node with two arms from the network
           guarantees sanity after operation

        survivor names the neighbour that takes the removed node's place and
        so inherits the endpoint slot old_tag held, as exadis' merge_nodes does
        (core/exadis/src/network.cpp). Left None, the joining segment's
        direction comes from old_tag's arm order instead. See
        segment_as_stored.

        position, which needs a survivor, is where the survivor ends up, as
        in exadis' merge_nodes_position. With plast_acc the area swept is
        accounted as old_tag and the survivor moving there together, or as
        old_tag collapsing onto the survivor, or onto its first neighbour,
        when no position is given.
        """
        if not self.has_node(old_tag):
            raise ValueError("remove_two_arm_node: Node %s does not exist" % str(old_tag))
        if self.out_degree(old_tag) != 2:
            raise ValueError("remove_two_arm_node: Node %s does not have two arms" % str(old_tag))
        if position is not None and survivor is None:
            raise ValueError("remove_two_arm_node: position needs a survivor")

        tag1, tag2 = self.neighbors_tags(old_tag)
        if plast_acc is not None:
            plast_acc.add_moves(self, self._collapse_moves(
                old_tag, survivor if survivor is not None else tag1, position))
        end_nodes_connected = self.has_segment(tag1, tag2)

        if survivor is None:
            prev_link_attr = self.segments((old_tag, tag2))
            new_link_attr = DisEdge(tag2, tag1, prev_link_attr.burg_vec_from(tag2).copy(), prev_link_attr.plane_normal.copy())
            self._combine_edge(tag1, tag2, new_link_attr)
            self._remove_node(old_tag)
            return

        if survivor not in (tag1, tag2):
            raise ValueError("remove_two_arm_node: survivor %s is not a "
                             "neighbour of %s" % (str(survivor), str(old_tag)))
        other = tag1 if survivor == tag2 else tag2
        arm = self.segments((old_tag, other))
        old_first = next(first for nbr, _, first
                         in self._arms_with_direction(old_tag) if nbr == other)
        # the survivor stands where old_tag stood, inheriting its slot and its
        # view of this arm, so the Burgers vector is read from whichever node
        # becomes first
        first, second = (survivor, other) if old_first else (other, survivor)
        b_from = arm.burg_vec_from(old_tag if first == survivor else other)
        new_link_attr = DisEdge(first, second, b_from.copy(),
                                arm.plane_normal.copy())
        self._combine_edge(first, second, new_link_attr)
        self._remove_node(old_tag)

        # remove neighbor nodes if they become orphaned
        self.remove_empty_arms(tag1)
        self.remove_empty_arms(tag2)
        # the removal can orphan the survivor and take it with it
        if position is not None and self.has_node(survivor):
            self.nodes(survivor).R = np.array(position, dtype=float)

    def _collapse_moves(self, dead: Tag, target: Tag, position=None) -> dict:
        """_collapse_moves: {tag: (old, new)} for dead merging into target

        Both end at position when it is given; otherwise target stays and
        dead moves onto it. The form PlasticStrainAccum.add_moves takes.
        """
        end = self.nodes(target).R if position is None else position
        moves = {dead: (self.nodes(dead).R, end)}
        if position is not None:
            moves[target] = (self.nodes(target).R, position)
        return moves

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

    def seg_vector(self, tag1: Tag, tag2: Tag) -> np.ndarray:
        """seg_vector: vector from tag1 to tag2, under the minimum image
        """
        R1 = self.nodes(tag1).R
        return self.cell.closest_image(Rref=R1, R=self.nodes(tag2).R) - R1

    def seg_length(self, tag1: Tag, tag2: Tag) -> float:
        """seg_length: distance between two connected nodes, under the minimum image
        """
        return float(np.linalg.norm(self.seg_vector(tag1, tag2)))

    def arm_vectors(self, tag: Tag) -> list:
        """arm_vectors: one entry per arm of a node

        Each entry is (nbr_tag, vec, burg_vec, plane_normal), with vec the
        minimum-image vector from this node to that neighbour and burg_vec
        oriented as leaving this node. plane_normal is None when the edge
        carries none.

        Collects in one place what several callers otherwise open-code:
        the neighbour walk, the periodic reduction, and getting the
        Burgers vector the right way round.
        """
        arms = []
        for nbr_tag, edge in self.neighbor_segments_dict(tag).items():
            arms.append((nbr_tag,
                         self.seg_vector(tag, nbr_tag),
                         edge.burg_vec_from(tag).copy(),
                         getattr(edge, "plane_normal", None)))
        return arms

    def merge_node(self, tag1: Tag, tag2: Tag, position: np.ndarray=None,
                   ignore_constraints: bool=False, defer_purge: bool=False,
                   plast_acc=None):
        """merge_node: merge two nodes into one
           guarantees sanity after operation
           return mergedTag (tag1 or tag2) if merge is successful, None otherwise
                  and status (MERGE_NOT_PERMITTED, MERGE_NODE_ORPHANED or MERGE_NODE_SUCCESS)

        Remove any links between node1 and node2
        If merged node has double-links to any neighbor, combine them into one (or zero) link

        position, when given, is where the surviving node is placed. Without
        it the survivor stays where it was, which is the previous and still
        the default behaviour. Corresponds to exadis SerialDisNet::merge_nodes
        (no position) and merge_nodes_position (with one), in
        core/exadis/src/network.cpp.

        ignore_constraints merges tag1 into tag2 even when both are pinned,
        which otherwise returns MERGE_NOT_PERMITTED. The caller then owns
        the decision, having already chosen which of the two survives.

        defer_purge keeps the dead node and any annihilated arm in the
        network instead of deleting them, leaving both for a later
        purge_network. This is exadis behaviour; ParaDiS deletes eagerly
        here, in MergeNode and RemoveDoubleLinks. The two give different
        results, because a pass that runs before the purge sees the
        annihilated arms as ordinary connectivity. Collision selects it
        with its purge_at_end argument.

        plast_acc accounts the area swept as both nodes reach the merged
        position, measured on the arms before they are combined. It is kept
        when the merge orphans the node: the arms did sweep that area.
        """
        node1Deletable = self.nodes(tag1).constraint != DisNode.Constraints.PINNED_NODE
        node2Deletable = self.nodes(tag2).constraint != DisNode.Constraints.PINNED_NODE

        if ignore_constraints:
            targetNode, deadNode = tag2, tag1
        elif node1Deletable:
            targetNode, deadNode = tag2, tag1
        elif node2Deletable:
            targetNode, deadNode = tag1, tag2
        else:
            mergedTag = None
            status = 'MERGE_NOT_PERMITTED'
            return mergedTag, status

        if plast_acc is not None:
            plast_acc.add_moves(self, self._collapse_moves(deadNode, targetNode, position))

        self._merge_arms(deadNode, targetNode, defer_purge)

        if defer_purge:
            # the dead node and the annihilated arms stay until purge_network
            mergedTag = targetNode
            status = 'MERGE_NODE_SUCCESS'
            if position is not None:
                self.nodes(mergedTag).R = np.array(position, dtype=float)
            return mergedTag, status

        self._remove_node(deadNode)

        self.remove_empty_arms(targetNode)

        if self.has_node(targetNode) and self.out_degree(targetNode) > 0:
            mergedTag = targetNode
            status = 'MERGE_NODE_SUCCESS'
            if position is not None:
                self.nodes(mergedTag).R = np.array(position, dtype=float)
        else:
            mergedTag = None
            status = 'MERGE_NODE_ORPHANED'

        return mergedTag, status
    
    def _merge_arms(self, deadNode: Tag, targetNode: Tag, defer_purge: bool=False) -> None:
        """merge_arms: move the dead node's arms onto the target node
           user is not supposed to call this low level function, does not guarantee sanity

        With defer_purge the dead node keeps every arm it cannot hand over,
        zeroed rather than deleted, which is what exadis
        SerialDisNet::merge_nodes_position does. See merge_node.
        """
        if self.has_segment(targetNode, deadNode):
            if defer_purge:
                self._zero_edge(targetNode, deadNode)
            else:
                self._remove_edge(targetNode, deadNode)

        for nbr_tag, link_attr, dead_first in list(self._arms_with_direction(deadNode)):
            if nbr_tag == targetNode or nbr_tag == deadNode:
                continue
            shared = self.has_segment(targetNode, nbr_tag)
            new_link_attr = DisEdge(nbr_tag, targetNode, link_attr.burg_vec_from(nbr_tag).copy(), link_attr.plane_normal.copy())
            # The transferred arm keeps its stored direction, the survivor
            # taking the slot the dead node held, as exadis
            # merge_nodes_position does. Rebuilding it as survivor -> neighbour
            # reverses every arm that pointed into the dead node. See
            # segment_as_stored.
            if dead_first:
                self._combine_edge(targetNode, nbr_tag, new_link_attr)
            else:
                self._combine_edge(nbr_tag, targetNode, new_link_attr)

            # To do: reset seg forces

            if not defer_purge:
                continue
            if shared:
                # the two arms have been summed onto the target's own arm,
                # so this one is annihilated but stays until purge_network
                self._zero_edge(deadNode, nbr_tag)
            else:
                self._remove_edge(deadNode, nbr_tag)

    def _arms_with_direction(self, tag: Tag):
        """_arms_with_direction: (nbr_tag, attr, tag_is_first) for each arm

        tag_is_first says whether this node is the segment's stored first
        endpoint, which neighbor_segments_dict does not carry.
        """
        for edge in self.tags_to_nodes[tag].edges():
            tag_is_first = edge.source.tag == tag
            nbr_tag = edge.target.tag if tag_is_first else edge.source.tag
            yield nbr_tag, edge.attr, tag_is_first

    def purge_network(self, tol: float=1e-8) -> None:
        """purge_network: drop annihilated segments and the nodes they orphan
           guarantees sanity after operation

        Corresponds to exadis SerialDisNet::purge_network in
        core/exadis/src/network.cpp. Only needed after merges made with
        defer_purge; eager merges leave nothing to collect.
        """
        for (tag1, tag2), attr in list(self.all_segments_mapping()):
            if np.max(np.abs(attr.burg_vec)) < tol:
                self._remove_edge(tag1, tag2)

        for tag in list(self.all_nodes_tags()):
            if self.out_degree(tag) == 0:
                self._remove_node(tag)

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

    def split_node(self, tag: Tag, pos1: np.ndarray, pos2: np.ndarray, nbrs_to_split: list, eps_b = 1e-3,
                   plast_acc=None) -> (Tag, Tag):
        """split_node: split a node into two nodes with neighbors split according to partition
           guarantees sanity after operation

        inputs:
           tag: tag of the node to be split
           nbrs_to_split: list of tags of neighbors to be split off to the new node
           plast_acc: collects the area swept as both products leave tag's
              position, including any segment created between them, which
              exadis' execute_split counts as zero

        outputs:
           split_node1: tag of the node to which all unselected arms are connected
           split_node2: tag of the node to which all selected arms are connected
        """
        pos_old = self.nodes(tag).R.copy()
        split_node1 = tag
        split_node2 = self.get_new_tag()
        self._add_node(split_node2, DisNode(R=pos2.copy()))

        self.nodes(split_node1).R = pos1.copy()

        bv = np.zeros(3)
        # read before the loop starts removing arms
        tag_first = {nbr: first for nbr, _, first in self._arms_with_direction(tag)}
        for nbr in nbrs_to_split:
            if not self.has_segment(tag, nbr):
                raise ValueError("split_node: Node %s and %s are not connected" % (str(tag), str(nbr)))

            link_attr = self.segments((tag, nbr))
            # A transferred arm keeps its stored direction, the new node
            # taking the slot the split node held, as exadis split_node does.
            # See segment_as_stored.
            if tag_first[nbr]:
                first, second, b_source = split_node2, nbr, tag
            else:
                first, second, b_source = nbr, split_node2, nbr
            new_link_attr = DisEdge(first, second,
                                    link_attr.burg_vec_from(b_source).copy(),
                                    link_attr.plane_normal.copy())
            self._add_edge(first, second, new_link_attr)
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

            self._add_edge(split_node1, split_node2, DisEdge(split_node1, split_node2, burg_vec=bv.copy(), plane_normal=pn.copy()))

        if plast_acc is not None:
            plast_acc.add_moves(self, {split_node1: (pos_old, pos1),
                                       split_node2: (pos_old, pos2)})
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

    @staticmethod
    def convert_nodevel_dict_to_array(state: dict) -> dict:
        nodevel_dict = state["vel_dict"]
        state["nodevels"] = np.array([v for v in nodevel_dict.values()])
        state["nodeveltags"] = np.array([ [domainID, index] for domainID, index in nodevel_dict.keys()])
        return state

    @staticmethod
    def convert_nodevel_array_to_dict(state: dict) -> dict:
        nodevels = state["nodevels"]
        nodeveltags = state["nodeveltags"]
        nodevel_dict = {(nodeveltags[i][0],nodeveltags[i][1]): nodevels[i] for i in range(len(nodeveltags))}
        state["vel_dict"] = nodevel_dict
        return state

    # Reading, writing and generating networks, as ExaDisNet does
    # (core/exadis/python/pyexadis_base.py), so that a script can build the same network with
    # either code: the readers and writers are ports of exadis' C++ ones, and the generators of
    # the Python functions of core/exadis/python/pyexadis_utils.py, with the same random numbers.

    def _set_network(self, cell, nodes, segs) -> None:
        """_set_network: replace the network by one of nodes [x,y,z,constraint] and segments
           [node1,node2,burg,plane], as ExaDisNet(cell, nodes, segs) does, without a sanity check
        """
        self.clear_graph()
        self.cell = cell
        for i, node in enumerate(nodes):
            self._add_node((0, i), DisNode(R=np.array(node[:3], dtype=float),
                                           constraint=int(node[3])))
        for seg in segs:
            tag, nbr_tag = (0, int(seg[0])), (0, int(seg[1]))
            self._add_edge(tag, nbr_tag, DisEdge(tag, nbr_tag, burg_vec=np.array(seg[2:5], dtype=float),
                                                 plane_normal=np.array(seg[5:8], dtype=float)))

    def read_paradis(self, datafile: str, verbose: bool=True):
        """read_paradis: read the network of a ParaDiS data file, and return self

        A port of exadis' read_paradis (core/exadis/src/generate.cpp). The cell spans the file's
        minCoordinates to maxCoordinates, periodic in all three directions. The nodes keep the
        file's order and constraints, with the tags (0, i). Each segment is the arm of its
        lower-numbered node, in the file's order. As in exadis, a Burgers vector that is not
        conserved at a node gives a warning, not an error.
        """
        if verbose: print("Reading ParaDiS configuration")
        with open(datafile, encoding="utf-8") as f:
            lines = f.read().splitlines()
        if verbose: print(" Input file: %s" % datafile)

        bounds, node_count, start = {}, -1, None
        for i, line in enumerate(lines):
            for name in ("minCoordinates", "maxCoordinates"):
                if line.startswith(name + " = ["):
                    words = " ".join(lines[i+1:i+5]).replace("]", " ").split()
                    bounds[name] = np.array([float(w) for w in words[:3]])
            if line.startswith("nodeCount"):
                node_count = int(line.split("=")[1])
            if line.startswith("nodalData ="):
                start = i + 1
                break
        if len(bounds) < 2 or start is None or node_count < 0:
            raise ValueError("read_paradis: invalid ParaDiS file %s" % datafile)

        words = " ".join(line for line in lines[start:] if not line.startswith("#"))
        words = words.replace(",", " ").split()
        index, nodes, arms, k = {}, [], [], 0
        while k + 7 <= len(words):
            tag = (int(words[k]), int(words[k+1]))
            pos = [float(w) for w in words[k+2:k+5]]
            narms, flag = int(words[k+5]), int(words[k+6])
            k += 7
            if tag not in index:
                index[tag] = len(nodes)
                nodes.append(pos + [flag])
            for _ in range(narms):
                nbr = (int(words[k]), int(words[k+1]))
                burg = np.array([float(w) for w in words[k+2:k+5]])
                plane = np.array([float(w) for w in words[k+5:k+8]])
                arms.append((tag, nbr, burg, plane))
                k += 8
        if node_count != len(nodes):
            raise ValueError("read_paradis: node inconsistency found in ParaDiS file")

        segs, first = [], {}
        for i, (tag, nbr, burg, plane) in enumerate(arms):
            if tag not in index or nbr not in index:
                raise ValueError("read_paradis: invalid segment found in ParaDiS file")
            first.setdefault((index[tag], index[nbr]), i)
        for tag, nbr, burg, plane in arms:
            n1, n2 = index[tag], index[nbr]
            back = first.get((n2, n1))
            if back is None or np.sum((burg + arms[back][2])**2) > 1e-5:
                raise ValueError("read_paradis: segment inconsistency found in ParaDiS file")
            if n1 < n2:
                segs.append([n1, n2, *burg, *plane])

        cell = Cell(h=np.diag(bounds["maxCoordinates"] - bounds["minCoordinates"]),
                    origin=bounds["minCoordinates"], is_periodic=[True, True, True])
        self._set_network(cell, nodes, segs)
        if verbose: print(" nodes: %d, segments: %d" % (self.num_nodes(), self.num_segments()))

        bsum = np.zeros((len(nodes), 3))
        for n1, n2, *rest in segs:
            bsum[int(n1)] += rest[:3]
            bsum[int(n2)] -= rest[:3]
        nb = int(np.sum(np.sum(bsum**2, axis=1) > 1e-5))
        if nb > 0: print(" Warning: Burgers vector is not conserved for %d node(s)" % nb)
        elif verbose: print(" Burgers vector is conserved for all nodes")
        return self
    read_data = read_paradis

    def write_data(self, datafile: str, verbose: bool=True) -> None:
        """write_data: write the network in ParaDiS format, version 4

        A port of exadis' SerialDisNet::write_data (core/exadis/src/network.cpp), which gives
        the same file for the same network: the nodes folded into the cell, each with its arms
        in the order of the segments, the Burgers vector taken from the node.
        """
        if verbose:
            print("Writing configuration in legacy data format")
            print(" Output file: %s" % datafile)
        if not self.cell.is_orthogonal():
            raise ValueError("write_data: volume must be orthorombic")
        data = self.export_data()
        positions, constraints = data["nodes"]["positions"], data["nodes"]["constraints"]
        nodeids, burgers, planes = (data["segs"]["nodeids"], data["segs"]["burgers"],
                                    data["segs"]["planes"])
        arms = [[] for _ in range(len(positions))]
        for s, (n1, n2) in enumerate(nodeids):
            arms[n1].append((n2, burgers[s], planes[s]))
            arms[n2].append((n1, -burgers[s], planes[s]))
        bmin = np.asarray(self.cell.origin, dtype=float)
        bmax = bmin + np.diag(self.cell.h)
        domain = 0
        with open(datafile, "w", encoding="utf-8") as f:
            f.write("dataFileVersion = %d\n" % 4)
            f.write("numFileSegments = %d\n" % 1)
            f.write("minCoordinates = [\n %f\n %f\n %f\n ]\n" % tuple(bmin))
            f.write("maxCoordinates = [\n %f\n %f\n %f\n ]\n" % tuple(bmax))
            f.write("nodeCount = %d\n" % len(positions))
            f.write("dataDecompType = 2\n")
            f.write("dataDecompGeometry = [\n 1\n 1\n 1\n ]\n\n")
            f.write("#\n#  END OF DATA FILE PARAMETERS\n#\n\n")
            f.write("domainDecomposition = \n")
            f.write("# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds\n")
            f.write("  %d  %f  %f  %f  %f  %f  %f\n" % (domain, *bmin, *bmax))
            f.write("nodalData =\n")
            f.write("# Primary lines: node_tag, x, y, z, num_arms, constraint\n")
            f.write("# Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n")
            for i, pos in enumerate(positions):
                pos = self.cell.fold(np.array(pos))
                f.write("%d, %4d %16.12f %16.12f %16.12f %4d %4d\n"
                        % (domain, i, *pos, len(arms[i]), int(constraints[i][0])))
                for k, b, p in arms[i]:
                    f.write("%10d, %4d %16.12f %16.12f %16.12f\n" % (domain, k, *b))
                    f.write("%33.12f %16.12f %16.12f\n" % tuple(p))

    def write_vtk(self, vtkfile: str, segprops: dict={}, pbc_wrap: bool=True) -> None:
        """write_vtk: write the network in vtk format

        A port of write_vtk in core/exadis/python/pyexadis_utils.py: the cell's corners, then
        each segment from its first node, its far end the nearest image; with pbc_wrap, a
        segment that leaves the cell is cut at its faces.
        """
        data = self.export_data()
        cell = self.cell
        cell_origin, cell_center, h = (np.array(cell.origin, dtype=float), np.array(cell.center()),
                                       np.array(cell.h, dtype=float))
        c = cell_origin + np.array([np.zeros(3), h[0], h[1], h[2], h[0]+h[1],
                                    h[0]+h[2], h[1]+h[2], h[0]+h[1]+h[2]])
        rn = data["nodes"]["positions"]
        segsnid = data["segs"]["nodeids"]
        r1 = np.array(cell.closest_image(Rref=cell_center, R=rn[segsnid[:,0]]))
        r2 = np.array(cell.closest_image(Rref=r1, R=rn[segsnid[:,1]]))
        b = data["segs"]["burgers"]
        p = data["segs"]["planes"]
        segprops = dict(segprops)

        is_periodic = np.array(cell.is_periodic)
        if np.all(is_periodic == 0): pbc_wrap = False
        if pbc_wrap:
            eps = 1e-10
            hinv = np.linalg.inv(h)

            def outside_box(p):
                s = np.matmul(hinv, p - cell_origin)
                return np.any(((s < -eps)|(s > 1.0+eps))&(is_periodic))

            def facet_intersection_position(r1, r2):
                s1 = np.matmul(hinv, r1 - cell_origin)
                s2 = np.matmul(hinv, r2 - cell_origin)
                t = s2 - s1
                t0 = -(s1 - 0.0) / (t + eps)
                t1 = -(s1 - 1.0) / (t + eps)
                s = np.hstack((t0, t1))
                s[s < eps] = 1.0
                facet = np.argmin(s)
                if s[facet] < 1.0:
                    pos = np.matmul(h, s1 + s[facet]*t) + cell_origin
                else:
                    facet = -1
                    pos = r2
                return pos, facet

            segsid, rsegs = [], []
            for i in range(segsnid.shape[0]):
                n1, n2 = segsnid[i]
                r1 = np.array(cell.closest_image(Rref=cell_center, R=rn[n1]))
                r2 = np.array(cell.closest_image(Rref=r1, R=rn[n2]))
                out = outside_box(r2)
                while out:
                    pos, facet = facet_intersection_position(r1, r2)
                    if facet < 0: break
                    segsid.append(i)
                    rsegs.append([r1, pos])
                    r1 = pos + (1-2*np.floor(facet/3)) * (1.0-2*eps)*h[:,facet%3]
                    r2 = np.array(cell.closest_image(Rref=r1, R=rn[n2]))
                    out = outside_box(r2)
                segsid.append(i)
                rsegs.append([r1, r2])
            segsid = np.array(segsid).astype(int)
            nsegs = segsid.shape[0]
            rsegs = np.array(rsegs).reshape(-1,3)
            b = b[segsid]
            p = p[segsid]
            for k, v in segprops.items():
                segprops[k] = v[segsid]
        else:
            nsegs = segsnid.shape[0]
            rsegs = np.hstack((r1, r2)).reshape(-1,3)

        with open(vtkfile, 'w', encoding="utf-8") as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Configuration exported from OpenDiS\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS %d FLOAT\n" % (c.shape[0]+2*nsegs))
            np.savetxt(f, c, fmt='%f')
            np.savetxt(f, rsegs, fmt='%f')
            f.write("CELLS %d %d\n" % (1+nsegs, 9+3*nsegs))
            f.write("8 0 1 4 2 3 5 7 6\n")
            nid = np.hstack((2*np.ones((nsegs,1)), np.arange(2*nsegs).reshape(-1,2)+8))
            np.savetxt(f, nid, fmt='%d')
            f.write("CELL_TYPES %d\n" % (1+nsegs))
            f.write("12\n")
            np.savetxt(f, 4*np.ones(nsegs), fmt='%d')
            f.write("CELL_DATA %d\n" % (1+nsegs))
            f.write("VECTORS Burgers FLOAT\n")
            f.write("%f %f %f\n" % tuple(np.zeros(3)))
            np.savetxt(f, b, fmt='%f')
            f.write("VECTORS Planes FLOAT\n")
            f.write("%f %f %f\n" % tuple(np.zeros(3)))
            np.savetxt(f, p, fmt='%f')
            for k, v in segprops.items():
                vals = np.atleast_2d(v.T).T
                if vals.shape[0] != nsegs:
                    raise ValueError('segprop value must the same size as the number of segments')
                f.write("SCALARS %s FLOAT %d\n" % (str(k), vals.shape[1]))
                f.write("LOOKUP_TABLE default\n")
                np.savetxt(f, np.vstack((np.zeros(vals.shape[1]), vals)), fmt='%f')

    def generate_line_config(self, crystal, Lbox, num_lines, theta=None, maxseg=-1, Rorient=None,
                             seed=-1, verbose=True):
        """generate_line_config: make the network a set of straight, infinite lines, and return self

        A port of generate_line_config in core/exadis/python/pyexadis_utils.py, with the same
        random numbers: for the same arguments it gives the network that ExaDisNet's gives.
        The lines cycle through the 12 slip systems of the crystal, 'fcc' or 'bcc', with the
        sign of their Burgers vectors alternating every 12 lines (dipoles). theta is the list of
        character angles to pick from (degrees); None picks, for each slip system, the angle
        that makes the lines' lengths closest across the systems. Lbox is the box's edge or a
        Cell; seed > 0 seeds numpy's random numbers.
        """
        if verbose: print('generate_line_config()')

        if crystal in ['BCC', 'bcc']:
            # the 12 <111>{110} slip systems
            b = np.array([
                [-1.,1.,1.], [1.,1.,1.], [-1.,-1.,1.], [1.,-1.,1.],
                [-1.,1.,1.], [1.,1.,1.], [-1.,-1.,1.], [1.,-1.,1.],
                [-1.,1.,1.], [1.,1.,1.], [-1.,-1.,1.], [1.,-1.,1.]
            ])
            n = np.array([
                [0.,-1.,1.], [0.,-1.,1.], [0.,1.,1.], [0.,1.,1.],
                [1.,0.,1.], [-1.,0.,1.], [1.,0.,1.], [-1.,0.,1.],
                [1.,1.,0.], [-1.,1.,0.], [-1.,1.,0.], [1.,1.,0.]
            ])
        elif crystal in ['FCC', 'fcc']:
            # the 12 <110>{111} slip systems
            b = np.array([
                [0.,1.,-1.], [1.,0.,-1.], [1.,-1.,0.],
                [0.,1.,-1.], [1.,0.,1.], [1.,1.,0.],
                [0.,1.,1.], [1.,0.,-1.], [1.,1.,0.],
                [0.,1.,1.], [1.,0.,1.], [1.,-1.,0.]
            ])
            n = np.array([
                [1.,1.,1.], [1.,1.,1.], [1.,1.,1.],
                [-1.,1.,1.], [-1.,1.,1.], [-1.,1.,1.],
                [1.,-1.,1.], [1.,-1.,1.], [1.,-1.,1.],
                [1.,1.,-1.], [1.,1.,-1.], [1.,1.,-1.]
            ])
        else:
            raise ValueError('Error: unsupported crystal type = %s in generate_line_config()' % crystal)

        nsys = b.shape[0]
        b = b / np.linalg.norm(b, axis=1)[:,None]
        n = n / np.linalg.norm(n, axis=1)[:,None]
        if Rorient is not None:
            Rorient = np.array(Rorient)
            Rorient = Rorient / np.linalg.norm(Rorient, axis=1)[:,None]
            b = np.matmul(b, Rorient.T)
            n = np.matmul(n, Rorient.T)

        cell = _periodic_cell(Lbox)
        Lmax = np.max(np.linalg.norm(cell.h, axis=0))

        if theta is None:
            # the character angle of each slip system that makes its line length the closest
            # to the longest of the shortest lines across the systems
            ntheta = 19
            theta = 90.0/(ntheta-1)*np.arange(ntheta)
            theta_minlength = np.zeros((nsys, ntheta))
            for isys in range(nsys):
                burg, plane = b[isys], n[isys]
                c = np.array(cell.center())
                for t in range(ntheta):
                    nodes, segs = [], []
                    length = _insert_infinite_line(cell, nodes, segs, burg, plane, c,
                                                   theta=theta[t], maxseg=maxseg, trial=True)
                    theta_minlength[isys,t] = length
            theta_minlength = np.ma.masked_less(theta_minlength, 0.0)
            minlength = theta_minlength.min(axis=1).filled(-1.0)
            maxlength = np.max(minlength)
            if maxlength > 10*Lmax or np.min(minlength) < 0.0:
                raise ValueError('Error: cannot find appropriate line to insert')
            theta_sys = np.argmin(np.abs(theta_minlength-maxlength), axis=1)
            theta_sys = theta[theta_sys][:,None]
        else:
            theta_sys = np.tile(np.array(theta), (nsys, 1))

        if seed > 0: np.random.seed(seed)
        pos = np.random.rand(num_lines, 3)
        pos = np.array(cell.origin) + np.matmul(pos, np.array(cell.h).T)
        ithe = np.random.randint(0, theta_sys.shape[1], num_lines)
        nodes, segs = [], []
        for i in range(num_lines):
            isys = i % nsys
            burg, plane = b[isys], n[isys]
            idip = np.floor(i/nsys).astype(int) % 2  # alternate sign to create dipoles
            lsign = 1-2*idip
            edir = np.cross(plane, burg)
            edir = edir / np.linalg.norm(edir)
            theta = theta_sys[isys,ithe[i-idip*nsys]]
            ldir = np.cos(theta*np.pi/180.0)*burg + np.sin(theta*np.pi/180.0)*edir
            if verbose: print(' insert dislocation: b = %.3f %.3f %.3f, n = %.3f %.3f %.3f, theta = %.1f deg' % (*burg, *plane, theta))
            nodes, segs = _insert_infinite_line(cell, nodes, segs, burg, plane, pos[i],
                                                linedir=lsign*ldir, maxseg=maxseg)
        self._set_network(cell, nodes, segs)
        return self

    def generate_prismatic_config(self, crystal, Lbox, num_loops, radius, maxseg=-1, Rorient=None,
                                  seed=1234, uniform=False):
        """generate_prismatic_config: make the network a set of prismatic loops, and return self

        A port of generate_prismatic_config in core/exadis/python/pyexadis_utils.py, with the
        same random numbers: for the same arguments it gives the network that ExaDisNet's gives.
        The loops cycle through the crystal's Burgers vectors, the 6 1/2<110> of 'fcc' or the 4
        1/2<111> of 'bcc'. radius is one radius, or [min, max]; uniform spreads the loops' centres
        on a grid, shuffled and jittered; Lbox is the box's edge or a Cell; seed > 0 seeds
        numpy's random numbers.
        """
        if crystal in ['BCC', 'bcc']:
            b = np.array([[1.,1.,1.],[-1.,1.,1.],[1.,-1.,1.],[1.,1.,-1.]])
        elif crystal in ['FCC', 'fcc']:
            b = np.array([[1.,1.,0.],[-1.,1.,0.],[1.,0.,1.],[-1.,0.,1.],[0.,1.,1.],[0.,-1.,1.]])
        else:
            raise ValueError('Error: unsupported crystal type = %s in generate_prismatic_config()' % crystal)
        nburg = b.shape[0]
        b = b / np.linalg.norm(b, axis=1)[:,None]

        cell = _periodic_cell(Lbox)
        if seed > 0: np.random.seed(seed)
        if uniform:
            ngrid = np.ceil((1.0*num_loops)**(1.0/3.0))
            H = 1.0/ngrid
            x = 0.5*H + H*np.arange(ngrid)
            x, y, z = np.meshgrid(x, x, x)
            p = np.random.permutation(len(x.flatten()))
            x, y, z = x.flatten()[p], y.flatten()[p], z.flatten()[p]
            pos = np.vstack((x, y, z)).T + 0.5*H*(np.random.rand(len(x), 3)-0.5)
        else:
            pos = np.random.rand(num_loops, 3)
        pos = np.array(cell.origin) + np.matmul(pos, np.array(cell.h).T)
        if isinstance(radius, list):
            R = np.random.uniform(radius[0], radius[1], size=(num_loops,))
        else:
            R = radius*np.ones(num_loops)

        nodes, segs = [], []
        for i in range(num_loops):
            burg = b[i % nburg]
            nodes, segs = _insert_prismatic_loop(crystal, cell, nodes, segs, burg,
                                                 R[i], pos[i], maxseg, Rorient)
        self._set_network(cell, nodes, segs)
        return self


def _periodic_cell(Lbox) -> Cell:
    """The cell of generate_line_config and generate_prismatic_config: a Cell as it is, or a
    periodic box from the origin, cubic for one number and orthorhombic for three, as
    exadis' Cell(Lbox) is."""
    if isinstance(Lbox, Cell):
        return Lbox
    edges = np.array(Lbox, dtype=float) * np.ones(3)
    return Cell(h=np.diag(edges), origin=np.zeros(3), is_periodic=[True, True, True])


def _insert_infinite_line(cell, nodes, segs, burg, plane, origin, theta=0.0, linedir=None,
                          maxseg=-1, trial=False):
    """Add an infinite line to the lists of nodes and segments, as insert_infinite_line in
    core/exadis/python/pyexadis_utils.py does; with trial, only return its length."""
    plane = plane / np.linalg.norm(plane)
    if np.abs(np.dot(burg, plane)) >= 1e-5:
        print('Warning: Burgers vector and plane normal are not orthogonal')
    if not linedir is None:
        ldir = np.array(linedir)
        ldir = ldir / np.linalg.norm(ldir)
    else:
        b = burg / np.linalg.norm(burg)
        y = np.cross(plane, b)
        y = y / np.linalg.norm(y)
        ldir = np.cos(theta*np.pi/180.0)*b+np.sin(theta*np.pi/180.0)*y

    h = np.array(cell.h)
    Lmin = np.min(np.linalg.norm(h, axis=0))
    seglength = 0.15*Lmin
    if maxseg > 0:
        seglength = np.min([seglength, maxseg])

    meet = 0
    maxnodes = 1000
    numnodes = 0
    origin = np.array(origin)
    p = 1.0*origin
    originpbc = 1.0*origin
    while ((~meet) & (numnodes < maxnodes)):
        p += seglength*ldir
        pp = np.asarray(cell.closest_image(Rref=origin, R=p))
        dist = np.linalg.norm(pp-origin)
        if ((numnodes > 0) & (dist < seglength)):
            originpbc = np.asarray(cell.closest_image(Rref=p, R=origin))
            meet = 1
        numnodes += 1

    if numnodes == maxnodes:
        if trial:
            return -1.0
        print('Warning: infinite line is too long, aborting')
        return nodes, segs
    if trial:
        return np.linalg.norm(originpbc-origin)
    istart = len(nodes)
    for i in range(numnodes):
        p = origin + 1.0*i/numnodes*(originpbc-origin)
        nodes.append(np.concatenate((p, [DisNode.Constraints.UNCONSTRAINED])))
    for i in range(numnodes):
        segs.append(np.concatenate(([istart+i, istart+(i+1)%numnodes], burg, plane)))
    return nodes, segs


def _insert_prismatic_loop(crystal, cell, nodes, segs, burg, radius, center, maxseg=-1,
                           Rorient=None):
    """Add a prismatic loop to the lists of nodes and segments, as insert_prismatic_loop in
    core/exadis/python/pyexadis_utils.py does: a hexagon on {110} planes in bcc, a rhombus on
    {111} planes in fcc."""
    b = -1.0*burg
    if crystal in ['BCC', 'bcc']:
        b0 = 1.0/np.sqrt(3.0)*np.array([[1.,1.,1.],[-1.,1.,1.],[1.,-1.,1.],[1.,1.,-1.]])
        bcol = np.abs(np.abs(np.dot(b0, b))-1.0)
        ib = bcol.argmin()
        if bcol[ib] > 1e-5:
            raise ValueError('BCC Burgers vector must be of the 1/2<111> type in insert_prismatic_loop()')
        Nsides = 6
        e = np.array([[-2.0*b[0],b[1],b[2]],[-b[0],-b[1],2.0*b[2]],
                      [b[0],-2.0*b[1],b[2]],[2.0*b[0],-b[1],-b[2]],
                      [b[0],b[1],-2.0*b[2]],[-b[0],2.0*b[1],-b[2]]])
        n = np.cross(b, e[(np.arange(6)+1)%6]-e[np.arange(6)])
        e = e / np.linalg.norm(e, axis=1)[:,None]
    elif crystal in ['FCC', 'fcc']:
        Nsides = 4
        b0 = 1.0/np.sqrt(2.0)*np.array([[0,1,1],[0,-1,1],[1,0,1],[-1,0,1],[1,1,0],[-1,1,0]])
        n01 = np.array([[-1,-1,1],[-1,1,1],[-1,1,1],[1,1,1],[-1,1,1],[1,1,1]])
        n02 = np.array([[1,-1,1],[1,1,1],[-1,-1,1],[1,-1,1],[-1,1,-1],[1,1,-1]])
        bcol = np.abs(np.abs(np.dot(b0, b))-1.0)
        ib = bcol.argmin()
        if bcol[ib] > 1e-5:
            raise ValueError('FCC Burgers vector must be of the 1/2<110> type in insert_prismatic_loop()')
        p1 = n01[ib] / np.linalg.norm(n01[ib])
        p2 = n02[ib] / np.linalg.norm(n02[ib])
        l1 = np.cross(p1, b)
        l1 = l1 / np.linalg.norm(l1)
        l2 = np.cross(p2, b)
        l2 = l2 / np.linalg.norm(l2)
        e = np.array([-0.5*l1-0.5*l2, +0.5*l1-0.5*l2, +0.5*l1+0.5*l2, -0.5*l1+0.5*l2])
        n = np.array([p1, p2, p1, p2])
    else:
        raise ValueError('Error: unsupported crystal type = %s in insert_prismatic_loop()' % crystal)

    n = n / np.linalg.norm(n, axis=1)[:,None]
    if Rorient is not None:
        Rorient = np.array(Rorient)
        Rorient = Rorient / np.linalg.norm(Rorient, axis=1)[:,None]
        b = np.matmul(b, Rorient.T)
        e = np.matmul(e, Rorient.T)
        n = np.matmul(n, Rorient.T)

    istart = len(nodes)
    Nnodes = 0
    for i in range(Nsides):
        l = radius*(e[(i+1)%Nsides]-e[i])
        Nseg = int(np.ceil(np.linalg.norm(l)/maxseg)) if maxseg > 0 else 1
        for j in range(Nseg):
            p = radius*e[i]+1.0*j/Nseg*l+center
            nodes.append(np.concatenate((p, [DisNode.Constraints.UNCONSTRAINED])))
            n1 = istart+Nnodes
            n2 = istart if (i == Nsides-1 and j == Nseg-1) else n1+1
            segs.append(np.concatenate(([n1, n2], b, n[i])))
            Nnodes += 1
    return nodes, segs
