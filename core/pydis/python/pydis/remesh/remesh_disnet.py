"""@package docstring
Remesh_DisNet: class for defining Remesh functions

Provide remesh functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet, DisNode
from framework.disnet_manager import DisNetManager

class Remesh:
    """Remesh: class for remeshing dislocation network

    """
    def __init__(self, state: dict={}, remesh_rule: str='LengthBased') -> None:
        self.remesh_rule = remesh_rule
        self.maxseg = state.get("maxseg", None)
        self.minseg = state.get("minseg", None)

        self.Remesh_Functions = {
            'LengthBased': self.Remesh_LengthBased,
            'RemeshRule_2_ParaDiS': self.RemeshRule_2_ParaDiS }
        
    def Remesh(self, DM: DisNetManager, state: dict) -> None:
        """Remesh: remesh dislocation network according to remesh_rule
        """
        G = DM.get_disnet(DisNet)
        self.Remesh_Functions[self.remesh_rule](G)
        return state

    def Remesh_LengthBased(self, G: DisNet) -> None:
        """Remesh_LengthBased: remesh dislocation network according to segment length
        """
        # mesh coarsen
        nodes_to_remove = []
        segs_data_with_positions = G.get_segs_data_with_positions()
        Nseg = segs_data_with_positions["nodeids"].shape[0]
        source_tags = segs_data_with_positions["tag1"]
        target_tags = segs_data_with_positions["tag2"]
        R1 = segs_data_with_positions["R1"]
        R2 = segs_data_with_positions["R2"]
        for i in range(Nseg):
            tag1, tag2 = tuple(source_tags[i]), tuple(target_tags[i])
            node1, node2 = G.nodes(tag1), G.nodes(tag2)
            #R1, R2 = node1.R, node2.R
            r1 = R1[i,:].copy()
            r2 = R2[i,:].copy()
            # apply PBC
            r2 = G.cell.closest_image(Rref=r1, R=r2)
            L = np.linalg.norm(r2-r1)
            if (L < self.minseg):
                if G.out_degree(tag1) == 2 and node1.constraint != DisNode.Constraints.PINNED_NODE:
                    nodes_to_remove.append(tag1)
                elif G.out_degree(tag2) == 2 and node2.constraint != DisNode.Constraints.PINNED_NODE:
                    nodes_to_remove.append(tag2)
        for tag in set(nodes_to_remove):
            if G.has_node(tag):
                G.remove_two_arm_node(tag)

        if not G.is_sane():
            raise ValueError("Remesh_LengthBased: sanity check failed 1")

        # mesh refine
        all_segments_list = list(G.all_segments())
        for tag1, tag2 in all_segments_list:
            node1, node2 = G.nodes(tag1), G.nodes(tag2)
            r1, r2 = node1.R, node2.R
            # apply PBC
            r2 = G.cell.closest_image(Rref=r1, R=r2)
            L = np.linalg.norm(r2-r1)
            if (L > self.maxseg) and ((node1.constraint != DisNode.Constraints.PINNED_NODE) or (node2.constraint != DisNode.Constraints.PINNED_NODE)):
                # insert new node on segment
                new_tag = G.get_new_tag()
                r = (r1 + r2)/2.0
                G.insert_node_between(tag1, tag2, new_tag, r)
                if not G.is_sane():
                    raise ValueError("Remesh_LengthBased: sanity check failed 1a")

        if not G.is_sane():
            raise ValueError("Remesh_LengthBased: sanity check failed 2")

    def RemeshRule_2_ParaDiS(self, G: DisNet) -> None:
        """RemeshRule_2_ParaDiS: using RemeshRule_2 of ParaDiS
        """
        raise NotImplementedError("RemeshRule_2_ParaDiS: not implemented yet")
