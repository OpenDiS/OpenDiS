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
        for segment in G.seg_list():
            tag1, tag2 = segment["edge"][0], segment["edge"][1]
            node1, node2 = G.nodes(tag1), G.nodes(tag2)
            R1, R2 = node1.R, node2.R
            # apply PBC
            R2 = G.cell.closest_image(Rref=R1, R=R2)
            L = np.linalg.norm(R2-R1)
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
        for segment in G.seg_list():
            tag1, tag2 = segment["edge"][0], segment["edge"][1]
            node1, node2 = G.nodes(tag1), G.nodes(tag2)
            R1, R2 = node1.R, node2.R
            # apply PBC
            R2 = G.cell.closest_image(Rref=R1, R=R2)
            L = np.linalg.norm(R2-R1)
            if (L > self.maxseg) and ((node1.constraint != DisNode.Constraints.PINNED_NODE) or (node2.constraint != DisNode.Constraints.PINNED_NODE)):
                # insert new node on segment
                new_tag = G.get_new_tag()
                # To do: apply PBC here
                R = (R1 + R2)/2.0
                G.insert_node_between(tag1, tag2, new_tag, R)
                if not G.is_sane():
                    raise ValueError("Remesh_LengthBased: sanity check failed 1a")

        if not G.is_sane():
            raise ValueError("Remesh_LengthBased: sanity check failed 2")

    def RemeshRule_2_ParaDiS(self, G: DisNet) -> None:
        """RemeshRule_2_ParaDiS: using RemeshRule_2 of ParaDiS
        """
        raise NotImplementedError("RemeshRule_2_ParaDiS: not implemented yet")
