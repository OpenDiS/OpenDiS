"""@package docstring
Surface_Topology: class for enforcing topology rules on the surface

Provide topology handling functions for nodes on the surfaces
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from copy import deepcopy
import itertools

from pydis.disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager

class Surface_Topology:
    """Surface_Topology: class for selecting and handling multi node splitting
    """
    def __init__(self, state: dict={}, enforce_mode: str='CutSeg', **kwargs) -> None:
        self.enforce_mode = enforce_mode
        self.Handle_Functions = {
            'CutSeg': self.Handle_CutSeg }
        
    def Handle(self, DM: DisNetManager, state: dict) -> dict:
        """Handle: handle topology according to split_mode
        """
        G = DM.get_disnet(DisNet)
        self.Handle_Functions[self.enforce_mode](G, state)
        return state

    def Handle_CutSeg(self, G: DisNet, state: dict) -> None:
        """Handle_MaxDiss: split_multi_nodes
        """
        #print("Surface_Topology: Handle_CutSeg")
        #state = Topology.init_topology_exemptions(G, state)
        #state = Topology.split_multi_nodes(G, state, self.force, self.mobility)
        print("Surface_Topology: Handle_CutSeg")
        return state
