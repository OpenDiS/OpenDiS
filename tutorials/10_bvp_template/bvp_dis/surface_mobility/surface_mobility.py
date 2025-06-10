"""@package docstring
Surface_MobilityLaw: class for defining Mobility Laws for surface nodes

Provide mobility law functions for surface nodes
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from pydis.disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager
from framework.mobility_base import MobilityLaw_Base

from typing import Tuple


class Surface_MobilityLaw(MobilityLaw_Base):
    """Surface_MobilityLaw: class for mobility laws of surface nodes
    """
    def __init__(self, state: dict={}, mobility_law: str='Relax', vmax: float=1e9) -> None:
        self.mobility_law = mobility_law
        self.mob = state.get("mob", 1.0)
        self.vmax = vmax

        self.NodeMobility_Functions = {
            'project': self.NodeMobility_Project
        }
        
    def Mobility(self, DM: DisNetManager, state: dict) -> dict:
        """Mobility: calculate node velocity according to mobility law function
        """
        #G = DM.get_disnet(DisNet)
        #if "nodeforces" in state and "nodeforcetags" in state:
        #    DisNet.convert_nodeforce_array_to_dict(state)

        #nodeforce_dict = state["nodeforce_dict"]
        #vel_dict = nodeforce_dict.copy()
        #for tag in G.all_nodes_tags():
        #    f = vel_dict[tag].copy()
        #    vel_dict[tag] = self.NodeMobility_Functions[self.mobility_law](G, tag, f)
        #state["vel_dict"] = vel_dict

        # prepare nodeforces and nodeforce_tags arrays for compatibility with exadis
        #state = DisNet.convert_nodevel_dict_to_array(state)
        print("Surface_MobilityLaw: Mobility")
        return state
        
    def OneNodeMobility(self, DM: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
        """OneNodeMobility: compute and return the mobility of one node specified by its tag
        """
        G = DM.get_disnet(DisNet)
        v = self.NodeMobility_Functions[self.mobility_law](G, tag, f)
        # update velocity dictionary if needed
        #if update_state:
        #    if "nodevels" in state and "nodeveltags" in state:
        #        nodeveltags = state["nodeveltags"]
        #        ind = np.where((nodeveltags[:,0]==tag[0])&(nodeveltags[:,1]==tag[1]))[0]
        #        if ind.size == 1:
        #            state["nodevels"][ind[0]] = v
        #        else:
        #            state["nodevels"] = np.vstack((state["nodevels"], v))
        #            state["nodeveltags"] = np.vstack((state["nodeveltags"], tag))
        #    else:
        #        state["nodevels"] = np.array([v])
        #        state["nodeveltags"] = np.array([tag])
        #return v
        return None
    
    def NodeMobility_Project(self, G: DisNet, tag: Tag, f: np.array) -> np.array:
        """NodeMobility_Project: project surface node velocity to the surface
        """
        print("NodeMobility_Project: project surface node velocity to surface")
        #vel = 1.0*f
        #for tag in G.all_nodes_tags():
        #    if G.nodes(tag).constraint == DisNode.Constraints.SURFACE_NODE:
        #        vel = np.zeros(3)
        return None
        
