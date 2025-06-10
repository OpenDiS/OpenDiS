"""@package docstring
MobilityLaw: class for defining Mobility Laws with special treatment for surface nodes

Provide mobility law functions for all nodes including surface nodes
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from pydis.disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager
from framework.mobility_base import MobilityLaw_Base

from typing import Tuple


class MobilityLaw(MobilityLaw_Base):
    """MobilityLaw: class for mobility laws of surface nodes
    """
    def __init__(self, state: dict={}, mobility_bulk=None) -> None:
        self.mobility_bulk = mobility_bulk
        
    def OneNodeMobility_Project(self, DM: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
        """NodeMobility_Project: project surface node velocity to the surface
        """
        G = DM.get_disnet(DisNet)
        vel = np.zeros([0.0, 0.0, 0.0])
        return vel

    def Mobility(self, DM: DisNetManager, state: dict) -> dict:
        """Mobility: compute all nodal velocities and store them in the state dictionary
        """
        state = self.mobility_bulk.Mobility(DM, state)
        print("Mobility: project surface node velocity to surface")
        #nodeforce_dict = state["nodeforce_dict"]
        #vel_dict = nodeforce_dict.copy()
        #G = DM.get_disnet(DisNet)
        #for tag in G.all_nodes_tags():
        #    f = vel_dict[tag].copy()
        #    vel_dict[tag] = self.OneNodeMobility_Project(DM, tag, f)
        #state["vel_dict"] = vel_dict
        return state

    def OneNodeMobility(self, DM: DisNetManager, state: dict, tag: Tag, f: np.array, update_state: bool=True) -> np.array:
        """OneNodeMobility: compute and return the mobility of one node specified by its tag
        """
        v = mobility_bulk.OneNodeMobility(DM, state, tag, f, update_state)
        #v = self.OneNodeMobility_Project(DM, state, tag, f, update_state)
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
        return v
    
        
