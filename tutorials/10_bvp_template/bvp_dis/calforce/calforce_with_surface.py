"""@package docstring
CalImageForce: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from pydis.disnet import DisNet, Tag
from framework.disnet_manager import DisNetManager
from framework.calforce_base import CalForce_Base

class CalForce(CalForce_Base):
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}, calforce_bulk=None) -> None:
        self.calforce_bulk = calforce_bulk
        pass

    def AddImageForce(self, DM: DisNetManager, state: dict) -> dict:
        """AddImageForce: add image force from image stress on nodes

        """
        print("CalImageForce: AddImageForce")
        return state

    def NodeForce(self, DM: DisNetManager, state: dict, pre_compute: bool=True) -> dict:
        """NodeForce: compute all nodal forces and store them in the state dictionary
        """
        state = self.calforce_bulk.NodeForce(DM, state, pre_compute)

        state = self.AddImageForce(DM, state)
        return state

    def PreCompute(self, DM: DisNetManager, state: dict) -> dict:
        """PreCompute: pre-compute some data for force calculation
        """
        return state

    def OneNodeForce(self, DM: DisNetManager, state: dict, tag: Tag, update_state: bool=True) -> np.array:
        """OneNodeForce: compute and return the force on one node specified by its tag
        """
        state = self.calforce_bulk.OneNodeForce(DM, state, tag, update_state)

        #ToDo: Add Image force contribution to OneNodeForce
        return state
