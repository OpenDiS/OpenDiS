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

class CalImageForce():
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}) -> None:
        pass

    def AddImageForce(self, DM: DisNetManager, state: dict) -> dict:
        """AddImageForce: add image force from image stress on nodes

        """
        print("CalImageForce: AddImageForce")
        return state
