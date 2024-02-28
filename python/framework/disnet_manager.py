"""@package docstring
DisNetManager: class for managing multiple implementations of dislocation network

Implements synchronization between different implementations of DisNet
"""

import numpy as np
import networkx as nx
from typing import Tuple
from copy import deepcopy
from enum import IntEnum
import itertools

from base_classes.disnet_base import DisNet_BASE

class DisNetManager:
    """Class for managing multiple implementations of dislocation network

    Implements synchronization between different implementations of DisNet
    """
    def __init__(self, disnet_dict: dict={}):
        self.disnet_dict = disnet_dict

    def add_disnet(self, disnet):
        """Add DisNet object of disnet_type
        """
        self.disnet_dict[type(disnet)] = disnet

    def get_disnet(self, disnet_type):
        """Get DisNet object of disnet_type
        """
        return self.disnet_dict[disnet_type]
    
    def is_active(self):
        """Return the type of DisNet that is active
        """
        return True
