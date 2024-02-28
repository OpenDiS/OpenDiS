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
        self._active_type = None
        self._last_active_type = None

    def add_disnet(self, disnet):
        """Add DisNet object of disnet_type
        """
        self.disnet_dict[type(disnet)] = disnet

    def synchronize_disnet(self, disnet_src, disnet_des):
        """Synchronize DisNet between disnet_src and disnet_des
        """
        if disnet_src == disnet_des:
            raise ValueError("synchronize_disnet: disnet_src and disnet_des are the same")

        if disnet_src in self.disnet_dict:
            G_src = self.disnet_dict[disnet_src]
        else:
            raise ValueError("synchronize_disnet: disnet_src not found")

        if disnet_des in self.disnet_dict:
            G_des = self.disnet_dict[disnet_des]
        else:
            # call constructor to create default DisNet object
            G_des = disnet_des()
            self.add_disnet(G_des)

        G_des.import_data(G_src.export_data())

        self._last_active_type = disnet_des

    def get_disnet(self, disnet_type):
        """Get DisNet object of disnet_type
        """
        self._active_type = disnet_type
        if self._last_active_type is not None and self._last_active_type != self._active_type:
            self.synchronize_disnet(self._last_active_type, self._active_type)
        return self.disnet_dict[disnet_type]
    
    def get_active_type(self):
        """Return the type of DisNet that is active
        """
        return self._active_type
