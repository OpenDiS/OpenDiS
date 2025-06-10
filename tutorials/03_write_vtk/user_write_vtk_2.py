import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../core/exadis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

from framework.disnet_manager import DisNetManager
from pydis import SimulateNetwork


class My_SimulateNetwork(SimulateNetwork):
    """My_SimulateNetwork: User-defined simulation driver
    """
    def step_write_files(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                from pyexadis_utils import write_vtk
                write_vtk(DM, f'{self.write_dir}/config.{state["istep"]}.vtk')
