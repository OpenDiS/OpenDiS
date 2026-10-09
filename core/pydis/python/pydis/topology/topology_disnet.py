"""@package docstring
Topology_DisNet: class for splitting multi nodes

Provide topology handling functions given a DisNet object. Each split mode
lives in its own module and is a plain function of
(G, tag, state, force, mobility, ...), so a mode can be called and tested
without constructing a Topology.

Two modes, both selecting the partition of a multi-arm node that releases the
most power, differing in where that power is measured:

'MaxDiss'
    Measures it with the two new nodes left on top of each other. Kept as it
    was, so results produced with it stay reproducible. Needs nothing from
    state.

'Serial'
    Matches the behavior of the 'TopologySerial' topology model of ExaDiS, and
    follows ParaDiS SplitMultiNodes. Moves the two new nodes apart before
    measuring. Reads rann, minseg and a from state.
"""

from functools import partial

from ..disnet import DisNet
from framework.disnet_manager import DisNetManager
from .topology_ops import (build_split_list, init_topology_exemptions,
                           exempt_from_collisions,
                           split_node_and_update_forces, split_multi_nodes)
from .topology_maxdiss import Topology_MaxDiss
from .topology_serial import Topology_Serial, TopologyParams


class Topology:
    """Topology: class for selecting and handling multi node splitting
    """
    def __init__(self, state: dict={}, split_mode: str='MaxDiss', **kwargs) -> None:
        self.split_mode = split_mode
        self.force = kwargs.get('force')
        if not self.force.__module__.split('.')[0] in ['pydis', 'pyexadis_base']:
            print("Topology: force.__module__ = ", self.force.__module__)
            raise ValueError("Topology: force must come compatible modules")
        self.mobility = kwargs.get('mobility')
        if not self.mobility.__module__.split('.')[0] in ['pydis']:
            print("Topology: mobility.__module__ = ", self.mobility.__module__)
            raise ValueError("Topology: mobility must come compatible modules")

        self.Handle_Functions = {
            'MaxDiss': self.Handle_MaxDiss,
            'Serial': self.Handle_Serial }

        if split_mode not in self.Handle_Functions:
            raise ValueError("Topology: unknown split_mode '%s', expected one of %s"
                             % (split_mode, sorted(self.Handle_Functions)))

        # only 'Serial' needs these, and building them here reports a missing
        # one at construction rather than part way into a run
        self.params = TopologyParams(state) if split_mode == 'Serial' else None

    # the low-level operations are kept reachable through the class, since
    # callers outside this package refer to them that way
    build_split_list = staticmethod(build_split_list)
    init_topology_exemptions = staticmethod(init_topology_exemptions)
    exempt_from_collisions = staticmethod(exempt_from_collisions)
    split_node_and_update_forces = staticmethod(split_node_and_update_forces)
    split_multi_nodes = staticmethod(split_multi_nodes)

    def Handle(self, DM: DisNetManager, state: dict) -> dict:
        """Handle: handle topology according to split_mode
        """
        if self.params is not None:
            # dt is written by the time integrator, so it is only available
            # from the second step onwards
            self.params.dt = state.get("dt", self.params.dt)
        G = DM.get_disnet(DisNet)
        self.Handle_Functions[self.split_mode](G, state)
        return state

    def Handle_MaxDiss(self, G: DisNet, state: dict) -> None:
        """Handle_MaxDiss: split_multi_nodes with the 'MaxDiss' criterion
        """
        state = init_topology_exemptions(G, state)
        state = split_multi_nodes(
            G, state, partial(Topology_MaxDiss,
                              force=self.force, mobility=self.mobility))
        return state

    def Handle_Serial(self, G: DisNet, state: dict) -> None:
        """Handle_Serial: split_multi_nodes with the 'Serial' criterion
        """
        state = init_topology_exemptions(G, state)
        state = split_multi_nodes(
            G, state, partial(Topology_Serial,
                              force=self.force, mobility=self.mobility,
                              params=self.params))
        return state
