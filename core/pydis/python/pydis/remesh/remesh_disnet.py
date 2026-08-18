"""@package docstring
Remesh_DisNet: class for defining Remesh functions

Provide remesh functions given a DisNet object. Each rule lives in its own
module and is a plain function of (network, params), so a rule can be called
and tested without constructing a Remesh.
"""

from ..disnet import DisNet
from framework.disnet_manager import DisNetManager
from .remesh_lengthbased import Remesh_LengthBased


class RemeshParams:
    """RemeshParams: the lengths the remesh rules work to"""

    def __init__(self, state: dict) -> None:
        self.maxseg = state.get("maxseg", None)
        self.minseg = state.get("minseg", None)


class Remesh:
    """Remesh: class for remeshing dislocation network

    """
    def __init__(self, state: dict={}, remesh_rule: str='LengthBased') -> None:
        self.remesh_rule = remesh_rule
        self.params = RemeshParams(state)

        self.Remesh_Functions = {
            'LengthBased': Remesh_LengthBased }

        if remesh_rule not in self.Remesh_Functions:
            raise ValueError("Remesh: unknown remesh_rule '%s', expected one of %s"
                             % (remesh_rule, sorted(self.Remesh_Functions)))

    # kept so callers that read them still work
    @property
    def maxseg(self):
        return self.params.maxseg

    @property
    def minseg(self):
        return self.params.minseg

    def Remesh(self, DM: DisNetManager, state: dict) -> None:
        """Remesh: remesh dislocation network according to remesh_rule
        """
        G = DM.get_disnet(DisNet)
        self.Remesh_Functions[self.remesh_rule](G, self.params)
        return state
