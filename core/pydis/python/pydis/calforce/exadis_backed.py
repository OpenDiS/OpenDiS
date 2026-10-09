"""ExadisBackedForce: let PyDiS's Topology use ExaDiS's own force module,
with arm order preserved through the conversion

Topology.__init__ (topology_disnet.py) already accepts a bare ExaDiS
CalForce directly for its force= argument (force.__module__ check allows
both 'pydis' and 'pyexadis_base'), and PYDIS_FORCE_SOURCE='exadis'
(test_topology_mode_pydis_exadis.py) has used that directly since it was
added. This wrapper exists because that plain usage is not enough:
CalForce.OneNodeForce(DM, state, tag, ...) does N.get_disnet(ExaDisNet),
which -- when DM wraps a plain PyDiS DisNet, as Topology's trial
evaluation always does -- reaches ExaDiS' own compute_node_force through a
fresh ExaDisNet rebuilt from DisNet.export_data(). That export's segment
order is the DisNet graph's global edge-insertion history, not any node's
own arm order or ExaDiS' own live conn history (see
framework.arm_order's module docstring), and OneNodeForce's underlying
pair sum (ExaDiS' own neighbor/cell-list construction, which candidate
pairs get summed in what order) is sensitive to that array order across
the whole network -- unlike NodeMobility's Lsum, which is a two-term sum
for a 2-armed node and therefore order-insensitive on its own (see
mobility.exadis_backed's docstring for why *that* module still needed
fixing). Measured directly on test5_topology_mode step 194: even with
PYDIS_FORCE_SOURCE='exadis' and PYDIS_MOBILITY_SOURCE='exadis' both set,
comparing PyDiS' own resulting split direction against ExaDiS' true
authoritative one (backed out of its own recorded final position, not a
reconstruction) still showed a ~2.15e-14 relative gap -- this wrapper is
the fix for that remaining gap, built the same way
mobility.exadis_backed.ExadisBackedMobility already is.

Same two-tier construction as that module, and for the same reason: if
state carries the before-trial snapshot (state['_exadis_backed_before']),
framework.arm_order.rebase_export_data(G, before) is used in place of
export_data_with_arms_first, closing the ~2.15e-14 gap above to ~2.5e-15,
a further ~9x -- see rebase_export_data's own docstring for why placing
each tag's own arms first is not, on its own, enough to fix a whole-
network sum. Falls back to export_data_with_arms_first when no before
snapshot is available.
"""


class ExadisBackedForce:
    """ExadisBackedForce: wraps an ExaDiS CalForce for PyDiS's Topology,
    registering an arm-order-preserving ExaDisNet on DM before calling
    through, rather than letting DisNetManager's generic conversion build
    one blind to arm order. See this module's docstring for why.
    """

    def __init__(self, force):
        self.force = force

    def _register_ordered_net(self, DM, state):
        from ..disnet import DisNet
        from pyexadis_base import ExaDisNet
        from framework.arm_order import rebase_export_data, export_data_with_arms_first
        G = DM.get_disnet(DisNet)
        before = state.get('_exadis_backed_before')
        data = (rebase_export_data(G, before) if before is not None
               else export_data_with_arms_first(G, list(G.all_nodes_tags())))
        exadis_net = ExaDisNet()
        exadis_net.import_data(data)
        DM.add_disnet(exadis_net)

    def OneNodeForce(self, DM, state, tag, update_state=True, match_global=False):
        self._register_ordered_net(DM, state)
        return self.force.OneNodeForce(DM, state, tag, update_state=update_state,
                                       match_global=match_global)

    def NodeForce(self, DM, state, pre_compute=True):
        self._register_ordered_net(DM, state)
        return self.force.NodeForce(DM, state, pre_compute)
