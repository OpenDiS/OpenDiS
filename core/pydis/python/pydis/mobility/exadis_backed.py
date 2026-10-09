"""ExadisBackedMobility: let PyDiS's Topology use ExaDiS's own mobility law

CalForce already supports this (Topology.__init__ accepts a force from
either 'pydis' or 'pyexadis_base', examples/02_frank_read_src/
test_frank_read_src_pydis_exadis.py demonstrates it); mobility does not --
Topology.__init__ (topology_disnet.py) requires mobility.__module__ to start
with 'pydis'. This module exists so mobility can be shared the same way,
without loosening that check: the wrapper's own module is 'pydis', so it
passes, while its actual computation is ExaDiS's own compiled MobilityGlide.

Why this is needed at all, not just force sharing: plan_topology_mode.md
section 10 measured test5_topology_mode's last position residual (down to
~1.6e-13 after arm-order fixes) to a node's mobility-scaled velocity itself,
not its force -- pydis's own numpy NodeMobility_SimpleGlide and ExaDiS's own
MobilityGlide::node_velocity compute the same formula, but numpy and
Kokkos/C++ do not round every step of it identically. Sharing the force
alone (PYDIS_FORCE_SOURCE='exadis') cannot close that; sharing the mobility
computation too does, since then both sides call the identical compiled
code for the entire force-to-velocity path.

Bridges through MobilityLaw.Mobility(DM, state), not OneNodeMobility: the
latter has a confirmed upstream bug (core/exadis/python/pyexadis_base.py,
`return f` instead of `return v`, reported to the ExaDiS developer
2026-08-25) that silently returns the input force unchanged. Mobility()
has no such bug -- it writes the correct value into
state["nodevels"]/state["nodeveltags"] and its return value is state
itself, not a value that could go stale the same way.

Still not the whole story: MobilityLaw.Mobility(DM, state) does
N.get_disnet(ExaDisNet), which -- when DM wraps a plain PyDiS DisNet, as
Topology's trial evaluation always does -- reaches ExaDiS' own
node_velocity through a fresh ExaDisNet rebuilt from DisNet.export_data(),
whose segment order does not preserve any node's own arm order (see
framework.arm_order's module docstring). ExaDiS' node_velocity sums
directly over conn[i].node[j] for j in range(nconn), so this rebuild can
still scramble the very sum PYDIS_MOBILITY_SOURCE='exadis' exists to share
bit-for-bit, even though the compiled code computing it is now identical
on both sides. Measured directly on test5_topology_mode step 194: sharing
the mobility computation alone, without also fixing this, left the
velocity as far from ExaDiS' own authoritative value as PyDiS's own numpy
formula did. Mobility() below builds that ExaDisNet itself before calling
through, rather than letting DisNetManager's generic conversion build one
blind to arm order.

Two ways to build it, in order of preference: if state carries the
before-trial ExaDisNet.export_data() snapshot (state['_exadis_backed_before'],
set by a caller that has it -- test5_topology_mode's harness sets it right
before Topology.Handle, since it already builds one fresh per step for its
own comparison), framework.arm_order.rebase_export_data(G, before) is used;
otherwise, framework.arm_order.export_data_with_arms_first(G, tags) is the
fallback. rebase_export_data is the stronger fix -- measured directly on
step 194's real trial, comparing PyDiS' resulting split direction against
ExaDiS' own true one, export_data_with_arms_first alone left a 2.15e-14
relative gap despite every node's own arm order already being right;
rebase_export_data closes it to 2.5e-15, a further ~9x. See
rebase_export_data's own docstring for why: it is not enough for a node's
own arms to be in the right relative order, since a whole-network pair sum
(what this feeds, through the force this mobility is scaled from) is
sensitive to every node's arms, not just the one being evaluated, and
export_data_with_arms_first only fixes that for whichever tags it is told
to place first. Not every caller has a before snapshot on hand, though
(most never trigger a topology trial where "before" has a distinct
meaning), so the fallback stays.
"""

from ..disnet import DisNet


class ExadisBackedMobility:
    """ExadisBackedMobility: wraps an ExaDiS MobilityLaw for PyDiS's Topology

    mobility is a pyexadis_base.MobilityLaw instance (or anything with a
    matching Mobility(DM, state) method). Mobility() below calls it, then
    converts its array-convention output (state["nodevels"]/
    state["nodeveltags"]) into the dict convention (state["vel_dict"])
    PyDiS's own topology code reads, via DisNet.convert_nodevel_array_to_dict
    -- the same bridge NodeMobility_SimpleGlide's own callers already rely
    on for the reverse direction.
    """

    def __init__(self, mobility):
        self.mobility = mobility

    def Mobility(self, DM, state):
        # PyDiS's own MobilityLaw.Mobility() bridges state["nodeforces"]/
        # nodeforcetags into state["nodeforce_dict"] as a side effect
        # (mobility_disnet.py), which evaluate_trial_split reads from
        # afterwards; ExaDiS' own Mobility() only ever touches arrays, so
        # that bridge has to happen here instead, in both directions.
        if "nodeforces" in state and "nodeforcetags" in state:
            state = DisNet.convert_nodeforce_array_to_dict(state)

        # Build the ExaDisNet ourselves, arm order intact, and hand it to
        # DM directly -- see this module's docstring for why the generic
        # DisNetManager conversion (N.get_disnet(ExaDisNet), which
        # self.mobility.Mobility would otherwise trigger) is not enough.
        from pyexadis_base import ExaDisNet
        from framework.arm_order import rebase_export_data, export_data_with_arms_first
        G = DM.get_disnet(DisNet)
        before = state.get('_exadis_backed_before')
        data = (rebase_export_data(G, before) if before is not None
               else export_data_with_arms_first(G, list(G.all_nodes_tags())))
        exadis_net = ExaDisNet()
        exadis_net.import_data(data)
        DM.add_disnet(exadis_net)

        state = self.mobility.Mobility(DM, state)
        state = DisNet.convert_nodevel_array_to_dict(state)
        return state
