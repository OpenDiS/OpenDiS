"""@package docstring
Topology_DisNet: class for splitting multi nodes

Provide topology handling functions given a DisNet object
"""

import numpy as np
from copy import deepcopy
import itertools

from ..disnet import DisNet, DisNode, Tag
from framework.disnet_manager import DisNetManager

class Topology:
    """Topology: class for selecting and handling multi node splitting
    """
    def __init__(self, state: dict={}, split_mode: str='MaxDiss', **kwargs) -> None:
        self.split_mode = split_mode
        self.force = kwargs.get('force')
        self.mobility = kwargs.get('mobility')
        self.Handle_Functions = {
            'MaxDiss': self.Handle_MaxDiss }
        
    def Handle(self, DM: DisNetManager, state: dict) -> None:
        """Handle: handle topology according to split_mode
        """
        G = DM.get_disnet(DisNet)
        self.Handle_Functions[self.split_mode](G, state)
        return state

    @staticmethod
    def build_split_list(n: int) -> list:
        """build_split_list: build a list of length n with True's and False's
           there must be at least two True's and at least 2 False's
        """
        indices = list(range(n))
        all_bool_list = list(itertools.product([True, False], repeat=n))
        # the first element must be selected to avoid double counting
        selected_list = [list(itertools.compress(indices,item)) for item in all_bool_list if item[0] and sum(item) >= 2 and sum(item) <= n-2]
        return selected_list

    @staticmethod
    def init_topology_exemptions(G) -> None:
        """init_topology_exemptions: initialize the topology exemptions
        """
        for node in G.all_nodes_dict().values():
            node.flag &= ~(DisNode.Flags.NO_COLLISIONS | DisNode.Flags.NO_MESH_COARSEN)
        return

    @staticmethod
    def trial_split_multi_node(G, tag: Tag, state: dict, force, mobility, power_th=1e-3) -> dict:
        """trial_split_multi_node: try to split multi-arm node in different ways
            and select the way that maximizes the power dissipation
        """
        # To do: implement this step
        """ ParaDiS SplitMultiNode.c:
            Temporary kludge!  If any of the arms of the multinode is less
            than the length <shortSeg>, don't do the split.  This is an
            attempt to prevent node splits that would leave extremely
            small segments which would end up oscillating with very high
            velocities (and hence significantly impacting the timestep)
        """

        vel_dict = state["vel_dict"]
        nodeforce_dict = state["nodeforce_dict"]
        segforce_dict = state["segforce_dict"]

        state_trial = deepcopy(state)

        n_degree = G.out_degree(tag)
        nbrs = list(G.neighbors(tag))
        nbr_idx_list = Topology.build_split_list(n_degree)

        power0 = np.dot(nodeforce_dict[tag], vel_dict[tag])
        #print("trial_split_multi_node (%s): power0 = %e"%(tag, power0))
        #print("edges = %s"%(str(G.edges(tag))))

        pos0 = G.nodes(tag).R
        n_splits = len(nbr_idx_list)
        power_diss = np.zeros(n_splits)
        for k in range(n_splits):
            nbrs_to_split = [nbrs[i] for i in nbr_idx_list[k]]

            # make a copy of the network G to make trial splits
            G_trial = G.copy()

            segforce_trial_dict = deepcopy(segforce_dict)

            # attempt to split node
            split_node1, split_node2 = G_trial.split_node(tag, pos0.copy(), pos0.copy(), nbrs_to_split)
            #print("split_node1 = %s, split_node2 = %s"%(str(split_node1), str(split_node2)))

            # modify the segforce_trial_dict to reflect the trial split
            for segment in segforce_dict:
                tag1, tag2 = segment
                if tag1 == tag or tag2 == tag:
                    if not G_trial.has_segment(tag1, tag2):
                        #print("remove segment (%s, %s) from segforce_trial_dict" % (str(tag1), str(tag2)))
                        segforce_trial_dict.pop(segment)
                        if tag1 == tag and G_trial.has_segment(split_node2, tag2):
                            new_segment = (split_node2, tag2)
                        elif tag2 == tag and G_trial.has_segment(tag1, split_node2):
                            new_segment = (tag1, split_node2)
                        else:
                            raise ValueError("trial_split_multi_node: cannot find corresponding segment (%s, %s) in G_trial" % (str(tag1), str(tag2)))
                        
                        #print("add segment %s to segforce_trial_dict" % str(new_segment))
                        segforce_trial_dict[new_segment] = segforce_dict[segment]

            # calculate nodal forces and velocities for the trial split
            state_trial["segforce_dict"] = segforce_trial_dict
            state_trial = force.NodeForce_from_SegForce(G_trial, state_trial)
            DM_trial = DisNetManager(G_trial)
            mobility.Mobility(DM_trial, state_trial)
            nodeforce_dict_trial = state_trial["nodeforce_dict"]
            vel_dict_trial = state_trial["vel_dict"]

            power_diss[k] = np.dot(nodeforce_dict_trial[split_node1], vel_dict_trial[split_node1]) \
                          + np.dot(nodeforce_dict_trial[split_node2], vel_dict_trial[split_node2])

            #print("trial_split_multi_node (%s): power_diss[%d] = %e"%(str(tag), k, power_diss[k]))

        # To do: adjust power_th to be consistent with ParaDiS
        if np.max(power_diss) - power0 > power_th:
            # select the split that leads to the maximum power dissipation
            k_sel = np.argmax(power_diss)
            do_split = True
            #print("trial_split_multi_node (%s): power0 = %e power_diss[%d] = %e"%(str(tag), power0, k_sel, power_diss[k_sel]))
        else:
            do_split = False

        # To do: perhaps we should just record which split is selected
        #        and do the actual split outside of this function
        if do_split:
            nbrs_to_split = [nbrs[i] for i in nbr_idx_list[k_sel]]
            split_node1, split_node2 = G.split_node(tag, pos0.copy(), pos0.copy(), nbrs_to_split)

            # To do: remove repeated code
            segforce_trial_dict = deepcopy(segforce_dict)
            # modify the segforce_trial_dict to reflect the trial split
            for segment in segforce_dict:
                tag1, tag2 = segment
                if tag1 == tag or tag2 == tag:
                    if not G.has_segment(tag1, tag2):
                        #print("remove segment (%s, %s) from segforce_trial_dict" % (str(tag1), str(tag2)))
                        segforce_trial_dict.pop(segment)
                        if tag1 == tag and G.has_segment(split_node2, tag2):
                            new_segment = (split_node2, tag2)
                        elif tag2 == tag and G.has_segment(tag1, split_node2):
                            new_segment = (tag1, split_node2)
                        else:
                            raise ValueError("trial_split_multi_node: cannot find corresponding segment (%s, %s) in G" % (str(tag1), str(tag2)))

                        #print("add segment %s to segforce_trial_dict" % str(new_segment))
                        segforce_trial_dict[new_segment] = segforce_dict[segment]

            state = force.NodeForce_from_SegForce(G, state)
            state = mobility.Mobility(DisNetManager(G), state)

            # Mark both nodes involved in the split as 'exempt' from subsequent collisions this time step
            G.nodes(split_node1).flag |= DisNode.Flags.NO_COLLISIONS
            G.nodes(split_node2).flag |= DisNode.Flags.NO_COLLISIONS

            # Update segforce_dict to segforce_trial_dict
            state["segforce_dict"] = segforce_trial_dict
        #else:
        #    return segforce_dict
        return state

    @staticmethod
    def split_multi_nodes(G, state: dict, force, mobility, max_degree=15) -> None:
        """split_multi_nodes: examines all nodes with at least four arms and decides
           if the node should be split and some of the node's arms moved to a new node.
           guarantees sanity after operation

           This function calls the lower level split_node() function.
        """
        nodes = list(G.all_nodes())
        for tag in nodes:
            n_degree = G.out_degree(tag)

            if n_degree < 4:
                continue
            elif n_degree > max_degree:
                raise ValueError("split_multi_node: Node %s has more than %d arms" % (str(tag), n_degree))

            state = Topology.trial_split_multi_node(G, tag, state, force, mobility)

        return state

    def Handle_MaxDiss(self, G: DisNet, state: dict) -> None:
        """Handle_MaxDiss: split_multi_nodes
        """
        Topology.init_topology_exemptions(G)
        state = Topology.split_multi_nodes(G, state, self.force, self.mobility)
        return state
