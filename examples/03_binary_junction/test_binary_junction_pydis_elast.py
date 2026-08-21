import numpy as np
import sys, os

from pathlib import Path
# this script lives in examples/03_binary_junction/, so the repository root is
# 2 levels up
opendis_root = Path(__file__).resolve().parents[2]
pydis_paths = [str(opendis_root / p)
               for p in ['python', 'lib', 'core/pydis/python']]
[sys.path.append(path) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from framework.simulation_setup import check_cutoff_maxseg
from framework.simulation_setup import remesh_initial_config
from framework.testing import write_ref_npz
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

# how far line 1's two pinned endpoints are shifted along the line's own
# direction, as a fraction of z0; see the comment in init_two_disl_lines
EPS_ARM_ASYMMETRY = 1.0e-7


def init_two_disl_lines(z0=1.0, box_length=8.0,
                        b1=np.array([-1.0, 1.0, 1.0]),
                        b2=np.array([1.0, -1.0, 1.0]),
                        pbc=False, maxseg=None):
    '''Generate an initial configuration for two dislocation lines.
    '''
    print("init_two_disl_lines: z0 = %f" % (z0))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc, pbc, pbc])

    PINNED = DisNode.Constraints.PINNED_NODE
    FREE   = DisNode.Constraints.UNCONSTRAINED
    rn = np.array([[0.0, -z0, -z0, PINNED],
                   [0.0,  0.0, 0.0, FREE],
                   [0.0,  z0,  z0,  PINNED],
                   [-z0,  0.0, -z0, PINNED],
                   [0.0,  0.0, 0.0, FREE],
                   [z0,   0.0, z0,  PINNED]])
    rn[:, 0:3] += cell.center()

    # Break the exact point-inversion symmetry of the two lines about the box
    # centre (node 0 is otherwise -node 2, and node 3 is otherwise -node 5)
    # by shifting line 1's two pinned endpoints (nodes 0 and 2) along the
    # line's own direction, by the same amount and in the same sense, so
    # node 1 (the free centre node, which is what actually collides and
    # splits) does not move at all. This shortens one of line 1's two arms
    # and lengthens the other where they meet at node 1 -- a length
    # asymmetry local to the arms directly involved in the split, rather
    # than a transverse perturbation of a distant point. That distinction
    # was measured to matter (.plan/2026-08-21/debug_topology_stage2.md):
    # a transverse perturbation of a far node only reaches the split
    # decision through the long-range elastic term, which Ec=2.8e10 (see
    # the note on CalForce below) now dwarfs, so it could be swept over five
    # orders of magnitude with no effect at all on which side the split
    # decision fell on. An arm-length asymmetry instead changes the self
    # force ParaDiS's SelfForceIsotropic computes per segment
    # (external/paradis/src/NodeForce.c:2628, the 'S' term, an explicit
    # function of segment length), which is local to whichever arm carries
    # it and is not swamped by Ec the same way. Applied identically in the
    # exadis version of this example.
    line1_dir = rn[2, :3] - rn[1, :3]
    line1_dir = line1_dir / np.linalg.norm(line1_dir)
    rn[0, :3] += EPS_ARM_ASYMMETRY*z0*line1_dir
    rn[2, :3] += EPS_ARM_ASYMMETRY*z0*line1_dir

    xi1, xi2 = rn[2, :3] - rn[1, :3], rn[5, :3] - rn[4, :3]
    n1, n2 = np.cross(b1, xi1), np.cross(b2, xi2)
    n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0, :] = np.concatenate(([0, 1], b1, n1))
    links[1, :] = np.concatenate(([1, 2], b1, n1))
    links[2, :] = np.concatenate(([3, 4], b2, n2))
    links[3, :] = np.concatenate(([4, 5], b2, n2))

    # Condition the configuration so that no segment exceeds maxseg. Each arm
    # here starts at z0*sqrt(2) long, several times maxseg, and the neighbor
    # search relies on the cutoff + maxseg bound holding from the first force
    # evaluation, before any remesh has run. Neither arm is pinned at both
    # ends, so LengthBased remesh would split them on its own after one step;
    # conditioning up front simply means the bound is never violated. Inserted
    # nodes are PINNED only when both parents are, which is never the case
    # here, so the new nodes are free and the geometry is unchanged.
    if maxseg is not None:
        rn, links = remesh_initial_config(rn, links, maxseg, cell.h,
                                          [pbc]*3)

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))


def main(plot=True, max_step=200, dt=1.0e-9, print_freq=10, write_freq=10):
    global net, sim, state

    # Lengths are in units of the Burgers vector, and this case is set up at
    # the same scale as 01_loop and 02_frank_read_src rather than at the
    # box_length=8 of test_binary_junction_pydis.py. At that scale the core
    # radius is a third of a segment, so log(L/a) is order 1 and the elastic
    # self-interaction is degenerate with the core; here a = 1 against
    # segments of order 40 separates the two.
    Lbox = 1000.0
    z0 = 0.125*Lbox
    state = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0,
             "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}

    # cutoff + maxseg <= d/3; see the note on CalForce below
    cutoff = 0.25*Lbox
    check_cutoff_maxseg(Lbox*np.eye(3), cutoff, state["maxseg"])

    net = init_two_disl_lines(z0=z0, box_length=Lbox, pbc=False,
                              maxseg=state["maxseg"])
    nbrlist = CellList(cell=net.cell, n_div=[4, 4, 4])
    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])

    if plot:
        try:
            vis = VisualizeNetwork(bounds=bounds)
        except:
            print("")
            print("Failed to create VisualizeNetwork object")
            print("Try run with option  --no-plot")
            print("")
            raise
    else:
        vis = None

    # Elasticity_SBA is the full segment-segment elastic interaction, in
    # contrast to the LineTension mode of test_binary_junction_pydis.py. It is
    # what makes this case a junction problem rather than a line-tension one:
    # the two lines attract or repel according to their Burgers vectors, and
    # whether a junction forms and how far it zips is an elastic question that
    # line tension alone cannot answer.
    #
    # The cutoff is stated explicitly and matched by the companion exadis run
    # so the two truncate the segment-segment sum identically. Note that the
    # two lines span a cube of side 2*z0 = 250, so the most distant pairs in
    # this configuration sit at ~433 apart and ARE truncated at cutoff = 250.
    # That is a deliberate choice for comparability rather than a claim that
    # the far field is negligible; both codes drop the same pairs.
    #
    # cutoff + maxseg <= d/3 is checked above. ExaDiS bins segments by
    # mid-point and clamps the bin count up to a minimum of 3, and past that
    # limit the +-1 bin scan compared mid-points with a periodic shift
    # quantized to the bin grid rather than the true minimum image, silently
    # dropping pairs. Upstream fixed that in exadis commit f8f9861 by falling
    # back to the true minimum image whenever the clamp is active, so the
    # bound is no longer load-bearing for correctness; it is kept because it
    # documents the regime and costs nothing.
    # Ec=2.8e10, not the ParaDiS default mu/(4*pi)*log(a/0.1) (~2.9317e10 for
    # this case's mu and a). The mirror-symmetric tie at the first split (both
    # candidate split directions equally valid by symmetry, decided only by
    # ~1e-13-level roundoff) turned out not to have a clean fix: sweeping Ec
    # with an identical literal value on both codes found agreement only in
    # the narrow band ~2.8e10-2.9e10, with mismatches immediately outside it
    # on both sides -- including at the exact default value. That is a
    # fragile, coincidental overlap of each code's own roundoff, not a
    # physical threshold, so this value is not expected to generalize; see
    # the companion exadis script for the same value.
    calforce  = CalForce(force_mode='Elasticity_SBA', state=state,
                         Ec=2.8e10, cutoff=cutoff)
    # vmax is raised above its 1e9 default for the same reason as in
    # 01_loop/test_disl_loop_pydis_elast.py: the default would clamp nodes
    # partway through and make this disagree with the exadis run, whose GLIDE
    # mobility applies no velocity cap.
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state,
                            vmax=1.0e15)
    timeint   = TimeIntegration(integrator='EulerForward', dt=dt, state=state)

    # 'Serial' rather than 'MaxDiss', to match the 'TopologySerial' model the
    # companion exadis run uses. The two lines start with a node each at the box
    # centre, so the first collision merges them into a node with 4 arms, and
    # splitting that node into two 3-arm nodes joined by a junction segment IS
    # the physics this case exists to show. 'MaxDiss' measures the power
    # released with the two trial nodes still on top of each other, which cannot
    # see the energy released by pulling the two lines apart, and so leaves the
    # 4-arm node alone. 'Serial' moves them apart first, as ParaDiS does, and
    # reads rann, minseg and a from state.
    #
    # Topology cannot simply be disabled: the Proximity collision handler reads
    # the nodeflag_dict that only Topology.init_topology_exemptions creates, so
    # topology=None fails with KeyError: 'nodeflag_dict'.
    topology  = Topology(split_mode='Serial', state=state,
                         force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Proximity', state=state,
                          nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    # No applied stress. The junction forms from the mutual elastic attraction
    # of the two lines and from their line tension alone, which is what makes
    # this a clean test of the elastic interaction: there is no external
    # driving force to mask an error in it.
    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, topology=topology,
                          collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=max_step,
                          loading_mode="stress",
                          applied_stress=np.array(
                              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                          print_freq=print_freq, plot_freq=10,
                          plot_pause_seconds=0.1,
                          write_freq=write_freq, write_dir='output',
                          save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plot', dest='plot', action='store_false',
                        default=True)
    parser.add_argument('--max-step', dest='max_step', type=int, default=200,
                        help='steps to run; the junction forms at step 2 and '
                             'grows for the rest of the run')
    parser.add_argument('--dt', dest='dt', type=float, default=1.0e-9)
    parser.add_argument('--print-freq', dest='print_freq', type=int,
                        default=10,
                        help='steps between progress lines')
    parser.add_argument('--write-freq', dest='write_freq', type=int,
                        default=10,
                        help='steps between intermediate configuration '
                             'dumps under output/')
    parser.add_argument('--write-ref', dest='write_ref', action='store_true',
                        default=False,
                        help='also save the final configuration as a '
                             'reference .npz under output/')
    args = parser.parse_args()

    main(plot=args.plot, max_step=args.max_step, dt=args.dt,
         print_freq=args.print_freq, write_freq=args.write_freq)

    # explore the network after simulation
    G = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/binary_junction_pydis_elast_final.json')

    if args.write_ref:
        write_ref_npz(net, 'output/binary_junction_elast_ref.npz',
                      source='pydis')
