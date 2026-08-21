import numpy as np
import sys, os

from pathlib import Path
# this script lives in examples/02_frank_read_src/, so the repository root is 2
# levels up
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

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0,
                             burg_vec=np.array([1.0,0.0,0.0]),
                             pbc=False, maxseg=None):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    PINNED = DisNode.Constraints.PINNED_NODE
    FREE   = DisNode.Constraints.UNCONSTRAINED
    rn = np.array([[0.0, -arm_length/2.0, 0.0,         PINNED],
                   [0.0,  0.0,            0.0,         FREE],
                   [0.0,  arm_length/2.0, 0.0,         PINNED],
                   [0.0,  arm_length/2.0, -arm_length, PINNED],
                   [0.0, -arm_length/2.0, -arm_length, PINNED]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    # Condition the configuration so that no segment exceeds maxseg. The back
    # edges of the source run between pinned corners and LengthBased remesh
    # never splits those, so without this they would stay at their initial
    # length for the whole run and break the cutoff + maxseg bound the neighbor
    # search relies on. Inserted nodes are PINNED here, leaving geometry and
    # boundary condition unchanged.
    if maxseg is not None:
        rn, links = remesh_initial_config(rn, links, maxseg, cell.h, [pbc]*3)

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))

def main(plot=True, max_step=200, print_freq=10, write_freq=10):
    global net, sim, state

    Lbox = 1000.0
    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
             "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}
    cutoff = 0.25*Lbox
    check_cutoff_maxseg(Lbox*np.eye(3), cutoff, state["maxseg"])

    net = init_frank_read_src_loop(box_length=Lbox,
                                   arm_length=0.125*Lbox, pbc=True,
                                   maxseg=state["maxseg"])
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    if plot:
        try:
            vis = VisualizeNetwork()
        except:
            print("")
            print("Failed to create VisualizeNetwork object")
            print("Try run with option  --no-plot")
            print("")
            raise
    else:
        vis = None

    # Elasticity_SBA includes the full segment-segment elastic interaction
    # (plus the i==j self term, regularized by the core radius state["a"]), in
    # contrast to the LineTension mode used in test_frank_read_src_pydis.py.
    # Note: this is an O(Nseg^2) double loop in python, so it is much slower
    # than LineTension. The cutoff is stated explicitly and matched by the
    # companion exadis run, so the two truncate the segment-segment interaction
    # identically and can be compared.
    #
    # It must satisfy  cutoff + maxseg <= d/3  (d = cell width), checked above.
    # ExaDiS bins segments by mid-point and clamps the bin count up to a
    # minimum of 3; past that limit the bins are smaller than the search radius
    # and the +-1 bin scan compares mid-points with a periodic shift quantized
    # to the bin grid instead of the true minimum image, silently dropping
    # pairs. The companion requirement is that no segment exceeds maxseg, which
    # the remesh in the init function above guarantees. LengthBased remesh
    # does not, since it never splits a segment pinned at both ends.
    # Ec=0.0 switches off the core term Elasticity_SBA now adds
    # (selfforcevec_LineTension, matching ParaDiS SelfForceIsotropic(coreOnly=0)).
    # It is off so this run stays comparable with the exadis one, which is given
    # Ec=0.0 for the same reason. Without it pydis would use the ParaDiS default,
    # mu/(4*pi)*log(a/0.1), and the two codes would no longer agree.
    calforce  = CalForce(force_mode='Elasticity_SBA', state=state,
                         Ec=0.0, cutoff=cutoff)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8,
                                state=state)
    # KNOWN LIMITATION: this example only runs to step 269.
    # Topology(split_mode='MaxDiss') calls OneNodeForce for any node with 4 or
    # more arms, and OneNodeForce is not implemented for the Elasticity_* force
    # modes. The bowing source stays free of such nodes until the expanding
    # loop collides with itself; the first 4-arm node appears at step 270, and
    # the run then stops with NotImplementedError: OneNodeForce_Elasticity_SBA
    # not implemented yet Hence the default max_step=200 below. Topology cannot
    # simply be disabled as a workaround: the Proximity collision handler reads
    # the nodeflag_dict that only Topology.init_topology_exemptions creates, so
    # topology=None fails earlier still with KeyError: 'nodeflag_dict'. To do:
    # implement OneNodeForce_Elasticity_SBA (pydis/calforce/calforce_disnet.py)
    # and raise max_step here. This test case should not be considered closed
    # until then.
    topology  = Topology(split_mode='MaxDiss', state=state,
                         force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Proximity', state=state,
                          nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, topology=topology,
                          collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=max_step,
                          loading_mode="stress",
                          applied_stress=np.array(
                              [0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=print_freq, plot_freq=10,
                          plot_pause_seconds=0.01,
                          write_freq=write_freq, write_dir='output',
                          save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plot', dest='plot', action='store_false',
                        default=True)
    # max_step > 269 currently fails, see the note on Topology in main()
    parser.add_argument('--max-step', dest='max_step', type=int, default=200)
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

    main(plot=args.plot, max_step=args.max_step,
         print_freq=args.print_freq, write_freq=args.write_freq)

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_pydis_elast_final.json')

    if args.write_ref:
        write_ref_npz(net, 'output/frank_read_src_elast_ref.npz',
                      source='pydis')
