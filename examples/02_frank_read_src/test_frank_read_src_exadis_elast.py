import numpy as np
import sys, os

# Import pyexadis
from pathlib import Path
# this script lives in examples/02_frank_read_src/, so the repository root is 2
# levels up
opendis_root = Path(__file__).resolve().parents[2]
pyexadis_paths = [str(opendis_root / p)
                  for p in ['python', 'lib', 'core/pydis/python',
                            'core/exadis/python']]
[sys.path.append(path) for path in pyexadis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

try:
    import pyexadis
    from framework.disnet_manager import DisNetManager
    from framework.simulation_setup import check_cutoff_maxseg
    from framework.simulation_setup import remesh_initial_config
    from framework.testing import write_ref_npz
    from pyexadis_base import ExaDisNet, NodeConstraints
    from pyexadis_base import SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
    from pyexadis_base import Collision, Topology, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

# Offset of the free middle node along the source line, in b. It breaks the
# mirror symmetry of the initial configuration, which is what makes the first
# self-collision an exact tie; see init_frank_read_src_loop. Both example
# scripts must use the same value or the comparison is meaningless.
#
# The value is empirical and does not generalize. Measured first step at which
# pydis and exadis part, and the residual there:
#     0.0    step 260, 1.6200e+01      0.5   step 250, 7.4949e-05
#     0.001  step 260, 2.2247e+00      2.0   step 260, 4.3081e-05
#     0.05   step 263, 5.6169e+01      5.0   step 256, 5.1477e-05
# 0.05 is the only value tried that carries both codes through the first two
# collisions at round-off (9.17e-08 at step 260, 1.11e-07 at 262). That is not
# a physical threshold, it is an overlap of each code's own round-off, so do
# not expect it to survive a compiler, platform or code change on either side.
MID_OFFSET = 0.05

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0,
                             burg_vec=np.array([1.0,0.0,0.0]),
                             pbc=False, maxseg=None, mid_offset=0.0):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    PINNED = NodeConstraints.PINNED_NODE
    FREE   = NodeConstraints.UNCONSTRAINED
    rn = np.array([[0.0, -arm_length/2.0, 0.0,         PINNED],
                   [0.0,  0.0,            0.0,         FREE],
                   [0.0,  arm_length/2.0, 0.0,         PINNED],
                   [0.0,  arm_length/2.0, -arm_length, PINNED],
                   [0.0, -arm_length/2.0, -arm_length, PINNED]])
    # mid_offset slides the free middle node along the line, breaking the
    # mirror symmetry of the source about its own mid-plane. Without it the two
    # halves of the expanding loop reach that plane simultaneously and the
    # first collision is an exact tie between several equally valid segment
    # pairs, which PyDiS and ExaDiS then resolve differently.
    rn[1,1] += mid_offset
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    # Condition the configuration so that no segment exceeds maxseg. The back
    # edges run between pinned corners and LengthBased remesh never splits
    # those. Inserted nodes are PINNED, leaving geometry and boundary condition
    # unchanged.
    if maxseg is not None:
        rn, links = remesh_initial_config(rn, links, maxseg,
                                          box_length*np.eye(3), [pbc]*3)

    return DisNetManager(ExaDisNet(cell, rn, links))

def main(plot=True, force_mode='DDD_FFT_MODEL', max_step=200,
         print_freq=10, write_freq=10):
    global net, sim, state

    Lbox = 1000.0
    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
             "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}
    cutoff = 0.25*Lbox
    check_cutoff_maxseg(Lbox*np.eye(3), cutoff, state["maxseg"])

    net = init_frank_read_src_loop(box_length=Lbox,
                                   arm_length=0.125*Lbox, pbc=True,
                                   maxseg=state["maxseg"],
                                   mid_offset=MID_OFFSET)

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

    # Full elastic interactions, in contrast to the LineTension mode used in
    # test_frank_read_src_exadis.py:
    # - DDD_FFT_MODEL: short-range segment-segment pairs + long-range FFT
    # (handles PBC) - CUTOFF_MODEL:  segment-segment pairs truncated at
    # 'cutoff' (no long-range part)
    #
    # CUTOFF_MODEL is the mode to use when comparing against
    # test_frank_read_src_pydis_elast.py: pydis' Elasticity_SBA applies the
    # minimum image convention (cell.closest_image) and sums no periodic images
    # beyond that, which is what a truncated pair sum does. After 200 steps the
    # two reach a max bow-out of 447.3 and 447.1 respectively, while
    # DDD_FFT_MODEL reaches 474.0 because it adds the long-range image
    # contribution pydis omits. That is a real physical difference between the
    # two force models, not a discrepancy.
    #
    # Ec=0.0 disables the core energy contribution of FORCE_CORE_SELF_PKEXT,
    # which would otherwise default to mu/(4*pi)*log(a/0.1). This is done only
    # so that this run can be compared directly against
    # test_frank_read_src_pydis_elast.py: pydis' Elasticity_SBA sums the
    # segment-segment forces (including the i==j self term regularized by 'a')
    # but adds no core term, so it has no Ec equivalent. With Ec=0.0 the two
    # codes give identical nodal forces; with the exadis default they agree
    # only while the arms are collinear, and the extra core line tension holds
    # the source near its equilibrium bow-out instead of letting it expand (the
    # applied stress here sits right at the Frank-Read threshold mu*b/L = 4e8
    # Pa). To do: add the Ec core term to pydis' Elasticity_SBA, as in the
    # legacy ParaDiS code, and then drop Ec=0.0 from this example.
    #
    # The cutoff matches test_frank_read_src_pydis_elast.py so the two runs
    # truncate identically. It must satisfy  cutoff + maxseg <= d/3  (checked
    # in main above): ExaDiS bins segments by mid-point and clamps the bin
    # count up to a minimum of 3, and past that limit the +-1 bin scan compares
    # mid-points using a periodic shift quantized to the bin grid rather than
    # the true minimum image, silently dropping pairs whose real separation is
    # inside the cutoff. The companion requirement, that no segment is longer
    # than maxseg, is met by the remesh in the init function, since LengthBased
    # remesh never splits a segment pinned at both ends.
    Ec = 0.0
    if force_mode == 'DDD_FFT_MODEL':
        calforce = CalForce(force_mode='DDD_FFT_MODEL', state=state,
                            Ec=Ec, Ngrid=32, cell=net.cell)
    elif force_mode == 'CUTOFF_MODEL':
        calforce = CalForce(force_mode='CUTOFF_MODEL', state=state,
                            Ec=Ec, cutoff=cutoff)
    else:
        raise ValueError('Unsupported force_mode %s for this example'
                         % force_mode)

    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8,
                                state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    # 'TopologySerial' to match the companion pydis run's split_mode='Serial'.
    # Topology was None here while pydis could not split a multi-arm node under
    # an Elasticity_* force mode. The expanding loop first touches itself at
    # step 256, and from there on a run without splitting is not comparable.
    topology  = Topology(topology_mode='TopologySerial', state=state,
                         force=calforce, mobility=mobility)
    # coarsen_mode=0 is the segment-centric coarsening branch: a segment
    # shorter than minseg has its two endpoints merged to their mid-point. It
    # is the only branch pydis implements, and the only one verified against
    # exadis (tests/unit_tests/test3_remesh_rule, bitwise over 500 steps).
    # pyexadis_base.Remesh defaults to coarsen_mode=1, which is node-centric
    # instead: a two-arm node with either arm shorter than minseg is merged
    # into its nearer neighbour, at that neighbour's own position. Leaving the
    # default in place compares pydis against an algorithm it does not have,
    # and the two answers differ by half the short segment.
    remesh    = Remesh(remesh_rule='LengthBased', state=state, coarsen_mode=0)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, collision=collision,
                          topology=topology, remesh=remesh, vis=vis,
                          state=state, max_step=max_step,
                          loading_mode='stress',
                          applied_stress=np.array(
                              [0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=print_freq, plot_freq=10,
                          plot_pause_seconds=0.01,
                          write_freq=write_freq, write_dir='output')
    sim.run(net, state)


if __name__ == "__main__":
    pyexadis.initialize()

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plot', dest='plot', action='store_false',
                        default=True)
    parser.add_argument('--force-mode', dest='force_mode', type=str,
                        default='DDD_FFT_MODEL',
                        choices=['DDD_FFT_MODEL', 'CUTOFF_MODEL'])
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

    main(plot=args.plot, force_mode=args.force_mode,
         max_step=args.max_step, print_freq=args.print_freq,
         write_freq=args.write_freq)

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_exadis_elast_final.json')

    if args.write_ref:
        write_ref_npz(net, 'output/frank_read_src_elast_ref.npz',
                      source='exadis')

    if not sys.flags.interactive:
        pyexadis.finalize()
