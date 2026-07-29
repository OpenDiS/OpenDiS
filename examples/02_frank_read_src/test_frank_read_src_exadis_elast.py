import numpy as np
import sys, os

# Import pyexadis
pyexadis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../core/exadis/python/']
[sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

try:
    import pyexadis
    from framework.disnet_manager import DisNetManager
    from pyexadis_base import ExaDisNet, NodeConstraints, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         NodeConstraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         NodeConstraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         NodeConstraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, NodeConstraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, NodeConstraints.PINNED_NODE]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(ExaDisNet(cell, rn, links))

def main(plot=True, force_mode='DDD_FFT_MODEL', max_step=200):
    global net, sim, state

    Lbox = 1000.0
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)

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

    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}

    # Full elastic interactions, in contrast to the LineTension mode used in
    # test_frank_read_src_exadis.py:
    #  - DDD_FFT_MODEL: short-range segment-segment pairs + long-range FFT (handles PBC)
    #  - CUTOFF_MODEL:  segment-segment pairs truncated at 'cutoff' (no long-range part)
    if force_mode == 'DDD_FFT_MODEL':
        calforce = CalForce(force_mode='DDD_FFT_MODEL', state=state, Ngrid=32, cell=net.cell)
    elif force_mode == 'CUTOFF_MODEL':
        calforce = CalForce(force_mode='CUTOFF_MODEL', state=state, cutoff=0.5*Lbox)
    else:
        raise ValueError('Unsupported force_mode %s for this example' % force_mode)

    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = None
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          state=state, max_step=max_step, loading_mode='stress',
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output')
    sim.run(net, state)


if __name__ == "__main__":
    pyexadis.initialize()

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plot', dest='plot', action='store_false', default=True)
    parser.add_argument('--force-mode', dest='force_mode', type=str, default='DDD_FFT_MODEL',
                        choices=['DDD_FFT_MODEL', 'CUTOFF_MODEL'])
    parser.add_argument('--max-step', dest='max_step', type=int, default=200)
    args = parser.parse_args()

    main(plot=args.plot, force_mode=args.force_mode, max_step=args.max_step)

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_exadis_elast_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
