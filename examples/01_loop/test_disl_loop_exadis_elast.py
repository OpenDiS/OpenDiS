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

def init_circular_loop(radius=100.0, N=20, box_length=1000.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=True):
    '''Generate a circular dislocation loop lying in the z = box_length/2 plane
    '''
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    theta = np.arange(N)*2.0*np.pi/N
    rn = np.vstack([radius*np.cos(theta), radius*np.sin(theta), np.zeros_like(theta)]).T
    rn += cell.center()
    rn = np.hstack([rn, np.full((N,1), NodeConstraints.UNCONSTRAINED)])

    # The loop is planar and burg_vec lies in its plane, so every segment glides on the
    # same (0,0,1) plane. Note this cannot be obtained as cross(burg_vec, t) the way it is
    # in 02_frank_read_src: that expression vanishes at the two screw orientations.
    links = np.zeros((N, 8))
    for i in range(N):
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, [0.0, 0.0, 1.0]))

    return DisNetManager(ExaDisNet(cell, rn, links))

def main(plot=True, force_mode='DDD_FFT_MODEL', max_step=200, dt=1.0e-9):
    global net, sim, state

    Lbox = 1000.0
    net = init_circular_loop(radius=0.1*Lbox, box_length=Lbox)

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

    # Geometry and discretization follow 02_frank_read_src (box 1000, core radius a = 1,
    # maxseg/minseg of that order) rather than test_disl_loop_exadis.py, which uses a loop of
    # radius 1 with a = 0.01 in a box of 10. That scale only works with the softened
    # Ec = 1.0e6 line tension used there: with real elasticity the self-force of a radius-1
    # loop is ~2e10 and no attainable applied stress can balance it (the critical stress
    # would be ~0.4*mu), so the loop collapses and annihilates almost immediately. At radius
    # 100 the critical stress is ~6.5e8 Pa, which is what the applied stress below is set
    # against. The material constants mu and nu are unchanged from test_disl_loop_exadis.py.
    state = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 1.0, "maxseg": 60.0, "minseg": 20.0, "rann": 3.0}

    # Full elastic interactions, in contrast to the LineTension mode used in
    # test_disl_loop_exadis.py:
    #  - DDD_FFT_MODEL: short-range segment-segment pairs + long-range FFT (handles PBC)
    #  - CUTOFF_MODEL:  segment-segment pairs truncated at 'cutoff' (no long-range part)
    #
    # CUTOFF_MODEL is the mode to use when comparing against test_disl_loop_pydis_elast.py:
    # pydis' Elasticity_SBA applies the minimum image convention (cell.closest_image) and
    # sums no periodic images beyond that, which is what a truncated pair sum does. After 200
    # steps the two agree to 3 decimals (mean loop radius 131.576 in both). DDD_FFT_MODEL
    # adds the long-range image contribution pydis omits and expands faster (143.0); that is
    # a real physical difference between the two force models, not a discrepancy.
    #
    # Ec=0.0 disables the core energy contribution of FORCE_CORE_SELF_PKEXT, which would
    # otherwise default to mu/(4*pi)*log(a/0.1). As in 02_frank_read_src, this is done only
    # so that this run can be compared directly against test_disl_loop_pydis_elast.py:
    # pydis' Elasticity_SBA sums the segment-segment forces (including the i==j self term
    # regularized by 'a') but adds no core term, so it has no Ec equivalent. It is also not
    # the Ec=1.0e6 of test_disl_loop_exadis.py, which is a softened line tension standing in
    # for the elastic self-force that is now computed explicitly.
    # To do: add the Ec core term to pydis' Elasticity_SBA, as in the legacy ParaDiS code,
    # and then drop Ec=0.0 from this example.
    Ec = 0.0
    if force_mode == 'DDD_FFT_MODEL':
        calforce = CalForce(force_mode='DDD_FFT_MODEL', state=state, Ec=Ec, Ngrid=32, cell=net.cell)
    elif force_mode == 'CUTOFF_MODEL':
        calforce = CalForce(force_mode='CUTOFF_MODEL', state=state, Ec=Ec, cutoff=0.5*Lbox)
    else:
        raise ValueError('Unsupported force_mode %s for this example' % force_mode)

    # SimpleGlide matches 02_frank_read_src and test_disl_loop_exadis.py. Unlike pydis'
    # SimpleGlide it applies no vmax cap, which is why the pydis counterpart has to raise
    # vmax explicitly for the two runs to agree.
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=dt, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = None
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          state=state, max_step=max_step, loading_mode='stress',
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -1.0e9, 0.0]),
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
    parser.add_argument('--dt', dest='dt', type=float, default=1.0e-9)
    args = parser.parse_args()

    main(plot=args.plot, force_mode=args.force_mode, max_step=args.max_step, dt=args.dt)

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/disl_loop_exadis_elast_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
