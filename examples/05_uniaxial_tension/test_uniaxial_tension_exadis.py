import numpy as np
import sys, os

# Import pyexadis
pyexadis_paths = ['../../python', '../../lib', '../../core/pydis/python',
                  '../../core/exadis/python/']
[sys.path.append(os.path.abspath(path)) for path in pyexadis_paths
 if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

try:
    import pyexadis
    from framework.disnet_manager import DisNetManager
    from pyexadis_base import ExaDisNet, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
    from pyexadis_base import Collision, Topology, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

def main():
    global net, sim, state

    state = {"crystal": 'bcc', "burgmag": 2.85e-10, "mu": 55.0e9, "nu": 0.339,
             "a": 6.0, "maxseg": 1000.0, "minseg": 200.0, "rtol": 1.0,
             "rann": 1.0, "nextdt": 1e-12}

    Lbox = 10000.0
    G = ExaDisNet()
    G.generate_prismatic_config(state["crystal"], Lbox, 4, 0.21*Lbox,
                                state["maxseg"], uniform=True, seed=123)
    net = DisNetManager(G)

    vis = VisualizeNetwork()

    calforce  = CalForce(force_mode='LINE_TENSION_MODEL', state=state)
    mobility  = MobilityLaw(mobility_law='BCC_0B', state=state, Medge=2600.0,
                            Mscrew=2600.0, Mclimb=1e-4, vmax=3400.0)
    timeint   = TimeIntegration(integrator='Trapezoid', multi=2, state=state,
                                force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(topology_mode='TopologyParallel', state=state,
                         force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, topology=topology,
                          collision=collision, remesh=remesh, vis=vis,
                          loading_mode='strain_rate', erate=3e3,
                          edir=np.array([0.,0.,1.]), max_step=500,
                          burgmag=state["burgmag"], state=state,
                          print_freq=1, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output')
    sim.run(net, state)


if __name__ == "__main__":
    pyexadis.initialize()

    main()

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/uniaxial_tension_exadis_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
