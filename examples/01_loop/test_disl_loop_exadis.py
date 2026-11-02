import numpy as np
import sys, os

pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)

try:
    import pyexadis
    # from framework.disnet_manager import DisNetManager
    from pyexadis_base import ExaDisNet, NodeConstraints, DisNetManager, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

def init_circular_loop(box_length, radius, N=20, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    theta = np.arange(N)*2.0*np.pi/N
    rn    = np.vstack([radius*np.cos(theta) + box_length/2, radius*np.sin(theta) + box_length/2, np.zeros_like(theta) + box_length/2]).T
    links = np.zeros((N, 5))
    for i in range(N):
        links[i,:] = np.array([i, (i+1)%N, burg_vec[0], burg_vec[1], burg_vec[2]])
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    return DisNetManager(ExaDisNet(cell, rn, links))
    
def main():
    
    pyexadis.initialize()

    Lbox = 1000.0
    N = init_circular_loop(Lbox, radius=0.1*Lbox)
    
    state = {
        "burgmag": 3e-10,
        "mu": 50e9,
        "nu": 0.3,
        "a": 1.0,
        "maxseg": 0.05*Lbox,
        "minseg": 0.01*Lbox,
        "rann": 0.5,
    }

    vis = VisualizeNetwork()
    
    calforce  = CalForce(force_mode='LineTension', state=state)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1e-8, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = None #Topology(split_mode='TopologySerial', state=state, force=calforce, mobility=mobility)
    remesh    = None #Remesh(remesh_rule='LengthBased', state=state)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, 
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          max_step=200, loading_mode='stress', state=state,
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -200e6, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=10, write_dir='test_loop_results')
    sim.run(N, state)
    
    pyexadis.finalize()


if __name__ == "__main__":
    main()
