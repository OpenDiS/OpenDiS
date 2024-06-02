import numpy as np
import sys, os
# from pydis.disnet import DisNode, Cell, CellList, DisNet
# from pydis.topology.topology_disnet import Topology
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)

# [sys.path.append(os.path.abspath(path)) for path in pyexadis_path if not path in sys.path]
# np.set_printoptions(threshold=20, edgeitems=5)

try:
    import pyexadis
    # from framework.disnet_manager import DisNetManager
    from pyexadis_base import ExaDisNet, NodeConstraints, DisNetManager, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

# for demonstrating how to convert data between pydis and pyexadis
# from pydis import DisNet

def init_circular_loop(arm_length=1.0, box_length=10, burg_vec=np.array([1.0,0.0,0.0]),pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    radius=1.0; nlinks = 20
    print("init_circular_loop: radius = %f, N = %d" % (radius, nlinks))
    theta = np.arange(nlinks)*2.0*np.pi/nlinks
    rn    = np.vstack([radius*np.cos(theta) + box_length/2, radius*np.sin(theta) + box_length/2, np.zeros_like(theta) + box_length/2]).T
    links = np.zeros((nlinks, 5))
    for i in range(nlinks):
        links[i,:] = np.array([i, (i+1)%nlinks, burg_vec[0], burg_vec[1], burg_vec[2]])
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    return ExaDisNet(cell, rn, links)
    
def main():
    
    pyexadis.initialize()
    
    G = init_circular_loop(pbc=True)
    N = DisNetManager({type(G): G})

    vis = VisualizeNetwork()
    
    params = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.02}
    
    calforce  = CalForce(force_mode='LineTension', params=params, Ec=1.0e6)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', params=params)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-9, params=params)
    collision = Collision(collision_mode='Retroactive', params=params)
    topology  = None
    remesh    = None#Remesh(remesh_rule='LengthBased', params=params)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, 
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          max_step=100, loading_mode='stress',
                          applied_stress=np.array([1.0e6, 1.0e6, 0.0, 0.0, 0.0, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=10, write_dir='output')
    sim.run(N)
    
    pyexadis.finalize()


if __name__ == "__main__":
    main()
