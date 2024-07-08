import os, sys
import numpy as np

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


'''
Example of a function to generate an initial
Frank-Read source configuration
'''
def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    center = np.array(cell.center())
    
    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         NodeConstraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         NodeConstraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         NodeConstraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, NodeConstraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, NodeConstraints.PINNED_NODE]])
    rn[:,0:3] += center
    
    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(ExaDisNet(cell, rn, links))
    

'''
Example of a script to perform a simple Frank-Read source
simulation using the pyexadis binding to ExaDiS
'''
def main():
    global net, sim, state
    
    Lbox = 1000.0
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)

    vis = VisualizeNetwork()
    
    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}
    
    calforce  = CalForce(force_mode='LineTension', state=state)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = None
    remesh    = Remesh(remesh_rule='LengthBased', state=state)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, 
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          state=state, max_step=200, loading_mode='stress',
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=10, write_dir='output')
    sim.run(net, state)


if __name__ == "__main__":
    pyexadis.initialize()

    main()

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_exadis_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
