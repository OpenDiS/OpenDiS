import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_circular_loop(box_length, radius, N=20, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    theta = np.arange(N)*2.0*np.pi/N
    rn    = np.vstack([radius*np.cos(theta)+box_length/2, radius*np.sin(theta)+box_length/2, np.zeros_like(theta)+box_length/2]).T
    links = np.zeros((N, 8))
    for i in range(N):
        links[i,:] = np.array([i, (i+1)%N, burg_vec[0], burg_vec[1], burg_vec[2], 0, 0, 0])
    cell = Cell(h=box_length*np.eye(3), is_periodic=[True,True,True])
    N = DisNetManager(DisNet(cell=cell, rn=rn, links=links))
    return N

def main():
    global N, sim

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

    bounds = np.array([-0.5*np.diag(N.cell.h), 0.5*np.diag(N.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)
    # for debugging purposes
    #vis.plot_disnet(N)
    
    calforce  = CalForce(force_mode='LineTension', state=state)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1e-8, state=state)
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce, mobility=mobility)
    collision = None #Collision(collision_mode='Proximity', state=state)
    remesh    = None #Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          max_step=200, loading_mode="stress", state=state,
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -200e6, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='test_loop_results')
    
    sim.run(N, state)
    
    return N.is_sane()


if __name__ == "__main__":
    main()

    graph = N.get_disnet()
