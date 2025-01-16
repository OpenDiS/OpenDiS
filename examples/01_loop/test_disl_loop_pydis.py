import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_circular_loop(radius=1.0, N=20, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    box_length = 10.0
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
    
    N = init_circular_loop()

    bounds = np.array([-0.5*np.diag(N.cell.h), 0.5*np.diag(N.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)
    # for debugging purposes
    #vis.plot_disnet(N)

    state = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.02}
    
    calforce  = CalForce(force_mode='LineTension', state=state, Ec=1.0e6)
    mobility  = MobilityLaw(mobility_law='Relax', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce, mobility=mobility)
    collision = None #Collision(collision_mode='Proximity', state=state)
    remesh    = None #Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          max_step=200, loading_mode="stress", state=state,
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -2.0e6, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='test_loop_results')
    
    sim.run(N, state)
    
    return N.is_sane()


if __name__ == "__main__":
    main()

    graph = N.get_disnet()
