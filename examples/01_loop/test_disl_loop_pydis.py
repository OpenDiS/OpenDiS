import numpy as np
import sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib'),os.path.abspath('../../core/pydis/python')])

from framework.disnet_manager import DisNetManager
from pydis.disnet import DisNode, Cell, CellList, DisNet
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.topology.topology_disnet import Topology
from pydis.collision.collision_disnet import Collision
from pydis.remesh.remesh_disnet import Remesh
from pydis.visualize.vis_disnet import VisualizeNetwork
from pydis.simulate.sim_disnet import SimulateNetwork

def init_circular_loop(radius=1.0, N=20, burg_vec=np.array([1.0,0.0,0.0]),pbc=False):
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    theta = np.arange(N)*2.0*np.pi/N
    rn    = np.vstack([radius*np.cos(theta), radius*np.sin(theta), np.zeros_like(theta)]).T
    links = np.zeros((N, 5))
    for i in range(N):
        links[i,:] = np.array([i, (i+1)%N, burg_vec[0], burg_vec[1], burg_vec[2]])
    cell = Cell(h=8.0*np.eye(3), is_periodic=[True,True,True])
    cell_list = CellList(cell=cell, n_div=[8,8,8])
    G = DisNetManager(DisNet(cell=cell, cell_list=cell_list, rn=rn, links=links))
    return G

def main():
    global G, sim
    G = init_circular_loop()

    bounds = np.array([-0.5*np.diag(G.cell.h), 0.5*np.diag(G.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)
    # for debugging purposes
    #vis.plot_disnet(G)

    params = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.0316}
    calforce = CalForce(force_mode='LineTension', params=params, Ec=1.0e6)
    applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -2.0e6, 0.0])
    nodeforce_dict = calforce.NodeForce(G, applied_stress=applied_stress)
    # lt_force_array = np.array([nodeforce_dict[node] for node in G.nodes])
    node_force = [array for array in nodeforce_dict[0].values()] # force at nodes
    seg_force = [array for array in nodeforce_dict[1].values()]
    # print(np.array(a1), np.array(a2))
    # for debugging purposes
    # nodeforce_dict = calforce.NodeForce(G)
    # lt_force_array = np.array([nodeforce_dict[node] for node in G.nodes])
    # np.savetxt("force.dat", lt_force_array)
    # print("save force to 'force.dat'")

    mobility  = MobilityLaw(mobility_law='Relax')
    graph = G.get_disnet()
    node_vel = mobility.Mobility_Relax(graph, nodeforce_dict=node_force)
    timeint   = TimeIntegration(integrator='EulerForward')
    topology  = Topology(split_mode='MaxDiss')
    collision = None #Collision(collision_mode='Proximity')
    remesh    = None #Remesh(remesh_rule='LengthBased')

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          max_step=200, loading_mode="stress",
                          applied_stress=applied_stress,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='output')
    
    sim.run(G)
    
    return G.is_sane()


if __name__ == "__main__":
    main()

    graph = G.get_disnet()
