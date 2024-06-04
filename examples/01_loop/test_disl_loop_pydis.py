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
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))
    cell = Cell(h=8.0*np.eye(3), is_periodic=[True,True,True])
    cell_list = CellList(cell=cell, n_div=[8,8,8])
    G = DisNetManager(DisNet(cell=cell, cell_list=cell_list, rn=rn, links=links))
    return G

def main():
    global net, sim
    net = init_circular_loop(pbc=True)

    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    # "burgmag" not yet used in pydis, make parameters (Ec) more physical
    params = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.0316}
    calforce = CalForce(force_mode='LineTension', params=params, Ec=1.0e6)
    mobility  = MobilityLaw(mobility_law='GNN', params=params)#SimpleGlide
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-9, params=params)
    topology  = Topology(split_mode='MaxDiss', params=params)
    collision = None#Collision(collision_mode='Proximity', params=params)
    remesh    = None#Remesh(remesh_rule='LengthBased', params=params)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          max_step=100, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='output')
    sim.run(net)

    return net.is_sane()


if __name__ == "__main__":
    main()

    G = net.get_disnet()