import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    cell_list = CellList(cell=cell, n_div=[8,8,8])

    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE]])
    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(disnet=DisNet(cell=cell, cell_list=cell_list, rn=rn, links=links))

def main():
    global net, sim
    net = init_frank_read_src_loop(pbc=True)

    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    # "burgmag" not yet used in pydis, make parameters (Ec) more physical
    params = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.0316}
    calforce = CalForce(force_mode='LineTension', params=params, Ec=1.0e6)
    mobility  = MobilityLaw(mobility_law='SimpleGlide', params=params)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, params=params)
    topology  = Topology(split_mode='MaxDiss', params=params)
    collision = Collision(collision_mode='Proximity', params=params)
    remesh    = Remesh(remesh_rule='LengthBased', params=params)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          max_step=200, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e6, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='output')
    sim.run(net)

    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G  = net.get_disnet()
