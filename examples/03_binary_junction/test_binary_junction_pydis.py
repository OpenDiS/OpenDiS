import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def init_two_disl_lines(z0=1.0, box_length=8.0, b1=np.array([-1.0,1.0,1.0]), b2=np.array([1.0,-1.0,1.0]), pbc=False):
    '''Generate an initial configuration for two dislocation lines.
    '''
    print("init_two_disl_lines: z0 = %f" % (z0))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    cell_list = CellList(cell=cell, n_div=[4,4,4])
    rn    = np.array([[0.0, -z0, -z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  z0,  z0,  DisNode.Constraints.PINNED_NODE],
                      [-z0,  0.0,-z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [ z0,  0.0, z0,  DisNode.Constraints.PINNED_NODE]])
    
    xi1, xi2 = rn[2,:3] - rn[1,:3], rn[5,:3] - rn[4,:3]
    n1, n2 = np.cross(b1, xi1), np.cross(b2, xi2)
    n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0, :] = np.concatenate(([0, 1], b1, n1))
    links[1, :] = np.concatenate(([1, 2], b1, n1))
    links[2, :] = np.concatenate(([3, 4], b2, n2))
    links[3, :] = np.concatenate(([4, 5], b2, n2))

    G = DisNetManager(DisNet(cell=cell, cell_list=cell_list, rn=rn, links=links))

    return G

def main():
    global net, sim, state

    Lbox = 8; z0 = 1
    net = init_two_disl_lines(z0=z0, box_length=Lbox, pbc=False)

    bounds = np.array([-0.5*np.diag(net.cell.h), 0.5*np.diag(net.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)

    state = {
        "burgmag": 3e-10, 
        "mu": 160e9, 
        "nu": 0.31, 
        "a": 0.1, 
        "maxseg": 0.04 * Lbox,
        "minseg": 0.01 * Lbox,
        "rann": 0.0316
    }

    calforce = CalForce(force_mode='LineTension', state=state, Ec=1.0e6)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=state)
    timeint = TimeIntegration(integrator='EulerForward', dt=1.0e-9, state=state)
    topology = Topology(split_mode='MaxDiss', state=state)
    collision = Collision(collision_mode='Proximity', state=state)
    remesh = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=200, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.1,
                          write_freq=10, write_dir='output')
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/two_disl_lines_final.json')