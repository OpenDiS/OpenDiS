import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork
from pydis.util.network_util import insert_frank_read_src

def init_fcc_junction(junction_type, disloc_length, phi1, phi2):
    '''Generate two Frank-Read sources on two fcc slip systems that cross at the box center
    '''
    Lbox = 1.0*disloc_length
    cell = Cell(h=Lbox*np.eye(3), is_periodic=[1,1,1])
    nodes, segs = [], []

    b1 = 1.0/np.sqrt(2.0)*np.array([0.0, 1.0, 1.0])
    p1 = np.array([1.0, 1.0, -1.0])

    if junction_type == 'lomer':
        b2 = 1.0/np.sqrt(2.0)*np.array([1.0, 0.0,-1.0])
        p2 = np.array([1.0, 1.0, 1.0])
    elif junction_type == 'glissile':
        b2 = 1.0/np.sqrt(2.0)*np.array([1.0, 0.0,-1.0])
        p2 = np.array([1.0,-1.0, 1.0])
    elif junction_type == 'hirth':
        b2 = 1.0/np.sqrt(2.0)*np.array([0.0,-1.0, 1.0])
        p2 = np.array([1.0, 1.0, 1.0])
    elif junction_type == 'colinear':
        b2 = 1.0/np.sqrt(2.0)*np.array([0.0,-1.0,-1.0])
        p2 = np.array([1.0,-1.0, 1.0])
    else:
        raise ValueError('Unknown junction type %s'%junction_type)

    linter = np.cross(p1, p2)
    linter = linter / np.linalg.norm(linter)
    y1 = np.cross(linter, p1)
    y1 = y1 / np.linalg.norm(y1)
    ldir1 = np.cos(phi1*np.pi/180.0)*linter+np.sin(phi1*np.pi/180.0)*y1
    y2 = np.cross(linter, p2)
    y2 = y2 / np.linalg.norm(y2)
    ldir2 = np.cos(phi2*np.pi/180.0)*linter+np.sin(phi2*np.pi/180.0)*y2

    center = 0.5*Lbox * np.ones(3)
    delta = 0.01 * disloc_length * np.array([1.0, 1.0, 0.0])

    nodes, segs = insert_frank_read_src(cell, nodes, segs, b1, p1, disloc_length, center+delta, linedir=ldir1)
    nodes, segs = insert_frank_read_src(cell, nodes, segs, b2, p2, disloc_length, center-delta, linedir=ldir2)

    return DisNetManager(DisNet(cell=cell, rn=np.array(nodes), links=np.array(segs)))

def main():
    global net, sim, state

    junction_type = 'lomer'
    #junction_type = 'glissile'
    #junction_type = 'hirth'
    #junction_type = 'colinear'

    disloc_length = 1000.0 # length of dislocation lines
    phi1 = 30.0 # angle (deg) of the first dislocation with intersection direction
    phi2 = 30.0 # angle (deg) of the second dislocation with intersection direction
    net = init_fcc_junction(junction_type, disloc_length, phi1, phi2)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    state = {
        "crystal": 'fcc',
        "burgmag": 2.55e-10,
        "mu": 54.6e9,
        "nu": 0.324,
        "a": 6.0,
        "maxseg": 0.10*disloc_length,
        "minseg": 0.02*disloc_length,
        "rann": 10.0,
        "rtol": 5.0,
        "nextdt": 1e-12,
        "maxdt": 1e-10,
        "use_glide_planes": True,
        # ExaDiS' defaults, so that the two codes run the same simulation
        "splitMultiNodeAlpha": 1e-3,
        "coarsen_mode": 1,
    }

    calforce  = CalForce(force_mode='Elasticity_SBA', state=state, cutoff=300.0)
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint   = TimeIntegration(integrator='Trapezoid', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state, nbrlist=nbrlist)
    topology  = Topology(split_mode='Serial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          state=state, max_step=200, loading_mode='stress',
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                          print_freq=1, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=10, write_dir='output')
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/fcc_junctions_pydis_final.json')
