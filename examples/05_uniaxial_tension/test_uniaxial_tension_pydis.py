import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths
 if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNet, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork

def main():
    global net, sim, state

    state = {"crystal": 'bcc', "burgmag": 2.85e-10, "mu": 55.0e9, "nu": 0.339,
             "a": 6.0, "maxseg": 1000.0, "minseg": 200.0, "rtol": 1.0,
             "rann": 1.0, "nextdt": 1e-12}

    Lbox = 10000.0
    G = DisNet()
    G.generate_prismatic_config(state["crystal"], Lbox, 4, 0.21*Lbox,
                                state["maxseg"], uniform=True, seed=123)
    net = DisNetManager(G)

    vis = VisualizeNetwork()

    calforce  = CalForce(force_mode='LineTension', state=state)
    mobility  = MobilityLaw(mobility_law='BCC_0B', state=state, Medge=2600.0,
                            Mscrew=2600.0, Mclimb=1e-4, vmax=3400.0)
    timeint   = TimeIntegration(integrator='Trapezoid', multi=2, state=state,
                                force=calforce, mobility=mobility)
    nbrlist   = CellList(cell=net.cell, n_div=[8,8,8])
    collision = Collision(collision_mode='Retroactive', state=state,
                          nbrlist=nbrlist)
    topology  = Topology(split_mode='Serial', state=state,
                         force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                          timeint=timeint, topology=topology,
                          collision=collision, remesh=remesh, vis=vis,
                          loading_mode='strain_rate', erate=3e3,
                          edir=np.array([0.,0.,1.]), max_step=50,
                          burgmag=state["burgmag"], state=state,
                          print_freq=1, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output', save_state=False)
    sim.run(net, state)
    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/uniaxial_tension_pydis_final.json')
