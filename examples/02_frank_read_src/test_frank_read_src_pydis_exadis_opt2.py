import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../core/exadis/python/']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce as PyDiS_CalForce, MobilityLaw as PyDiS_MobilityLaw
from pydis import TimeIntegration as PyDiS_TimeIntegration, Topology as PyDiS_Topology
from pydis import Collision as PyDiS_Collision, Remesh as PyDiS_Remesh
from pydis import VisualizeNetwork, SimulateNetwork

try:
    import pyexadis
    from pyexadis_base import ExaDisNet
    from pyexadis_base import CalForce as ExaDiS_CalForce, MobilityLaw as ExaDiS_MobilityLaw, Topology as ExaDiS_Topology
    from pyexadis_base import TimeIntegration as ExaDiS_TimeIntegration, Collision as ExaDiS_Collision, Remesh as ExaDiS_Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

# import the function to set up the frank-read source
from test_frank_read_src_pydis_exadis import init_frank_read_src_loop

def main():
    global net, sim, state

    Lbox = 1000.0
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    state = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0}

    pydis_calforce   = PyDiS_CalForce(force_mode='LineTension', state=state)
    exadis_calforce  = ExaDiS_CalForce(force_mode='LineTension', state=state)

    pydis_mobility   = PyDiS_MobilityLaw(mobility_law='SimpleGlide', state=state)
    exadis_mobility  = ExaDiS_MobilityLaw(mobility_law='SimpleGlide', state=state)

    pydis_timeint    = PyDiS_TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    exadis_timeint   = ExaDiS_TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state, force=pydis_calforce, mobility=exadis_mobility)

    pydis_topology   = PyDiS_Topology(split_mode='MaxDiss', state=state, force=exadis_calforce, mobility=pydis_mobility)
    exadis_topology  = ExaDiS_Topology(topology_mode='TopologySerial', state=state, force=pydis_calforce, mobility=exadis_mobility)

    pydis_collision  = PyDiS_Collision(collision_mode='Proximity', state=state, nbrlist=nbrlist)
    exadis_collision = ExaDiS_Collision(collision_mode='Retroactive', state=state)

    pydis_remesh     = PyDiS_Remesh(remesh_rule='LengthBased', state=state)
    exadis_remesh    = ExaDiS_Remesh(remesh_rule='LengthBased', state=state)


    sim = SimulateNetwork(calforce=exadis_calforce, mobility=pydis_mobility, timeint=exadis_timeint,
                          topology=exadis_topology, collision=exadis_collision, remesh=exadis_remesh, vis=vis,
                          state=state, max_step=200, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.01,
                          write_freq=10, write_dir='output', save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    pyexadis.initialize()
    main()

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_pydis_exadis_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
