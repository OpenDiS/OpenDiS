import os, sys
import numpy as np

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetwork, NodeConstraints, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
    from pyexadis_utils import insert_frank_read_src
except ImportError:
    raise ImportError('Cannot import pyexadis')
    
from user_stress import CalForceUserStress


if __name__ == "__main__":
    
    pyexadis.initialize()
    
    state = {
        "crystal": 'fcc',
        "burgmag": 2.55e-10,
        "mu": 82.0e9,
        "nu": 0.3,
        "a": 2.0,
        "maxseg": 400.0,
        "minseg": 80.0,
        "rtol": 1.0,
        "rann": 1.0,
        "nextdt": 5e-13,
    }
    
    # Insert a Frank-Read source (dislocation pinned at both ends)
    Lbox = 8000.0
    Lsrc = 0.3*Lbox
    cell = pyexadis.Cell(Lbox)
    nodes, segs = [], []
    burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.])
    plane = np.array([-1.,1.,1.])
    center = cell.center()
    nodes, segs = insert_frank_read_src(cell, nodes, segs, burg, plane, Lsrc, center, theta=0.0)
    G = ExaDisNet(cell, nodes, segs)
    N = DisNetManager(G)
    
    vis = VisualizeNetwork()
    
    # Define user-defined CalForce that adds a user-defined stress to a base force model
    calforce_base = CalForce(force_mode='LINE_TENSION_MODEL', state=state)
    calforce  = CalForceUserStress(calforce_base=calforce_base, state=state)
    
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=5.0e4, Mscrew=5.0e4, vmax=4000.0)
    timeint   = TimeIntegration(integrator='RKF', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(topology_mode='TopologySerial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    applied_stress = np.zeros(6)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, collision=collision, 
                          topology=topology, remesh=remesh, vis=vis,
                          loading_mode='stress', applied_stress=applied_stress,
                          max_step=1000, burgmag=state["burgmag"], state=state,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=100, write_dir='output_user_stress')
    
    sim.run(N, state)
    

    if not sys.flags.interactive:
        pyexadis.finalize()
