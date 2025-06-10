import os, sys
import numpy as np

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh, CrossSlip
except ImportError:
    raise ImportError('Cannot import pyexadis')

from user_rate_tensor import SimulationDriver


def example_fcc_Cu_15um_1e3_rate_tensor():
    """example_fcc_Cu_15um_1e3_rate_tensor
    """
    global net, state, sim
    
    pyexadis.initialize()
    
    state = {
        "crystal": 'fcc',
        "burgmag": 2.55e-10,
        "mu": 54.6e9,
        "nu": 0.324,
        "a": 6.0,
        "maxseg": 2000.0,
        "minseg": 300.0,
        "rtol": 10.0,
        "rann": 10.0,
        "nextdt": 1e-10,
        "maxdt": 1e-9,
    }
    
    output_dir = 'output_fcc_Cu_15um_rate_tensor'
    
    # Initial configuration
    Lbox = np.round(15e-6 / state["burgmag"]) # 15um
    G = ExaDisNet().generate_prismatic_config(state["crystal"], Lbox, num_loops=12, radius=0.2*Lbox, maxseg=state["maxseg"])
    net = DisNetManager(G)
    
    #vis = None
    vis = VisualizeNetwork()
    
    calforce  = CalForce(force_mode='CUTOFF_MODEL', state=state, cutoff=0.2*Lbox)
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint   = TimeIntegration(integrator='RKF', multi=10, state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(topology_mode='TopologyParallel', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)
    cross_slip = None
    
    strain_rate_tensor = 1e3*np.array([1.,1.,-1.,0.,0.,0.]) # xx,yy,zz,yz,xz,xy
    
    sim = SimulationDriver(calforce=calforce, mobility=mobility, timeint=timeint, collision=collision, 
                           topology=topology, remesh=remesh, cross_slip=cross_slip, vis=vis,
                           loading_mode='strain_rate_tensor', strain_rate_tensor=strain_rate_tensor,
                           max_step=100, burgmag=state["burgmag"], state=state,
                           print_freq=1, plot_freq=10, plot_pause_seconds=0.0001,
                           write_freq=100, write_dir=output_dir)
    sim.run(net, state)
    
    
    if not sys.flags.interactive:
        pyexadis.finalize()
    

if __name__ == "__main__":
    example_fcc_Cu_15um_1e3_rate_tensor()
