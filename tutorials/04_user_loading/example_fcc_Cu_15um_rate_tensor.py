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


class SimulationDriver(SimulateNetwork):
    """ User-defined simulation driver to impose a strain-rate tensor loading
    """
    
    def __init__(self, *args, **kwargs) -> None:
        super(SimulationDriver, self).__init__(*args, **kwargs)
        
        self.strain_rate_tensor = kwargs.get("strain_rate_tensor")
        state = kwargs.get("state")
        self.MU, self.NU = state["mu"], state["nu"]
        self.LA = 2*self.MU*self.NU/(1-2*self.NU)
    
    # Override step_update_response() function to apply strain-rate tensor loading
    def step_update_response(self, N: DisNetManager, state: dict):
        """step_update_response: update applied stress and rotation if needed
        """
        if self.loading_mode == 'strain_rate_tensor':
            
            # get values of plastic strain, plastic spin, and density computed internally in exadis
            dEp, dWp, state["density"] = N.get_disnet(ExaDisNet).net.get_plastic_strain()
            dEp = np.array(dEp).ravel()[[0,4,8,5,2,1]] # xx,yy,zz,yz,xz,xy
            dWp = np.array(dWp).ravel()[[5,2,1]] # yz,xz,xy
            state["dEp"] = dEp
            state["dWp"] = dWp
            
            # update strain and stress states based on strain rate tensor
            dE = self.strain_rate_tensor * state["dt"]
            dEe = dE - dEp # elastic strain
            dstress = self.LA*np.sum(dEe[0:3])*np.array([1,1,1,0,0,0]) + 2*self.MU*dEe
            
            # increment stress and strain tensors
            state["applied_stress"] += dstress
            state["Etot"] += dE
            
            def von_mises(T):
                S = np.array(T[[0,5,4,5,1,3,4,3,2]]).reshape(3,3)
                Sdev = S - np.trace(S)/3.0*np.eye(3)
                return np.sqrt(3.0/2.0*np.dot(Sdev.ravel(), Sdev.ravel()))
            
            # store strain and stress values used for the output (e.g. von Mises)
            state["strain"] = von_mises(state["Etot"])
            state["stress"] = von_mises(state["applied_stress"])
            
        else:
            # call base class function
            super().step_update_response(N, state)
            
        return state


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
    
    calforce  = CalForce(force_mode='SUBCYCLING_MODEL', state=state, Ngrid=64, cell=net.cell)
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint   = TimeIntegration(integrator='Subcycling', rgroups=[0.0, 100.0, 600.0, 1600.0], state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(topology_mode='TopologyParallel', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)
    cross_slip = None
    
    strain_rate_tensor = 1e3*np.array([1.,1.,1.,0.,0.,0.]) # xx,yy,zz,yz,xz,xy
    
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
