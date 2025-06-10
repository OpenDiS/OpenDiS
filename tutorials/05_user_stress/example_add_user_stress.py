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


"""
User-defined additional stress field for the simulation.
As an example, here the stress component xx (in units of Pa) is
defined as a toy function of the z-position and simulation time.
The function must return an array of size (size(x),6), where R is
the array of query positions, and the 6 stress components are in
Voigt order (xx,yy,zz,yz,xz,xy).
"""
def UserStress(R, cell, state):
    Lz = cell.h[2][2]
    R = np.atleast_2d(R)
    z = R[:,2] # z-coordinate
    t = state["time"] if "time" in state else 0.0 # time
    
    # toy function for the Sxx component as function of z and t
    Sxx = (1.0 + np.cos(z/Lz*2.0*np.pi) + np.sin(t/1e-9)) * 5e7 # Pa
    
    user_stress = np.zeros((R.shape[0],6)) # xx,yy,zz,yz,xz,xy in Pa
    user_stress[:,0] = Sxx
    return user_stress


"""
User-defined CalForce function to add a spatially and temporally
varying user-defined stress contribution to a base CalForce model.
"""
class CalForceUserStress:
    def __init__(self, state: dict, calforce_base: CalForce, **kwargs) -> None:
        self.calforce_base = calforce_base
        self.add_user_stress = 1
        
    def PreCompute(self, N: DisNetManager, state: dict) -> dict:
        self.calforce_base.PreCompute(N, state)
        return state
    
    def NodeForce(self, N: DisNetManager, state: dict, pre_compute=True) -> dict:
        
        # Compute base forces
        self.calforce_base.NodeForce(N, state, pre_compute)
        
        # Add user-defined stress contribution
        if self.add_user_stress:
            self.AddUserStress(N, state)
        
        return state
    
    def AddUserStress(self, N: DisNetManager, state: dict) -> dict:
        
        G = N.get_disnet(ExaDisNet)
        cell = G.cell
        rn = G.get_nodes_data()["positions"]
        segs = G.get_segs_data()
        segsnid = segs["nodeids"] # segment connectivity to nodes
        
        # segments end-nodes positions and Burgers vectors
        r1 = np.array(cell.closest_image(Rref=np.array(cell.center()), R=rn[segsnid[:,0]]))
        r2 = np.array(cell.closest_image(Rref=r1, R=rn[segsnid[:,1]]))
        r_ij = r2-r1
        b_ij = segs["burgers"]
        
        # compute user stress at segment mid-positions
        Rseg = 0.5*(r1+r2)
        user_stress = UserStress(Rseg, cell, state)
        
        # add PK force from user stress to current force
        sig_user = user_stress[:,[0,5,4,5,1,3,4,3,2]].reshape(-1,3,3)
        sigb = np.einsum('kij,kj->ki', sig_user, b_ij)
        fpkuser = 0.5 * np.cross(sigb, r_ij)
        f = G.get_forces() # current nodal forces
        np.add.at(f, segsnid[:,0], fpkuser)
        np.add.at(f, segsnid[:,1], fpkuser)
        
        # store new forces in state dictionnary
        state["nodeforces"] = f
        
        return state
    
    def OneNodeForce(self, N: DisNetManager, state: dict, tag, update_state=True) -> np.array:
        
        if self.add_user_stress:
            # here let's just hack the code by temporarily incrementing the global 
            # applied_stress by the user_stress value at the node position.
            # this should be fine as we are only computing the force on that node
            
            applied_stress0 = 1.0*state["applied_stress"] # save applied_stress
            
            # find node position from node tag
            G = N.get_disnet(ExaDisNet)
            tags = G.get_tags()
            ind = np.where((tags[:,0]==tag[0])&(tags[:,1]==tag[1]))[0]
            if ind.size != 1:
                raise ValueError("Cannot find node tag (%d,%d) in OneNodeForce" % tuple(tag))
            pos = G.get_positions()[ind]
            
            time = state["time"] if "time" in state else 0.0
            user_stress = UserStress(G.cell, pos, time)
            state["applied_stress"] += user_stress.squeeze() # increment applied_stress by user_stress
            f = self.calforce_base.OneNodeForce(N, state, tag, update_state) # compute node force
            
            state["applied_stress"] = applied_stress0 # reset applied_stress
            
        else:
            f = self.calforce_base.OneNodeForce(N, state, tag, update_state) # compute node force
        
        return f


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
