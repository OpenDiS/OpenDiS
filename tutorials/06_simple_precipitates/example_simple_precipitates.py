import os, sys
import numpy as np

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetworkPerf, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh, CrossSlip
except ImportError:
    raise ImportError('Cannot import pyexadis')


def init_line_config(Lbox, burg, plane, theta, maxseg):
    """init_line_config
    Initialize an infinite straight line along the y-direction
    with glide plane normal oriented along the z-direction.
    Returns the line configuration and appropriate crystal 
    orientation matrix.
    """
    from pyexadis_utils import insert_infinite_line
    
    # Crystal orientation
    plane = plane / np.linalg.norm(plane)
    b = burg / np.linalg.norm(burg)
    y = np.cross(plane, b)
    y = y / np.linalg.norm(y)
    linedir = np.cos(theta*np.pi/180.0)*b+np.sin(theta*np.pi/180.0)*y
    x = np.cross(linedir, plane)
    x = x / np.linalg.norm(x)
    Rorient = np.array([x, linedir, plane])
    
    burg = np.matmul(Rorient, burg)
    plane = np.matmul(Rorient, plane)
    linedir = -np.matmul(Rorient, linedir)
    
    # Create the cell and insert the line
    cell = pyexadis.Cell(Lbox)
    nodes, segs = [], []
    origin = np.array(cell.origin) + 0.5*np.sum(np.array(cell.h), axis=0)
    nodes, segs = insert_infinite_line(cell, nodes, segs, burg, plane, origin, linedir=linedir, maxseg=maxseg)
    
    N = DisNetManager(ExaDisNet(cell, nodes, segs))
    return N, Rorient


class Precipitate:
    """Precipitate
    Class to hold precipitate object
    """
    def __init__(self, pos: np.array, radius: float):
        self.pos = np.array(pos)
        self.radius = radius


import matplotlib.pyplot as plt
from pyexadis_base import VisualizeNetwork
class VisualizeNetworkPrecips(VisualizeNetwork):
    """VisualizeNetworkPrecips
    User-defined visualizer to plot dislocations and precipitates
    """
    def __init__(self, *args, **kwargs) -> None:
        super(VisualizeNetworkPrecips, self).__init__(*args, **kwargs)
        self.precips = kwargs.get("precips")
        
    def plot_sphere(self, ax, pos, radius):
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 20)
        x = pos[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = pos[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = pos[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, zorder=3.5, alpha=1.0)
    
    def plot_disnet(self, N: DisNetManager, state: dict={},
                    plot_nodes=True, plot_segs=True, plot_cell=True, trim=False,
                    fig=None, ax=None, block=False, pause_seconds=0.01):
        # plot dislocations using base class
        fig, ax = super().plot_disnet(N, state, plot_nodes, plot_segs, plot_cell,
                                      trim, fig, ax, block, pause_seconds=-1)
        # plot precipitates
        ax.computed_zorder = False
        for p in self.precips:
            self.plot_sphere(ax, p.pos, p.radius)  
        # display plot
        plt.draw()
        plt.show(block=block)
        plt.pause(pause_seconds)


from pyexadis_base import SimulateNetwork, NodeConstraints
class SimulationDriver(SimulateNetwork):
    """SimulationDriver
    User-defined simulation driver that adds the precipitate-dislocation interaction
    """
    def __init__(self, *args, **kwargs) -> None:
        super(SimulationDriver, self).__init__(*args, **kwargs)
        self.precips = kwargs.get("precips")
    
    def precipitates_intersection(self, N: DisNetManager, state: dict, xold):
        """precipitates_intersection:
        Adjust node positions of segments entering a precipitate
        """
        # Detect nodes that are inside precipitates and move
        # them to the precipitate surface
        data = N.export_data()
        pos = data["nodes"]["positions"]
        cons = data["nodes"]["constraints"]
        
        def sphere_intersect(ro, rd, ce, ra):
            oc = ro - ce
            b = np.einsum('ij,ij->i', oc, rd)
            qc = oc - b[:,None]*rd
            h = ra**2 - np.einsum('ij,ij->i', qc, qc)
            x, y = -1.0, -1.0
            ind = (h >= 0.0)
            h[ind] = np.sqrt(h[ind])
            x, y = -b-h, -b+h
            x[~ind], y[~ind] = -1, -1    
            return x, y
        
        for p in self.precips:
            dpos = xold-pos
            dposn = np.linalg.norm(dpos, axis=1)
            ind = np.argwhere((np.sum((pos-p.pos)**2, axis=1) < p.radius**2) & (dposn > 1e-10)).ravel()
            d = dpos[ind]/dposn[ind,None]
            x, y = sphere_intersect(pos[ind], d, p.pos, p.radius)
            ind = ind[x<0]
            pos[ind] += y[x<0,None]*d
            #cons[ind] = NodeConstraints.PINNED_NODE
        
        # Commit changes
        N.get_disnet(ExaDisNet).import_data(data)
    
    # Override step_begin to store beginning-of-step node positions
    def step_begin(self, N: DisNetManager, state: dict):
        """step_begin: invoked at the begining of each time step
        """
        self.xold = N.export_data()["nodes"]["positions"].copy()
    
    # Override step_post_integrate to handle dislocation-precipitate intersections
    def step_post_integrate(self, N: DisNetManager, state: dict):
        """step_post_integrate: invoked after time-integration of each time step
        """
        self.precipitates_intersection(N, state, self.xold)
    

def simple_precipitates():
    """
    Main simulation script to simulate dislocation interaction
    with unshearable precipitates
    """
    
    pyexadis.initialize()
    
    state = {
        "crystal": 'fcc',
        "burgmag": 2.50e-10,
        "mu": 50.0e9,
        "nu": 0.3,
        "a": 1.0,
        "maxseg": 2000.0,
        "minseg": 200.0,
        "rtol": 1.0,
        "rann": 2.0,
        "nextdt": 1e-11,
        "maxdt": 1e-10,
    }
    
    # Initial dislocation configuration
    L0 = 50000.0
    burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.])
    plane = np.array([1.,-1.,1.])
    N, state["Rorient"] = init_line_config(L0, burg, plane, 90.0, state["maxseg"])
    
    # Precipitates configuration
    if 0:
        # single precipitate
        precips = [Precipitate(N.cell.center()+L0*np.array([0.05,0.0,0.0]), 2000.0)]
    else:
        # distribution of precipitates
        Nprecip = 10 # number of precipitates
        Rprecip = [500., 2000.] # radius range
        np.random.seed(1234)
        ppos = L0*np.random.rand(Nprecip, 3)
        ppos[:,2] = 0.5*L0
        rp = Rprecip[0]+(Rprecip[1]-Rprecip[0])*np.random.rand(Nprecip)
        precips = [Precipitate(p, r) for p, r in zip(ppos, rp)]
    
    #vis = None
    vis = VisualizeNetworkPrecips(precips=precips)
    
    applied_stress = np.array([0.0, 0.0, 0.0, 0.0, 1e8, 0.0]) # xx,yy,zz,yz,xz,xy in Pa
    
    calforce   = CalForce(force_mode='CUTOFF_MODEL', state=state, cutoff=0.1*L0)
    mobility   = MobilityLaw(mobility_law='FCC_0', state=state, Medge=60000.0, Mscrew=60000.0, vmax=3400.0)
    timeint    = TimeIntegration(integrator='Trapezoid', state=state, force=calforce, mobility=mobility)
    collision  = Collision(collision_mode='Retroactive', state=state)
    topology   = Topology(topology_mode='TopologySerial', state=state, force=calforce, mobility=mobility)
    remesh     = Remesh(remesh_rule='LengthBased', state=state)
    cross_slip = CrossSlip(cross_slip_mode='ForceBasedSerial', state=state, force=calforce)
    
    sim = SimulationDriver(
        calforce=calforce, mobility=mobility, timeint=timeint, collision=collision,
        topology=topology, remesh=remesh, cross_slip=cross_slip, vis=vis,
        loading_mode='stress', applied_stress=applied_stress,
        max_step=200, burgmag=state["burgmag"], state=state,
        print_freq=1, plot_freq=10, plot_pause_seconds=0.0001,
        write_freq=10, write_dir='output_precip',
        precips=precips
    )
    sim.run(N, state)
    
    plt.show()
    
    if not sys.flags.interactive:
        pyexadis.finalize()


if __name__ == "__main__":
    simple_precipitates()
