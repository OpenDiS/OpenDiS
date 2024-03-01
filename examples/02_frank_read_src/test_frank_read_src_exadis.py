import os, sys
import numpy as np
from typing import Tuple

# Import pydis
sys.path.extend([
    os.path.abspath('../../python'),
    os.path.abspath('../../lib'),
    os.path.abspath('../../core/pydis/python')
])
try:
    from framework.disnet_manager import DisNetManager
    from pydis.disnet import DisNet, Cell, DisNode
    from pydis.calforce.calforce_disnet import CalForce
    from pydis.mobility.mobility_disnet import MobilityLaw
    from pydis.timeint.timeint_disnet import TimeIntegration
    from pydis.collision.collision_disnet import Collision
    from pydis.remesh.remesh_disnet import Remesh
    from pydis.visualize.vis_disnet import VisualizeNetwork
except ImportError:
    raise ImportError('Cannot import pydis')

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet
    from pyexadis_base import CalForce as CalForce1
    from pyexadis_base import MobilityLaw as MobilityLaw1
    from pyexadis_base import TimeIntegration as TimeIntegration1
    from pyexadis_base import Collision as Collision1
    from pyexadis_base import Remesh as Remesh1
except ImportError:
    raise ImportError('Cannot import pyexadis')

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')


class SimulateNetwork:
    def __init__(self, calforce=None, mobility=None, timeint=None, 
                 collision=None, topology=None, remesh=None, vis=None,
                 dt0: float=1.0e-8,
                 max_step: int=10,
                 print_freq: int=None,
                 plot_freq: int=None,
                 plot_pause_seconds: float=None,
                 **kwargs) -> None:
        self.calforce = calforce
        self.mobility = mobility
        self.timeint = timeint
        self.collision = collision
        self.topology = topology
        self.remesh = remesh
        self.vis = vis
        self.dt0 = dt0
        self.max_step = max_step
        self.print_freq = print_freq
        self.plot_freq = plot_freq
        self.plot_pause_seconds = plot_pause_seconds
        
    def step(self, N: DisNetManager):
        """step: take a time step of DD simulation on DisNet G
        """
        nodeforce_dict, segforce_dict = self.calforce.NodeForce(N)

        vel_dict = self.mobility.Mobility(N, nodeforce_dict)

        # using a constant time step (for now)
        self.timeint.dt = self.dt0
        self.timeint.Update(N, vel_dict)

        if self.collision is not None:
            self.collision.HandleCol(N)
            
        if self.topology is not None:
            self.topology.Handle(N)

        if self.remesh is not None:
            self.remesh.Remesh(N)
        
    def run(self, N: DisNetManager):
        
        import time
        t0 = time.perf_counter()
        
        if self.vis != None and self.plot_freq != None:
            try: 
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            G = N.get_disnet(DisNet)
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)

        for tstep in range(self.max_step):
            self.step(N)

            if self.print_freq != None:
                if tstep % self.print_freq == 0:
                    dt = self.timeint.dt if self.timeint else 0.0
                    print("step = %d dt = %e"%(tstep, dt))

            if self.vis != None and self.plot_freq != None:
                if tstep % self.plot_freq == 0:
                    G = N.get_disnet(DisNet)
                    self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False, pause_seconds=self.plot_pause_seconds)

        # plot final configuration
        if self.vis != None:
            G = N.get_disnet(DisNet)
            self.vis.plot_disnet(G, fig=fig, ax=ax, trim=True, block=False)
            
        t1 = time.perf_counter()
        print('RUN TIME: %f sec' % (t1-t0))


def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    G = DisNet(cell=cell)
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

    G.add_nodes_links_from_list(rn, links)
    return G
    

def main():
    
    pyexadis.initialize()
    
    G = init_frank_read_src_loop(pbc=False)
    N = DisNetManager({type(G): G})

    bounds = np.array([-0.5*np.diag(G.cell.h), 0.5*np.diag(G.cell.h)])
    vis = VisualizeNetwork(bounds=bounds)
    
    calforce = CalForce1(mu=160e9, nu=0.31, a=0.01, Ec=1.0e6,
                        applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e6, 0.0]),
                        force_mode='LineTension')
    mobility  = MobilityLaw1(mobility_law='SimpleGlide')
    timeint   = TimeIntegration1(integrator='EulerForward')
    collision = Collision1(collision_mode='Proximity')
    topology  = None
    remesh    = Remesh1(remesh_rule='LengthBased', Lmin=0.1, Lmax=0.3)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, 
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          dt0 = 1.0e-8, max_step=200,
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001)
    sim.run(N)
    
    pyexadis.finalize()


if __name__ == "__main__":
    main()
