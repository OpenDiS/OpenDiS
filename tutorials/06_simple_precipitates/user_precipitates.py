import os, sys
import numpy as np

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, VisualizeNetwork
    from pyexadis_base import SimulateNetwork, NodeConstraints
except ImportError:
    raise ImportError('Cannot import pyexadis')

import matplotlib.pyplot as plt


class Precipitate:
    """Precipitate
    Class to hold precipitate object
    """
    def __init__(self, pos: np.array, radius: float):
        self.pos = np.array(pos)
        self.radius = radius


class My_VisualizeNetwork(VisualizeNetwork):
    """My_VisualizeNetwork
    User-defined visualizer to plot dislocations and precipitates
    """
    def __init__(self, *args, **kwargs) -> None:
        super(My_VisualizeNetwork, self).__init__(*args, **kwargs)
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


class My_SimulateNetwork(SimulateNetwork):
    """My_SimulateNetwork
    User-defined simulation driver that adds the precipitate-dislocation interaction
    """
    def __init__(self, *args, **kwargs) -> None:
        super(My_SimulateNetwork, self).__init__(*args, **kwargs)
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
        
