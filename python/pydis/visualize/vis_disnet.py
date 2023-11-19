"""@package docstring
Vis_DisNet: class for visualizing dislocation network

Provide plotting functions given a DisNet object
"""

import numpy as np
from ..disnet import DisNet

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

class VisualizeNetwork:
    """VisualizeNetwork: class for plotting dislocation network

    To do: clean up code, improve ways to iterate over nodes and edges 
    """
    def __init__(self, bounds: np.ndarray, **kwargs) -> None:
        self.bounds = bounds
        pass

    def plot_disnet(self, G: DisNet,
                    plot_links=True, trim=False,
                    fig=None, ax=None, block=False, pause_seconds=0.01):
        if fig==None:
            try: fig = plt.figure(figsize=(8,8))
            except NameError: print('plt not defined'); return
        if ax==None:
            try: ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return

        # to do: extend to non-cubic box
        L = self.bounds[1][0] - self.bounds[0][0]

        rn = G.pos_array()
        p_link = np.empty((0,6))

        plt.cla()
        if plot_links:
            for my_tag, node_attr in G.nodes().items():
                my_coords = node_attr['R']
                for arm in G.edges()(my_tag):
                    nbr_tag = arm[1]
                    if my_tag < nbr_tag:
                        r_link = np.zeros((2,3))
                        nbr_coords = G.nodes()[nbr_tag]['R']
                        r_link[0,:] = my_coords
                        # to do: extend to non-cubic box
                        r_link[1,:] = nbr_coords
                        #r_link[1,:] = self.pbc_position_L(my_coords, nbr_coords, L)
                        if (not trim) or np.max(np.absolute(r_link)) <= L/2:
                            p_link = np.append(p_link, [r_link[0,:], r_link[1,:]])

        ls = p_link.reshape((-1,2,3))
        lc = Line3DCollection(ls, linewidths=0.5, colors='b')
        ax.add_collection(lc)

        ax.scatter(rn[:,0], rn[:,1], rn[:,2], c='r', s=4)
        ax.set_xlim(self.bounds[0][0], self.bounds[1][0])
        ax.set_ylim(self.bounds[0][1], self.bounds[1][1])
        ax.set_zlim(self.bounds[0][2], self.bounds[1][2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        try: ax.set_box_aspect([1,1,1])
        except AttributeError: print('ax.set_box_aspect does not work')

        plt.draw()
        plt.show(block=block)
        plt.pause(pause_seconds)

        return fig, ax
