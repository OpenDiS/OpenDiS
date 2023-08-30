from abc import ABC, abstractmethod
import networkx as nx
import numpy as np

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

from compute_stress_force_analytic_paradis import compute_segseg_force_vec
from compute_stress_force_analytic_python  import python_segseg_force_vec
from compute_stress_analytic_paradis       import compute_seg_stress

try:
    opendis_lib = __import__('opendis')
    found_opendis = True
except ImportError:
    found_opendis = False
    print('---------------------------------------')
    print(' cannot find opendis.py                ')
    print(' to create: cd src ; make libopendis.so')
    print('            cd python ; make           ')
    print('---------------------------------------')

# DisNetwork_ABSTRACT is an abstract (meta) class that contains methods
#  specific for dislocation networks, built upon general network methods
class DisNetwork_ABSTRACT(ABC):
    def __init__(self, data=None, **attr):
         super(DisNetwork_ABSTRACT, self).__init__(data, **attr)

    @abstractmethod
    def nodes(self):
        pass

    @abstractmethod
    def edges(self):
        pass

    @abstractmethod
    def has_node(self, tag):
        pass

    @abstractmethod
    def add_node(self, tag, **attr):
        pass

    @abstractmethod
    def remove_node(self, tag):
        pass

    @abstractmethod
    def has_edge(self, node1, node2):
        pass

    @abstractmethod
    def add_edge(self, node1, node2, **attr):
        pass

    @abstractmethod
    def remove_edge(self, node1, node2):
        pass

    def pos_array(self):
        return np.array([self.nodes[node]['R'] for node in self.nodes])
    
    def seg_list(self):
        # construct segment list (each link appear once: node1 < node2)
        segments = []
        for edge in self.edges():
            node1 = edge[0]
            node2 = edge[1]
            if node1 < node2:
                segments.append({"edge":edge, "burg_vec":self.edges[edge]["burg_vec"], "R1":self.nodes[node1]["R"], "R2":self.nodes[node2]["R"]})
        return segments

    def insert_node(self, node1, node2, tag, pos=None):
        # insert a new node on the link connecting node1 and node2
        if not self.has_edge(node1, node2) or not self.has_edge(node2, node1):
            raise ValueError('insert_node: no link exists between '+str(node1)+' and '+str(node2)) 

        self.add_node(tag, R = pos)
        link12_attr = self.edges[(node1, node2)]
        self.add_edge(node1, tag, **link12_attr)
        self.add_edge(tag, node2, **link12_attr)
        link21_attr = self.edges[(node2, node1)]
        self.add_edge(node2, tag, **link21_attr)
        self.add_edge(tag, node1, **link21_attr)
        self.remove_edge(node1, node2)
        self.remove_edge(node2, node1)

    def remove_two_arm_node(self, node):
        # remove two-arm node
        if not self.has_node(node):
            raise ValueError('remove_two_arm_node: node '+str(node)+' does not exist')
        if self.out_degree(node) != 2:
            raise ValueError('remove_two_arm_node: node '+str(node)+' does not have 2 arms')
        
        node1 = list(self.neighbors(node))[0]
        node2 = list(self.neighbors(node))[1]
        link12_attr = self.edges[(node, node2)]
        self.add_edge(node1, node2, **link12_attr)
        link21_attr = self.edges[(node, node1)]
        self.add_edge(node2, node1, **link21_attr)
        self.remove_edge(node , node1)
        self.remove_edge(node1, node )
        self.remove_edge(node , node2)
        self.remove_edge(node2, node )
        self.remove_node(node)

    def merge_node(self, node1, node2):
        # merge two nodes, node1 and node2, into node1
        # redirect all links to/from node2 into node1

        # remove any links between node1 and node2

        # if node1 has double-links to any neighbor, combine them into one (or zero) link
        pass

    def split_node(self, node, partition):
        # split node into two nodes, with neighbors split according to partition
        pass



class DDParam():
    def __init__(self, mu=1.0, nu=0.3, a=None, Ec=None, bounds=None, boundary_cond=None, force_mode=None,
                       remesh_rule=None, collision_mode=None, mobility_law=None, time_integrator=None):
        self.mu = mu
        self.nu = nu
        self.a  = a
        self.Ec = Ec
        self.bounds = bounds
        self.boundary_cond = boundary_cond
        self.force_mode = force_mode
        self.remesh_rule = remesh_rule
        self.collision_mode = collision_mode
        self.mobility_law = mobility_law
        self.time_integrator = time_integrator

        self.load_type       = None
        self.applied_stress  = np.zeros(6)
        self.dt0 = None
        self.max_step = None
        self.print_freq = None
        self.plot_freq = None
        self.plot_pause_seconds = 0.01

def voigt_vector_to_tensor(voigt_vector):
    return np.array([[voigt_vector[0], voigt_vector[5], voigt_vector[4]],
                     [voigt_vector[5], voigt_vector[1], voigt_vector[3]],
                     [voigt_vector[4], voigt_vector[3], voigt_vector[2]]])

def pkforcevec(sigext, segments):
    nseg = len(segments)
    fpk = np.zeros((nseg, 3))
    for idx, segment in enumerate(segments):
        sigb = sigext @ segment["burg_vec"]
        dR = segment["R2"] - segment["R1"]
        fpk[idx] = np.cross(sigb, dR)
    return fpk

def selfforcevec_LineTension(MU, NU, Ec, segments):
    # to do: vectorize the calculations
    nseg = len(segments)
    fs0 = np.zeros((nseg, 3))
    fs1 = np.zeros((nseg, 3))
    omninv = 1.0/(1.0-NU)
    for idx, segment in enumerate(segments):
        dR = segment["R2"] - segment["R1"]
        L = np.linalg.norm(dR)
        t = dR / L
        bs = np.dot(segment["burg_vec"], t)
        bs2 = bs*bs
        bev = segment["burg_vec"] - bs*t
        be2 = np.sum(bev*bev)
        Score = 2.0*NU*omninv*Ec*bs
        LTcore = (bs2+be2*omninv)*Ec
        fs1[idx] = Score*bev - LTcore*t
    fs0 = -fs1
    return fs0, fs1

class DDSim():
    # DisNetwork is an extension of DiGraph for which each node has attributes such as R, F, V
    #  and each edge has attributes such as burg_vec and plane_normal
    class DisNetwork (nx.DiGraph, DisNetwork_ABSTRACT):
        pass

    def __init__(self, param=None):
        self.disnet = self.DisNetwork()
        self.param = DDParam()
        self.NodeForce_Functions = {'LineTension': self.NodeForce_LineTension}
        self.Remesh_Functions = {'LengthBased': self.Remesh_LengthBased}
        self.Collision_Functions = {'Proximity': self.Collision_Proximity}
        self.MobilityLaw_Functions = {'FCC0': self.MobilityLaw_FCC0}
        self.TimeIntegration_Functions = {'EulerForward': self.TimeIntegration_EulerForward}

    def pbc_position_L(self, r1, r2, L):
        ds = (r2-r1)/L
        ds = ds - ds.round()
        return r1 + ds*L

    def NodeForce(self):
        self.NodeForce_Functions[self.param.force_mode]()

    def Remesh(self):
        self.Remesh_Functions[self.param.remesh_rule]()

    def Collision(self):
        self.Collision_Functions[self.param.collision_mode]()

    def MobilityLaw(self):
        self.MobilityLaw_Functions[self.param.mobility_law]()

    def TimeIntegration(self):
        self.TimeIntegration_Functions[self.param.time_integrator]()

    def Step(self):
        self.NodeForce()
        self.MobilityLaw()
        self.TimeIntegration()
        self.Collision()

    def Run(self):
        if self.param.plot_freq != None:
            try: 
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return
            # plot initial configuration
            self.plot_disnet(fig=fig, ax=ax, trim=True, block=False)

        for tstep in range(self.param.max_step):
            self.Step()

            if self.param.print_freq != None:
                if tstep % self.param.print_freq == 0:
                    print("step = %d dt = %e"%(tstep, self.param.dt))

            if self.param.plot_freq != None:
                if tstep % self.param.plot_freq == 0:
                    self.plot_disnet(fig=fig, ax=ax, trim=True, block=False, pause_seconds=self.param.plot_pause_seconds)

        # plot final configuration
        self.plot_disnet(fig=fig, ax=ax, trim=True, block=False)


    def NodeForce_LineTension(self):
        # pk force from external stress and line tension forces only (from DDLab/src/segforcevec.m)
        # Note: PBC not handled
        #print('NodeForce_LineTension')
        self.disnet.segments = self.disnet.seg_list()
        sigext = voigt_vector_to_tensor(self.param.applied_stress)
        fpk = pkforcevec(sigext, self.disnet.segments)
        fs0, fs1 = selfforcevec_LineTension(self.param.mu, self.param.nu, self.param.Ec, self.disnet.segments)
        self.disnet.fseg = np.hstack((fpk + fs0, fpk + fs1))

        for node in self.disnet.nodes:
            self.disnet.nodes[node]["F"] = np.array([0.0,0.0,0.0])
        for idx, segment in enumerate(self.disnet.segments):
            node1 = segment["edge"][0]
            node2 = segment["edge"][1]
            self.disnet.nodes[node1]["F"] += self.disnet.fseg[idx, 0:3]
            self.disnet.nodes[node2]["F"] += self.disnet.fseg[idx, 3:6]

    def Remesh_LengthBased(self):
        # remesh based on segment length
        #print('Remesh_LengthBased')
        # to be implemented ...
        pass

    def Collision_Proximity(self):
        # use current node position to detect collision
        #print('Collision_Proximity')
        # to be implemented ...
        pass

    def MobilityLaw_FCC0(self):
        # compute velocity on all nodes
        #print('MobilityLaw_FCC0')
        # to be properly implemented (a place holder for now) - steepest descent
        for node in self.disnet.nodes:
            self.disnet.nodes[node]["V"] = self.disnet.nodes[node]["F"]

    def TimeIntegration_EulerForward(self):
        # advance nodes by one time step using Eurler-Forward method
        #print('TimeIntegration_EulerForward')
        self.param.dt = self.param.dt0
        for node in self.disnet.nodes:
            self.disnet.nodes[node]["R"] += self.disnet.nodes[node]["V"] * self.param.dt

    def LinkForce(self):
        # compute forces on every link by calling
        # opendis_lib.SegSegForce(...)
        # shall we move it to an extension ?
        pass

    def plot_disnet(self, plot_links=True, trim=False, fig=None, ax=None, block=False, pause_seconds=0.01):
        if fig==None:
            try: fig = plt.figure(figsize=(8,8))
            except NameError: print('plt not defined'); return
        if ax==None:
            try: ax = plt.axes(projection='3d')
            except NameError: print('plt not defined'); return

        # to do: extend to non-cubic box
        L = self.param.bounds[1][0] - self.param.bounds[0][0]

        rn = self.disnet.pos_array()
        p_link = np.empty((0,6))

        plt.cla()
        if plot_links:
            for my_tag in list(self.disnet.nodes()):
                my_coords = self.disnet.nodes[my_tag]['R']
                for arm in self.disnet.out_edges(my_tag):
                    nbr_tag = arm[1]
                    if my_tag < nbr_tag:
                        r_link = np.zeros((2,3))
                        nbr_coords = self.disnet.nodes[nbr_tag]['R']
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
        ax.set_xlim(self.param.bounds[0][0], self.param.bounds[1][0])
        ax.set_ylim(self.param.bounds[0][1], self.param.bounds[1][1])
        ax.set_zlim(self.param.bounds[0][2], self.param.bounds[1][2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_box_aspect([1,1,1])

        plt.draw()
        plt.show(block=block)
        plt.pause(pause_seconds)

        return fig, ax
