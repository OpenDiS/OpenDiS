import numpy as np
import sys, os, itertools, urllib.request

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork
from pydis.util.network_util import insert_infinite_line

def init_infinite_line(Lbox, burg, plane, theta, seed, maxseg):
    '''Generate a straight, infinite dislocation line through the periodic box, at theta
    degrees from its Burgers vector, through a random point
    '''
    cell = Cell(h=Lbox*np.eye(3), origin=-0.5*Lbox*np.ones(3), is_periodic=[1,1,1])
    origin = Lbox*(np.random.default_rng(seed).random(3) - 0.5)
    # segments just shorter than maxseg, so that remeshing does not split them
    nodes, segs = insert_infinite_line(cell, [], [], burg, plane, origin, theta=theta,
                                       maxseg=0.95*maxseg)
    return DisNetManager(DisNet(cell=cell, rn=np.array(nodes), links=np.array(segs)))

# The internal stress of a dislocation network, frozen, on a 128^3 grid of its periodic box
stress_url = 'https://www.dropbox.com/scl/fi/q8025nk18zrovpfupo1vt/rs0084.stress.npz?rlkey=z1ertzr41fn68c6sjzptrysyq&st=5kvnt6qi&dl=1'
stress_file = 'rs0084.stress.npz'

def load_stress_field(path, url):
    '''The rows x y z Sxx Syy Szz Syz Sxz Sxy of a stress field file, and its grid size;
    the file is downloaded first if it is not there
    '''
    if not os.path.exists(path):
        print('Downloading %s' % url)
        part = '%s.%d.part' % (path, os.getpid())
        urllib.request.urlretrieve(url, part)
        os.replace(part, path)
    with np.load(path) as f:
        return f['data'], f['ngrid']

class StressField:
    '''A stress field given at the centres of a grid of the periodic cell, interpolated
    linearly between them
    '''
    def __init__(self, data, ngrid, cell):
        self.n = np.array(ngrid)
        self.S = data[:,3:9].reshape(*ngrid, 6) # Sxx Syy Szz Syz Sxz Sxy, x slowest, z fastest
        self.L = np.diag(np.array(cell.h))
        self.h = self.L/self.n
        self.origin = np.array(cell.origin)

    def __call__(self, R):
        '''The stress (Pa, xx yy zz yz xz xy) at positions R (b), in any periodic image'''
        u = (np.atleast_2d(R) - self.origin)/self.h - 0.5
        i = np.floor(u).astype(int)
        f = u - i
        S = np.zeros((u.shape[0], 6))
        for corner in itertools.product((0, 1), repeat=3):
            w = np.prod(np.where(corner, f, 1.0-f), axis=1)
            j = (i + corner) % self.n
            S += w[:,None]*self.S[j[:,0], j[:,1], j[:,2]]
        return S

class CalForceUserStress:
    '''The forces of a base force module, plus the Peach-Koehler force of a user stress
    field, as in user_stress.py; the stress is integrated along each segment
    '''
    def __init__(self, calforce_base, user_stress, n_gauss=3):
        self.calforce_base = calforce_base
        self.user_stress = user_stress
        x, w = np.polynomial.legendre.leggauss(n_gauss)
        self.s, self.w = 0.5*(x+1.0), 0.5*w # Gauss points and weights along a segment

    def PreCompute(self, N, state):
        self.calforce_base.PreCompute(N, state)
        return state

    def NodeForce(self, N, state, pre_compute=True):
        self.calforce_base.NodeForce(N, state, pre_compute)
        data = N.export_data()
        rn, tags = data["nodes"]["positions"], data["nodes"]["tags"]
        nid, burg = data["segs"]["nodeids"], data["segs"]["burgers"]
        L = self.user_stress.L
        r1 = rn[nid[:,0]]
        rij = rn[nid[:,1]] - r1
        rij -= L*np.round(rij/L) # across the periodic boundaries
        f = np.zeros_like(rn)
        for s, w in zip(self.s, self.w):
            S = self.user_stress(r1 + s*rij)
            sigb = np.einsum('kij,kj->ki', S[:,[0,5,4,5,1,3,4,3,2]].reshape(-1,3,3), burg)
            fpk = w*np.cross(sigb, rij)
            np.add.at(f, nid[:,0], (1.0-s)*fpk) # shared by the two nodes linearly
            np.add.at(f, nid[:,1], s*fpk)
        rows = {tuple(tag): i for i, tag in enumerate(state["nodeforcetags"])}
        nodeforces = np.array(state["nodeforces"], dtype=float)
        nodeforces[[rows[tuple(tag)] for tag in tags]] += f
        state["nodeforces"] = nodeforces
        return state

    def OneNodeForce(self, N, state, tag, update_state=True):
        # the user stress at the node, added to the applied stress for this node only
        data = N.export_data()
        tags = data["nodes"]["tags"]
        i = np.where((tags[:,0] == tag[0]) & (tags[:,1] == tag[1]))[0]
        applied_stress = state["applied_stress"]
        state["applied_stress"] = applied_stress + self.user_stress(data["nodes"]["positions"][i])[0]
        try:
            return self.calforce_base.OneNodeForce(N, state, tag, update_state)
        finally:
            state["applied_stress"] = applied_stress

def line_length(N):
    '''The total length of the network's segments (b), across the periodic boundaries'''
    data = N.export_data()
    rn, nid = data["nodes"]["positions"], data["segs"]["nodeids"]
    L = np.diag(np.array(data["cell"]["h"]))
    rij = rn[nid[:,1]] - rn[nid[:,0]]
    return np.linalg.norm(rij - L*np.round(rij/L), axis=1).sum()

class SimulateNetworkDisp(SimulateNetwork):
    '''SimulateNetwork that also writes the line's mean displacement against time into
    displacement.dat, with stress_strain_dens.dat: the area that the line has swept on
    its slip plane, from the plastic strain, over its initial length
    '''
    def __init__(self, burg, plane, **kwargs):
        super().__init__(**kwargs)
        self.burg, self.plane = burg, plane

    def run(self, N, state):
        self.length0, self.volume = line_length(N), abs(np.linalg.det(np.array(N.cell.h)))
        self.Ep, self.records = np.zeros(6), [[0, 0.0, 0.0]]
        return super().run(N, state)

    def step_update_response(self, N, state):
        super().step_update_response(N, state)
        self.Ep += state["dEp"]
        gamma = 2.0*self.burg @ self.Ep[[0,5,4,5,1,3,4,3,2]].reshape(3,3) @ self.plane
        area = gamma*self.volume # swept area (b^2)
        self.records.append([len(self.records), state["time"], area/self.length0*state["burgmag"]])
        return state

    def write_results(self):
        super().write_results()
        np.savetxt(os.path.join(self.write_dir, 'displacement.dat'), self.records,
                   fmt='%d %e %e', header='Step Time(s) Displacement(m)')

def main():
    global net, sim, state

    Lbox = 58824.0 # the stress field's periodic box (b)
    burg = 1.0/np.sqrt(2.0)*np.array([-1.0, 0.0, 1.0]) # the slip system B4
    plane = 1.0/np.sqrt(3.0)*np.array([1.0, 1.0, 1.0])
    theta = 0.0 # the line's character angle (deg): 0 screw, 90 edge
    seed = 0 # where the line starts
    net = init_infinite_line(Lbox, burg, plane, theta, seed, maxseg=0.017*Lbox)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    state = {
        "crystal": 'fcc',
        "burgmag": 2.55e-10,
        "mu": 54.6e9,
        "nu": 0.324,
        "a": 6.0,
        "maxseg": 0.017*Lbox,
        "minseg": 0.005*Lbox,
        "rann": 10.0,
        "rtol": 10.0,
        "nextdt": 1e-11,
        "maxdt": 1e-9,
        "use_glide_planes": True,
        # ExaDiS' defaults, so that the two codes run the same simulation
        "splitMultiNodeAlpha": 1e-3,
        "coarsen_mode": 1,
    }

    data, ngrid = load_stress_field(stress_file, stress_url)
    calforce_base = CalForce(force_mode='LineTension', state=state)
    calforce  = CalForceUserStress(calforce_base, StressField(data, ngrid, net.cell))
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint   = TimeIntegration(integrator='Trapezoid', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state, nbrlist=nbrlist)
    topology  = Topology(split_mode='Serial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    tau = 3.0e6 # the resolved shear stress on the slip system (Pa)
    S = np.outer(burg, plane) + np.outer(plane, burg)
    applied_stress = tau*S[[0,1,2,1,0,0],[0,1,2,2,2,1]] # xx yy zz yz xz xy

    sim = SimulateNetworkDisp(burg, plane, calforce=calforce, mobility=mobility,
                              timeint=timeint, topology=topology,
                              collision=collision, remesh=remesh, vis=vis,
                              state=state, max_step=300, loading_mode='stress',
                              applied_stress=applied_stress,
                              print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                              write_freq=10, write_dir='output')
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/user_stress_pydis_final.json')
