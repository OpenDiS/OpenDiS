import os, sys, itertools, urllib.request
import numpy as np

# Import pyexadis
pyexadis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../core/exadis/python/']
[sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

try:
    import pyexadis
    from framework.disnet_manager import DisNetManager
    from pyexadis_base import ExaDisNet, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
    from pyexadis_utils import insert_infinite_line
except ImportError:
    raise ImportError('Cannot import pyexadis')


# The internal stress of a dislocation network (ParaDiS restart rs0084), frozen,
# on a 128^3 grid of its periodic box. The file is downloaded once.
stress_url = 'https://www.dropbox.com/scl/fi/q8025nk18zrovpfupo1vt/rs0084.stress.npz?rlkey=z1ertzr41fn68c6sjzptrysyq&st=5kvnt6qi&dl=1'
stress_file = 'rs0084.stress.npz'
if not os.path.exists(stress_file):
    urllib.request.urlretrieve(stress_url, stress_file)
with np.load(stress_file) as npz:
    ngrid = npz['ngrid']
    stress_grid = npz['data'][:,3:9].reshape(*ngrid, 6) # xx,yy,zz,yz,xz,xy in Pa


"""
User-defined additional stress field for the simulation.
Here it is the internal stress above, interpolated linearly
between the grid points, at the centres of the grid's cells.
The function must return an array of size (size(x),6), where R is
the array of query positions, and the 6 stress components are in
Voigt order (xx,yy,zz,yz,xz,xy).
"""
def UserStress(R, cell, state):
    h = np.diag(np.array(cell.h))/ngrid # grid spacing
    u = (np.atleast_2d(R) - np.array(cell.origin))/h - 0.5
    i = np.floor(u).astype(int)
    f = u - i

    # weighted sum over the 8 grid points around each position
    user_stress = np.zeros((u.shape[0],6)) # xx,yy,zz,yz,xz,xy in Pa
    for corner in itertools.product((0, 1), repeat=3):
        w = np.prod(np.where(corner, f, 1.0-f), axis=1)
        j = (i + corner) % ngrid # periodic grid
        user_stress += w[:,None]*stress_grid[j[:,0], j[:,1], j[:,2]]
    return user_stress


"""
User-defined CalForce function to add a spatially varying
user-defined stress contribution to a base CalForce model,
as in user_stress.py.
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
        f = state["nodeforces"] # current nodal forces, in the order of the nodes
        np.add.at(f, segsnid[:,0], fpkuser)
        np.add.at(f, segsnid[:,1], fpkuser)

        # store new forces in state dictionary
        state["nodeforces"] = f

        return state

    def OneNodeForce(self, N: DisNetManager, state: dict, tag, update_state=True) -> np.array:

        if self.add_user_stress:
            # here let's just hack the code by temporarily incrementing the global
            # applied_stress by the user_stress value at the node position.
            # this should be fine as we are only computing the force on that node

            applied_stress0 = state["applied_stress"] # save applied_stress

            # find node position from node tag
            G = N.get_disnet(ExaDisNet)
            tags = G.get_tags()
            ind = np.where((tags[:,0]==tag[0])&(tags[:,1]==tag[1]))[0]
            if ind.size != 1:
                raise ValueError("Cannot find node tag (%d,%d) in OneNodeForce" % tuple(tag))
            pos = G.get_positions()[ind]

            user_stress = UserStress(pos, G.cell, state)
            state["applied_stress"] = applied_stress0 + user_stress.squeeze() # increment applied_stress by user_stress
            f = self.calforce_base.OneNodeForce(N, state, tag, update_state) # compute node force

            state["applied_stress"] = applied_stress0 # reset applied_stress

        else:
            f = self.calforce_base.OneNodeForce(N, state, tag, update_state) # compute node force

        return f


def main():
    global net, sim, state

    # A straight screw dislocation on the slip system B4, through a random point of the box
    Lbox = 58824.0 # the stress field's periodic box (b)
    cell = pyexadis.Cell(h=Lbox*np.eye(3), origin=-0.5*Lbox*np.ones(3), is_periodic=[1,1,1])
    burg = 1.0/np.sqrt(2.0)*np.array([-1.0, 0.0, 1.0])
    plane = 1.0/np.sqrt(3.0)*np.array([1.0, 1.0, 1.0])
    origin = Lbox*(np.random.default_rng(0).random(3) - 0.5)
    maxseg = 0.017*Lbox
    # segments just shorter than maxseg, so that remeshing does not split them
    nodes, segs = insert_infinite_line(cell, [], [], burg, plane, origin, theta=0.0,
                                       maxseg=0.95*maxseg)
    net = DisNetManager(ExaDisNet(cell, nodes, segs))

    vis = VisualizeNetwork()

    state = {
        "crystal": 'fcc',
        "burgmag": 2.55e-10,
        "mu": 54.6e9,
        "nu": 0.324,
        "a": 6.0,
        "maxseg": maxseg,
        "minseg": 0.005*Lbox,
        "rann": 10.0,
        "rtol": 10.0,
        "nextdt": 1e-11,
        "maxdt": 1e-9,
        "use_glide_planes": True,
    }

    # Define user-defined CalForce that adds a user-defined stress to a base force model
    calforce_base = CalForce(force_mode='CUTOFF_MODEL', state=state, cutoff=5000.0)
    calforce  = CalForceUserStress(calforce_base=calforce_base, state=state)

    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint   = TimeIntegration(integrator='Trapezoid', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(topology_mode='TopologySerial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    tau = 5.0e6 # the resolved shear stress on the slip system (Pa)
    S = np.outer(burg, plane) + np.outer(plane, burg)
    applied_stress = tau*S[[0,1,2,1,0,0],[0,1,2,2,2,1]] # xx yy zz yz xz xy

    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          state=state, max_step=3000, loading_mode='stress',
                          applied_stress=applied_stress,
                          print_freq=30, plot_freq=30, plot_pause_seconds=0.0001,
                          write_freq=30, write_dir='output')
    sim.run(net, state)


if __name__ == "__main__":
    pyexadis.initialize()

    main()

    # explore the network after simulation
    G  = net.get_disnet(ExaDisNet)

    os.makedirs('output', exist_ok=True)
    net.write_json('output/user_stress_exadis_final.json')

    if not sys.flags.interactive:
        pyexadis.finalize()
