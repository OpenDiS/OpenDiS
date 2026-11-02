import os, sys
import numpy as np

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork
from pydis.util.network_util import insert_infinite_line


def init_line_config(Lbox, burg, plane, theta, maxseg):
    """An infinite line along y on a glide plane normal to z, and the crystal's orientation"""
    plane = plane / np.linalg.norm(plane)
    b = burg / np.linalg.norm(burg)
    y = np.cross(plane, b)
    y = y / np.linalg.norm(y)
    linedir = np.cos(theta*np.pi/180.0)*b + np.sin(theta*np.pi/180.0)*y
    x = np.cross(linedir, plane)
    Rorient = np.array([x / np.linalg.norm(x), linedir, plane])
    cell = Cell(h=Lbox*np.eye(3), is_periodic=[1,1,1])
    origin = np.array(cell.origin) + 0.5*np.sum(np.array(cell.h), axis=0)
    nodes, segs = insert_infinite_line(cell, [], [], Rorient @ burg, Rorient @ plane, origin,
                                       linedir=-(Rorient @ linedir), maxseg=maxseg)
    return DisNetManager(DisNet(cell=cell, rn=np.array(nodes), links=np.array(segs))), Rorient


class Precipitate:
    """A spherical precipitate: its centre and radius (b)"""
    def __init__(self, pos, radius):
        self.pos, self.radius = np.array(pos), radius


class My_SimulateNetwork(SimulateNetwork):
    """SimulateNetwork that moves nodes entering a precipitate back to its surface"""
    def __init__(self, *args, precips=(), **kwargs):
        super().__init__(*args, **kwargs)
        self.precips = precips

    def step_begin(self, N, state):
        self.xold = N.export_data()["nodes"]["positions"].copy()

    def step_post_integrate(self, N, state):
        data = N.export_data()
        pos = data["nodes"]["positions"]
        cell = Cell(**data["cell"])
        back = np.array(cell.closest_image(Rref=pos, R=self.xold)) - pos
        moved = np.linalg.norm(back, axis=1) > 1e-10
        for p in self.precips:
            r = np.array(cell.closest_image(Rref=p.pos, R=pos)) - p.pos
            inside = np.flatnonzero(moved & (np.einsum('ij,ij->i', r, r) < p.radius**2))
            if inside.size == 0:
                continue
            r, d = r[inside], back[inside] / np.linalg.norm(back[inside], axis=1)[:, None]
            old = np.array(cell.closest_image(Rref=p.pos, R=self.xold[inside])) - p.pos
            born = np.einsum('ij,ij->i', old, old) < p.radius**2  # already inside: go out radially
            d[born] = r[born] / np.linalg.norm(r[born], axis=1)[:, None]
            rd = np.einsum('ij,ij->i', r, d)
            q = r - rd[:, None]*d
            pos[inside] += (np.sqrt(p.radius**2 - np.einsum('ij,ij->i', q, q)) - rd)[:, None]*d
        N.get_disnet(DisNet).import_data(data)


def main():
    global net, sim, state

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
        "use_glide_planes": True,
        # ExaDiS' defaults, so that the two codes run the same simulation
        "splitMultiNodeAlpha": 1e-3,
        "coarsen_mode": 1,
    }

    L0 = 50000.0
    burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.])
    plane = np.array([1.,-1.,1.])
    net, state["Rorient"] = init_line_config(L0, burg, plane, 90.0, state["maxseg"])
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    # precipitates on the glide plane
    Nprecip = 10
    Rprecip = [500., 2000.] # radius range
    np.random.seed(1234)
    ppos = L0*np.random.rand(Nprecip, 3)
    ppos[:,2] = 0.5*L0
    rp = Rprecip[0]+(Rprecip[1]-Rprecip[0])*np.random.rand(Nprecip)
    precips = [Precipitate(p, r) for p, r in zip(ppos, rp)]

    vis = VisualizeNetwork()

    calforce  = CalForce(force_mode='Elasticity_SBA', state=state, cutoff=0.1*L0)
    mobility  = MobilityLaw(mobility_law='FCC_0', state=state, Medge=60000.0, Mscrew=60000.0, vmax=3400.0)
    timeint   = TimeIntegration(integrator='Trapezoid', state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode='Retroactive', state=state, nbrlist=nbrlist)
    topology  = Topology(split_mode='Serial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    applied_stress = np.array([0.0, 0.0, 0.0, 0.0, 1e8, 0.0]) # xx,yy,zz,yz,xz,xy in Pa

    sim = My_SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                             topology=topology, collision=collision, remesh=remesh, vis=vis,
                             state=state, max_step=200, loading_mode='stress',
                             applied_stress=applied_stress, burgmag=state["burgmag"],
                             print_freq=1, plot_freq=10, plot_pause_seconds=0.0001,
                             write_freq=10, write_dir='output', precips=precips)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    main()

    # explore the network after simulation
    G = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/precipitates_pydis_final.json')
