"""common.codegen: the PyDiS or ExaDiS script of a spec, step by step and in full.

The full script is the tutorial's script, such as
``examples/02_frank_read_src/test_frank_read_src_pydis.py`` or ``..._exadis.py``, with the
spec's values. For an example's own values it runs the same simulation as the tutorial's
script, which ``tests/test_codegen.py`` checks. It differs from that script in four ways: it
finds PyDiS and ExaDiS from ``$OPENDIS_DIR`` (else from the repository the app runs in), so
that it runs from any folder; each step's code is one block, so a few lines come in another
order; an optional value (Ec, mob, a Burgers vector) appears only when it is not the default;
and ExaDiS's script also writes its snapshots as PyDiS's does, ``output/disnet_<step>.json``,
which the app plays.

The user stress's script also defines what its tutorials in ``tutorials/05_user_stress/``
define: the reader of its frozen field and the force module that adds it (step 3), and a
driver that writes the line's displacement (step 7). It downloads the field from the spec's
link once, into that folder, a file for each link.

A script that builds a network and runs no simulation is one option of
``tutorials/01_initial_configurations/initial_configurations.py``, with either code: the
network that it builds is the tutorial's, which ``tests/test_codegen.py`` checks too. It writes
its files into ``output/``, the network as its final one.
"""
import hashlib
import textwrap
from dataclasses import replace

import numpy as np

from common.paths import OPENDIS_DIR, PYDIS_PATHS
from common.spec import (FCCJunction, FrankReadSource, InfiniteLines, LoopSources, ParaDiSFile,
                         SimpleSource, Spec, TwoLines, UserStressLine)

EXADIS_PATHS = PYDIS_PATHS + ("core/exadis/python",)  # as in the ExaDiS tutorials' scripts
# The folder where the user stress's scripts keep a copy of each link's field, from the
# repository's root
STRESS_DIR = "tutorials/05_user_stress"


def num(x: float) -> str:
    """A float as Python code: plain from 0.001 to 1e5, else as d.ddde±N, the shortest exact one."""
    x = float(x)
    if x == 0.0 or 1e-3 <= abs(x) < 1e5:
        return repr(x)
    for digits in range(1, 18):
        text = f"{x:.{digits}e}"
        if float(text) == x:
            break
    mantissa, exponent = text.split("e")
    mantissa = mantissa.rstrip("0")
    return f"{mantissa}0e{int(exponent)}" if mantissa.endswith(".") else f"{mantissa}e{int(exponent)}"


def array(values: tuple[float, ...]) -> str:
    """A NumPy array of floats, as the tutorials write one."""
    return f"np.array([{', '.join(num(v) for v in values)}])"


_PATHS = '''\
import numpy as np
import sys, os

# The root of the OpenDiS repository: $OPENDIS_DIR, else the one that the app ran from
opendis_dir = os.environ.get('OPENDIS_DIR', {opendis_dir!r})
{name} = [os.path.join(opendis_dir, path) for path in {paths!r}]
[sys.path.append(os.path.abspath(path)) for path in {name}
 if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

'''
IMPORTS = {
    "pydis": '''\
from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork
''',
    "exadis": '''\
import pyexadis
from framework.disnet_manager import DisNetManager
from pyexadis_base import ExaDisNet, NodeConstraints
from pyexadis_base import SimulateNetwork, VisualizeNetwork
from pyexadis_base import CalForce, MobilityLaw, TimeIntegration
from pyexadis_base import Collision, Topology, Remesh
''',
}

# The tutorials' functions that make the initial configuration, word for word, by code and kind
INIT_FUNCTIONS = {
    ("pydis", FrankReadSource.kind): '''\
def init_frank_read_src_loop(arm_length=1.0, box_length=8.0,
                             burg_vec=np.array([1.0,0.0,0.0]),
                             pbc=False):
    \'\'\'Generate an initial Frank-Read source configuration
    \'\'\'
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    PINNED = DisNode.Constraints.PINNED_NODE
    FREE   = DisNode.Constraints.UNCONSTRAINED
    rn = np.array([[0.0, -arm_length/2.0, 0.0,         PINNED],
                   [0.0,  0.0,            0.0,         FREE],
                   [0.0,  arm_length/2.0, 0.0,         PINNED],
                   [0.0,  arm_length/2.0, -arm_length, PINNED],
                   [0.0, -arm_length/2.0, -arm_length, PINNED]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))
''',
    ("exadis", FrankReadSource.kind): '''\
def init_frank_read_src_loop(arm_length=1.0, box_length=8.0,
                             burg_vec=np.array([1.0,0.0,0.0]),
                             pbc=False):
    \'\'\'Generate an initial Frank-Read source configuration
    \'\'\'
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    PINNED = NodeConstraints.PINNED_NODE
    FREE   = NodeConstraints.UNCONSTRAINED
    rn = np.array([[0.0, -arm_length/2.0, 0.0,         PINNED],
                   [0.0,  0.0,            0.0,         FREE],
                   [0.0,  arm_length/2.0, 0.0,         PINNED],
                   [0.0,  arm_length/2.0, -arm_length, PINNED],
                   [0.0, -arm_length/2.0, -arm_length, PINNED]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(ExaDisNet(cell, rn, links))
''',
    ("pydis", TwoLines.kind): '''\
def init_two_disl_lines(z0=1.0, box_length=8.0, b1=np.array([-1.0,1.0,1.0]), b2=np.array([1.0,-1.0,1.0]), pbc=False):
    \'\'\'Generate an initial configuration for two dislocation lines.
    \'\'\'
    print("init_two_disl_lines: z0 = %f" % (z0))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    rn    = np.array([[0.0, -z0, -z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  z0,  z0,  DisNode.Constraints.PINNED_NODE],
                      [-z0,  0.0,-z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [ z0,  0.0, z0,  DisNode.Constraints.PINNED_NODE]])
    rn[:,0:3] += cell.center()

    xi1, xi2 = rn[2,:3] - rn[1,:3], rn[5,:3] - rn[4,:3]
    n1, n2 = np.cross(b1, xi1), np.cross(b2, xi2)
    n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0, :] = np.concatenate(([0, 1], b1, n1))
    links[1, :] = np.concatenate(([1, 2], b1, n1))
    links[2, :] = np.concatenate(([3, 4], b2, n2))
    links[3, :] = np.concatenate(([4, 5], b2, n2))

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))
''',
    ("exadis", TwoLines.kind): '''\
def init_two_disl_lines(z0=1.0, box_length=8.0, b1=np.array([-1.0,1.0,1.0]), b2=np.array([1.0,-1.0,1.0]), pbc=False):
    \'\'\'Generate an initial configuration for two dislocation lines.
    \'\'\'
    print("init_two_disl_lines: z0 = %f" % (z0))
    cell = pyexadis.Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
    center = np.array(cell.center())

    rn    = np.array([[0.0, -z0, -z0,  NodeConstraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, NodeConstraints.UNCONSTRAINED],
                      [0.0,  z0,  z0,  NodeConstraints.PINNED_NODE],
                      [-z0,  0.0,-z0,  NodeConstraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, NodeConstraints.UNCONSTRAINED],
                      [ z0,  0.0, z0,  NodeConstraints.PINNED_NODE]])
    rn[:,0:3] += center

    xi1, xi2 = rn[2,:3] - rn[1,:3], rn[5,:3] - rn[4,:3]
    n1, n2 = np.cross(b1, xi1), np.cross(b2, xi2)
    n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    links = np.zeros((4, 8))
    links[0, :] = np.concatenate(([0, 1], b1, n1))
    links[1, :] = np.concatenate(([1, 2], b1, n1))
    links[2, :] = np.concatenate(([3, 4], b2, n2))
    links[3, :] = np.concatenate(([4, 5], b2, n2))

    return DisNetManager(ExaDisNet(cell, rn, links))
''',
}
# The fcc junction's function, the same in both codes but for its cell and its network
_FCC_JUNCTION = '''\
def init_fcc_junction(junction_type, disloc_length, phi1, phi2):
    \'\'\'Generate two Frank-Read sources on two fcc slip systems that cross at the box center
    \'\'\'
    Lbox = 1.0*disloc_length
    cell = {cell}(h=Lbox*np.eye(3), is_periodic=[1,1,1])
    nodes, segs = [], []

    b1 = 1.0/np.sqrt(2.0)*np.array([0.0, 1.0, 1.0])
    p1 = np.array([1.0, 1.0, -1.0])

    if junction_type == 'lomer':
        b2 = 1.0/np.sqrt(2.0)*np.array([1.0, 0.0,-1.0])
        p2 = np.array([1.0, 1.0, 1.0])
    elif junction_type == 'glissile':
        b2 = 1.0/np.sqrt(2.0)*np.array([1.0, 0.0,-1.0])
        p2 = np.array([1.0,-1.0, 1.0])
    elif junction_type == 'hirth':
        b2 = 1.0/np.sqrt(2.0)*np.array([0.0,-1.0, 1.0])
        p2 = np.array([1.0, 1.0, 1.0])
    elif junction_type == 'colinear':
        b2 = 1.0/np.sqrt(2.0)*np.array([0.0,-1.0,-1.0])
        p2 = np.array([1.0,-1.0, 1.0])
    else:
        raise ValueError('Unknown junction type %s'%junction_type)

    linter = np.cross(p1, p2)
    linter = linter / np.linalg.norm(linter)
    y1 = np.cross(linter, p1)
    y1 = y1 / np.linalg.norm(y1)
    ldir1 = np.cos(phi1*np.pi/180.0)*linter+np.sin(phi1*np.pi/180.0)*y1
    y2 = np.cross(linter, p2)
    y2 = y2 / np.linalg.norm(y2)
    ldir2 = np.cos(phi2*np.pi/180.0)*linter+np.sin(phi2*np.pi/180.0)*y2

    center = 0.5*Lbox * np.ones(3)
    delta = 0.01 * disloc_length * np.array([1.0, 1.0, 0.0])

    nodes, segs = insert_frank_read_src(cell, nodes, segs, b1, p1, disloc_length, center+delta, linedir=ldir1)
    nodes, segs = insert_frank_read_src(cell, nodes, segs, b2, p2, disloc_length, center-delta, linedir=ldir2)

    return DisNetManager({network})
'''
INIT_FUNCTIONS.update({
    ("pydis", FCCJunction.kind): _FCC_JUNCTION.format(
        cell="Cell", network="DisNet(cell=cell, rn=np.array(nodes), links=np.array(segs))"),
    ("exadis", FCCJunction.kind): _FCC_JUNCTION.format(
        cell="pyexadis.Cell", network="ExaDisNet(cell, nodes, segs)"),
})
# The user stress's line, the same in both codes but for its cell and its network
_INFINITE_LINE = '''\
def init_infinite_line(Lbox, burg, plane, theta, seed, maxseg):
    \'\'\'Generate a straight, infinite dislocation line through the periodic box, at theta
    degrees from its Burgers vector, through a random point
    \'\'\'
    cell = {cell}(h=Lbox*np.eye(3), origin=-0.5*Lbox*np.ones(3), is_periodic=[1,1,1])
    origin = Lbox*(np.random.default_rng(seed).random(3) - 0.5)
    # segments just shorter than maxseg, so that remeshing does not split them
    nodes, segs = insert_infinite_line(cell, [], [], burg, plane, origin, theta=theta,
                                       maxseg=0.95*maxseg)
    return DisNetManager({network})
'''
INIT_FUNCTIONS.update({
    ("pydis", UserStressLine.kind): _INFINITE_LINE.format(
        cell="Cell", network="DisNet(cell=cell, rn=np.array(nodes), links=np.array(segs))"),
    ("exadis", UserStressLine.kind): _INFINITE_LINE.format(
        cell="pyexadis.Cell", network="ExaDisNet(cell, nodes, segs)"),
})
# What a kind's function imports besides the code's own modules
KIND_IMPORTS = {
    ("pydis", FCCJunction.kind): "from pydis.util.network_util import insert_frank_read_src\n",
    ("exadis", FCCJunction.kind): "from pyexadis_utils import insert_frank_read_src\n",
    ("pydis", UserStressLine.kind): "import hashlib, itertools, urllib.request\n"
                                    "from pydis.util.network_util import insert_infinite_line\n",
    ("exadis", UserStressLine.kind): "import hashlib, itertools, urllib.request\n"
                                     "from pyexadis_utils import insert_infinite_line\n",
}
# The user stress's field, and the force module that adds it to a base module's forces, as in
# tutorials/05_user_stress/user_stress.py; after the lines that give its link and file
USER_STRESS = '''\
def load_stress_field(path, url, max_mb=300):
    \'\'\'The rows x y z Sxx Syy Szz Syz Sxz Sxy of a stress field file, and its grid size;
    the file, of at most max_mb, is downloaded first if it is not there
    \'\'\'
    if not os.path.exists(path):
        print('Downloading %s' % url)
        part = '%s.%d.part' % (path, os.getpid())
        def check_size(blocks, block_size, size):
            if max(blocks*block_size, size) > max_mb*1e6:
                raise ValueError('the file at %s is larger than %d MB' % (url, max_mb))
        try:
            urllib.request.urlretrieve(url, part, reporthook=check_size)
            os.replace(part, path)
        finally:
            if os.path.exists(part):
                os.remove(part)
    with np.load(path) as f:
        if f['columns'].tolist() != ['x', 'y', 'z', 'Sxx', 'Syy', 'Szz', 'Syz', 'Sxz', 'Sxy']:
            raise ValueError('the columns of %s are not x y z Sxx Syy Szz Syz Sxz Sxy' % path)
        return f['data'], f['ngrid']

class StressField:
    \'\'\'A stress field given at the centres of a grid of the periodic cell, interpolated
    linearly between them
    \'\'\'
    def __init__(self, data, ngrid, cell):
        self.n = np.array(ngrid)
        self.S = data[:,3:9].reshape(*ngrid, 6) # Sxx Syy Szz Syz Sxz Sxy, x slowest, z fastest
        self.L = np.diag(np.array(cell.h))
        self.h = self.L/self.n
        self.origin = np.array(cell.origin)
        # the grid's first points along z, y and x, and its last, at the centres of its cells
        ijk = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], self.n - 1])
        rows = ijk @ [self.n[1]*self.n[2], self.n[2], 1]
        if not np.allclose(data[rows,:3], self.origin + (ijk + 0.5)*self.h, rtol=0.0,
                           atol=0.01*self.h.min()):
            raise ValueError('the stress field is not on the grid of the cell, x slowest and '
                             'z fastest: its first point is %s, its last %s'
                             % (data[0,:3], data[-1,:3]))

    def __call__(self, R):
        \'\'\'The stress (Pa, xx yy zz yz xz xy) at positions R (b), in any periodic image\'\'\'
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
    \'\'\'The forces of a base force module, plus the Peach-Koehler force of a user stress
    field, as in user_stress.py; the stress is integrated along each segment
    \'\'\'
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
'''
# The user stress's driver, which also writes the line's displacement against time
DISP_DRIVER = '''\
def line_length(N):
    \'\'\'The total length of the network's segments (b), across the periodic boundaries\'\'\'
    data = N.export_data()
    rn, nid = data["nodes"]["positions"], data["segs"]["nodeids"]
    L = np.diag(np.array(data["cell"]["h"]))
    rij = rn[nid[:,1]] - rn[nid[:,0]]
    return np.linalg.norm(rij - L*np.round(rij/L), axis=1).sum()

class SimulateNetworkDisp({base}):
    \'\'\'{base} that also writes the line's mean displacement against time into
    displacement.dat, with stress_strain_dens.dat: the area that the line has swept on
    its slip plane, from the plastic strain, over its initial length
    \'\'\'
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
'''
# ExaDiS's SimulateNetwork writes its snapshots as config.<step>.data; this one also writes
# them as PyDiS's does, for the app's views
JSON_DRIVER = '''\
class SimulateNetworkJSON(SimulateNetwork):
    \'\'\'ExaDiS's SimulateNetwork, which also writes each snapshot as
    output/disnet_<step>.json, as PyDiS's does
    \'\'\'
    def run(self, N, state):
        if self.write_freq is not None:
            N.write_json(os.path.join(self.write_dir, 'disnet_0.json'))
        return super().run(N, state)

    def step_write_files(self, N, state):
        super().step_write_files(N, state)
        if self.write_freq is not None and state["istep"] % self.write_freq == 0:
            N.write_json(os.path.join(self.write_dir, 'disnet_%d.json' % state["istep"]))
'''
INIT_FUNCTIONS.update({(engine, LoopSources.kind): "" for engine in ("pydis", "exadis")})
DEFAULT_BURGERS = {"burg_vec": (1.0, 0.0, 0.0), "b1": (-1.0, 1.0, 1.0), "b2": (1.0, -1.0, 1.0)}


def _optional_burgers(config: FrankReadSource | TwoLines, names: tuple[str, ...]) -> str:
    """The Burgers vector arguments that differ from the function's defaults."""
    return "".join(f", {name}={array(getattr(config, name))}" for name in names
                   if getattr(config, name) != DEFAULT_BURGERS[name])


def _config(spec: Spec) -> str:
    """Step 1: the box and the initial configuration."""
    c = spec.config
    if isinstance(c, FrankReadSource):
        return (f"Lbox = {num(c.box_length)}\n"
                f"net = init_frank_read_src_loop(box_length=Lbox,\n"
                f"                               arm_length={num(c.arm_ratio)}*Lbox, "
                f"pbc={c.pbc}{_optional_burgers(c, ('burg_vec',))})\n")
    if isinstance(c, LoopSources):  # with step 2's crystal and maxseg
        m = spec.material
        return (f"Lbox = {num(c.box_length)}\n"
                f"G = {'ExaDisNet' if spec.engine == 'exadis' else 'DisNet'}()\n"
                f"G.generate_prismatic_config({m.crystal!r}, Lbox, {c.num_loops}, "
                f"{num(c.radius_ratio)}*Lbox,\n"
                f"                            {num(m.maxseg_ratio)}*Lbox, uniform={c.uniform}, "
                f"seed={c.seed})\n"
                f"net = DisNetManager(G)\n")
    if isinstance(c, FCCJunction):  # the lines are as long as the box
        return (f"Lbox = {num(c.box_length)}\n"
                f"net = init_fcc_junction(junction_type={c.junction_type!r}, disloc_length=Lbox,\n"
                f"                        phi1={num(c.phi1)}, phi2={num(c.phi2)})\n")
    if isinstance(c, UserStressLine):  # with step 2's maxseg
        return (f"Lbox = {num(c.box_length)} # the stress field's periodic box (b)\n"
                f"burg = 1.0/np.sqrt(2.0)*np.array([-1.0, 0.0, 1.0]) # the slip system B4\n"
                f"plane = 1.0/np.sqrt(3.0)*np.array([1.0, 1.0, 1.0])\n"
                f"theta = {num(c.theta)} # the line's character angle (deg): 0 screw, 90 edge\n"
                f"seed = {c.seed} # where the line starts\n"
                f"net = init_infinite_line(Lbox, burg, plane, theta, seed, "
                f"maxseg={num(spec.material.maxseg_ratio)}*Lbox)\n")
    return (f"Lbox = {num(c.box_length)}\n"
            f"z0 = {num(c.z0)}\n"
            f"net = init_two_disl_lines(z0=z0, box_length=Lbox, "
            f"pbc={c.pbc}{_optional_burgers(c, ('b1', 'b2'))})\n")


def _material(spec: Spec) -> str:
    """Step 2: the state dictionary."""
    m = spec.material
    items = [f'"{name}": {value!r}' if isinstance(value, (str, bool)) else f'"{name}": {num(value)}'
             for name, value in (("crystal", m.crystal), ("rtol", m.rtol), ("nextdt", m.nextdt),
                                 ("maxdt", m.maxdt), ("use_glide_planes", m.use_glide_planes))
             if value is not None]
    extras = ",\n         ".join(", ".join(items[i:i + 3]) for i in range(0, len(items), 3))
    return (f'state = {{"burgmag": {num(m.burgmag)}, "mu": {num(m.mu)}, "nu": {num(m.nu)}, '
            f'"a": {num(m.a)},\n'
            f'         "maxseg": {num(m.maxseg_ratio)}*Lbox, "minseg": {num(m.minseg_ratio)}*Lbox, '
            f'"rann": {num(m.rann)}' + (f",\n         {extras}" if extras else "") + "}\n")


def _forces(spec: Spec) -> str:
    """Step 3: the force and the mobility modules."""
    f = spec.forces
    ec = "" if f.Ec is None else f", Ec={num(f.Ec)}"
    cutoff = "" if f.cutoff is None else f", cutoff={num(f.cutoff)}"
    calforce = f"calforce  = CalForce(force_mode={f.force_mode!r}, state=state{ec}{cutoff})\n"
    if f.stress_field:  # the user stress: its field added to a base module's forces
        calforce = ("data, ngrid = load_stress_field(stress_file, stress_url)\n"
                    f"calforce_base = CalForce(force_mode={f.force_mode!r}, state=state{ec}{cutoff})\n"
                    "calforce  = CalForceUserStress(calforce_base, StressField(data, ngrid, net.cell))\n")
    if f.mobility_law in ("BCC_0B", "FCC_0"):
        climb = f", Mclimb={num(f.Mclimb)}" if f.mobility_law == "BCC_0B" else ""
        vmax = "" if f.vmax is None else f", vmax={num(f.vmax)}"
        return (calforce +
                f"mobility  = MobilityLaw(mobility_law={f.mobility_law!r}, state=state, "
                f"Medge={num(f.Medge)},\n"
                f"                        Mscrew={num(f.Mscrew)}{climb}{vmax})\n")
    if spec.engine == "exadis":  # ExaDiS's mobility takes mob as an argument
        mob = "" if f.mob is None else f", mob={num(f.mob)}"
        return (calforce +
                f"mobility  = MobilityLaw(mobility_law={f.mobility_law!r}, state=state{mob})\n")
    mob = "" if f.mob is None else f'state["mob"] = {num(f.mob)}\n'
    return (calforce + mob +
            f"mobility  = MobilityLaw(mobility_law={f.mobility_law!r}, state=state)\n")


def _timeint(spec: Spec) -> str:
    """Step 4: the time integrator."""
    t = spec.timeint
    if t.integrator == "Trapezoid":  # its time step adapts, from state["nextdt"]
        multi = f"multi={t.multi}, " if t.multi else ""
        return (f"timeint   = TimeIntegration(integrator='Trapezoid', {multi}state=state,\n"
                f"                            force=calforce, mobility=mobility)\n")
    return (f"timeint   = TimeIntegration(integrator={t.integrator!r}, dt={num(t.dt)},\n"
            f"                            state=state)\n")


def _topology(spec: Spec) -> str:
    """Step 5: the cell list (PyDiS), and the topology, collision and remesh modules."""
    t = spec.topology
    if spec.engine == "exadis":
        topology = ("None\n" if t.topology_mode == "None" else
                    f"Topology(topology_mode={t.topology_mode!r}, state=state,\n"
                    f"                     force=calforce, mobility=mobility)\n")
        return (f"collision = Collision(collision_mode={t.collision_mode!r}, state=state)\n"
                f"topology  = {topology}"
                f"remesh    = Remesh(remesh_rule={t.remesh_rule!r}, state=state)\n")
    n_div = ",".join(str(n) for n in t.n_div)
    matched = "".join(f'state["{name}"] = {value}\n' for name, value in  # as ExaDiS's defaults
                      (("splitMultiNodeAlpha", None if t.split_alpha is None else num(t.split_alpha)),
                       ("coarsen_mode", t.coarsen_mode)) if value is not None)
    return (matched +
            f"nbrlist   = CellList(cell=net.cell, n_div=[{n_div}])\n"
            f"topology  = Topology(split_mode={t.split_mode!r}, state=state,\n"
            f"                     force=calforce, mobility=mobility)\n"
            f"collision = Collision(collision_mode={t.collision_mode!r}, state=state,\n"
            f"                      nbrlist=nbrlist)\n"
            f"remesh    = Remesh(remesh_rule={t.remesh_rule!r}, state=state)\n")


def _loading(spec: Spec) -> str:
    """Step 6: the applied stress, or the strain rate and its direction."""
    loading = spec.loading
    if loading.loading_mode == "strain_rate":
        return f"erate = {num(loading.erate)}\nedir = {array(loading.edir)}\n"
    if loading.tau is not None:  # the user stress: a shear stress on its line's slip system
        return (f"tau = {num(loading.tau)} # the resolved shear stress on the slip system (Pa)\n"
                "S = np.outer(burg, plane) + np.outer(plane, burg)\n"
                "applied_stress = tau*S[[0,1,2,1,0,0],[0,1,2,2,2,1]] # xx yy zz yz xz xy\n")
    return f"applied_stress = {array(loading.applied_stress)}\n"


def _run_user_stress(spec: Spec) -> str:
    """Step 7 of the user stress: its driver also writes the line's displacement."""
    r, pad = spec.run, " " * len("sim = SimulateNetworkDisp(")
    modules = (f"timeint=timeint, collision=collision,\n{pad}topology=topology, remesh=remesh, "
               f"vis=vis,\n" if spec.engine == "exadis" else
               f"timeint=timeint, topology=topology,\n{pad}collision=collision, remesh=remesh, "
               f"vis=vis,\n")
    return (f"vis = VisualizeNetwork()\n"
            f"sim = SimulateNetworkDisp(burg, plane, calforce=calforce, mobility=mobility,\n"
            f"{pad}{modules}"
            f"{pad}state=state, max_step={r.max_step}, loading_mode='stress',\n"
            f"{pad}applied_stress=applied_stress,\n"
            f"{pad}print_freq={r.print_freq}, plot_freq={r.plot_freq}, "
            f"plot_pause_seconds={num(r.plot_pause_seconds)},\n"
            f"{pad}write_freq={r.write_freq}, write_dir='output'"
            f"{'' if spec.engine == 'exadis' else ', save_state=False'})\n"
            f"sim.run(net, state)\n")


def _run(spec: Spec) -> str:
    """Step 7: the simulation, and its run."""
    if isinstance(spec.config, UserStressLine):
        return _run_user_stress(spec)
    r, loading_mode = spec.run, spec.loading.loading_mode
    load = ('erate=erate, edir=edir, burgmag=state["burgmag"],' if loading_mode == "strain_rate"
            else "applied_stress=applied_stress,")
    if spec.engine == "exadis":
        return (f"vis = VisualizeNetwork()\n"
                f"sim = SimulateNetworkJSON(calforce=calforce, mobility=mobility,\n"
                f"                          timeint=timeint, collision=collision,\n"
                f"                          topology=topology, remesh=remesh, vis=vis,\n"
                f"                          state=state, max_step={r.max_step}, "
                f"loading_mode={loading_mode!r},\n"
                f"                          {load}\n"
                f"                          print_freq={r.print_freq}, plot_freq={r.plot_freq}, "
                f"plot_pause_seconds={num(r.plot_pause_seconds)},\n"
                f"                          write_freq={r.write_freq}, write_dir='output')\n"
                f"sim.run(net, state)\n")
    return (f"vis = VisualizeNetwork()\n"
            f"sim = SimulateNetwork(calforce=calforce, mobility=mobility,\n"
            f"                      timeint=timeint, topology=topology,\n"
            f"                      collision=collision, remesh=remesh, vis=vis,\n"
            f"                      state=state, max_step={r.max_step}, "
            f"loading_mode={loading_mode!r},\n"
            f"                      {load}\n"
            f"                      print_freq={r.print_freq}, plot_freq={r.plot_freq}, "
            f"plot_pause_seconds={num(r.plot_pause_seconds)},\n"
            f"                      write_freq={r.write_freq}, write_dir='output', "
            f"save_state=False)\n"
            f"sim.run(net, state)\n")


STEP_CODE = {"config": _config, "material": _material, "forces": _forces, "timeint": _timeint,
             "topology": _topology, "loading": _loading, "run": _run}


# ── A script that builds a network and runs no simulation ──

BUILD_IMPORTS = {
    "pydis": '''\
from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, VisualizeNetwork
''',
    "exadis": '''\
import pyexadis
from framework.disnet_manager import DisNetManager
from pyexadis_base import ExaDisNet, NodeConstraints, VisualizeNetwork
''',
}
# The tutorial's first line for each network (its options 1 and 2 for the simple source)
BUILD_TITLES = {"simple_source": "Create simple dislocation source using {code}",
                "paradis_file": "Read initial configuration from ParaDiS file",
                "line_config": "Create initial configuration of infinite lines",
                "prismatic_config": "Create initial configuration of prismatic loops"}


def _simple_source(c: SimpleSource, exadis: bool) -> str:
    """The tutorial's option 1 (PyDiS) or 2 (ExaDiS), with the spec's values."""
    if exadis:
        cell, free, pinned, nodes, segs, net = (
            "pyexadis.Cell", "NodeConstraints.UNCONSTRAINED", "NodeConstraints.PINNED_NODE",
            "nodes", "segs", "ExaDisNet(cell, nodes, segs)")
    else:
        cell, free, pinned, nodes, segs, net = (
            "Cell", "DisNode.Constraints.UNCONSTRAINED", "DisNode.Constraints.PINNED_NODE",
            "rn", "links", "DisNet(cell=cell, rn=rn, links=links)")
    pad = " " * len(f"{nodes} = np.array(")
    return (f"Ldis = {num(c.Ldis)} # dislocation length\n"
            f"Lbox = 2*Ldis # simulation box size\n"
            f"burg = {array(c.burg)} # Burgers vector\n"
            f"burg = burg/np.linalg.norm(burg) # of unit length, 1/sqrt(2)*[1,1,0] for [1,1,0]\n"
            f"plane = {array(c.plane)} # plane normal\n"
            f"\n"
            f"# Simulation cell object\n"
            f"cell = {cell}(h=Lbox*np.eye(3), is_periodic=[True,True,True])\n"
            f"\n"
            f"# List of nodes, nodes = [x,y,z,constraint]\n"
            f"linevec = Ldis*burg # line vector\n"
            f"{nodes} = np.array([[*(-0.5*linevec), {pinned}],\n"
            f"{pad}[*( 0.0*linevec), {free}],\n"
            f"{pad}[*( 0.5*linevec), {pinned}]])\n"
            f"{nodes}[:,0:3] += cell.center() # translate line to the center of the box\n"
            f"\n"
            f"# List of segments, {segs} = [node1,node2,burg,plane]\n"
            f"{segs} = np.array([[0, 1, *burg, *plane],\n"
            f"{' ' * len(f'{segs} = np.array(')}[1, 2, *burg, *plane]])\n"
            f"\n"
            f"G = {net}\n")


def _build_config(spec: Spec) -> str:
    """Step 1 of a build: the network, made as the tutorial makes it."""
    c, exadis = spec.config, spec.engine == "exadis"
    net = "ExaDisNet" if exadis else "DisNet"
    title = BUILD_TITLES[c.kind].format(code="pyexadis" if exadis else "pydis")
    if isinstance(c, SimpleSource):
        code = _simple_source(c, exadis)
    elif isinstance(c, ParaDiSFile):
        code = f"G = {net}().read_paradis(os.path.join(opendis_dir, {c.datafile!r}))\n"
    elif isinstance(c, InfiniteLines):
        code = (f"G = {net}().generate_line_config(\n"
                f"    crystal={c.crystal!r}, Lbox={num(c.Lbox)}, num_lines={c.num_lines}, "
                f"theta={num(c.theta)}, maxseg={num(c.maxseg)}, seed={c.seed}\n"
                f")\n")
    else:
        code = (f"G = {net}().generate_prismatic_config(\n"
                f"    crystal={c.crystal!r}, Lbox={num(c.Lbox)}, num_loops={c.num_loops}, "
                f"radius={num(c.radius)}, maxseg={num(c.maxseg)}, seed={c.seed}, "
                f"uniform={c.uniform}\n"
                f")\n")
    return f"print({title!r})\n\n{code}net = DisNetManager(G)\n"


def _build_run(spec: Spec) -> str:
    """Step 2 of a build: the plot, and the infinite lines' files, as the tutorial writes them."""
    code = "VisualizeNetwork().plot_disnet(net, block=True)\n"
    if isinstance(spec.config, InfiniteLines):
        code += "\nos.makedirs('output', exist_ok=True)\n"
        if spec.engine == "exadis":
            code += ("from pyexadis_utils import write_data, write_vtk\n"
                     "write_data(net, 'output/infinite_lines.data')\n"
                     "write_vtk(net, 'output/infinite_lines.vtk')\n")
        else:
            code += ("G.write_data('output/infinite_lines.data')\n"
                     "G.write_vtk('output/infinite_lines.vtk')\n")
    return code


def _main_guard(spec: Spec, final_json: str) -> str:
    """The end of the script: the run, then the final network written to ``output/``."""
    exadis = spec.engine == "exadis"
    return ('if __name__ == "__main__":\n'
            + ("    pyexadis.initialize()\n" if exadis else "")
            + "    main()\n"
            "\n"
            + ("    # explore the network after simulation\n" if spec.config.simulates
               else "    # explore the network\n")
            + ("    G  = net.get_disnet(ExaDisNet)\n" if exadis else "    G  = net.get_disnet()\n")
            + "\n"
            "    os.makedirs('output', exist_ok=True)\n"
            f"    net.write_json({'output/' + final_json!r})\n"
            + ("\n    if not sys.flags.interactive:\n        pyexadis.finalize()\n" if exadis
               else ""))


def header(spec: Spec) -> str:
    """The imports, with PyDiS (and ExaDiS) found from ``$OPENDIS_DIR``."""
    name, paths = (("pyexadis_paths", EXADIS_PATHS) if spec.engine == "exadis"
                   else ("pydis_paths", PYDIS_PATHS))
    imports = IMPORTS if spec.config.simulates else BUILD_IMPORTS
    return (_PATHS.format(opendis_dir=str(OPENDIS_DIR), name=name, paths=list(paths))
            + imports[spec.engine] + KIND_IMPORTS.get((spec.engine, spec.config.kind), ""))


def init_function(spec: Spec) -> str:
    """The tutorial's function that makes the spec's initial configuration, in its code."""
    return INIT_FUNCTIONS[spec.engine, spec.config.kind]


def stress_file(url: str) -> str:
    """Where the user stress's script keeps its copy of the link's field, from the root."""
    return f"{STRESS_DIR}/stress_{hashlib.sha256(url.encode()).hexdigest()}.npz"


def line_origin(box: float, seed: int) -> np.ndarray:
    """The point that the user stress's script starts its line through (init_infinite_line)."""
    return box * (np.random.default_rng(seed).random(3) - 0.5)


def user_stress_code(url: str) -> str:
    """The user stress's field, its link and its copy's file, and the force module that adds it."""
    return ("# The internal stress of a dislocation network, frozen, on a grid of its periodic\n"
            "# box: the file at this link, which the script downloads once, a copy for each link\n"
            f"stress_url = {url!r}\n"
            f"stress_file = os.path.join(opendis_dir, {STRESS_DIR!r},\n"
            "                           'stress_%s.npz' % "
            "hashlib.sha256(stress_url.encode()).hexdigest())\n\n" + USER_STRESS)


def _drivers(spec: Spec) -> str:
    """The script's own drivers: ExaDiS's that also writes JSON snapshots, and the user
    stress's that also writes the line's displacement."""
    exadis = spec.engine == "exadis"
    code = JSON_DRIVER + "\n" if exadis else ""
    if isinstance(spec.config, UserStressLine):
        code += DISP_DRIVER.format(base="SimulateNetworkJSON" if exadis else "SimulateNetwork")
        code += "\n"
    return code


def step_code(spec: Spec, step: str, final_json: str) -> str:
    """The code panel's "this step": the step's block, as it reads outside ``main()``."""
    if not spec.config.simulates:
        return (_build_config(spec) if step == "config"
                else _build_run(spec) + "\n" + _main_guard(spec, final_json))
    code = STEP_CODE[step](spec)
    if step == "config":
        return init_function(spec) + "\n" + code
    if step == "forces" and spec.forces.stress_field:
        return user_stress_code(spec.forces.stress_url) + "\n" + code
    if step == "run":
        return _drivers(spec) + code + "\n" + _main_guard(spec, final_json)
    return code


def full_script(spec: Spec, final_json: str) -> str:
    """The script the job runs, and that the code panel offers for download."""
    if not spec.config.simulates:
        blocks = "\n".join(textwrap.indent(code, "    ")
                           for code in (_build_config(spec), _build_run(spec)))
        return (header(spec) + "\n"
                "def main():\n"
                "    global net\n"
                "\n"
                f"{blocks}\n\n"
                f"{_main_guard(spec, final_json)}")
    blocks = "\n".join(textwrap.indent(STEP_CODE[step](spec), "    ") for step in STEP_CODE)
    exadis = spec.engine == "exadis"
    return (header(spec) + "\n" + init_function(spec) + "\n"
            + (user_stress_code(spec.forces.stress_url) + "\n" if spec.forces.stress_field else "")
            + _drivers(spec)
            + "def main():\n"
            "    global net, sim, state\n"
            "\n"
            f"{blocks}"
            + ("" if exadis else "\n    return net.is_sane()\n")
            + "\n\n"
            f"{_main_guard(spec, final_json)}")


def init_code(spec: Spec) -> str:
    """Code that builds ``net``, the initial configuration, for the app to draw it.

    PyDiS builds it whichever code runs the script: the same values give the same network.
    The app puts PyDiS on ``sys.path`` itself, so this has no path lines.
    """
    if not spec.config.simulates:
        return ("import os\n"
                "import numpy as np\n"
                + BUILD_IMPORTS["pydis"]
                + f"opendis_dir = {str(OPENDIS_DIR)!r}\n\n"
                + _build_config(replace(spec, engine="pydis")))
    return ("import numpy as np\n"
            "from framework.disnet_manager import DisNetManager\n"
            "from pydis import DisNode, DisNet, Cell\n"
            + KIND_IMPORTS.get(("pydis", spec.config.kind), "") + "\n"
            + INIT_FUNCTIONS["pydis", spec.config.kind] + "\n" + _config(replace(spec, engine="pydis")))
