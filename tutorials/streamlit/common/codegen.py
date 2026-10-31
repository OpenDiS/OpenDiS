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

A script that builds a network and runs no simulation is one option of
``tutorials/01_initial_configurations/initial_configurations.py``, with either code: the
network that it builds is the tutorial's, which ``tests/test_codegen.py`` checks too. It writes
its files into ``output/``, the network as its final one.
"""
import textwrap
from dataclasses import replace

from common.paths import OPENDIS_DIR, PYDIS_PATHS
from common.spec import (FrankReadSource, InfiniteLines, LoopSources, ParaDiSFile,
                         SimpleSource, Spec, TwoLines)

EXADIS_PATHS = PYDIS_PATHS + ("core/exadis/python",)  # as in the ExaDiS tutorials' scripts


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
}
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
    return (f"Lbox = {num(c.box_length)}\n"
            f"z0 = {num(c.z0)}\n"
            f"net = init_two_disl_lines(z0=z0, box_length=Lbox, "
            f"pbc={c.pbc}{_optional_burgers(c, ('b1', 'b2'))})\n")


def _material(spec: Spec) -> str:
    """Step 2: the state dictionary."""
    m = spec.material
    extras = ", ".join(f'"{name}": {value!r}' if isinstance(value, str) else f'"{name}": {num(value)}'
                       for name, value in (("crystal", m.crystal), ("rtol", m.rtol),
                                           ("nextdt", m.nextdt)) if value is not None)
    return (f'state = {{"burgmag": {num(m.burgmag)}, "mu": {num(m.mu)}, "nu": {num(m.nu)}, '
            f'"a": {num(m.a)},\n'
            f'         "maxseg": {num(m.maxseg_ratio)}*Lbox, "minseg": {num(m.minseg_ratio)}*Lbox, '
            f'"rann": {num(m.rann)}' + (f",\n         {extras}" if extras else "") + "}\n")


def _forces(spec: Spec) -> str:
    """Step 3: the force and the mobility modules."""
    f = spec.forces
    ec = "" if f.Ec is None else f", Ec={num(f.Ec)}"
    if f.mobility_law == "BCC_0B":
        vmax = "" if f.vmax is None else f", vmax={num(f.vmax)}"
        return (f"calforce  = CalForce(force_mode={f.force_mode!r}, state=state{ec})\n"
                f"mobility  = MobilityLaw(mobility_law='BCC_0B', state=state, Medge={num(f.Medge)},\n"
                f"                        Mscrew={num(f.Mscrew)}, Mclimb={num(f.Mclimb)}{vmax})\n")
    if spec.engine == "exadis":  # ExaDiS's mobility takes mob as an argument
        mob = "" if f.mob is None else f", mob={num(f.mob)}"
        return (f"calforce  = CalForce(force_mode={f.force_mode!r}, state=state{ec})\n"
                f"mobility  = MobilityLaw(mobility_law={f.mobility_law!r}, state=state{mob})\n")
    mob = "" if f.mob is None else f'state["mob"] = {num(f.mob)}\n'
    return (f"calforce  = CalForce(force_mode={f.force_mode!r}, state=state{ec})\n"
            f"{mob}"
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
    return (f"nbrlist   = CellList(cell=net.cell, n_div=[{n_div}])\n"
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
    return f"applied_stress = {array(loading.applied_stress)}\n"


def _run(spec: Spec) -> str:
    """Step 7: the simulation, and its run."""
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
            + imports[spec.engine])


def init_function(spec: Spec) -> str:
    """The tutorial's function that makes the spec's initial configuration, in its code."""
    return INIT_FUNCTIONS[spec.engine, spec.config.kind]


def step_code(spec: Spec, step: str, final_json: str) -> str:
    """The code panel's "this step": the step's block, as it reads outside ``main()``."""
    if not spec.config.simulates:
        return (_build_config(spec) if step == "config"
                else _build_run(spec) + "\n" + _main_guard(spec, final_json))
    code = STEP_CODE[step](spec)
    if step == "config":
        return init_function(spec) + "\n" + code
    if step == "run":
        driver = JSON_DRIVER + "\n" if spec.engine == "exadis" else ""
        return driver + code + "\n" + _main_guard(spec, final_json)
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
            + (JSON_DRIVER + "\n" if exadis else "")
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
            "from pydis import DisNode, DisNet, Cell\n\n"
            + INIT_FUNCTIONS["pydis", spec.config.kind] + "\n" + _config(replace(spec, engine="pydis")))
