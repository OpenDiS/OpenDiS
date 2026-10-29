"""common.codegen: the PyDiS script of a spec, step by step and in full.

The full script is the tutorial's script, ``examples/02_frank_read_src/test_frank_read_src_pydis.py``
or ``examples/03_binary_junction/test_binary_junction_pydis.py``, with the spec's values. For an
example's own values it runs the same simulation as the tutorial's script, which
``tests/test_codegen.py`` checks. It differs from that script in three ways: it finds PyDiS
from ``$OPENDIS_DIR`` (else from the repository the app runs in), so that it runs from any
folder; each step's code is one block, so a few lines come in another order; and an optional
value (Ec, mob, a Burgers vector) appears only when it is not the default.
"""
import textwrap

from common.paths import OPENDIS_DIR, PYDIS_PATHS
from common.spec import FrankReadSource, Spec, TwoLines


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


HEADER = '''\
import numpy as np
import sys, os

# The root of the OpenDiS repository: $OPENDIS_DIR, else the one that the app ran from
opendis_dir = os.environ.get('OPENDIS_DIR', {opendis_dir!r})
pydis_paths = [os.path.join(opendis_dir, path) for path in {pydis_paths!r}]
[sys.path.append(os.path.abspath(path)) for path in pydis_paths
 if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork, SimulateNetwork
'''

# The tutorials' functions that make the initial configuration, word for word
INIT_FUNCTIONS = {
    FrankReadSource.kind: '''\
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
    TwoLines.kind: '''\
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
    return (f"Lbox = {num(c.box_length)}\n"
            f"z0 = {num(c.z0)}\n"
            f"net = init_two_disl_lines(z0=z0, box_length=Lbox, "
            f"pbc={c.pbc}{_optional_burgers(c, ('b1', 'b2'))})\n")


def _material(spec: Spec) -> str:
    """Step 2: the state dictionary."""
    m = spec.material
    return (f'state = {{"burgmag": {num(m.burgmag)}, "mu": {num(m.mu)}, "nu": {num(m.nu)}, '
            f'"a": {num(m.a)},\n'
            f'         "maxseg": {num(m.maxseg_ratio)}*Lbox, "minseg": {num(m.minseg_ratio)}*Lbox, '
            f'"rann": {num(m.rann)}}}\n')


def _forces(spec: Spec) -> str:
    """Step 3: the force and the mobility modules."""
    f = spec.forces
    ec = "" if f.Ec is None else f", Ec={num(f.Ec)}"
    mob = "" if f.mob is None else f'state["mob"] = {num(f.mob)}\n'
    return (f"calforce  = CalForce(force_mode={f.force_mode!r}, state=state{ec})\n"
            f"{mob}"
            f"mobility  = MobilityLaw(mobility_law={f.mobility_law!r}, state=state)\n")


def _timeint(spec: Spec) -> str:
    """Step 4: the time integrator."""
    t = spec.timeint
    return (f"timeint   = TimeIntegration(integrator={t.integrator!r}, dt={num(t.dt)},\n"
            f"                            state=state)\n")


def _topology(spec: Spec) -> str:
    """Step 5: the cell list, and the topology, collision and remesh modules."""
    t = spec.topology
    n_div = ",".join(str(n) for n in t.n_div)
    return (f"nbrlist   = CellList(cell=net.cell, n_div=[{n_div}])\n"
            f"topology  = Topology(split_mode={t.split_mode!r}, state=state,\n"
            f"                     force=calforce, mobility=mobility)\n"
            f"collision = Collision(collision_mode={t.collision_mode!r}, state=state,\n"
            f"                      nbrlist=nbrlist)\n"
            f"remesh    = Remesh(remesh_rule={t.remesh_rule!r}, state=state)\n")


def _loading(spec: Spec) -> str:
    """Step 6: the applied stress."""
    return f"applied_stress = {array(spec.loading.applied_stress)}\n"


def _run(spec: Spec) -> str:
    """Step 7: the simulation, and its run."""
    r, loading_mode = spec.run, spec.loading.loading_mode
    return (f"vis = VisualizeNetwork()\n"
            f"sim = SimulateNetwork(calforce=calforce, mobility=mobility,\n"
            f"                      timeint=timeint, topology=topology,\n"
            f"                      collision=collision, remesh=remesh, vis=vis,\n"
            f"                      state=state, max_step={r.max_step}, "
            f"loading_mode={loading_mode!r},\n"
            f"                      applied_stress=applied_stress,\n"
            f"                      print_freq={r.print_freq}, plot_freq={r.plot_freq}, "
            f"plot_pause_seconds={num(r.plot_pause_seconds)},\n"
            f"                      write_freq={r.write_freq}, write_dir='output', "
            f"save_state=False)\n"
            f"sim.run(net, state)\n")


STEP_CODE = {"config": _config, "material": _material, "forces": _forces, "timeint": _timeint,
             "topology": _topology, "loading": _loading, "run": _run}


def _main_guard(final_json: str) -> str:
    """The end of the script: the run, then the final network written to ``output/``."""
    return ('if __name__ == "__main__":\n'
            "    main()\n"
            "\n"
            "    # explore the network after simulation\n"
            "    G  = net.get_disnet()\n"
            "\n"
            "    os.makedirs('output', exist_ok=True)\n"
            f"    net.write_json({'output/' + final_json!r})\n")


def header() -> str:
    """The imports, with PyDiS found from ``$OPENDIS_DIR``."""
    return HEADER.format(opendis_dir=str(OPENDIS_DIR), pydis_paths=list(PYDIS_PATHS))


def step_code(spec: Spec, step: str, final_json: str) -> str:
    """The code panel's "this step": the step's block, as it reads outside ``main()``."""
    code = STEP_CODE[step](spec)
    if step == "config":
        return INIT_FUNCTIONS[spec.config.kind] + "\n" + code
    if step == "run":
        return code + "\n" + _main_guard(final_json)
    return code


def full_script(spec: Spec, final_json: str) -> str:
    """The script the job runs, and that the code panel offers for download."""
    blocks = "\n".join(textwrap.indent(STEP_CODE[step](spec), "    ") for step in STEP_CODE)
    return (header() + "\n" + INIT_FUNCTIONS[spec.config.kind] + "\n"
            "def main():\n"
            "    global net, sim, state\n"
            "\n"
            f"{blocks}\n"
            "    return net.is_sane()\n"
            "\n\n"
            f"{_main_guard(final_json)}")


def init_code(spec: Spec) -> str:
    """Code that builds ``net``, the initial configuration, for the app to draw it.

    The app puts PyDiS on ``sys.path`` itself, so this has no path lines.
    """
    return ("import numpy as np\n"
            "from framework.disnet_manager import DisNetManager\n"
            "from pydis import DisNode, DisNet, Cell\n\n"
            + INIT_FUNCTIONS[spec.config.kind] + "\n" + _config(spec))
