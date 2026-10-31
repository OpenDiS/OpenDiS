"""common.steps: each step's selection widgets, which return the spec with the step's part changed.

The widget keys start with the context's prefix, which changes whenever code replaces the spec
(a reset), so that the widgets start again from the new values. The values are checked by
``Spec.validate``, whose messages the app shows below the widgets. The choices are those of the
code that runs the example, PyDiS or ExaDiS (``CHOICES`` in ``common/spec.py``).
"""
import math
from dataclasses import dataclass, replace
from typing import Any

import streamlit as st

from common.catalog import Example
from common.spec import (BUILT, CHOICES, CRYSTALS, ENGINE_NAMES, MAX_BUILT, MAX_STEP,
                         PARADIS_FILES, STRESS_LABELS, FrankReadSource, InfiniteLines,
                         LoopSources, ParaDiSFile, PrismaticLoops, Run, SimpleSource, Spec)

STRESS_STEP = 1.0e7  # Pa, the +/- of an applied stress component
TUTORIAL_BY = {"pydis": "by pure Python", "exadis": "by Python calling ExaDiS"}  # the pages' titles
# What each choice does, for the help of its select box
HELP = {
    "LineTension": "each segment pulls its nodes along the line, as in the tutorials",
    "Elasticity_SBA": "the elastic interaction of all the segments (slower)",
    "SimpleGlide": "the nodes move in their glide planes",
    "Relax": "they move along their forces",
    "None": "no topology module, as in ExaDiS's tutorial",
    "MaxDiss": "splits each multi-node the way that dissipates the most energy, the two new "
               "nodes on top of each other",
    "Serial": "the same, with the two new nodes moved apart, as in ParaDiS and ExaDiS's "
              "TopologySerial",
    "TopologySerial": "splits each multi-node the way that dissipates the most energy",
    "TopologyParallel": "the same, in parallel",
    "BCC_0B": "bcc's edge, screw and climb mobilities, and a smooth speed cap",
    "EulerForward": "a fixed time step",
    "Trapezoid": "a time step that adapts, from nextdt, to keep the error below rtol (step 2)",
    "stress": "a constant applied stress",
    "strain_rate": "uniaxial loading at a constant strain rate along edir; the stress follows "
                   "from the plastic strain",
}


def _help(choices: tuple[str, ...]) -> str:
    """The help of a select box: what each of its choices does."""
    return " ".join(f"{choice}: {HELP[choice]}." for choice in choices)


@dataclass(frozen=True)
class StepContext:
    """What a step's widgets get from the app."""

    key: str          # widget-key prefix; it changes when code replaces the spec
    example: Example


def _step(value: float) -> float:
    """The +/- of a number input: a tenth of the value's power of ten (0.1 for zero)."""
    return 10.0 ** (math.floor(math.log10(abs(value))) - 1) if value else 0.1


def _number(label: str, value: float, key: str, help_text: str | None = None,
            step: float | None = None) -> float:
    """A number input for a float, shown with ``%g``."""
    return float(st.number_input(label, value=float(value), step=step or _step(value),
                                 format="%g", key=key, help=help_text))


def _integer(label: str, value: int, key: str, help_text: str | None = None,
             max_value: int | None = None) -> int:
    """A number input for a whole number, 1 or more."""
    return int(st.number_input(label, value=int(value), min_value=1, max_value=max_value,
                               step=1, key=key, help=help_text))


def _vector(label: str, value: tuple, key: str, kind: type = float,
            help_text: str | None = None) -> tuple:
    """A text input of numbers separated by commas; the old value if the text is not that."""
    text = st.text_input(label, value=", ".join(f"{v:g}" for v in value), key=key,
                         help=help_text)
    try:
        numbers = tuple(kind(float(part)) for part in text.split(","))
    except ValueError:
        numbers = ()
    if len(numbers) != len(value) or (kind is int and any(
            float(part) != int(float(part)) for part in text.split(","))):
        st.error(f"{label}: give {len(value)} {'whole ' if kind is int else ''}numbers, "
                 "separated by commas")
        return value
    return numbers


def _select(label: str, options: tuple[str, ...], value: str, key: str,
            help_text: str | None = None) -> str:
    """A select box, on the spec's value."""
    return st.selectbox(label, options, index=options.index(value) if value in options else 0,
                        key=key, help=help_text, disabled=len(options) == 1)


def _config(spec: Spec, ctx: StepContext) -> Spec:
    """Step 1: the example's text, the box, and the initial configuration's own values."""
    example, c, k = ctx.example, spec.config, ctx.key
    with st.expander("About this example", expanded=True):
        st.markdown(example.body)
        if example.tutorial_url:
            title = example.tutorial_title or f"{example.label} {TUTORIAL_BY[example.engine]}"
            st.markdown(f"From the tutorial [{title}]({example.tutorial_url}).")
    if not c.simulates:
        return spec.with_part("config", _network(c, k))
    box = _number("Box length, Lbox (b)", c.box_length, f"{k}/box_length",
                  "The edge of the cubic simulation cell, in units of the Burgers vector's length.")
    if isinstance(c, LoopSources):
        st.caption("The box is periodic in all three directions.")
        count = _integer("Number of loops", c.num_loops, f"{k}/num_loops",
                         "They cycle through the crystal's Burgers vectors (step 2's crystal).",
                         max_value=MAX_BUILT)
        radius = _number("Radius / Lbox", c.radius_ratio, f"{k}/radius_ratio",
                         "The loops' radius, as a fraction of the box.", 0.01)
        st.caption(f"radius = {radius * box:g} b; the segments are at most step 2's maxseg")
        seed = _integer("Seed of the random positions", c.seed, f"{k}/seed")
        uniform = st.checkbox("Spread the loops evenly", c.uniform, key=f"{k}/uniform",
                              help="Their centres on a grid, shuffled and jittered, rather than "
                                   "at random.")
        return spec.with_part("config", replace(c, box_length=box, num_loops=count,
                                                radius_ratio=radius, seed=seed, uniform=uniform))
    pbc = st.checkbox("Periodic boundary conditions", c.pbc, key=f"{k}/pbc")
    if isinstance(c, FrankReadSource):
        arm = _number("Arm length / Lbox", c.arm_ratio, f"{k}/arm_ratio",
                      "The loop is a square of this side; its top arm is the source.", 0.005)
        st.caption(f"arm length = {arm * box:g} b")
        burg = _vector("Burgers vector b", c.burg_vec, f"{k}/burg_vec",
                       help_text="The loop is in the x = 0 plane, so b along x makes it "
                                 "prismatic, as in the tutorial.")
        return spec.with_part("config", replace(c, box_length=box, pbc=pbc, arm_ratio=arm,
                                                burg_vec=burg))
    z0 = _number("z0 (b)", c.z0, f"{k}/z0",
                 "The lines' pinned ends are at ±z0 from the cell's centre.")
    b1 = _vector("b1, the first line's Burgers vector", c.b1, f"{k}/b1")
    b2 = _vector("b2, the second line's Burgers vector", c.b2, f"{k}/b2")
    return spec.with_part("config", replace(c, box_length=box, pbc=pbc, z0=z0, b1=b1, b2=b2))


def _network(c, k: str):
    """Step 1 of a script that builds a network: which of the tutorial's networks, and its
    values. Another network starts from the tutorial's values."""
    labels = tuple(cls.label for cls in BUILT)
    label = _select("Network", labels, c.label, f"{k}/network",
                    "The tutorial's ways to build a network; each runs with PyDiS or ExaDiS.")
    cls = BUILT[labels.index(label)]
    c = c if isinstance(c, cls) else cls()
    k = f"{k}/{c.kind}"  # the networks share names, not values
    if isinstance(c, SimpleSource):
        ldis = _number("Ldis, the dislocation's length (b)", c.Ldis, f"{k}/Ldis",
                       "The box is a cube twice as long, periodic in all three directions.")
        st.caption(f"Lbox = 2 Ldis = {2 * ldis:g} b")
        burg = _vector("Burgers vector b, the line's direction", c.burg, f"{k}/burg",
                       help_text="The script makes it of unit length. The dislocation lies "
                                 "along it: a screw dislocation.")
        plane = _vector("Plane normal n", c.plane, f"{k}/plane",
                        help_text="Normal to b: the dislocation's glide plane.")
        return SimpleSource(ldis, burg, plane)
    if isinstance(c, ParaDiSFile):
        datafile = _select("Data file", PARADIS_FILES, c.datafile, f"{k}/datafile",
                           "The repository's ParaDiS data files, from its root.")
        return ParaDiSFile(datafile)
    crystal = _select("Crystal", CRYSTALS, c.crystal, f"{k}/crystal",
                      "fcc: the 12 <110>{111} slip systems. bcc: the 12 <111>{110}.")
    box = _number("Lbox, the box's edge (b)", c.Lbox, f"{k}/Lbox",
                  "A cube, periodic in all three directions.")
    if isinstance(c, InfiniteLines):
        count = _integer("Number of lines", c.num_lines, f"{k}/num_lines",
                         "They cycle through the slip systems, their signs alternating every "
                         "12 lines (dipoles).", max_value=MAX_BUILT)
        theta = _number("theta, the lines' character angle (degrees)", c.theta, f"{k}/theta",
                        "0: screw. 90: edge.", step=5.0)
        maxseg = _number("maxseg, the longest segment (b)", c.maxseg, f"{k}/maxseg")
        seed = _integer("Seed of the random positions", c.seed, f"{k}/seed")
        return InfiniteLines(crystal, box, count, theta, maxseg, seed)
    count = _integer("Number of loops", c.num_loops, f"{k}/num_loops",
                     "They cycle through the crystal's Burgers vectors.", max_value=MAX_BUILT)
    radius = _number("Radius of the loops (b)", c.radius, f"{k}/radius")
    maxseg = _number("maxseg, the longest segment (b)", c.maxseg, f"{k}/maxseg")
    seed = _integer("Seed of the random positions", c.seed, f"{k}/seed")
    uniform = st.checkbox("Spread the loops evenly", c.uniform, key=f"{k}/uniform",
                          help="Their centres on a grid, shuffled and jittered, rather than "
                               "at random.")
    return PrismaticLoops(crystal, box, count, radius, maxseg, seed, uniform)


def _material(spec: Spec, ctx: StepContext) -> Spec:
    """Step 2: the material and the lengths, the entries of ``state``."""
    m, k, box = spec.material, ctx.key, spec.config.box_length
    st.markdown("**Material**")
    burgmag = _number("burgmag, the length of b (m)", m.burgmag, f"{k}/burgmag")
    mu = _number("mu, the shear modulus (Pa)", m.mu, f"{k}/mu")
    nu = _number("nu, Poisson's ratio", m.nu, f"{k}/nu", step=0.01)
    a = _number("a, the core radius (b)", m.a, f"{k}/a")
    st.markdown("**Lengths**")
    maxseg = _number("maxseg / Lbox", m.maxseg_ratio, f"{k}/maxseg_ratio",
                     "Remeshing splits the segments longer than maxseg.", 0.005)
    minseg = _number("minseg / Lbox", m.minseg_ratio, f"{k}/minseg_ratio",
                     "Remeshing removes the nodes of segments shorter than minseg.", 0.001)
    rann = _number("rann, the annihilation distance (b)", m.rann, f"{k}/rann",
                   "Segments closer than rann collide.")
    st.caption(f"maxseg = {maxseg * box:g} b, minseg = {minseg * box:g} b")
    extra = {}
    if m.crystal is not None:
        extra["crystal"] = _select("Crystal", CRYSTALS, m.crystal, f"{k}/crystal",
                                   "state[\"crystal\"], which BCC_0B and the loops need.")
    if m.rtol is not None or m.nextdt is not None:
        st.markdown("**Time step control** (Trapezoid)")
    if m.rtol is not None:
        extra["rtol"] = _number("rtol, the error tolerance (b)", m.rtol, f"{k}/rtol",
                                "The Trapezoid integrator halves the time step until the "
                                "positions' error is below rtol.")
    if m.nextdt is not None:
        extra["nextdt"] = _number("nextdt, the first time step (s)", m.nextdt, f"{k}/nextdt")
    return spec.with_part("material", replace(m, burgmag=burgmag, mu=mu, nu=nu, a=a,
                                              maxseg_ratio=maxseg, minseg_ratio=minseg,
                                              rann=rann, **extra))


def _forces(spec: Spec, ctx: StepContext) -> Spec:
    """Step 3: the force model and the mobility law, with their optional values."""
    f, m, k = spec.forces, spec.material, ctx.key
    choices, name = CHOICES[spec.engine], ENGINE_NAMES[spec.engine]
    mode = _select("Force", choices["force_mode"], f.force_mode, f"{k}/force_mode",
                   _help(choices["force_mode"]))
    default_ec = m.mu / (4.0 * math.pi) * math.log(m.a / 0.1)
    if st.checkbox("Set the core energy Ec", f.Ec is not None, key=f"{k}/own_ec"):
        ec = _number("Ec (Pa)", f.Ec if f.Ec is not None else default_ec or 1.0, f"{k}/Ec")
    else:
        ec = None
        st.caption(f"Ec = mu/(4π) ln(a/0.1) = {default_ec:.4g} Pa, {name}'s default"
                   + ("; it is not positive, so set Ec" if default_ec <= 0 else ""))
    law = _select("Mobility law", choices["mobility_law"], f.mobility_law, f"{k}/mobility_law",
                  _help(choices["mobility_law"]))
    if law == "BCC_0B":
        medge = _number("Medge, the edge mobility (1/(Pa s))", f.Medge or 2600.0, f"{k}/Medge")
        mscrew = _number("Mscrew, the screw mobility (1/(Pa s))", f.Mscrew or 2600.0,
                         f"{k}/Mscrew")
        mclimb = _number("Mclimb, the climb mobility (1/(Pa s))", f.Mclimb or 1e-4,
                         f"{k}/Mclimb")
        vmax = _number("vmax, the speed cap (m/s)", f.vmax or 3400.0, f"{k}/vmax",
                       "The speed tends to vmax smoothly as the force grows.")
        return spec.with_part("forces", replace(f, force_mode=mode, Ec=ec, mobility_law=law,
                                                Medge=medge, Mscrew=mscrew, Mclimb=mclimb,
                                                vmax=vmax))
    if st.checkbox("Set the mobility mob", f.mob is not None, key=f"{k}/own_mob"):
        mob = _number("mob", f.mob if f.mob is not None else 1.0, f"{k}/mob")
    else:
        mob = None
        st.caption(f"mob = 1.0, {name}'s default")
    return spec.with_part("forces", replace(f, force_mode=mode, Ec=ec, mobility_law=law,
                                            mob=mob))


def _timeint(spec: Spec, ctx: StepContext) -> Spec:
    """Step 4: the integrator and its time step."""
    t, k = spec.timeint, ctx.key
    choices = CHOICES[spec.engine]["integrator"]
    integrator = _select("Integrator", choices, t.integrator, f"{k}/integrator", _help(choices))
    if integrator == "Trapezoid":
        multi = int(st.number_input("multi: Trapezoid steps in a row, up to maxdt", value=t.multi,
                                    min_value=0, step=1, key=f"{k}/multi",
                                    help="0 or 1: one step at a time."))
        st.caption("The time step adapts: it starts from nextdt, and the error stays below "
                   "rtol (both in step 2).")
        return spec.with_part("timeint", replace(t, integrator=integrator, multi=multi))
    dt = _number("dt, the time step (s)", t.dt, f"{k}/dt")
    return spec.with_part("timeint", replace(t, integrator=integrator, dt=dt))


def _topology(spec: Spec, ctx: StepContext) -> Spec:
    """Step 5: the topological operations, collisions, remeshing and the cell list."""
    t, k, choices = spec.topology, ctx.key, CHOICES[spec.engine]
    if spec.engine == "exadis":
        mode = _select("Topology", choices["topology_mode"], t.topology_mode,
                       f"{k}/topology_mode", _help(choices["topology_mode"]))
        collision = _select("Collision", choices["collision_mode"], t.collision_mode,
                            f"{k}/collision_mode", "Retroactive: segments closer than rann "
                            "collide, and those whose paths crossed during the step. In ExaDiS, "
                            "Proximity is the same.")
        remesh = _select("Remesh rule", choices["remesh_rule"], t.remesh_rule,
                         f"{k}/remesh_rule")
        return spec.with_part("topology", replace(t, topology_mode=mode,
                                                  collision_mode=collision, remesh_rule=remesh))
    split = _select("Split mode", choices["split_mode"], t.split_mode, f"{k}/split_mode",
                    _help(choices["split_mode"]))
    collision = _select("Collision", choices["collision_mode"], t.collision_mode,
                        f"{k}/collision_mode", "Proximity: segments closer than rann collide. "
                        "Retroactive: also the segments whose paths crossed during the step.")
    remesh = _select("Remesh rule", choices["remesh_rule"], t.remesh_rule, f"{k}/remesh_rule")
    n_div = _vector("Cell list divisions, n_div", t.n_div, f"{k}/n_div", kind=int,
                    help_text="The cell list's cells along x, y and z, for finding the "
                              "segments that may collide.")
    return spec.with_part("topology", replace(t, split_mode=split, collision_mode=collision,
                                              remesh_rule=remesh, n_div=n_div))


def _loading(spec: Spec, ctx: StepContext) -> Spec:
    """Step 6: the applied stress."""
    loading, k = spec.loading, ctx.key
    choices = CHOICES[spec.engine]["loading_mode"]
    mode = _select("Loading", choices, loading.loading_mode, f"{k}/loading_mode", _help(choices))
    if mode == "strain_rate":
        erate = _number("erate, the strain rate (1/s)", loading.erate or 1e3, f"{k}/erate")
        edir = _vector("edir, the loading direction", loading.edir, f"{k}/edir",
                       help_text="The script makes it of unit length.")
        return spec.with_part("loading", replace(loading, loading_mode=mode, erate=erate,
                                                 edir=edir))
    st.markdown("**Applied stress (Pa)**")
    stress = tuple(_number(f"σ_{label}", value, f"{k}/sigma_{label}", step=STRESS_STEP)
                   for label, value in zip(STRESS_LABELS, loading.applied_stress))
    return spec.with_part("loading", replace(loading, loading_mode=mode, applied_stress=stress))


def run_settings(spec: Spec, ctx: StepContext) -> Spec:
    """Step 7's values: how long the run is, and how often it prints, writes and plots."""
    r, k = spec.run, ctx.key
    max_step = _integer("max_step, the number of time steps", r.max_step, f"{k}/max_step",
                        max_value=MAX_STEP)
    print_freq = _integer("print_freq: print a line every", r.print_freq, f"{k}/print_freq",
                          "The output column shows these lines.")
    write_freq = _integer("write_freq: write a snapshot every", r.write_freq,
                          f"{k}/write_freq", "The snapshots are what the simulation view plays.")
    with st.expander("Plots, when the script runs outside the app"):
        plot_freq = _integer("plot_freq: plot every", r.plot_freq, f"{k}/plot_freq")
        pause = _number("plot_pause_seconds", r.plot_pause_seconds, f"{k}/plot_pause_seconds",
                        step=0.01)
    return spec.with_part("run", Run(max_step, print_freq, write_freq, plot_freq, pause))


RENDER: dict[str, Any] = {"config": _config, "material": _material, "forces": _forces,
                          "timeint": _timeint, "topology": _topology, "loading": _loading}


def render(step: str, spec: Spec, ctx: StepContext) -> Spec:
    """The step's widgets; returns the spec with the user's changes."""
    return RENDER[step](spec, ctx)
