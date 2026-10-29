"""common.spec: an example's selections, the spec, one part per step.

A spec holds every value that the generated script takes, in one frozen dataclass per step,
so that a step's widgets replace only their own part. ``validate`` lists what would make the
script fail, or the run meaningless. The module imports neither Streamlit nor PyDiS, so the
tests can use it.
"""
import json
import math
import types
from collections.abc import Mapping
from dataclasses import MISSING, asdict, dataclass, fields, replace
from typing import Any, ClassVar, Union, get_args, get_origin

import numpy as np

Vec3 = tuple[float, float, float]
Stress = tuple[float, float, float, float, float, float]

# (key, title) of each step, in order; a key is also the name of the spec's part for that step
STEPS = (
    ("config", "Initial configuration"),
    ("material", "Material and lengths"),
    ("forces", "Forces and mobility"),
    ("timeint", "Time integration"),
    ("topology", "Topology, collision and remesh"),
    ("loading", "Loading"),
    ("run", "Run"),
)
PARTS = tuple(key for key, _ in STEPS)

# The choices that PyDiS's modules accept, as the app offers them
FORCE_MODES = ("LineTension", "Elasticity_SBA")
MOBILITY_LAWS = ("SimpleGlide", "Relax")
INTEGRATORS = ("EulerForward",)
SPLIT_MODES = ("MaxDiss",)
COLLISION_MODES = ("Proximity", "Retroactive")
REMESH_RULES = ("LengthBased",)
LOADING_MODES = ("stress",)
STRESS_LABELS = ("xx", "yy", "zz", "yz", "xz", "xy")  # PyDiS's Voigt order
MAX_STEP = 100_000
MAX_N_DIV = 64


class InputError(ValueError):
    """An input file's problems, one message each."""

    def __init__(self, problems: list[str]) -> None:
        super().__init__("; ".join(problems))
        self.problems = problems


def _parallel(b: Vec3, xi: Vec3) -> bool:
    """True when the Burgers vector is parallel to the line direction (no glide plane)."""
    return float(np.linalg.norm(np.cross(b, xi))) <= 1e-12 * np.linalg.norm(b) * np.linalg.norm(xi)


def _burgers_problems(label: str, b: Vec3, lines: Mapping[str, Vec3]) -> list[str]:
    """The problems of a Burgers vector: zero, or parallel to a line it is on."""
    if not any(b):
        return [f"{label}: must not be zero"]
    return [f"{label}: must not be parallel to the {name}, or that segment has no glide plane"
            for name, xi in lines.items() if _parallel(b, xi)]


@dataclass(frozen=True)
class FrankReadSource:
    """Step 1 of a Frank-Read source: a prismatic loop whose top arm has a free middle node."""

    kind: ClassVar[str] = "frank_read_src_loop"
    box_length: float
    pbc: bool
    arm_ratio: float                  # arm_length = arm_ratio * Lbox
    burg_vec: Vec3 = (1.0, 0.0, 0.0)

    def problems(self) -> list[str]:
        """What is wrong with these values."""
        errors = []
        if not 0.0 < self.arm_ratio <= 0.5:
            errors.append("arm length / box length: must be in (0, 0.5], so that the loop "
                          "fits in the box")
        lines = {"arms along y": (0.0, 1.0, 0.0), "arms along z": (0.0, 0.0, 1.0)}
        return errors + _burgers_problems("Burgers vector", self.burg_vec, lines)


@dataclass(frozen=True)
class TwoLines:
    """Step 1 of a binary junction: two lines, pinned at their ends, that cross at the centre."""

    kind: ClassVar[str] = "two_disl_lines"
    box_length: float
    pbc: bool
    z0: float
    b1: Vec3 = (-1.0, 1.0, 1.0)
    b2: Vec3 = (1.0, -1.0, 1.0)

    def problems(self) -> list[str]:
        """What is wrong with these values."""
        errors = []
        if not 0.0 < self.z0 < 0.5 * self.box_length:
            errors.append("z0: must be in (0, box length / 2), so that the lines fit in the box")
        errors += _burgers_problems("b1", self.b1, {"first line": (0.0, 1.0, 1.0)})
        return errors + _burgers_problems("b2", self.b2, {"second line": (1.0, 0.0, 1.0)})


CONFIGS = {cls.kind: cls for cls in (FrankReadSource, TwoLines)}
Config = Union[FrankReadSource, TwoLines]


@dataclass(frozen=True)
class Material:
    """Step 2: the material and the lengths, the entries of the state dictionary."""

    burgmag: float       # m
    mu: float            # Pa
    nu: float
    a: float             # core radius, in b
    maxseg_ratio: float  # maxseg = maxseg_ratio * Lbox
    minseg_ratio: float  # minseg = minseg_ratio * Lbox
    rann: float          # annihilation distance, in b


@dataclass(frozen=True)
class Forces:
    """Step 3: the force model and the mobility law."""

    force_mode: str = "LineTension"
    Ec: float | None = None   # core energy; None is PyDiS's default, mu/(4 pi) ln(a/0.1)
    mobility_law: str = "SimpleGlide"
    mob: float | None = None  # mobility, state["mob"]; None is PyDiS's default, 1.0


@dataclass(frozen=True)
class TimeInt:
    """Step 4: the time integrator and its time step."""

    integrator: str = "EulerForward"
    dt: float = 1.0e-8  # s


@dataclass(frozen=True)
class Topology:
    """Step 5: the topological operations, collisions, remeshing and the cell list."""

    split_mode: str = "MaxDiss"
    collision_mode: str = "Proximity"
    remesh_rule: str = "LengthBased"
    n_div: tuple[int, int, int] = (8, 8, 8)


@dataclass(frozen=True)
class Loading:
    """Step 6: the applied stress, constant in time."""

    loading_mode: str = "stress"
    applied_stress: Stress = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)  # Pa, Voigt order


@dataclass(frozen=True)
class Run:
    """Step 7: how long the run is, and how often it prints, writes and plots."""

    max_step: int
    print_freq: int = 10
    write_freq: int = 10
    plot_freq: int = 10               # matplotlib, for a run outside the app
    plot_pause_seconds: float = 0.01


PART_CLASSES = {"material": Material, "forces": Forces, "timeint": TimeInt,
                "topology": Topology, "loading": Loading, "run": Run}


def _numbers(spec: "Spec") -> list[tuple[str, float]]:
    """Every number of the spec, named, for the finite check."""
    found = []
    for part in PARTS:
        for key, value in asdict(getattr(spec, part)).items():
            items = value if isinstance(value, (tuple, list)) else [value]
            found += [(f"{part}.{key}", v) for v in items
                      if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return found


@dataclass(frozen=True)
class Spec:
    """Every selection of an example; each part is one step's."""

    config: Config
    material: Material
    forces: Forces
    timeint: TimeInt
    topology: Topology
    loading: Loading
    run: Run

    def with_part(self, part: str, value: Any) -> "Spec":
        """The spec with one step's part replaced."""
        return replace(self, **{part: value})

    def changed_parts(self, other: "Spec") -> set[str]:
        """The parts (step keys) whose values differ from ``other``'s."""
        return {part for part in PARTS if getattr(self, part) != getattr(other, part)}

    def validate(self) -> list[str]:  # pylint: disable=too-many-branches
        """The problems with the spec; an empty list means it can run."""
        bad = [name for name, value in _numbers(self) if not math.isfinite(value)]
        if bad:
            return [f"{name}: must be a finite number" for name in bad]
        errors = []
        if self.config.box_length <= 0:
            errors.append("box length: must be positive")
        errors += self.config.problems()
        m = self.material
        for label, value in (("burgmag", m.burgmag), ("mu", m.mu), ("a", m.a),
                             ("rann", m.rann)):
            if value <= 0:
                errors.append(f"{label}: must be positive")
        if not -1.0 < m.nu < 0.5:
            errors.append("nu: must be in (-1, 0.5)")
        if not 0.0 < m.minseg_ratio < m.maxseg_ratio < 1.0:
            errors.append("segment lengths: need 0 < minseg < maxseg < box length")
        elif m.minseg_ratio >= 0.5 * m.maxseg_ratio:
            errors.append("minseg: must be less than half of maxseg, or remeshing undoes itself")
        f = self.forces
        if f.force_mode not in FORCE_MODES:
            errors.append(f"force mode: must be one of {', '.join(FORCE_MODES)}")
        if f.mobility_law not in MOBILITY_LAWS:
            errors.append(f"mobility law: must be one of {', '.join(MOBILITY_LAWS)}")
        for label, value in (("Ec", f.Ec), ("mob", f.mob)):
            if value is not None and value <= 0:
                errors.append(f"{label}: must be positive")
        if self.timeint.integrator not in INTEGRATORS:
            errors.append(f"integrator: must be one of {', '.join(INTEGRATORS)}")
        if self.timeint.dt <= 0:
            errors.append("dt: must be positive")
        t = self.topology
        for label, value, choices in (("split mode", t.split_mode, SPLIT_MODES),
                                      ("collision mode", t.collision_mode, COLLISION_MODES),
                                      ("remesh rule", t.remesh_rule, REMESH_RULES)):
            if value not in choices:
                errors.append(f"{label}: must be one of {', '.join(choices)}")
        if not all(1 <= n <= MAX_N_DIV for n in t.n_div):
            errors.append(f"cell list divisions: each must be in [1, {MAX_N_DIV}]")
        if self.loading.loading_mode not in LOADING_MODES:
            errors.append(f"loading mode: must be one of {', '.join(LOADING_MODES)}")
        r = self.run
        if not 1 <= r.max_step <= MAX_STEP:
            errors.append(f"max_step: must be in [1, {MAX_STEP}]")
        for label, value in (("print_freq", r.print_freq), ("write_freq", r.write_freq),
                             ("plot_freq", r.plot_freq)):
            if value < 1:
                errors.append(f"{label}: must be at least 1")
        if r.plot_pause_seconds < 0:
            errors.append("plot_pause_seconds: must not be negative")
        return errors

    def to_json(self) -> str:
        """The spec as JSON, for the run folder."""
        data = {part: asdict(getattr(self, part)) for part in PARTS}
        data["config"] = {"kind": self.config.kind, **data["config"]}
        return json.dumps(data, indent=1)

    @classmethod
    def from_json(cls, text: str) -> "Spec":
        """The spec from ``to_json`` text; raises ValueError, KeyError or TypeError if bad."""
        return cls.from_tables(json.loads(text))

    @classmethod
    def from_tables(cls, data: Mapping[str, Any]) -> "Spec":
        """The spec of an input file's tables, or of ``to_json``; raises ``InputError``."""
        problems: list[str] = []
        config_table = dict(data.get("config", {}))
        kind = config_table.pop("kind", None)
        parts: dict[str, Any] = {}
        if kind not in CONFIGS:
            problems.append(f"config.kind: must be one of {', '.join(CONFIGS)}")
        else:
            parts["config"] = _build(CONFIGS[kind], config_table, "config", problems)
        for part, part_class in PART_CLASSES.items():
            parts[part] = _build(part_class, data.get(part, {}), part, problems)
        if problems:
            raise InputError(problems)
        return cls(**parts)


def _convert(value: Any, kind: Any, where: str) -> Any:
    """A TOML or JSON value as the field's type; raises ValueError with a message."""
    if kind is bool:
        if not isinstance(value, bool):
            raise ValueError(f"{where}: must be true or false")
        return value
    if kind in (int, float):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{where}: must be a number")
        if kind is int and not float(value).is_integer():
            raise ValueError(f"{where}: must be a whole number")
        return kind(value)
    if kind is str:
        if not isinstance(value, str):
            raise ValueError(f"{where}: must be text")
        return value
    if get_origin(kind) is tuple:
        args = get_args(kind)
        if not isinstance(value, (list, tuple)) or len(value) != len(args):
            raise ValueError(f"{where}: must be a list of {len(args)} numbers")
        return tuple(_convert(v, a, where) for v, a in zip(value, args))
    if get_origin(kind) in (Union, types.UnionType):  # float | None
        inner = next(a for a in get_args(kind) if a is not types.NoneType)
        return None if value is None else _convert(value, inner, where)
    raise TypeError(f"{where}: unsupported type {kind}")


def _build(cls: Any, table: Any, part: str, problems: list[str]) -> Any:
    """One part from its table; adds its problems (unknown, missing or bad keys) to the list."""
    if not isinstance(table, Mapping):
        problems.append(f"{part}: must be a table")
        return None
    known = {f.name: f for f in fields(cls)}
    for key in table:
        if key not in known:
            problems.append(f"{part}.{key}: unknown key; the keys are {', '.join(known)}")
    values, before = {}, len(problems)
    for name, f in known.items():
        if name in table:
            try:
                values[name] = _convert(table[name], f.type, f"{part}.{name}")
            except ValueError as exc:
                problems.append(str(exc))
        elif f.default is MISSING and f.default_factory is MISSING:
            problems.append(f"{part}.{name}: missing")
    return cls(**values) if len(problems) == before else None
