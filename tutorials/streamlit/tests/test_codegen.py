"""tests.test_codegen: the app's scripts against the tutorials' own scripts (README.md, Tests)."""
import argparse
import filecmp
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
import types
import urllib.request
from dataclasses import replace
from pathlib import Path

import numpy as np

from common import catalog, codegen, job, network
from common.paths import OPENDIS_DIR, PYDIS_PATHS
from common.spec import BUILT, MAX_PRECIPITATES, STEPS, InfiniteLines, Spec
from common.testing import check, has, summary

SHORT_STEPS = 60
# The option of tutorials/01_initial_configurations/initial_configurations.py for each network
# and code: its options 3 to 5 are ExaDiS's only, so PyDiS's network is checked against them
OPTIONS = {("simple_source", "pydis"): "1", ("simple_source", "exadis"): "2",
           "paradis_file": "3", "line_config": "4", "prismatic_config": "5"}
# The tutorial, run with one option, and its plot replaced by a write of the network
TUTORIAL_RUN = """\
import os, runpy, sys
sys.argv = ["initial_configurations.py", {option!r}]
sys.path[:0] = {paths!r}


def capture(self, N, *args, **kwargs):
    os.makedirs("output", exist_ok=True)
    N.write_json("output/tutorial.json")


import pydis
pydis.VisualizeNetwork.plot_disnet = capture
if {exadis!r}:
    import pyexadis_base
    pyexadis_base.VisualizeNetwork.plot_disnet = capture
runpy.run_path("tutorial.py", run_name="__main__")
"""


def outputs(folder: Path) -> dict[str, object]:
    """The networks that a run wrote into ``output/``, by name: the JSON files, and ExaDiS's
    ``config.<step>.data`` files, as text."""
    found: dict[str, object] = {path.name: json.loads(path.read_text(encoding="utf-8"))
                                for path in sorted((folder / "output").glob("*.json"))}
    found.update((path.name, path.read_text(encoding="utf-8"))
                 for path in sorted((folder / "output").glob("*.data")))
    return found


def canon(data: dict) -> dict:
    """A network's JSON, as arrays of one shape whatever the code: PyDiS writes the cell's
    periodic flags as booleans, ExaDiS as numbers."""
    cell, nodes, segs = data["cell"], data["nodes"], data["segs"]
    found = {"h": cell["h"], "origin": cell["origin"],
             "is_periodic": [bool(p) for p in cell["is_periodic"]]}
    for group, size in ((nodes, len(nodes["tags"])), (segs, len(segs["nodeids"]))):
        found.update((key, np.asarray(value).reshape(size, -1).tolist())
                     for key, value in group.items())
    return found


def tutorial_paths(text: str) -> str:
    """The tutorial's text with its paths, relative to its own folder, made absolute."""
    return re.sub(r"'\.\./\.\./([\w/.\-]+?)/?'",
                  lambda match: repr(str(OPENDIS_DIR / match.group(1))), text)


def build_checks(example: catalog.Example, tmp: Path) -> list:
    """A script that builds a network, for each of the tutorial's networks: the app's job, and
    the tutorial's option of that network, give the same network and the same files; and the
    3-D window's own build of it is that network too."""
    results = []
    for cls in BUILT:
        spec = example.default.with_part("config", cls())
        name, final = f"{example.name}: {cls.label}", example.final_json
        script = codegen.full_script(spec, final)
        compile(script, example.script_name, "exec")
        folder = job.RunFolder.for_script(f"{example.name}_{cls.kind}", example.script_name,
                                          final, script, runs=tmp / "app")
        folder.prepare(script, spec.to_json())
        status = job.run(folder)
        option = OPTIONS.get((cls.kind, spec.engine)) or OPTIONS[cls.kind]
        tutorial = tmp / "tutorial" / f"{example.name}_{cls.kind}"
        tutorial.mkdir(parents=True)
        text = (OPENDIS_DIR / example.tutorial_script).read_text(encoding="utf-8")
        (tutorial / "tutorial.py").write_text(tutorial_paths(text), encoding="utf-8")
        (tutorial / "run.py").write_text(TUTORIAL_RUN.format(
            option=option, exadis=option != "1",
            paths=[str(OPENDIS_DIR / p) for p in codegen.EXADIS_PATHS]), encoding="utf-8")
        ran = subprocess.run([sys.executable, "-u", "run.py"], cwd=tutorial,
                             env={**os.environ, **job.JOB_ENV}, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, check=False)
        theirs = tutorial / "output" / "tutorial.json"
        ours = (canon(json.loads(folder.results.read_text(encoding="utf-8")))
                if folder.results.exists() else None)
        drawn = network.initial(spec)
        first = codegen.BUILD_TITLES[cls.kind].format(
            code="pyexadis" if spec.engine == "exadis" else "pydis")
        results += [
            check(f"{name}: the app's job ends done", status.ok, status.summary()),
            check(f"{name}: the tutorial's option {option} exits with 0", ran.returncode == 0,
                  ran.stdout[-300:]),
            check(f"{name}: the same network as the tutorial's option {option}",
                  ours is not None and theirs.exists()
                  and ours == canon(json.loads(theirs.read_text(encoding="utf-8"))),
                  (ours is not None, theirs.exists())),
            check(f"{name}: the 3-D window's own build is that network",
                  ours is not None and np.array_equal(drawn.positions, ours["positions"])
                  and np.array_equal(drawn.segments, ours["nodeids"]), len(drawn.positions)),
            check(f"{name}: the tutorial's first line", first in status.lines,
                  status.lines[:2]),
        ]
        if cls is InfiniteLines:
            same = [n for n in ("infinite_lines.data", "infinite_lines.vtk")
                    if (folder.output / n).exists() and filecmp.cmp(
                        folder.output / n, tutorial / n, shallow=False)]
            results.append(check(f"{name}: the same .data and .vtk files as the tutorial's",
                                 len(same) == 2, same))
    return results


def field_checks(tmp: Path) -> list:
    """The user stress's field, on a small grid: UserStress gives the grid's own values at its
    points, and the reader refuses other columns, a grid of another box or in another order,
    and a file over its size limit."""
    names = {"np": np, "os": os, "itertools": itertools, "urllib": urllib}
    exec(codegen.STRESS_FIELD, names)  # pylint: disable=exec-used
    load = names["load_stress_field"]
    box, n = 100.0, np.array([4, 5, 6])
    cell = types.SimpleNamespace(h=box * np.eye(3), origin=-0.5 * box * np.ones(3))

    def grid(order: str) -> np.ndarray:  # the rows of a grid, x slowest (C) or fastest (F)
        ijk = np.indices(n).reshape(3, -1, order=order).T
        return np.hstack([cell.origin + (ijk + 0.5) * box / n, np.arange(6.0 * len(ijk))
                          .reshape(-1, 6)])

    def refusal(call) -> str:
        try:
            call()
        except ValueError as exc:
            return str(exc)
        return "no error"

    data, columns = grid("C"), ["x", "y", "z", "Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"]
    np.savez(tmp / "field.npz", data=data, columns=columns, ngrid=n)
    np.savez(tmp / "columns.npz", data=data, columns=columns[::-1], ngrid=n)
    np.savez(tmp / "fastest.npz", data=grid("F"), columns=columns, ngrid=n)
    (tmp / "large.bin").write_bytes(bytes(2_000_000))
    names["ngrid"], names["stress_grid"] = load(str(tmp / "field.npz"), "unused", box)
    too_large = refusal(lambda: load(str(tmp / "copy.npz"), (tmp / "large.bin").as_uri(), box, 1))
    left = sorted(p.name for p in tmp.glob("copy.npz*"))
    return [
        check("user stress field: the grid's own values at its points", np.allclose(
            names["UserStress"](data[[0, 7, 119], :3], cell, None), data[[0, 7, 119], 3:])),
        check("user stress field: refuses other columns", "columns" in refusal(
            lambda: load(str(tmp / "columns.npz"), "unused", box))),
        check("user stress field: refuses a grid of another box", "not on the grid" in refusal(
            lambda: load(str(tmp / "field.npz"), "unused", 2.0 * box))),
        check("user stress field: refuses x fastest", "not on the grid" in refusal(
            lambda: load(str(tmp / "fastest.npz"), "unused", box))),
        check("user stress field: refuses a file over its limit, and leaves no part of it",
              "larger than 1 MB" in too_large and not left, (too_large, left)),
    ]


def step_lines(lines: list[str]) -> list[str]:
    """The ``step = ...`` lines, without their wall time."""
    return [re.sub(r", elapsed = .*$", "", line) for line in lines if line.startswith("step")]


def main() -> int:  # pylint: disable=too-many-locals
    """Run every check; ``--full`` runs the examples' full length."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="run the examples at full length")
    full = parser.parse_args().full
    examples, problems = catalog.scan()
    values = [0.0, 1.0, 0.3, 0.125, 0.0316, 1000.0, 3e-10, 50e9, 160e9, 1e6, 1e-8, 1e-9, -4e8,
              1e5, 1e-3, 99999.9, 1.234567890123e-7]
    inexact = [v for v in values if float(codegen.num(v)) != v]
    looks = [codegen.num(v) for v in (50e9, 3e-10, 1e-8, -4e8, 1e6)]
    results = [
        check("the input files have no problems", not problems, problems),
        check("fourteen examples, in order", [e.name for e in examples] ==
              ["frank_read_src", "frank_read_src_exadis", "binary_junction",
               "binary_junction_exadis", "binary_junction_fcc", "binary_junction_fcc_exadis",
               "uniaxial_tension", "uniaxial_tension_exadis", "user_stress",
               "user_stress_exadis", "precipitates", "precipitates_exadis",
               "initial_configurations", "initial_configurations_exadis"],
              [e.name for e in examples]),
        check("num: every number comes back exactly", not inexact, inexact),
        check("num: 5.0e10, 3.0e-10, 1.0e-8, -4.0e8, 1.0e6",
              looks == ["5.0e10", "3.0e-10", "1.0e-8", "-4.0e8", "1.0e6"], looks),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        results += field_checks(Path(tmp))
        for example in examples:
            if not example.default.config.simulates:
                results += build_checks(example, Path(tmp))
                continue
            spec, final, name = example.default, example.final_json, example.name
            script = codegen.full_script(spec, final)
            compile(script, example.script_name, "exec")
            missing = [step for step, _ in STEPS
                       if textwrap.indent(codegen.STEP_CODE[step](spec), "    ") not in script
                       or codegen.STEP_CODE[step](spec) not in codegen.step_code(spec, step, final)]
            function = codegen.init_function(spec)
            results += [
                check(f"{name}: the script compiles and holds each step's block", not missing,
                      missing),
                check(f"{name}: the script and step 1 hold the tutorial's function",
                      function in script and function in codegen.step_code(spec, "config", final)),
                check(f"{name}: the spec's JSON round trip", Spec.from_json(spec.to_json()) == spec),
            ]
            refused = {  # what: (part, its new value, a word of the message)
                "a zero time step": ("timeint", replace(spec.timeint, dt=0.0), "dt"),
                "minseg over half of maxseg": ("material", replace(
                    spec.material, minseg_ratio=0.75 * spec.material.maxseg_ratio), "minseg"),
                "an unknown collision mode": ("topology", replace(
                    spec.topology, collision_mode="Nope"), "collision"),
                "a negative box": ("config", replace(spec.config, box_length=-1.0), "box"),
                "a NaN stress": ("loading", replace(
                    spec.loading, applied_stress=(float("nan"),) + (0.0,) * 5), "finite"),
            }
            if spec.config.kind == "frank_read_src_loop":
                refused["b along the arms"] = ("config", replace(
                    spec.config, burg_vec=(0.0, 1.0, 0.0)), "parallel")
            if spec.config.kind == "precipitates":
                refused.update({
                    "radii, the largest first": ("config", replace(
                        spec.config, radius_min=3000.0), "radii"),
                    "too many precipitates": ("config", replace(
                        spec.config, num_precipitates=MAX_PRECIPITATES + 1), "precipitates"),
                    "a cross-slip mode off the menu": ("topology", replace(
                        spec.topology, cross_slip="Other"), "cross-slip")})
            if spec.forces.mobility_law == "FCC_0":
                refused["FCC_0 in bcc"] = ("material", replace(spec.material, crystal="bcc"),
                                           "fcc")
            if spec.forces.force_mode == "CUTOFF_MODEL":
                refused["CUTOFF_MODEL without a cutoff"] = ("forces", replace(
                    spec.forces, cutoff=None), "cutoff")
            if spec.config.kind == "user_stress_line":
                refused.update({
                    "a mixed line, which does not close": ("config", replace(
                        spec.config, theta=45.0), "theta"),
                    "a box other than the field's": ("config", replace(
                        spec.config, box_length=1000.0), "field"),
                    "a shear stress off the menu": ("loading", replace(
                        spec.loading, tau=2e6), "tau"),
                    "a field kind off the menu": ("forces", replace(
                        spec.forces, stress_kind="other"), "kind")})
                link = spec.forces.stress_url
                for what, url, word in (
                        ("a link off Dropbox", "https://example.com/field.npz", "https link"),
                        ("an http link", link.replace("https:", "http:"), "https link"),
                        ("a quote in the link", link + "'", "quote"),
                        ("Dropbox's page, dl=0", link.replace("dl=1", "dl=0"), "dl=1")):
                    refused[what] = ("forces", replace(spec.forces, stress_url=url), word)
                other = "https://dl.dropboxusercontent.com/s/abc/other.npz"
                text = codegen.full_script(spec.with_part("forces", replace(
                    spec.forces, stress_url=other)), final)
                results.append(check(f"{name}: another field's link goes into the script",
                                     f"stress_url = {other!r}\n" in text and link not in text))
                # user_stress.py's field, which varies in space and time: no file, 20 steps
                toy = spec.with_part("forces", replace(spec.forces, stress_kind="space_time"))
                toy = toy.with_part("run", replace(toy.run, max_step=20))
                text = codegen.full_script(toy, final)
                toy_folder = job.RunFolder.for_script(f"{name}_space_time", example.script_name,
                                                      final, text, runs=Path(tmp) / "app")
                toy_folder.prepare(text, toy.to_json())
                toy_status = job.run(toy_folder)
                results += [
                    check(f"{name}: user_stress.py's field goes into the script, with no file",
                          codegen.SPACE_TIME_STRESS in text and "stress_url" not in text
                          and not toy.validate(), toy.validate()),
                    check(f"{name}: user_stress.py's field runs 20 steps", toy_status.ok,
                          toy_status.summary()),
                ]
            for what, (part, value, word) in refused.items():
                messages = spec.with_part(part, value).validate()
                results.append(check(f"{name}: refuses {what}", has(messages, word), messages))

            # the app's job and the tutorial's own script, side by side
            if not full:
                spec = spec.with_part("run", replace(spec.run, max_step=SHORT_STEPS))
                script = codegen.full_script(spec, final)
            folder = job.RunFolder.for_script(name, example.script_name, final, script,
                                              runs=Path(tmp) / "app")
            folder.prepare(script, spec.to_json())
            status = job.run(folder)
            text = (OPENDIS_DIR / example.tutorial_script).read_text(encoding="utf-8")
            paths = codegen.EXADIS_PATHS if example.engine == "exadis" else PYDIS_PATHS
            text, count = re.subn(r"'\.\./\.\./([\w/]+?)/?'",
                                  lambda match: repr(str(OPENDIS_DIR / match.group(1))), text)
            if count not in (1, len(paths)):  # the precipitates' tutorial has one path
                raise AssertionError(f"the path line of {example.tutorial_script} has changed")
            text = re.sub(r"write_dir='[\w./-]+'", "write_dir='output'", text)  # its own folder
            # the app's length, which may differ from the tutorial's (the precipitates)
            text, count = re.subn(r"max_step=\d+", f"max_step={spec.run.max_step}", text)
            if count != 1:
                raise AssertionError(f"max_step moved in {example.tutorial_script}")
            tutorial = Path(tmp) / "tutorial" / name
            tutorial.mkdir(parents=True)
            (tutorial / "tutorial.py").write_text(text, encoding="utf-8")
            here = (OPENDIS_DIR / example.tutorial_script).parent
            for module in re.findall(r"^from (\w+) import", text, re.MULTILINE):
                local = here / f"{module}.py"
                if local.exists():  # the tutorial's own module, as user_precipitates.py
                    (tutorial / local.name).write_text(re.sub(
                        r"'\.\./\.\./([\w/]+?)/?'",
                        lambda match: repr(str(OPENDIS_DIR / match.group(1))),
                        local.read_text(encoding="utf-8")), encoding="utf-8")
            if spec.forces.stress_field:  # the field that the app's job read, in its folder
                field = OPENDIS_DIR / codegen.stress_file(spec.forces.stress_url)
                results.append(check(f"{name}: the app's job keeps its field where "
                                     "codegen.stress_file says", field.exists(), field))
                (tutorial / "rs0084.stress.npz").symlink_to(field)  # the tutorial's name
            ran = subprocess.run([sys.executable, "-u", "tutorial.py"], cwd=tutorial,
                                 env={**os.environ, **job.JOB_ENV}, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, check=False)
            lines = ran.stdout.splitlines()
            ours, theirs = outputs(folder.path), outputs(tutorial)
            differ = sorted(n for n in theirs if ours.get(n) != theirs[n])
            played = [n for n in ours if n.startswith("disnet_")]  # what the app's views read
            written = [n for n in theirs if n.startswith(("disnet_", "config."))]
            # ExaDiS warns of the pinned ends of the fcc junction's lines, in the tutorial too
            warnings = [line for line in status.lines if "Warning" in line and line not in lines]
            title = f"{name} ({spec.run.max_step} steps)"
            results += [
                check(f"{title}: the app's job ends done", status.ok, status.summary()),
                check(f"{title}: the tutorial's script exits with 0", ran.returncode == 0,
                      lines[-5:]),
                check(f"{title}: the job's output has no warnings but the tutorial's",
                      not warnings, warnings[:2]),
                check(f"{title}: the same {len(theirs)} networks, snapshots and final",
                      len(theirs) > 2 and not differ, f"these differ: {differ[:5]}"),
                check(f"{title}: a JSON snapshot for each of the tutorial's snapshots",
                      len(played) == len(written) > 2, (len(played), len(written))),
                check(f"{title}: the same step lines", step_lines(status.lines) == step_lines(lines)
                      and bool(step_lines(lines)), (step_lines(status.lines)[:2],
                                                    step_lines(lines)[:2])),
            ]
    results.append(check("no job is left running", not job.running()))
    return summary(results)


if __name__ == "__main__":
    sys.exit(main())
