"""tests.test_codegen: the app's scripts against the tutorials' own scripts (README.md, Tests)."""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
from dataclasses import replace
from pathlib import Path

from common import catalog, codegen, job
from common.paths import OPENDIS_DIR, PYDIS_PATHS
from common.spec import STEPS, Spec
from common.testing import check, has, summary

SHORT_STEPS = 60


def outputs(folder: Path) -> dict[str, object]:
    """The networks that a run wrote into ``output/``, by name: the JSON files, and ExaDiS's
    ``config.<step>.data`` files, as text."""
    found: dict[str, object] = {path.name: json.loads(path.read_text(encoding="utf-8"))
                                for path in sorted((folder / "output").glob("*.json"))}
    found.update((path.name, path.read_text(encoding="utf-8"))
                 for path in sorted((folder / "output").glob("*.data")))
    return found


def step_lines(lines: list[str]) -> list[str]:
    """The ``step = ...`` lines, without their wall time."""
    return [re.sub(r", elapsed = .*$", "", line) for line in lines if line.startswith("step")]


def main() -> int:  # pylint: disable=too-many-locals
    """Run every check; ``--full`` runs the tutorials' full length."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="run the full tutorials")
    full = parser.parse_args().full
    examples, problems = catalog.scan()
    values = [0.0, 1.0, 0.3, 0.125, 0.0316, 1000.0, 3e-10, 50e9, 160e9, 1e6, 1e-8, 1e-9, -4e8,
              1e5, 1e-3, 99999.9, 1.234567890123e-7]
    inexact = [v for v in values if float(codegen.num(v)) != v]
    looks = [codegen.num(v) for v in (50e9, 3e-10, 1e-8, -4e8, 1e6)]
    results = [
        check("the input files have no problems", not problems, problems),
        check("three examples, in order", [e.name for e in examples] ==
              ["frank_read_src", "frank_read_src_exadis", "binary_junction"],
              [e.name for e in examples]),
        check("num: every number comes back exactly", not inexact, inexact),
        check("num: 5.0e10, 3.0e-10, 1.0e-8, -4.0e8, 1.0e6",
              looks == ["5.0e10", "3.0e-10", "1.0e-8", "-4.0e8", "1.0e6"], looks),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        for example in examples:
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
            if count != len(paths):
                raise AssertionError(f"the path line of {example.tutorial_script} has changed")
            if not full:
                text, count = re.subn(r"max_step=\d+", f"max_step={SHORT_STEPS}", text)
                if count != 1:
                    raise AssertionError(f"max_step moved in {example.tutorial_script}")
            tutorial = Path(tmp) / "tutorial" / name
            tutorial.mkdir(parents=True)
            (tutorial / "tutorial.py").write_text(text, encoding="utf-8")
            ran = subprocess.run([sys.executable, "-u", "tutorial.py"], cwd=tutorial,
                                 env={**os.environ, **job.JOB_ENV}, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, check=False)
            lines = ran.stdout.splitlines()
            ours, theirs = outputs(folder.path), outputs(tutorial)
            differ = sorted(n for n in theirs if ours.get(n) != theirs[n])
            played = [n for n in ours if n.startswith("disnet_")]  # what the app's views read
            written = [n for n in theirs if n.startswith(("disnet_", "config."))]
            warnings = [line for line in status.lines if "Warning" in line]
            title = f"{name} ({spec.run.max_step} steps)"
            results += [
                check(f"{title}: the app's job ends done", status.ok, status.summary()),
                check(f"{title}: the tutorial's script exits with 0", ran.returncode == 0,
                      lines[-5:]),
                check(f"{title}: the job's output has no warnings", not warnings, warnings[:2]),
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
