"""tests.test_build: PyDiS's library and ExaDiS's module, built as the app builds them."""
import filecmp
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

from common import catalog, codegen, paths
from common.job import JOB_ENV
from common.testing import check, summary

# The binary junction's ExaDiS tutorial, as it is but for its plots, with the module built here
EXADIS_RUN = """\
import sys
sys.path[:0] = {paths!r}
import pyexadis
import pyexadis_base
pyexadis_base.VisualizeNetwork = lambda: None
import test_binary_junction_exadis as tutorial
pyexadis.initialize()
tutorial.main()
print("module:", pyexadis.__file__)
pyexadis.finalize()
"""


def short_run(example: catalog.Example, root: Path, folder: Path) -> str:
    """A 20-step run of the example's own values, with PyDiS from ``root``; returns its output."""
    spec = example.default.with_part("run", replace(example.default.run, max_step=20))
    folder.mkdir(parents=True)
    (folder / example.script_name).write_text(codegen.full_script(spec, example.final_json))
    done = subprocess.run([sys.executable, example.script_name], cwd=folder, capture_output=True,
                          text=True, timeout=600, check=False,
                          env={**os.environ, **JOB_ENV, "OPENDIS_DIR": str(root)})
    return done.stdout + done.stderr


def exadis_checks(tmp: Path) -> list:
    """ExaDiS's module, built as the app builds it (on macOS, without OpenMP, which Apple's
    clang lacks), into a temporary folder; then the binary junction's ExaDiS tutorial with it."""
    target, repo = tmp / "pyexadis", paths.OPENDIS_DIR
    options = paths.EXADIS_OPTIONS + (("-DKokkos_ENABLE_OPENMP=Off",)
                                      if sys.platform == "darwin" else ())
    try:
        note = paths.build_exadis(repo, target, options, jobs=os.cpu_count() or 1)
    except (OSError, subprocess.SubprocessError) as exc:
        note = f"{exc}\n{(getattr(exc, 'stderr', None) or '')[-1000:]}"
    script = EXADIS_RUN.format(paths=[str(target)] + [str(repo / part) for part in (
        *paths.PYDIS_PATHS, "core/exadis/python", "examples/03_binary_junction")])
    (tmp / "exadis_run").mkdir()
    done = subprocess.run([sys.executable, "-c", script], cwd=tmp / "exadis_run",
                          capture_output=True, text=True, timeout=600, check=False,
                          env={**os.environ, **JOB_ENV, "OMP_NUM_THREADS": "1"})
    output = done.stdout + done.stderr
    return [
        check("ExaDiS's module is built", (target / paths.EXADIS_MODULE).is_file(), note),
        check("the binary junction's ExaDiS tutorial runs with it",
              f"module: {target / paths.EXADIS_MODULE}" in output and "step = 200" in output,
              output[-300:]),
    ]


def main() -> int:
    """Run every check."""
    example = next(e for e in catalog.examples() if e.name == "binary_junction")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "opendis"  # the repository's Python code, without its library
        for part in ("python", "core/pydis/python"):
            shutil.copytree(paths.OPENDIS_DIR / part, root / part,
                            ignore=shutil.ignore_patterns("__pycache__"))
        note = paths.build_pydis_library(paths.OPENDIS_DIR, root / "lib")
        built = short_run(example, root, Path(tmp) / "built")
        results = [
            check("the library and its wrapper are built",
                  all((root / "lib" / name).is_file() for name in paths.LIBRARY), note),
            check("a run uses them", "pydis_lib not found" not in built
                  and "step = 10" in built, built[-300:]),
        ]
        if all((paths.OPENDIS_DIR / "lib" / name).is_file() for name in paths.LIBRARY):
            own = short_run(example, paths.OPENDIS_DIR, Path(tmp) / "own")
            outputs = sorted(p.name for p in (Path(tmp) / "own" / "output").glob("*.json"))
            same = [name for name in outputs if filecmp.cmp(
                Path(tmp) / "built" / "output" / name, Path(tmp) / "own" / "output" / name,
                shallow=False)]
            results.append(check("its run equals the run with the repository's own library, "
                                 "byte for byte (if that is a _repro build)",
                                 outputs and same == outputs,
                                 f"{len(same)} of {len(outputs)} files, {own[-200:]}"))
        results += exadis_checks(Path(tmp))
    return summary(results)


if __name__ == "__main__":
    sys.exit(main())
