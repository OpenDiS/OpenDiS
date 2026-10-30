"""common.paths: where the app, its examples, its runs and the OpenDiS repository are.

The app lives in ``tutorials/streamlit/`` of the OpenDiS repository, and runs PyDiS from that
repository: ``OPENDIS_DIR`` is its root, as in the tutorials' ``cd ${OPENDIS_DIR}``. A clone
that lacks PyDiS's graph submodule, or its compiled library, gets them here, as Streamlit
Community Cloud's does.
"""
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
EXAMPLES = APP / "examples"
ASSETS = APP / "assets"
RUNS = Path(os.environ.get("OPENDIS_APP_RUNS", APP / "runs"))
OPENDIS_DIR = Path(os.environ.get("OPENDIS_DIR", APP.parents[1])).resolve()
PYDIS_PATHS = ("python", "lib", "core/pydis/python")  # as in the tutorials' scripts
GRAPH = "core/pydis/python/pydis/graph"  # a git submodule, which PyDiS imports
CTYPESGEN = "core/pydis/python/ctypesgen"  # a git submodule, which writes the library's wrapper
LIBRARY = ("libpydis.so", "pydis_lib.py")  # PyDiS's C library in lib/, and its ctypes wrapper
C_FLAGS = ("-O3", "-fPIC", "-ffp-contract=off")  # contraction off, as in the _repro builds
_GRAPH_LOCK = threading.Lock()
_BUILD = {"started": False}  # the library's build: once per app process
_BUILD_LOCK = threading.Lock()


def add_pydis_paths() -> None:
    """Put PyDiS on ``sys.path``, for the app's own use of it (the initial configuration)."""
    for path in PYDIS_PATHS:
        full = str(OPENDIS_DIR / path)
        if full not in sys.path:
            sys.path.append(full)


def _fetch_submodule(path: str, marker: str) -> None:
    """Fetch a git submodule if the clone lacks its ``marker`` file, as Streamlit Community
    Cloud's clone lacks them all; a failure shows later, as the missing file."""
    if (OPENDIS_DIR / path / marker).is_file():
        return
    try:
        subprocess.run(["git", "submodule", "update", "--init", path], cwd=OPENDIS_DIR,
                       check=False, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        pass


def fetch_graph() -> None:
    """Fetch PyDiS's graph submodule if the clone lacks it, as Streamlit Community Cloud's does."""
    with _GRAPH_LOCK:
        _fetch_submodule(GRAPH, "graph.py")  # else PyDiS's import error names the module


def _cmake_set(path: Path, name: str) -> list[str]:
    """The words of ``set(NAME ...)`` in a CMakeLists.txt."""
    match = re.search(rf"set\(\s*{name}\s(.*?)\)", path.read_text(encoding="utf-8"), re.I | re.S)
    return match.group(1).split() if match else []


def build_pydis_library(repo: Path = OPENDIS_DIR, target: Path | None = None) -> str:
    """Build PyDiS's C library into ``target`` (``lib/`` of the repository) as CMake does, from
    the same sources and headers, with a _repro build's flags; returns what it built.

    CMake is not needed: gcc compiles each source of ``core/pydis/c``'s CMakeLists.txt files,
    with the portable log/atan shim for SegSegForce.c, and ctypesgen, a git submodule, writes
    the wrapper. Each file appears whole, the library before the wrapper that loads it.
    """
    target = target or repo / "lib"
    c_dir, py_dir = repo / "core/pydis/c", repo / "core/pydis/python"
    sources = [c_dir / part / name for part in ("calforce", "collision", "remesh", "util")
               for name in _cmake_set(c_dir / part / "CMakeLists.txt", "SOURCES")]
    headers = ([c_dir / "calforce" / name
                for name in _cmake_set(py_dir / "CMakeLists.txt", "CALFORCE_HEADER_FILES")]
               + [c_dir / "include" / name
                  for name in _cmake_set(py_dir / "CMakeLists.txt", "INCLUDE_HEADER_FILES")])
    with tempfile.TemporaryDirectory(prefix="opendis_pydis_build_") as tmp:
        work = Path(tmp)

        def compile_one(source: Path) -> Path:
            shim = (["-include", str(c_dir / "calforce" / "portable_math_shim.h")]
                    if source.name == "SegSegForce.c" else [])
            obj = work / f"{source.parent.name}_{source.stem}.o"
            subprocess.run(["gcc", "-c", *C_FLAGS, *shim, "-I", str(c_dir / "include"),
                            str(source), "-o", str(obj)], cwd=source.parent, check=True,
                           capture_output=True, text=True, timeout=600)
            return obj

        with concurrent.futures.ThreadPoolExecutor(min(8, os.cpu_count() or 1)) as pool:
            objects = list(pool.map(compile_one, sources))
        subprocess.run(["gcc", "-shared", "-o", str(work / LIBRARY[0]), *map(str, objects),
                        "-lm"], check=True, capture_output=True, text=True, timeout=600)
        subprocess.run([sys.executable, str(py_dir / "ctypesgen" / "run.py"), "--allow-gnu-c",
                        "--cpp=gcc -E", *map(str, headers), "-l", LIBRARY[0], "-o", LIBRARY[1]],
                       cwd=work, check=True, capture_output=True, text=True, timeout=600)
        target.mkdir(exist_ok=True)
        for name in LIBRARY:
            shutil.copy2(work / name, target / f".{name}.new")
            os.replace(target / f".{name}.new", target / name)
    return f"built {target / LIBRARY[0]} from {len(sources)} C files"


def build_pydis_library_once() -> None:
    """If the clone lacks PyDiS's C library, build it once, in the background.

    Until it is there, and if the build fails, PyDiS uses its Python versions, and its output
    starts with "pydis_lib not found". The outcome goes to the app's log (on Streamlit
    Community Cloud, under Manage app).
    """
    with _BUILD_LOCK:
        if _BUILD["started"] or all((OPENDIS_DIR / "lib" / name).is_file() for name in LIBRARY):
            return
        _BUILD["started"] = True
    if shutil.which("gcc") is None:
        print("[OpenDiS app] no gcc, so PyDiS uses its Python versions", file=sys.stderr)
        return

    def build() -> None:
        start = time.monotonic()
        _fetch_submodule(CTYPESGEN, "run.py")
        try:
            note = build_pydis_library()
            print(f"[OpenDiS app] {note} in {time.monotonic() - start:.1f} s", file=sys.stderr,
                  flush=True)
        except (OSError, subprocess.SubprocessError) as exc:
            detail = (getattr(exc, "stderr", None) or "")[-2000:]
            print(f"[OpenDiS app] PyDiS's C library was not built: {exc}\n{detail}",
                  file=sys.stderr, flush=True)

    threading.Thread(target=build, name="pydis-library-build", daemon=True).start()
