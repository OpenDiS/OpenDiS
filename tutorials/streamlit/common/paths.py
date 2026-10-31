"""common.paths: where the app, its examples, its runs and the OpenDiS repository are.

The app lives in ``tutorials/streamlit/`` of the OpenDiS repository, and runs PyDiS from that
repository: ``OPENDIS_DIR`` is its root, as in the tutorials' ``cd ${OPENDIS_DIR}``. A clone
that lacks PyDiS's graph submodule, PyDiS's compiled library or ExaDiS's Python module gets
them here, as Streamlit Community Cloud's does.
"""
import collections
import concurrent.futures
import os
import re
import resource
import shutil
import subprocess
import sys
import sysconfig
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
EXADIS = "core/exadis"  # a git submodule, with Kokkos and pybind11 as submodules of its own
EXADIS_MODULE = "pyexadis" + sysconfig.get_config_var("EXT_SUFFIX")  # in core/exadis/python
# ExaDiS's build: OpenDiS's Linux preset; no FFTW, which the Cloud lacks; and no link-time
# optimization, whose final link would compile all the files at once, in parallel
EXADIS_OPTIONS = ("-DSYS=ubuntu", "-DEXADIS_FFT=Off", "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF")
BUILD_THREAD = "opendis-build"
_PERCENT = re.compile(r"\[\s*(\d+)%\]")  # CMake's progress, as in "[ 45%] Building CXX object"
_OBJECT = re.compile(r"object \S*?([^/\s]+)\.o\s*$")  # ... and the file it compiles
_GRAPH_LOCK = threading.Lock()
_BUILD = {"started": False}  # the builds: once per app process
_BUILD_LOCK = threading.Lock()


def add_pydis_paths() -> None:
    """Put PyDiS on ``sys.path``, for the app's own use of it (the initial configuration)."""
    for path in PYDIS_PATHS:
        full = str(OPENDIS_DIR / path)
        if full not in sys.path:
            sys.path.append(full)


def _fetch_submodule(path: str, marker: str, *options: str) -> None:
    """Fetch a git submodule if the clone lacks its ``marker`` file, as Streamlit Community
    Cloud's clone lacks them all; a failure shows later, as the missing file."""
    if (OPENDIS_DIR / path / marker).is_file():
        return
    try:
        subprocess.run(["git", "submodule", "update", "--init", *options, path], cwd=OPENDIS_DIR,
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


def build_exadis(repo: Path = OPENDIS_DIR, target: Path | None = None,  # pylint: disable=too-many-locals
                 options: tuple[str, ...] = EXADIS_OPTIONS, jobs: int = 1,
                 progress: dict | None = None) -> str:
    """Build ExaDiS's Python module into ``target`` (``core/exadis/python`` of the repository)
    with OpenDiS's own CMake project, as its install instructions do; returns what it built.

    CMake comes from pip where the system has none (the root ``requirements.txt``). The build
    compiles ``jobs`` files at a time: ExaDiS's largest takes about 1 GiB, and the Cloud gives
    the app 3 GiB in all. ``progress`` gets the build's phase, CMake's percent and the file it
    compiles. The module appears whole.
    """
    progress = {} if progress is None else progress
    target = target or repo / EXADIS / "python"
    cmake = shutil.which("cmake") or shutil.which("cmake", path=str(Path(sys.executable).parent))
    if cmake is None:
        raise FileNotFoundError("no cmake")
    target.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="opendis_exadis_build_") as tmp:
        build, out = Path(tmp) / "build", Path(tmp) / "out"
        progress.update(phase="configuring with CMake")
        subprocess.run([cmake, "-S", str(repo), "-B", str(build), *options,
                        f"-DPYTHON_EXECUTABLE={sys.executable}",
                        f"-DPython3_EXECUTABLE={sys.executable}", f"-DPYEXADIS_OUTPUT_DIR={out}"],
                       check=True, capture_output=True, text=True, timeout=600)
        progress.update(phase="compiling", percent=0)
        command = [cmake, "--build", str(build), "--target", "pyexadis", "-j", str(jobs)]
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, errors="replace") as proc:
            tail: collections.deque = collections.deque(maxlen=40)
            for line in proc.stdout:
                tail.append(line)
                if match := _PERCENT.match(line):
                    progress["percent"] = int(match.group(1))
                if match := _OBJECT.search(line):
                    progress["file"] = match.group(1)
        if proc.returncode:
            raise subprocess.CalledProcessError(proc.returncode, command, output="".join(tail))
        shutil.copy2(out / EXADIS_MODULE, target / f".{EXADIS_MODULE}.new")
        os.replace(target / f".{EXADIS_MODULE}.new", target / EXADIS_MODULE)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss  # KiB, bytes on macOS
    peak /= 2**20 if sys.platform == "darwin" else 2**10
    return f"built {target / EXADIS_MODULE} (largest process: {peak:.0f} MiB)"


def _build_pydis(_progress: dict) -> str:
    _fetch_submodule(CTYPESGEN, "run.py")
    return build_pydis_library()


def _build_exadis(progress: dict) -> str:
    progress.update(phase="fetching its sources")
    _fetch_submodule(EXADIS, "kokkos/CMakeLists.txt", "--depth", "1", "--recursive")
    return build_exadis(progress=progress)


def _why(exc: Exception) -> str:
    """Why a build failed, in a few words: the command that failed, or the error."""
    if isinstance(exc, subprocess.CalledProcessError):
        return f"{Path(exc.cmd[0]).name} failed with exit status {exc.returncode}"
    return str(exc)


def _build(builds: list) -> None:
    """Run the builds in turn; each one's outcome goes to its progress and the app's log."""
    for what, build, progress in builds:
        start = time.monotonic()
        progress.update(phase="starting", start=start)
        print(f"[OpenDiS app] building {what}", file=sys.stderr, flush=True)
        try:
            note = build(progress)
            print(f"[OpenDiS app] {note} in {time.monotonic() - start:.1f} s", file=sys.stderr,
                  flush=True)
        except (OSError, subprocess.SubprocessError) as exc:
            detail = (getattr(exc, "stderr", None) or getattr(exc, "stdout", None) or "")[-2000:]
            progress["error"] = _why(exc)
            print(f"[OpenDiS app] {what} was not built: {exc}\n{detail}", file=sys.stderr,
                  flush=True)


def build_libraries_once() -> None:
    """Build what the clone lacks, once, in the background: PyDiS's C library, then, on Linux,
    ExaDiS's Python module, for the app's Python.

    Until PyDiS's library is there, and if its build fails, PyDiS uses its Python versions,
    and its output starts with "pydis_lib not found". Until ExaDiS's module is there, ``import
    pyexadis`` fails. The outcomes go to the app's log (on Streamlit Community Cloud, under
    Manage app). A build that runs on after Streamlit reloaded this module is not started again.
    """
    with _BUILD_LOCK:
        if _BUILD["started"] or any(t.name == BUILD_THREAD for t in threading.enumerate()):
            return
        _BUILD["started"] = True
    builds, exadis = [], None
    if not all((OPENDIS_DIR / "lib" / name).is_file() for name in LIBRARY):
        builds.append(("PyDiS's C library", _build_pydis, {}))
    if (sys.platform.startswith("linux")
            and not (OPENDIS_DIR / EXADIS / "python" / EXADIS_MODULE).is_file()):
        exadis = {"phase": "waiting for PyDiS's library" if builds else "starting"}
        builds.append(("ExaDiS's Python module", _build_exadis, exadis))
    _BUILD["exadis"] = exadis
    if builds:
        thread = threading.Thread(target=_build, args=(builds,), name=BUILD_THREAD, daemon=True)
        thread.exadis = exadis  # for exadis_status(), also after Streamlit reloads this module
        thread.start()


def exadis_status() -> dict:
    """ExaDiS's module for the app's Python: ``ready``; ``building``, with the build's phase,
    CMake's percent, the file it compiles and the time so far; ``failed``, with why; or
    ``missing``, where the app does not build it (on macOS)."""
    module = OPENDIS_DIR / EXADIS / "python" / EXADIS_MODULE
    if module.is_file():
        return {"state": "ready", "module": module.name}
    thread = next((t for t in threading.enumerate() if t.name == BUILD_THREAD), None)
    progress = getattr(thread, "exadis", None) or _BUILD.get("exadis") or {}
    if progress.get("error"):
        return {"state": "failed", "error": progress["error"]}
    if thread is not None and progress:
        return {"state": "building", "percent": 0, **progress,
                "elapsed_s": time.monotonic() - progress.get("start", time.monotonic())}
    return {"state": "missing"}
