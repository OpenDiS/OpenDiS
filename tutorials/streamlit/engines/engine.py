"""engines/engine.py: what each entry file of this folder runs, the app of ``app.py``.

Streamlit Community Cloud deploys one app per entry file. So each engine,
``opendis-NN.streamlit.app``, has its own, ``opendis_NN.py``, which calls ``main()`` (README.md).
"""
import os
import sys
from pathlib import Path

import streamlit as st

APP_DIR = str(Path(__file__).resolve().parents[1])  # tutorials/streamlit, where app.py is
RESET = Path(__file__).with_name("RESET")  # a push that changes its line resets every engine
SEEN = "OPENDIS_ENGINE_RESET"  # the line of RESET that this process saw last, in its environment


def first(folder: str, paths: str) -> str:
    """A search path, folders joined by ``os.pathsep``, with ``folder`` first and only once."""
    return os.pathsep.join([folder] + [p for p in paths.split(os.pathsep) if p and p != folder])


def reset_due() -> bool:
    """Whether RESET holds a line that this process has not seen since its first page run."""
    line = RESET.read_text(encoding="utf-8").strip() if RESET.is_file() else ""
    seen = os.environ.get(SEEN)
    os.environ[SEEN] = line
    return seen is not None and seen != line


def main() -> None:
    """Run the app's page as ``app.py`` does, after a reset if one is due.

    The app's folder goes first on sys.path, for the app's imports, and on PYTHONPATH:
    Streamlit reloads a changed module only if it lies in the entry file's folder or on
    PYTHONPATH. So a push updates an engine in place, as it updates app.py.
    """
    if APP_DIR not in sys.path:
        sys.path.insert(0, APP_DIR)
    os.environ["PYTHONPATH"] = first(APP_DIR, os.environ.get("PYTHONPATH", ""))
    if reset_due():
        from common import job  # pylint: disable=import-outside-toplevel
        deleted, kept = job.clear_runs()
        st.cache_data.clear()
        st.cache_resource.clear()
        print(f"[OpenDiS engine] reset: deleted the results of {deleted} run(s), kept {kept}",
              file=sys.stderr, flush=True)
    import app  # pylint: disable=import-outside-toplevel
    app.main()
