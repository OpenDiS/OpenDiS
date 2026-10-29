"""common.paths: where the app, its examples, its runs and the OpenDiS repository are.

The app lives in ``tutorials/streamlit/`` of the OpenDiS repository, and runs PyDiS from that
repository: ``OPENDIS_DIR`` is its root, as in the tutorials' ``cd ${OPENDIS_DIR}``.
"""
import os
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
EXAMPLES = APP / "examples"
ASSETS = APP / "assets"
RUNS = Path(os.environ.get("OPENDIS_APP_RUNS", APP / "runs"))
OPENDIS_DIR = Path(os.environ.get("OPENDIS_DIR", APP.parents[1])).resolve()
PYDIS_PATHS = ("python", "lib", "core/pydis/python")  # as in the tutorials' scripts


def add_pydis_paths() -> None:
    """Put PyDiS on ``sys.path``, for the app's own use of it (the initial configuration)."""
    for path in PYDIS_PATHS:
        full = str(OPENDIS_DIR / path)
        if full not in sys.path:
            sys.path.append(full)
