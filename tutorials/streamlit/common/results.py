"""common.results: what a run leaves in its folder, for the views and the downloads.

The job's script writes into ``output/`` a snapshot of the network every write_freq steps
(``disnet_<step>.json``) and, at its end, the final network. The history, the number of nodes
and the dislocation density against the time step, is read from the snapshots, which every
PyDiS writes. This module does not import Streamlit.
"""
import functools
import io
import zipfile
from pathlib import Path

import numpy as np

from common import job
from common.network import read_json, snapshots

MAX_FRAMES = 120  # an animation with more snapshots plays a subset, the first and last included


def frames(folder: job.RunFolder) -> list[tuple[str, Path]]:
    """The snapshots to play, labelled (``step 10``), then the final network (``final``)."""
    found = snapshots(folder.output)
    if len(found) > MAX_FRAMES:
        keep = np.unique(np.linspace(0, len(found) - 1, MAX_FRAMES).round().astype(int))
        found = [found[i] for i in keep]
    labelled = [(f"step {step}", path) for step, path in found]
    if folder.results.exists():
        labelled.append(("final", folder.results))
    return labelled


def history(folder: job.RunFolder, burgmag: float) -> tuple[list[int], list[int], list[float]]:
    """The snapshots' time steps, numbers of nodes and dislocation densities (m⁻²): the total
    line length over the cell's volume, both in units of b, over burgmag²."""
    steps, nodes, density = [], [], []
    for step, path in snapshots(folder.output):
        try:
            count, length, volume = _measure(path, path.stat().st_mtime)
        except (OSError, ValueError):  # being written
            continue
        steps.append(step)
        nodes.append(count)
        density.append(length / volume / burgmag**2)
    return steps, nodes, density


@functools.lru_cache(maxsize=4096)
def _measure(path: Path, _mtime: float) -> tuple[int, float, float]:
    """A snapshot's number of nodes, total line length and cell volume; ``_mtime`` tells a file
    rewritten since from the one read before."""
    net = read_json(path)
    r1, r2 = net.segment_ends()
    return (len(net.positions), float(np.linalg.norm(r2 - r1, axis=1).sum()),
            abs(float(np.linalg.det(net.h))))


def stress_strain(folder: job.RunFolder) -> tuple[np.ndarray, np.ndarray] | None:
    """The run's strain and stress (Pa), from output/stress_strain_dens.dat, which SimulateNetwork
    writes in the same format with either code; None without the file, or while it is written."""
    path = folder.output / "stress_strain_dens.dat"
    try:
        table = np.atleast_2d(np.loadtxt(path, comments="#"))
    except (OSError, ValueError):
        return None
    if table.shape[0] == 0 or table.shape[1] < 3:
        return None
    return table[:, 1], table[:, 2]


def output_zip(folder: job.RunFolder) -> bytes:
    """Everything in ``output/``, zipped, with the script, the spec and the log."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in [folder.script, folder.path / job.SPEC, folder.log]:
            if path.exists():
                archive.write(path, arcname=path.name)
        for path in sorted(folder.output.glob("*")):
            archive.write(path, arcname=f"{job.OUTPUT}/{path.name}")
    return buffer.getvalue()
