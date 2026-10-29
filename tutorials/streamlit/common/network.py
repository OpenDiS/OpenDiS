"""common.network: a dislocation network as arrays, from PyDiS's JSON or from PyDiS itself.

The views draw what this module gives: the node positions, the nodes' constraints, the segments
as pairs of node indices, and the cell. ``read_json`` reads what ``DisNetManager.write_json``
writes (a job's snapshots and its final network), without PyDiS. ``initial`` asks PyDiS, in the
app's process, for the network that the script's first lines build, from the same code.
The positions are mapped as PyDiS's own plot maps them: each node to its image nearest the
cell's centre, and each segment's far end to its image nearest the near end.
"""
import contextlib
import io
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from common import codegen
from common.paths import add_pydis_paths
from common.spec import Spec

PINNED = 7  # DisNode.Constraints.PINNED_NODE
SNAPSHOT = re.compile(r"disnet_(\d+)\.json")  # SimulateNetwork's write_freq files


@dataclass(frozen=True)
class Network:
    """A network's nodes, segments and cell, as arrays."""

    positions: np.ndarray    # (N, 3)
    constraints: np.ndarray  # (N,)
    segments: np.ndarray     # (M, 2), node indices
    h: np.ndarray            # (3, 3), the cell's edges as columns
    origin: np.ndarray       # (3,), the cell's lowest corner
    periodic: tuple[bool, bool, bool]

    @property
    def center(self) -> np.ndarray:
        """The cell's centre, as ``Cell.center``."""
        return self.origin + 0.5 * self.h @ np.ones(3)

    def closest_image(self, ref: np.ndarray, r: np.ndarray) -> np.ndarray:
        """``r`` mapped to its images nearest ``ref``, as ``Cell.closest_image`` maps it."""
        if not any(self.periodic):
            return np.array(r, dtype=float)
        off_diagonal = self.h - np.diag(np.diag(self.h))
        if not off_diagonal.any():  # an orthogonal cell: whole box edges, one axis at a time
            out = np.array(r, dtype=float)
            for k in range(3):
                if self.periodic[k]:
                    edge = self.h[k, k]
                    out[..., k] -= np.rint((out[..., k] - ref[..., k]) / edge) * edge
            return out
        ds = (np.linalg.inv(self.h) @ (np.asarray(r) - ref).T).T
        ds -= np.array(self.periodic) * np.round(ds)
        return (self.h @ ds.T).T + ref

    def nodes(self) -> np.ndarray:
        """Every node's image nearest the cell's centre."""
        return self.closest_image(self.center, self.positions)

    def segment_ends(self) -> tuple[np.ndarray, np.ndarray]:
        """Each segment's two ends: the first in the cell, the second its nearest image."""
        nodes = self.nodes()
        r1 = nodes[self.segments[:, 0]]
        return r1, self.closest_image(r1, nodes[self.segments[:, 1]])

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        """The lowest and highest corner of the cell's bounding box."""
        corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float)
        points = self.origin + corners @ self.h.T
        return points.min(axis=0), points.max(axis=0)


def from_data(data: dict) -> Network:
    """The network of PyDiS's export data, or of its JSON (lists instead of arrays)."""
    nodes, segs, cell = data["nodes"], data["segs"], data["cell"]
    segments = np.asarray(segs["nodeids"], dtype=int).reshape(-1, 2)
    return Network(positions=np.asarray(nodes["positions"], dtype=float).reshape(-1, 3),
                   constraints=np.asarray(nodes["constraints"], dtype=int).reshape(-1),
                   segments=segments,
                   h=np.asarray(cell["h"], dtype=float),
                   origin=np.asarray(cell.get("origin", np.zeros(3)), dtype=float),
                   periodic=tuple(bool(p) for p in cell["is_periodic"]))


def read_json(path: Path) -> Network:
    """The network of a file that ``DisNetManager.write_json`` wrote."""
    return from_data(json.loads(path.read_text(encoding="utf-8")))


def initial(spec: Spec) -> Network:
    """The initial configuration, built by PyDiS in this process from the script's own code.

    What the code prints (the tutorials' functions print a line) is not shown.
    """
    add_pydis_paths()
    namespace: dict = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(codegen.init_code(spec), "<initial configuration>", "exec"),  # pylint: disable=exec-used
             namespace)
    return from_data(namespace["net"].export_data())


def snapshots(output: Path) -> list[tuple[int, Path]]:
    """A run's snapshots, ``output/disnet_<step>.json``, by step."""
    found = [(int(m.group(1)), path) for path in output.glob("disnet_*.json")
             if (m := SNAPSHOT.fullmatch(path.name))]
    return sorted(found)
