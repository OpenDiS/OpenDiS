"""tests.test_views: the networks, a run's results and the figures (README.md, Tests)."""
import sys
import tempfile
from pathlib import Path

import numpy as np

from common import catalog, job, network, network3d, results, views
from common.paths import add_pydis_paths
from common.testing import check, summary


def main() -> int:  # pylint: disable=too-many-locals
    """Run every check."""
    found = {e.name: network.initial(e.default) for e in catalog.examples()}
    fr, bj = found["frank_read_src"], found["binary_junction"]
    checks = [
        check("the Frank-Read source: 5 nodes, 4 pinned, 5 segments",
              (len(fr.positions), int((fr.constraints == network.PINNED).sum()),
               len(fr.segments)) == (5, 4, 5)),
        check("the binary junction: 6 nodes, 4 pinned, 4 segments, not periodic",
              (len(bj.positions), int((bj.constraints == network.PINNED).sum()),
               len(bj.segments), any(bj.periodic)) == (6, 4, 4, False)),
        check("the source's arm is 125 long, its free node at the centre",
              np.isclose(np.ptp(fr.positions[:, 1]), 125.0)
              and np.allclose(fr.positions[fr.constraints != network.PINNED][0], 500.0)),
    ]

    # PyDiS's JSON read back without PyDiS: a closed triangle with b along z, whose first
    # segment crosses the periodic boundary at x = 10
    add_pydis_paths()
    from framework.disnet_manager import DisNetManager  # pylint: disable=import-outside-toplevel
    from pydis import Cell, DisNet, DisNode  # pylint: disable=import-outside-toplevel
    cell = Cell(h=10.0 * np.eye(3), is_periodic=[True, True, True])
    pinned, free = DisNode.Constraints.PINNED_NODE, DisNode.Constraints.UNCONSTRAINED
    rn = np.array([[9.5, 5.0, 5.0, pinned], [0.5, 5.0, 5.0, free], [0.5, 6.0, 5.0, free]])
    s = np.sqrt(0.5)
    links = np.array([[0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0],
                      [1, 2, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0],
                      [2, 0, 0.0, 0.0, 1.0, s, -s, 0.0]])
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        DisNetManager(DisNet(cell=cell, rn=rn, links=links)).write_json(str(tmp / "net.json"))
        net = network.read_json(tmp / "net.json")
        first = int(np.flatnonzero((net.segments == [0, 1]).all(axis=1) |
                                   (net.segments == [1, 0]).all(axis=1))[0])
        r1, r2 = (end[first] for end in net.segment_ends())
        rng = np.random.default_rng(1)
        points, refs = rng.uniform(-20, 30, (50, 3)), rng.uniform(0, 10, (50, 3))
        ours = net.closest_image(refs, points)
        theirs = np.array([cell.closest_image(Rref=r, R=p) for r, p in zip(refs, points)])
        checks += [
            check("the JSON reads back: 3 nodes, 1 pinned, 3 segments, periodic",
                  (len(net.positions), int((net.constraints == network.PINNED).sum()),
                   len(net.segments), all(net.periodic)) == (3, 1, 3, True)),
            check("a segment across the boundary is drawn as 1 long, not 9",
                  np.isclose(np.linalg.norm(r2 - r1), 1.0), np.linalg.norm(r2 - r1)),
            check("the periodic images equal PyDiS's Cell.closest_image",
                  np.array_equal(ours, theirs), np.abs(ours - theirs).max()),
        ]

        # the figures
        fig = views.network_figure(fr, "title")
        segments, free_nodes, pinned_nodes, box = fig.data
        eye = fig.layout.scene.camera.eye
        azimuth = np.degrees(np.arctan2(eye.y, eye.x))
        elevation = np.degrees(np.arctan2(eye.z, np.hypot(eye.x, eye.y)))
        history = views.history_figure([0, 10], [7, 9], [5e12, 6e12])
        tension = views.history_figure([0, 10], [7, 9], [5e12, 6e12], ([0.0, 1e-4], [0.0, 5e7]))
        payload = network3d.frames_payload([("step 0", fr), ("final", bj)])
        again = network3d.frames_payload([("step 0", fr), ("final", bj)])
        other = network3d.frames_payload([("step 0", fr)])
        junction = network3d.frames_payload([("initial configuration", bj)])
        checks += [
            check("both examples start in the same view: the same camera and aspect ratio",
                  all(other.layout["scene"][k] == junction.layout["scene"][k]
                      for k in ("camera", "aspectratio", "aspectmode"))),
            check("the 3-D window's payload: the frames in order, the first view, an id that "
                  "changes with them", [f["label"] for f in payload.frames] == ["step 0", "final"]
                  and payload.layout["scene"]["camera"]["eye"] == views.CAMERA["eye"]
                  and payload.data_id == again.data_id != other.data_id
                  and network3d.PLOTLY_JS.stat().st_size > 1_000_000),
            check("a network: 5 segments, 1 free node, 4 pinned, the cell's 12 edges",
                  (len(segments.x), len(free_nodes.x), len(pinned_nodes.x), len(box.x))
                  == (15, 1, 4, 36)),
            check("the camera is kept between redraws",
                  fig.layout.uirevision == "network" and fig.layout.scene.uirevision == "network"),
            check("the first view: x on the cube's left edge, y on its right, off the diagonal",
                  -90 < azimuth < -45 - 2 and 0 < elevation < 35.26 - 2
                  and np.isclose(np.linalg.norm([eye.x, eye.y, eye.z]), 1.15 * 1.25 * np.sqrt(3),
                                 atol=1e-3), (azimuth, elevation)),
            check("the history: nodes and density", len(history.data) == 2),
            check("a strain-rate run's history: also the stress (MPa) against the strain",
                  len(tension.data) == 3 and list(tension.data[2].y) == [0.0, 50.0]
                  and tension.layout.height > history.layout.height),
        ]

        # a run's results: 300 snapshots of the triangle, its final network, and one more
        # snapshot that is still being written
        folder = job.RunFolder.for_script("x", "s.py", "final.json", "pass\n", runs=tmp)
        folder.prepare("pass\n")
        folder.output.mkdir()
        text = (tmp / "net.json").read_text()
        for step in range(0, 3000, 10):
            (folder.output / f"disnet_{step}.json").write_text(text)
        (folder.output / "final.json").write_text(text)
        frames = results.frames(folder)
        (folder.output / "disnet_3000.json").write_text(text[:len(text) // 2])
        steps, nodes, density = results.history(folder, 3.0e-10)
        checks += [
            check("300 snapshots play as 120 frames, the first, the last and the final kept",
                  len(frames) == results.MAX_FRAMES + 1 and frames[0][0] == "step 0"
                  and frames[-2][0] == "step 2990" and frames[-1][0] == "final",
                  (len(frames), frames[:1], frames[-2:])),
            check("the history, from the whole snapshots: 3 nodes, and the line length 2 + √2 "
                  "over the volume 1000 and b²", steps == list(range(0, 3000, 10))
                  and set(nodes) == {3} and np.allclose(density, (2 + np.sqrt(2)) / 1e3 / 9e-20),
                  (steps[-1], set(nodes), density[0])),
            check("the zip holds output/", b"output/final.json" in results.output_zip(folder)),
        ]
    return summary(checks)


if __name__ == "__main__":
    sys.exit(main())
