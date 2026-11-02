"""tests.test_views: the networks, a run's results and the figures (README.md, Tests)."""
import sys
import tempfile
from pathlib import Path

import numpy as np

from common import catalog, codegen, job, network, network3d, results, views
from common.paths import add_pydis_paths
from common.spec import LINE_SLIP
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
        time = np.array([0.0, 1e-9, 2e-9, 3e-9, 4e-9])  # 100 m/s, then 50 m/s from 1 ns
        disp = np.array([0.0, 1e-7, 1.5e-7, 2e-7, 2.5e-7])
        motion = views.history_figure([0, 10], [7, 9], [5e12, 6e12], motion=(time, disp, 60.0))
        titles = [a.text for a in motion.layout.annotations]
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
            check("the user stress's history: also the displacement (µm) against the time (ns), "
                  "the velocity over the second half, and the drag velocity's line",
                  len(motion.data) == 5 and np.allclose(motion.data[2].y, disp * 1e6)
                  and np.allclose(motion.data[3].x, [2.0, 4.0])
                  and np.allclose(motion.data[3].y, [0.15, 0.25])
                  and np.allclose(motion.data[4].y, [0.0, 0.24]), list(motion.data[3].y)),
            check("... whose title gives the two velocities, 50 and 60 m/s", titles[2] ==
                  "displacement (µm) against time (ns): 50 m/s, against M b τ = 60 m/s", titles),
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
        (folder.output / "displacement.dat").write_text(
            "# Step Time(s) Displacement(m)\n0 0.000000e+00 0.000000e+00\n")
        one_row = results.displacement(folder)
        (folder.output / "displacement.dat").write_text(
            "# Step Time(s) Displacement(m)\n0 0.000000e+00 0.000000e+00\n"
            "1 1.000000e-11 1.036107e-09\n2 2.200000e-11 2.279167e-09\n")
        moved = results.displacement(folder)
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
            check("the user stress's displacement against the time, from output/displacement.dat; "
                  "none before its second row", one_row is None and moved is not None
                  and np.allclose(moved[0], [0.0, 1e-11, 2.2e-11])
                  and np.allclose(moved[1], [0.0, 1.036107e-9, 2.279167e-9]), moved),
        ]

    # the user stress's slip plane, over a uniform field whose resolved shear stress is 5 MPa
    spec = next(e.default for e in catalog.examples() if e.name == "user_stress")
    box, line = spec.config.box_length, found["user_stress"]
    b, n = (np.array(v) / np.linalg.norm(v) for v in LINE_SLIP)
    voigt = (5e6 * (np.outer(b, n) + np.outer(n, b)))[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]
    offset = views.plane_offset(box, LINE_SLIP, codegen.line_origin(box, spec.config.seed))
    surfaces = views.slip_plane(lambda r: np.tile(voigt, (len(r), 1)), box, LINE_SLIP, offset)
    period = box / np.sqrt(3)  # between the plane's periodic images

    def off_plane(r: np.ndarray) -> np.ndarray:
        d = r @ n - offset
        return np.abs(d - period * np.round(d / period))

    grids = [np.stack([np.array(s[axis], float) for axis in "xyz"], axis=-1) for s in surfaces]
    area = sum(0.5 * np.linalg.norm(np.cross(g[1:, 1:] - g[:-1, :-1], g[:-1, 1:] - g[1:, :-1]),
                                    axis=-1).sum() for g in grids)
    over = views.network_figure(line, "over the plane", on_plane=True)
    plain = network3d.frames_payload([("initial configuration", line)])
    payload = network3d.frames_payload([("initial configuration", line)], surfaces)
    checks += [
        check("the user stress's slip plane: three polygons in the cell, the colour bar on the "
              "first", [s["showscale"] for s in surfaces] == [True, False, False]),
        check("... on the plane through the line's start, to within the rounding to 1 b, and "
              "inside the cell", all(off_plane(g.reshape(-1, 3)).max() < 1.0
                                     and np.abs(g).max() <= 0.5 * box + 0.5 for g in grids)),
        check("... covering the plane's area in the cell, √3 L² (the rounding to 1 b gives the "
              "polygons' edges slivers of 3e-4 of it)",
              np.isclose(area, np.sqrt(3) * box**2, rtol=1e-3), area / (np.sqrt(3) * box**2)),
        check("... coloured by the resolved shear stress: 5 MPa in the uniform field",
              all(np.allclose(s["surfacecolor"], 5.0) for s in surfaces)),
        check("the user stress's initial line lies on that plane", off_plane(line.positions).max()
              < 1e-3, off_plane(line.positions).max()),
        check("over the plane, the network is black, in the plane's first view",
              over.data[0].line.color == views.DARK_COLOR == over.data[1].marker.color
              and over.layout.scene.camera.eye.to_plotly_json() == views.PLANE_CAMERA["eye"]),
        check("the payload keeps the plane apart from its frames, and its id changes with it",
              payload.background == surfaces and plain.background == []
              and payload.data_id != plain.data_id and len(payload.frames[0]["data"]) == 4),
    ]

    # user_stress.py's field on the plane: tau = n.sigma.b = -Sxx/sqrt(6), from z and the time
    toy = views.slip_plane(network.space_time_field(box), box, LINE_SLIP, offset,
                           views.SPACE_TIME_CLIP, views.SPACE_TIME_TITLE, views.SPACE_TIME_POINTS,
                           0)
    later = views.plane_colors(network.space_time_field(box, 0.5e-9 * np.pi), box, LINE_SLIP,
                               offset, views.SPACE_TIME_POINTS, 0)  # sin(t/1 ns) = 1
    z = np.concatenate([np.ravel(s["z"]) for s in toy])
    tau = [np.concatenate([np.ravel(c) for c in colours])
           for colours in ([s["surfacecolor"] for s in toy], later)]
    lines = ["step = 0, nodes = 87, dt = 1.0e-11, time = 1.0e-11, elapsed = 0.5 sec",
             "step = 10, nodes = 87, dt = 6.0e-11, time = 3.0e-10, elapsed = 1.0 sec"]
    moving = network3d.frames_payload([("step 0", line), ("step 10", line)], toy,
                                      colors=[[s["surfacecolor"] for s in toy], later])
    checks += [
        check("user_stress.py's field on the plane: -(1 + cos(2π z/L) + sin(t/1 ns)) 50/√6 MPa, "
              "at t = 0 and at π/2 ns, in whole MPa, under its own colour bar at +/- 60 MPa",
              all(np.allclose(t, -(1.0 + np.cos(2.0 * np.pi * z / box) + sine) * 50.0
                              / np.sqrt(6.0), atol=0.51) for t, sine in zip(tau, (0.0, 1.0)))
              and toy[0]["cmax"] == 60.0
              and toy[0]["colorbar"]["title"]["text"] == views.SPACE_TIME_TITLE,
              [(t.min(), t.max()) for t in tau]),
        check("the frames' times: a step's own line, else interpolated, and the final's apart",
              np.allclose(results.frame_times(lines, ["step 0", "step 5", "step 10", "final"],
                                              4e-10), [1e-11, 1.55e-10, 3e-10, 4e-10])),
        check("each frame carries the plane's colours at its time, on the plane's grids",
              [len(frame["colors"]) for frame in moving.frames] == [len(toy)] * 2
              and all(np.shape(c) == np.shape(s["surfacecolor"]) for c, s in zip(later, toy))),
    ]

    # the precipitates: a sphere of each, and of each image that reaches into the cell
    near = views.spheres(np.array([[100.0, 25000.0, 25000.0], [25000.0, 25000.0, 25000.0]]),
                         np.array([500.0, 1000.0]), 50000.0)
    xyz = np.stack([np.array(near[0][axis], float) for axis in "xyz"], axis=-1).reshape(-1, 3)
    checks.append(check("the precipitates: a sphere of each radius, and one more for the image "
                        "of the sphere at the cell's face", len(near) == 3 and np.allclose(
                            np.linalg.norm(xyz - [100.0, 25000.0, 25000.0], axis=1), 500.0,
                            atol=1.0) and np.array(near[1]["x"]).min() > 49000.0,
                        [np.array(s["x"]).mean().round() for s in near]))
    return summary(checks)


if __name__ == "__main__":
    sys.exit(main())
