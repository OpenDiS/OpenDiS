"""common.views: the Plotly figures of the 3-D window: a network, and a run's history.

A network is drawn as PyDiS's own plot draws it, blue segments and red nodes, with the pinned
nodes marked and the cell's edges in grey. The 3-D window draws these figures through
``common/network3d.py``, which keeps the camera that the user turned to. Under the user stress's
line, it draws the frozen field on the line's slip plane: the polygons where that plane and its
periodic images cut the cell, coloured by the field's resolved shear stress, as in the figure of
the field's own tutorial; the line is then black. This module does not import Streamlit.
"""
import itertools

import numpy as np
import plotly.colors
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from common.network import PINNED, Network

SEGMENT_COLOR, NODE_COLOR, PINNED_COLOR, BOX_COLOR = "#1f4fd1", "#d62728", "#1a2332", "#9aa8b8"
DARK_COLOR = "#111111"  # the segments and nodes over the user stress's slip plane
# That plane's colours: red for a positive resolved shear stress and blue for a negative one,
# clipped at +/- FIELD_CLIP MPa, as in the figure of the field's tutorial; its grid's points
# across the cell, along each of the plane's two directions
FIELD_SCALE = plotly.colors.make_colorscale(plotly.colors.diverging.RdBu_r)
FIELD_CLIP = 15.0
PLANE_POINTS = 96
_CORNERS = np.array(list(itertools.product((-0.5, 0.5), repeat=3)))  # of the cell, per box size
_EDGES = [(i, j) for i, j in itertools.combinations(range(8), 2) if bin(i ^ j).count("1") == 1]
HEIGHT = 560
HISTORY_HEIGHT = 500  # px: the history's two plots, with their labels, inside the panel
CURVE_HEIGHT = 220    # px more for its third, the stress against the strain
# The first view: Plotly's default camera turned 90° to the right (clockwise seen from above),
# so that the x axis runs along the cube's front-left edge and the y axis along its front-right
# edge, both growing from left to right; then, as a small drag of the mouse up and to the right
# would, 7.5° further round and 5° lower (azimuth -52.5°, elevation 30.3°), so that it does not
# look along the cube's diagonal; and 15% farther than Plotly's default, so that the cube's near
# corner and its tick labels fit in the window
CAMERA = {"eye": {"x": 1.309, "y": -1.706, "z": 1.255}, "up": {"x": 0.0, "y": 0.0, "z": 1.0},
          "center": {"x": 0.0, "y": 0.0, "z": 0.0}}
# The first view over the user stress's slip plane: 22.5° further round and 10° higher (azimuth
# -30°, elevation 40°), at the same distance, so that the (111) plane is seen at 58° from its
# normal, and not nearly edge on, at 79°
PLANE_CAMERA = {**CAMERA, "eye": {"x": 1.652, "y": -0.954, "z": 1.601}}


def _lines(r1: np.ndarray, r2: np.ndarray) -> tuple[list, list, list]:
    """Segments as one polyline, with a gap (None) after each."""
    xs, ys, zs = [], [], []
    for a, b in zip(r1, r2):
        xs += [a[0], b[0], None]
        ys += [a[1], b[1], None]
        zs += [a[2], b[2], None]
    return xs, ys, zs


def _box(net: Network) -> go.Scatter3d:
    """The cell's twelve edges."""
    corners = {(i, j, k): net.origin + net.h @ np.array([i, j, k], float)
               for i in (0, 1) for j in (0, 1) for k in (0, 1)}
    edges = [(a, b) for a in corners for b in corners
             if a < b and sum(abs(p - q) for p, q in zip(a, b)) == 1]
    xs, ys, zs = _lines(np.array([corners[a] for a, _ in edges]),
                        np.array([corners[b] for _, b in edges]))
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", name="cell", hoverinfo="skip",
                        line={"color": BOX_COLOR, "width": 2}, showlegend=False)


def _traces(net: Network, on_plane: bool = False) -> list[go.Scatter3d]:
    """The segments, the free nodes, the pinned nodes and the cell, in that order; on_plane draws
    the segments and the free nodes black, over the user stress's slip plane."""
    nodes = net.nodes()
    pinned = net.constraints == PINNED
    xs, ys, zs = _lines(*net.segment_ends())
    segment, node = (DARK_COLOR, DARK_COLOR) if on_plane else (SEGMENT_COLOR, NODE_COLOR)
    return [
        go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", name="segments", hoverinfo="skip",
                     line={"color": segment, "width": 5}),
        go.Scatter3d(x=nodes[~pinned, 0], y=nodes[~pinned, 1], z=nodes[~pinned, 2],
                     mode="markers", name="nodes", marker={"color": node, "size": 3},
                     hovertemplate="(%{x:.4g}, %{y:.4g}, %{z:.4g})<extra>node</extra>"),
        go.Scatter3d(x=nodes[pinned, 0], y=nodes[pinned, 1], z=nodes[pinned, 2],
                     mode="markers", name="pinned nodes",
                     marker={"color": PINNED_COLOR, "size": 5, "symbol": "diamond"},
                     hovertemplate="(%{x:.4g}, %{y:.4g}, %{z:.4g})<extra>pinned</extra>"),
        _box(net),
    ]


def _layout(fig: go.Figure, net: Network, title: str, camera: dict) -> None:
    """The scene: the cell's bounds, equal scales, the first view, the camera kept between
    redraws."""
    low, high = net.bounds()
    sizes = high - low
    axis = {"backgroundcolor": "#f7f9fc", "gridcolor": "#dde3ea", "zerolinecolor": "#dde3ea"}
    fig.update_layout(
        height=HEIGHT, margin={"l": 0, "r": 0, "t": 36, "b": 0}, uirevision="network",
        title={"text": title, "x": 0.5, "xanchor": "center", "font": {"size": 14}},
        legend={"orientation": "h", "x": 0.0, "y": 1.0, "bgcolor": "rgba(255,255,255,0.6)"},
        scene={"xaxis": {"range": [low[0], high[0]], "title": "x", **axis},
               "yaxis": {"range": [low[1], high[1]], "title": "y", **axis},
               "zaxis": {"range": [low[2], high[2]], "title": "z", **axis},
               "aspectmode": "manual", "uirevision": "network", "camera": camera,
               "aspectratio": dict(zip("xyz", sizes / sizes.max()))})


def describe(net: Network) -> str:
    """Nodes and segments, in words."""
    return f"{len(net.positions)} nodes, {len(net.segments)} segments"


def network_figure(net: Network, title: str, on_plane: bool = False) -> go.Figure:
    """One network; on_plane draws it black, from a view that shows the user stress's slip plane
    under it."""
    fig = go.Figure(data=_traces(net, on_plane))
    _layout(fig, net, title, PLANE_CAMERA if on_plane else CAMERA)
    return fig


def plane_offset(box: float, slip: tuple, point: np.ndarray) -> float:
    """The slip plane through the point: its offset (b) along its unit normal, for its periodic
    image nearest the cell's centre; slip is (b, n), n by its Miller indices."""
    n = np.asarray(slip[1], float)
    period = box / np.linalg.norm(n)  # between the plane's periodic images
    offset = float(np.asarray(point, float) @ n) / np.linalg.norm(n)
    return round(offset - period * round(offset / period), 3)


def _section(box: float, n: np.ndarray, e1: np.ndarray, e2: np.ndarray,
             c: float) -> np.ndarray | None:
    """The corners (along e1 and e2) of the polygon where the plane n·r = c cuts the cell, in
    order round it, counterclockwise; None if the plane misses the cell."""
    corners = box * _CORNERS
    s = corners @ n - c
    points = [corners[i] + s[i] / (s[i] - s[j]) * (corners[j] - corners[i])
              for i, j in _EDGES if s[i] * s[j] < 0.0]
    points += [corners[i] for i in range(8) if s[i] == 0.0]
    uv = np.unique((np.array(points) @ np.column_stack([e1, e2])).round(6), axis=0) \
        if points else np.zeros((0, 2))
    if len(uv) < 3:
        return None
    d = uv - uv.mean(axis=0)
    return uv[np.argsort(np.arctan2(d[:, 1], d[:, 0]))]


def _clamp(uv: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    """The points, with those outside the convex polygon (corners counterclockwise) moved to its
    nearest point: then a grid over the polygon's bounds covers the polygon, and nothing more."""
    a, e = polygon, np.roll(polygon, -1, axis=0) - polygon  # each edge, from a along e
    rel = uv[:, None, :] - a
    inside = (e[:, 0] * rel[..., 1] - e[:, 1] * rel[..., 0] >= 0.0).all(axis=1)
    nearest = a + np.clip((rel * e).sum(axis=2) / (e * e).sum(axis=1), 0.0, 1.0)[..., None] * e
    pick = ((uv[:, None, :] - nearest) ** 2).sum(axis=2).argmin(axis=1)
    return np.where(inside[:, None], uv, nearest[np.arange(len(uv)), pick])


def slip_plane(field, box: float, slip: tuple, offset: float) -> list[dict]:
    """The user stress's frozen field on its line's slip plane, as the Plotly surfaces (JSON)
    that the 3-D window draws under the networks. There is one for each polygon where the plane,
    at the offset (b) along its unit normal, or one of its periodic images cuts the cell: the
    line's pieces lie on them. Each is coloured by the field's resolved shear stress n·σ·b (MPa),
    to which the applied stress adds. The field takes positions (b) and gives the stresses (Pa:
    xx yy zz yz xz xy), as the script's StressField does; slip is (b, n)."""
    b, n = (np.asarray(v, float) / np.linalg.norm(v) for v in slip)
    e2 = np.cross(n, b)  # the plane's second direction, after b
    period = box / np.linalg.norm(slip[1])
    images = np.ceil(0.5 * box * np.abs(n).sum() / period)  # the cell's corners, along n
    spacing = box / PLANE_POINTS
    surfaces = []
    for c in offset + period * np.arange(-images, images + 1):
        polygon = _section(box, n, b, e2, c)
        if polygon is None:
            continue
        low, high = polygon.min(axis=0), polygon.max(axis=0)
        counts = np.maximum(np.ceil((high - low) / spacing).astype(int) + 1, 2)
        u, v = np.meshgrid(np.linspace(low[0], high[0], counts[0]),
                           np.linspace(low[1], high[1], counts[1]))
        uv = _clamp(np.column_stack([u.ravel(), v.ravel()]), polygon)
        r = c * n + uv[:, :1] * b + uv[:, 1:] * e2
        stress = np.asarray(field(r))[:, [0, 5, 4, 5, 1, 3, 4, 3, 2]].reshape(-1, 3, 3)
        tau = np.einsum("i,kij,j->k", n, stress, b) / 1e6

        def grid(values: np.ndarray, digits: int) -> list:  # rounded, so that the JSON is short
            values = values.reshape(u.shape).round(digits)
            return (values.astype(int) if digits == 0 else values).tolist()

        surfaces.append({
            "type": "surface", "name": "slip plane",
            "x": grid(r[:, 0], 0), "y": grid(r[:, 1], 0), "z": grid(r[:, 2], 0),
            "surfacecolor": grid(tau, 1), "colorscale": FIELD_SCALE,
            "cmin": -FIELD_CLIP, "cmax": FIELD_CLIP, "showscale": not surfaces,
            "colorbar": {"title": {"text": "internal resolved shear stress (MPa)",
                                   "side": "right"},
                         "thickness": 12, "len": 0.6, "x": 1.0, "xanchor": "right", "y": 0.45},
            "lighting": {"ambient": 1.0, "diffuse": 0.0, "specular": 0.0, "roughness": 1.0,
                         "fresnel": 0.0},
            # no hover: plotly.js reads a surface's arrays transposed for a hover template
            "hoverinfo": "skip",
        })
    return surfaces


def second_half(time: np.ndarray, disp: np.ndarray) -> tuple[float, float, float]:
    """The second half of a run: the time (s) and the displacement (m) at its start, which is
    interpolated, and the mean velocity over it (m/s)."""
    start = 0.5 * time[-1]
    at = float(np.interp(start, time, disp))
    return start, at, float((disp[-1] - at) / (time[-1] - start)) if start > 0 else float("nan")


def _motion_title(motion: tuple) -> str:
    """The title of the displacement's plot, with its velocity and the drag velocity."""
    time, disp, drag = motion
    title = f"displacement (µm) against time (ns): {second_half(time, disp)[2]:.4g} m/s"
    return title + ("" if drag is None else f", against M b τ = {drag:.4g} m/s")


def history_figure(steps: list[int], nodes: list[int], density: list[float],
                   curve: tuple | None = None, motion: tuple | None = None) -> go.Figure:
    """The number of nodes, and the dislocation density, against the time step; with a curve
    (strain, stress in Pa), also the stress against the strain; with motion (the time in s, the
    user stress's line's mean displacement in m, and the drag velocity M b τ in m/s or None),
    also the displacement against the time, with its velocity over the second half (dashed)
    and the drag velocity's line (dotted)."""
    if curve is None and motion is None:
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                            subplot_titles=("nodes", "dislocation density (m⁻²)"))
    else:
        third = "stress (MPa) against strain" if curve is not None else _motion_title(motion)
        fig = make_subplots(rows=3, cols=1, vertical_spacing=0.15,
                            subplot_titles=("nodes", "dislocation density (m⁻²)", third))
        fig.update_xaxes(matches="x2", showticklabels=False, row=1, col=1)
    fig.add_trace(go.Scatter(x=steps, y=nodes, mode="lines+markers", name="nodes",
                             line={"color": SEGMENT_COLOR}), row=1, col=1)
    fig.add_trace(go.Scatter(x=steps, y=density, mode="lines+markers", name="density",
                             line={"color": NODE_COLOR}), row=2, col=1)
    fig.update_xaxes(title_text="time step", row=2, col=1)
    if curve is not None:
        strain, stress = curve
        fig.add_trace(go.Scatter(x=list(strain), y=[s / 1e6 for s in stress], mode="lines",
                                 name="stress", line={"color": "#1a2332"}), row=3, col=1)
        fig.update_xaxes(title_text="strain", row=3, col=1)
    if motion is not None:
        time, disp, drag = motion
        start, at, _ = second_half(time, disp)
        end = float(time[-1])
        fig.add_trace(go.Scatter(x=list(time * 1e9), y=list(disp * 1e6), mode="lines",
                                 name="displacement", line={"color": "#1a2332"}), row=3, col=1)
        fig.add_trace(go.Scatter(x=[start * 1e9, end * 1e9], y=[at * 1e6, disp[-1] * 1e6],
                                 mode="lines", name="velocity over the second half",
                                 line={"color": NODE_COLOR, "dash": "dash"}), row=3, col=1)
        if drag is not None:
            fig.add_trace(go.Scatter(x=[0.0, end * 1e9], y=[0.0, drag * end * 1e6],
                                     mode="lines", name="M b τ",
                                     line={"color": BOX_COLOR, "dash": "dot"}), row=3, col=1)
        fig.update_xaxes(title_text="time (ns)", row=3, col=1)
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)  # the tick labels and the axis title stay in view
    fig.update_layout(height=HISTORY_HEIGHT + (0 if curve is None and motion is None
                                               else CURVE_HEIGHT),
                      margin={"l": 70, "r": 20, "t": 40, "b": 60},
                      showlegend=False, plot_bgcolor="#f7f9fc")
    return fig
