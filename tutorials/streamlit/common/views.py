"""common.views: the Plotly figures of the 3-D window: a network, and a run's history.

A network is drawn as PyDiS's own plot draws it, blue segments and red nodes, with the pinned
nodes marked and the cell's edges in grey. The 3-D window draws these figures through
``common/network3d.py``, which keeps the camera that the user turned to. This module does not
import Streamlit.
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from common.network import PINNED, Network

SEGMENT_COLOR, NODE_COLOR, PINNED_COLOR, BOX_COLOR = "#1f4fd1", "#d62728", "#1a2332", "#9aa8b8"
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


def _traces(net: Network) -> list[go.Scatter3d]:
    """The segments, the free nodes, the pinned nodes and the cell, in that order."""
    nodes = net.nodes()
    pinned = net.constraints == PINNED
    xs, ys, zs = _lines(*net.segment_ends())
    return [
        go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", name="segments", hoverinfo="skip",
                     line={"color": SEGMENT_COLOR, "width": 5}),
        go.Scatter3d(x=nodes[~pinned, 0], y=nodes[~pinned, 1], z=nodes[~pinned, 2],
                     mode="markers", name="nodes", marker={"color": NODE_COLOR, "size": 3},
                     hovertemplate="(%{x:.4g}, %{y:.4g}, %{z:.4g})<extra>node</extra>"),
        go.Scatter3d(x=nodes[pinned, 0], y=nodes[pinned, 1], z=nodes[pinned, 2],
                     mode="markers", name="pinned nodes",
                     marker={"color": PINNED_COLOR, "size": 5, "symbol": "diamond"},
                     hovertemplate="(%{x:.4g}, %{y:.4g}, %{z:.4g})<extra>pinned</extra>"),
        _box(net),
    ]


def _layout(fig: go.Figure, net: Network, title: str) -> None:
    """The scene: the cell's bounds, equal scales, the camera kept between redraws."""
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
               "aspectmode": "manual", "uirevision": "network", "camera": CAMERA,
               "aspectratio": dict(zip("xyz", sizes / sizes.max()))})


def describe(net: Network) -> str:
    """Nodes and segments, in words."""
    return f"{len(net.positions)} nodes, {len(net.segments)} segments"


def network_figure(net: Network, title: str) -> go.Figure:
    """One network."""
    fig = go.Figure(data=_traces(net))
    _layout(fig, net, title)
    return fig


def history_figure(steps: list[int], nodes: list[int], density: list[float],
                   curve: tuple | None = None) -> go.Figure:
    """The number of nodes, and the dislocation density, against the time step; with a curve
    (strain, stress in Pa), also the stress against the strain."""
    if curve is None:
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                            subplot_titles=("nodes", "dislocation density (m⁻²)"))
    else:
        fig = make_subplots(rows=3, cols=1, vertical_spacing=0.15,
                            subplot_titles=("nodes", "dislocation density (m⁻²)",
                                            "stress (MPa) against strain"))
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
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)  # the tick labels and the axis title stay in view
    fig.update_layout(height=HISTORY_HEIGHT + (0 if curve is None else CURVE_HEIGHT),
                      margin={"l": 70, "r": 20, "t": 40, "b": 60},
                      showlegend=False, plot_bgcolor="#f7f9fc")
    return fig
