"""common.network3d: the 3-D window, a small Streamlit component of the app's own.

``st.plotly_chart`` puts its camera back to the figure's whenever it draws a new figure, so a view
that the user turned to does not last. This component, ``components/network3d/``, is a Plotly
figure in a frame of its own, which keeps the user's camera as it draws other networks (one
network, or the snapshots of a run, which it replays in the browser). Its page loads plotly.js
from the same folder, where this module writes it from the Plotly package when it is imported;
that copy is not committed.
"""
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import plotly.offline
import streamlit.components.v1 as components

from common import views
from common.network import Network

FOLDER = Path(__file__).resolve().parents[1] / "components" / "network3d"
PLOTLY_JS = FOLDER / "plotly.min.js"
CONFIG = {"displaylogo": False, "responsive": True}


@dataclass(frozen=True)
class Payload:
    """What the component draws: the frames, the layout, and an id that changes with them."""

    frames: list
    layout: dict
    data_id: str


def _write_plotly_js() -> None:
    """Put the Plotly package's plotly.js next to the component's page, unless it is there."""
    source = plotly.offline.get_plotlyjs().encode("utf-8")
    if not PLOTLY_JS.exists() or PLOTLY_JS.stat().st_size != len(source):
        PLOTLY_JS.write_bytes(source)


_write_plotly_js()
_COMPONENT = components.declare_component("network3d", path=str(FOLDER))


def frames_payload(frames: list[tuple[str, Network]]) -> Payload:
    """The figures of the networks (label, network), as the component takes them."""
    figures = [views.network_figure(net, f"{label} · {views.describe(net)}")
               for label, net in frames]
    data = [{"label": label, "title": figure.layout.title.text,
             "data": json.loads(figure.to_json())["data"]}
            for (label, _), figure in zip(frames, figures)]
    layout = json.loads(figures[0].to_json())["layout"]
    data_id = hashlib.sha256(json.dumps([data, layout]).encode("utf-8")).hexdigest()[:16]
    return Payload(data, layout, data_id)


def show(payload: Payload, key: str, scope: str, start: str = "last", frame_ms: int = 300) -> None:
    """Draw the payload: one network, or a run's snapshots to replay, from the last (or first).

    The camera that the user turns to is kept for the ``scope`` (the example): every view of the
    same scope starts in it, and another scope starts in the first view.
    """
    _COMPONENT(frames=payload.frames, layout=payload.layout, config=CONFIG, start=start,
               frame_ms=frame_ms, data_id=payload.data_id, scope=scope, key=key, default=None)
