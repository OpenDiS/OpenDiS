"""common.network3d: the 3-D window, a small Streamlit component of the app's own.

``st.plotly_chart`` puts its camera back to the figure's whenever it draws a new figure, so a view
that the user turned to does not last. This component, ``components/network3d/``, is a Plotly
figure in a frame of its own, which keeps the user's camera as it draws other networks (one
network, or the snapshots of a run, which it replays in the browser). Its page loads plotly.js
from the same folder, where this module writes it from the Plotly package when it is imported;
that copy is not committed. During a run, the page fetches the job's newest snapshot from the
folder's ``live/``, where the app writes it, so that one figure keeps the camera all along. A
background (the user stress's field on its slip plane) is drawn under every network, and sent
once, not with each snapshot.
"""
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import plotly.offline
import streamlit.components.v1 as components

from common import views
from common.network import Network

FOLDER = Path(__file__).resolve().parents[1] / "components" / "network3d"
PLOTLY_JS = FOLDER / "plotly.min.js"
LIVE = FOLDER / "live"  # a running job's newest snapshot, <run_id>.json; not committed
LIVE_MS = 500           # how often the page fetches it
# No "Share chart..." button, which plotly.js 4 adds (it uploads the figure to Plotly Cloud)
CONFIG = {"displaylogo": False, "responsive": True, "showSendToCloud": False}


@dataclass(frozen=True)
class Payload:
    """What the component draws: the frames, the layout, and an id that changes with them; and
    the traces drawn under each frame's (the background)."""

    frames: list
    layout: dict
    data_id: str
    background: list = field(default_factory=list)


def _write_plotly_js() -> None:
    """Put the Plotly package's plotly.js next to the component's page, unless it is there."""
    source = plotly.offline.get_plotlyjs().encode("utf-8")
    if not PLOTLY_JS.exists() or PLOTLY_JS.stat().st_size != len(source):
        PLOTLY_JS.write_bytes(source)


_write_plotly_js()
_COMPONENT = components.declare_component("network3d", path=str(FOLDER))


def frames_payload(frames: list[tuple[str, Network]], background: list | None = None,
                   on_plane: bool = False) -> Payload:
    """The figures of the networks (label, network), as the component takes them, over the
    background's traces: the user stress's slip plane, or nothing; on_plane draws the networks
    as over that plane (``views.network_figure``), as they are over a background."""
    background = background or []
    figures = [views.network_figure(net, f"{label} · {views.describe(net)}",
                                    on_plane=on_plane or bool(background))
               for label, net in frames]
    data = [{"label": label, "title": figure.layout.title.text,
             "data": json.loads(figure.to_json())["data"]}
            for (label, _), figure in zip(frames, figures)]
    layout = json.loads(figures[0].to_json())["layout"]
    data_id = hashlib.sha256(json.dumps([data, layout, background]).encode("utf-8")).hexdigest()
    return Payload(data, layout, data_id[:16], background)


def show(payload: Payload, key: str, scope: str, start: str = "last", frame_ms: int = 300) -> None:
    """Draw the payload: one network, or a run's snapshots to replay, from the last (or first).

    The camera that the user turns to is kept for the ``scope`` (the example): every view of the
    same scope starts in it, and another scope starts in the first view.
    """
    _COMPONENT(frames=payload.frames, layout=payload.layout, background=payload.background,
               config=CONFIG, start=start, frame_ms=frame_ms, data_id=payload.data_id,
               scope=scope, key=key, default=None)


def write_live(run_id: str, label: str, net: Network, on_plane: bool = False) -> None:
    """Write a running job's newest snapshot where ``show_live``'s page fetches it; on_plane draws
    it as over the user stress's slip plane, which ``show_live`` draws under it."""
    payload = frames_payload([(label, net)], on_plane=on_plane)
    LIVE.mkdir(exist_ok=True)
    path = LIVE / f"{run_id}.json"
    part = path.with_suffix(".part")
    part.write_text(json.dumps({"frames": payload.frames, "layout": payload.layout,
                                "data_id": payload.data_id}), encoding="utf-8")
    part.replace(path)  # the page never reads half a file


def end_live(run_id: str) -> None:
    """Remove the run's newest snapshot when the run ends."""
    (LIVE / f"{run_id}.json").unlink(missing_ok=True)


def show_live(run_id: str, key: str, scope: str, background: list | None = None) -> None:
    """Draw a running job's snapshots as ``write_live`` writes them, over the background, with
    the camera kept as by ``show``: the page fetches the newest every ``LIVE_MS`` ms, so
    Streamlit draws it once."""
    _COMPONENT(live=f"live/{run_id}.json", live_ms=LIVE_MS, background=background or [],
               config=CONFIG, scope=scope, key=key, default=None)
