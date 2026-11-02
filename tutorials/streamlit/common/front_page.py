"""common.front_page: the tutorial's front page, a link to each engine with its users.

opendis.streamlit.app shows it in place of the app when the secret ``OPENDIS_FRONT_PAGE`` of its
deployment is true (``wanted()``). Each engine, ``opendis-NN.streamlit.app``, runs the app for
one user (``engines/README.md``). How many sessions an engine has comes from the gauge
``active_sessions`` that Streamlit publishes at ``/_stcore/metrics``, which needs no login.
Above the table, the page shows the announcements of ``engines/announcements.md``.
"""
import concurrent.futures
import html
import os
import re
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass

import streamlit as st

from common import style
from common.paths import APP, ASSETS
from engines.engine import IDLE

SWITCH = "OPENDIS_FRONT_PAGE"  # the deployment's secret, or an environment variable for tests
NEWS = APP / "engines" / "announcements.md"  # Markdown; an empty file shows no announcement
COUNT, ROWS = 50, 10  # the engines, opendis-00 to opendis-49, in 5 columns of 10
NAME = "opendis-{n:02d}"
# An engine's page and its gauges on the Cloud. The environment variables of the same names
# override them, as for engines on this computer: "http://localhost:87{n:02d}"
PAGE = os.environ.get("OPENDIS_ENGINE_PAGE", "https://{name}.streamlit.app")
METRICS = os.environ.get("OPENDIS_ENGINE_METRICS",
                         "https://{name}.streamlit.app/~/+/_stcore/metrics")
TIMEOUT_S = 5     # the longest wait for one engine's gauges
KEEP_S = 10       # the counts are read again after this long, for all visitors at once
REFRESH_S = 10    # the table refreshes itself this often...
STOP_S = 30 * 60  # ...until this long after the page opened; then a click refreshes it
OPENED = "front_page_opened"  # session key: when the table began to refresh itself
GAUGE = re.compile(r"^active_sessions\s+(\d+)(?:\.0*)?\s*$", re.MULTILINE)
CSS = """<style>
[data-testid="stLayoutWrapper"]:has(> .st-key-opendis_head), .st-key-opendis_head {
  position: fixed; top: 0; left: 0; width: 0; height: 0; z-index: 2147483000; }
table.opendis-engines { width: 100%; max-width: 59rem; margin: 0.2rem auto 0; table-layout: fixed;
  border-collapse: separate; border-spacing: 0.45rem 0.35rem; font-size: 0.85rem; }
.opendis-engines td { padding: 0.3rem 0.65rem; border-radius: 0.35rem; border: 1px solid #b8c5d4;
  color: #1a2332; white-space: nowrap; overflow: hidden; line-height: 1.4; }
.opendis-engines td a { font-size: 0.85rem !important; }
.opendis-engines a { color: #1a2332 !important; font-weight: 600; text-decoration: none; }
.opendis-engines a:hover { text-decoration: underline; }
.opendis-engines .state { float: right; margin-left: 0.5rem; font-size: 0.8rem !important; }
.opendis-engines .free { background: #dff3e4; }
.opendis-engines .busy { background: #fbe9c8; }
.opendis-engines .asleep { background: #eceff3; color: #5f6b78; }
.opendis-engines .missing { background: #f6f7f9; color: #9aa4b0; }
.opendis-engines .missing a { color: #9aa4b0 !important; }
.opendis-table { container-type: inline-size; }
.opendis-engines .prefix { font-size: inherit !important; }
/* A table narrower than 30.2rem has boxes too narrow for opendis-00 (4.79em, with the paddings
   and gaps): every box then shows only the number, 00 */
@container (max-width: 30.2rem) { .opendis-engines .prefix { display: none; } }
.st-key-opendis_news, .st-key-opendis_guide { max-width: 59rem; margin: 0 auto !important;
  padding: 0.7rem 1.1rem 0.95rem; border-radius: 0.4rem; }
.st-key-opendis_news { background: #fff6dc; border: 1px solid #e9c46a;
  border-left: 0.35rem solid #e0a800; }
.st-key-opendis_guide { background: #eef4fa; border: 1px solid #b8c5d4;
  border-left: 0.35rem solid #5b7fa6; }
.st-key-opendis_news p, .st-key-opendis_guide p { color: #1a2332; margin: 0 0 0.55rem !important; }
.st-key-opendis_news p:last-child, .st-key-opendis_guide p:last-child { margin: 0 !important; }
.st-key-opendis_guide ul { margin: 0 !important; padding-left: 1.3rem; }
.st-key-opendis_guide li { color: #1a2332; margin: 0.1rem 0 !important; }
.st-key-opendis_guide [data-testid="stMarkdownContainer"] p,
.st-key-opendis_guide [data-testid="stMarkdownContainer"] li,
.st-key-opendis_guide [data-testid="stMarkdownContainer"] p *,
.st-key-opendis_guide [data-testid="stMarkdownContainer"] li * { font-size: 0.95rem !important; }
.st-key-opendis_note { max-width: 59rem; margin: 0.4rem auto 0 !important; }
</style>"""
_LOCK = threading.Lock()  # one round of reads at a time, for all sessions of this process
_LAST = {"time": float("-inf"), "engines": []}  # the last round: when, and its engines
GUIDE = ("**Pick a free engine** by clicking on one of the links below. Each engine runs the "
         "OpenDiS app on a server of its own. "
         "Its count is the number of browser tabs open on it now.\n\n"
         "- Every engine is shared with anyone who opens it. Another user may then press "
         "*clear all results*, which wipes out your results too.\n"
         "- Closing the browser tab releases the engine.\n"
         f"- After {IDLE[1]:g} minutes of idling, the engine signs you off automatically, so "
         "that another user can have it.")
SIGNED_OFF = re.compile(r"opendis-\d\d")  # the name in the address after a sign-off


@dataclass(frozen=True)
class Engine:
    """One engine: its number, and its sessions when its gauges could be read."""
    n: int
    users: int | None = None
    missing: bool = False  # no such app on the Cloud: its gauges' address answered 404

    @property
    def name(self) -> str:
        """Its name, which is its subdomain on the Cloud, e.g. ``opendis-07``."""
        return NAME.format(n=self.n)

    def url(self, page: str | None = None) -> str:
        """Its page, by PAGE or another template of ``name`` and ``n``."""
        return (page or PAGE).format(name=self.name, n=self.n)

    @property
    def state(self) -> str:
        """``free``, ``busy`` (with users), ``asleep`` (or starting, or failing), ``missing``."""
        if self.missing:
            return "missing"
        if self.users is None:
            return "asleep"
        return "busy" if self.users else "free"

    @property
    def label(self) -> str:
        """What the page says of its state: ``free``, ``1 user``, ``3 users``, ``asleep``..."""
        if self.state == "busy":
            return f"{self.users} user" + ("s" if self.users > 1 else "")
        return {"missing": "not deployed"}.get(self.state, self.state)


def parse(text: str) -> int | None:
    """The number of active sessions in Streamlit's metrics, or None."""
    found = GAUGE.search(text)
    return int(found.group(1)) if found else None


def read(n: int, metrics: str | None = None, timeout: float = TIMEOUT_S) -> Engine:
    """Engine ``n``, with its sessions read from its gauges (by METRICS, or another template).
    An engine that does not answer with them is asleep, starting or failing: the Cloud then
    answers with a page of its own."""
    address = (metrics or METRICS).format(name=NAME.format(n=n), n=n)
    try:
        with urllib.request.urlopen(address, timeout=timeout) as reply:
            text = reply.read(1 << 16).decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        return Engine(n, missing=exc.code == 404)
    except (OSError, ValueError):  # no connection, a timeout, a bad address
        return Engine(n)
    return Engine(n, parse(text))


def read_all(metrics: str | None = None, timeout: float = TIMEOUT_S) -> list[Engine]:
    """Every engine, in order, with their gauges read at once."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=COUNT) as pool:
        return list(pool.map(lambda n: read(n, metrics, timeout), range(COUNT)))


def table_html(engines: list[Engine], page: str | None = None) -> str:
    """The table: 5 columns of 10 engines, opendis-00 to opendis-09 in the first. In a window
    too narrow for the names, the boxes show only the numbers (CSS)."""
    rows = []
    for row in range(ROWS):
        cells = []
        for engine in engines[row::ROWS]:
            prefix, number = engine.name.rsplit("-", 1)
            cells.append(
                f'<td class="{engine.state}"><a href="{html.escape(engine.url(page))}" '
                f'target="_blank" rel="noopener"><span class="prefix">{prefix}-</span>{number}</a>'
                f'<span class="state">{html.escape(engine.label)}</span></td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    table = f'<table class="opendis-engines">{"".join(rows)}</table>'
    return f'<div class="opendis-table">{table}</div>'  # the container of the CSS's @container


def wanted() -> bool:
    """Whether this deployment shows the front page: ``OPENDIS_FRONT_PAGE`` is true in the
    environment, or else in the secrets. A boolean secret does not reach the environment."""
    value = os.environ.get(SWITCH)
    if value is None and st.secrets.load_if_toml_exists():
        value = st.secrets.get(SWITCH)
    return str(value).strip().lower() in ("1", "true", "yes")


def recent() -> list[Engine]:
    """The engines as read at most KEEP_S ago: a session that comes during a round waits for it,
    so that the 50 engines are read once for all visitors."""
    with _LOCK:
        if time.monotonic() - _LAST["time"] > KEEP_S:
            _LAST.update(engines=read_all(), time=time.monotonic())
        return _LAST["engines"]


def news() -> str:
    """The announcements, in Markdown; nothing when their file is empty or missing."""
    return NEWS.read_text(encoding="utf-8").strip() if NEWS.is_file() else ""


def _body(live: bool) -> None:
    """The announcements and the guide, each in its box, then the table, which their fragment
    refreshes: a pushed announcement shows without a reload. After STOP_S, a full run stops
    the refresh."""
    text = news()
    if text:
        with st.container(key="opendis_news"):
            st.markdown(text)
    with st.container(key="opendis_guide"):
        st.markdown(GUIDE)
    st.markdown(table_html(recent()), unsafe_allow_html=True)
    if live and time.time() - st.session_state[OPENED] > STOP_S:
        st.rerun()


def _reopen() -> None:
    st.session_state[OPENED] = time.time()


def _signed_off() -> None:
    """A note when an idle engine sent the visitor here: ?signed_off=opendis-07&idle=60."""
    name, minutes = st.query_params.get("signed_off", ""), st.query_params.get("idle", "")
    if SIGNED_OFF.fullmatch(name):
        after = ""
        if minutes.replace(".", "", 1).isdigit():
            after = f"after {minutes} minute{'' if float(minutes) == 1 else 's'} "
        with st.container(key="opendis_note"):
            st.info(f"{name} signed you off {after}without activity. Open it again below, "
                    "if it is still free.")


def main() -> None:
    """The page: the app's banner, then the announcements, a short guide and the table."""
    st.set_page_config(page_title="OpenDiS tutorial",
                       page_icon=str(ASSETS / "opendis_favicon.png"), layout="wide")
    with st.container(key="opendis_head"):  # out of the page's flow: no gaps above the boxes
        st.markdown(style.css() + CSS, unsafe_allow_html=True)
        st.markdown(style.banner_html(style.BANNER_TITLE, ASSETS / "opendis_logo.png"),
                    unsafe_allow_html=True)
    _signed_off()
    live = time.time() - st.session_state.setdefault(OPENED, time.time()) < STOP_S
    st.fragment(_body, run_every=REFRESH_S if live else None)(live)
    if not live:
        _, middle, _ = st.columns([2, 1, 2])
        middle.button("↻ refresh the counts", width="stretch", on_click=_reopen)
