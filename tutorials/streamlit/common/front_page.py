"""common.front_page: the tutorial's front page, a link to each engine with its users.

opendis.streamlit.app shows it in place of the app when the secret ``OPENDIS_FRONT_PAGE`` of its
deployment is true (``wanted()``). Each engine, ``opendis-NN.streamlit.app``, runs the app for
one user (``engines/README.md``). How many sessions an engine has comes from the gauge
``active_sessions`` that Streamlit publishes at ``/_stcore/metrics``, which needs no login.
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
from common.paths import ASSETS

SWITCH = "OPENDIS_FRONT_PAGE"  # the deployment's secret, or an environment variable for tests
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
p.opendis-guide { max-width: 46rem; margin: 0.6rem auto 1.2rem !important; text-align: center; }
table.opendis-engines { margin: 0 auto; border-collapse: separate; border-spacing: 0.6rem 0.4rem; }
.opendis-engines td { min-width: 10.5rem; padding: 0.45rem 0.8rem; border-radius: 0.4rem;
  border: 1px solid #b8c5d4; color: #1a2332; white-space: nowrap; }
.opendis-engines a { color: #1a2332 !important; font-weight: 600; text-decoration: none; }
.opendis-engines a:hover { text-decoration: underline; }
.opendis-engines .state { float: right; margin-left: 1rem; font-size: 0.9em; }
.opendis-engines .free { background: #dff3e4; }
.opendis-engines .busy { background: #fbe9c8; }
.opendis-engines .asleep { background: #eceff3; color: #5f6b78; }
.opendis-engines .missing { background: #f6f7f9; color: #9aa4b0; }
.opendis-engines .missing a { color: #9aa4b0 !important; }
</style>"""
_LOCK = threading.Lock()  # one round of reads at a time, for all sessions of this process
_LAST = {"time": float("-inf"), "engines": []}  # the last round: when, and its engines
GUIDE = ("<b>Pick a free engine.</b> Each engine runs the OpenDiS app on a server of its own, "
         "for one user. Its count is the number of browser tabs open on it now. An engine "
         "opens in a new tab. An asleep engine wakes up when you open it: press its button, "
         "and wait a minute or two.")


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
    """The table: 5 columns of 10 engines, opendis-00 to opendis-09 in the first."""
    rows = []
    for row in range(ROWS):
        cells = []
        for engine in engines[row::ROWS]:
            cells.append(
                f'<td class="{engine.state}"><a href="{html.escape(engine.url(page))}" '
                f'target="_blank" rel="noopener">{engine.name}</a>'
                f'<span class="state">{html.escape(engine.label)}</span></td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return f'<table class="opendis-engines">{"".join(rows)}</table>'


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


def _table(live: bool) -> None:
    """The table, which its fragment refreshes. After STOP_S, a full run stops the refresh."""
    st.markdown(table_html(recent()), unsafe_allow_html=True)
    if live and time.time() - st.session_state[OPENED] > STOP_S:
        st.rerun()


def _reopen() -> None:
    st.session_state[OPENED] = time.time()


def main() -> None:
    """The page: the banner, a short guide, then the table of engines."""
    st.set_page_config(page_title="OpenDiS tutorial",
                       page_icon=str(ASSETS / "opendis_favicon.png"), layout="wide")
    st.markdown(style.css() + CSS, unsafe_allow_html=True)
    st.markdown(style.banner_html(style.BANNER_TITLE, ASSETS / "opendis_logo.png"),
                unsafe_allow_html=True)
    st.markdown(f'<p class="opendis-guide">{GUIDE}</p>', unsafe_allow_html=True)
    live = time.time() - st.session_state.setdefault(OPENED, time.time()) < STOP_S
    st.fragment(_table, run_every=REFRESH_S if live else None)(live)
    if not live:
        _, middle, _ = st.columns([2, 1, 2])
        middle.button("↻ refresh the counts", width="stretch", on_click=_reopen)
