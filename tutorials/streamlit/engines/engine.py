"""engines/engine.py: what each entry file of this folder runs, the app of ``app.py``.

Streamlit Community Cloud deploys one app per entry file. So each engine,
``opendis-NN.streamlit.app``, has its own, ``opendis_NN.py``, which calls ``main()`` (README.md).
The right end of the engine's banner shows its number, and its users at that moment. A session
without a click for an hour is signed off: its tab goes back to the front page, which frees the
engine. Use that stays in the browser, as a turn of the 3-D window, closes the warning first.
"""
import json
import math
import os
import sys
import time
from pathlib import Path

import streamlit as st

APP_DIR = str(Path(__file__).resolve().parents[1])  # tutorials/streamlit, where app.py is
RESET = Path(__file__).with_name("RESET")  # a push that changes its line resets every engine
SEEN = "OPENDIS_ENGINE_RESET"  # the line of RESET that this process saw last, in its environment
FRONT = os.environ.get("OPENDIS_FRONT_URL", "https://opendis.streamlit.app")  # the front page
# Minutes without a click before a session is signed off (0: never), and of the warning before
# it. An environment variable, or a secret of the deployment, of the same name overrides each.
IDLE = ("OPENDIS_IDLE_MINUTES", 60.0)
WARNING = ("OPENDIS_IDLE_WARNING_MINUTES", 5.0)
WATCH_S = 10  # the banner's count of users and the idle check refresh this often
TOUCHED = "engine_touched"  # session key: when the last page run began, which a click causes
CSS = """<style>
[data-testid="stLayoutWrapper"]:has(> .st-key-opendis_engine), .st-key-opendis_engine {
  position: fixed; top: 0; left: 0; width: 0; height: 0; z-index: 2147483002; }
.opendis-users { position: fixed; top: 0; right: 1rem; height: var(--opendis-banner-h, 2.9rem);
  display: flex; align-items: center; gap: 0.3rem; z-index: 2147483002; color: #4a5a6e;
  font-size: 0.85rem; }
.opendis-users .num { font-family: "Source Code Pro", ui-monospace, Menlo, Consolas, monospace;
  font-size: 1.05rem; font-weight: 600; color: #1a2332; padding: 0 0.3rem; border-radius: 0.3rem; }
.opendis-users .num.shared { background: #fbe9c8; box-shadow: 0 0 0 1px #e0a800; }
.opendis-users .sep { margin: 0 0.3rem; color: #8a96a3; }
.st-key-opendis_idle { position: fixed; top: calc(var(--opendis-banner-h, 2.9rem) + 0.6rem);
  left: 50%; transform: translateX(-50%); width: min(36rem, 92vw); z-index: 2147483002;
  background: #ffffff; border-radius: 0.5rem; padding: 0.6rem 0.8rem 0.8rem;
  box-shadow: 0 4px 18px rgba(26, 35, 50, 0.25); }
</style>"""
# Use of the page that stays in the browser, as a turn of the 3-D window (in its own frame of the
# same origin) or a zoom of a chart: presses, keys and the wheel anywhere in the tab. Not pointer
# moves, which the browser may send by itself when the page changes under a still cursor.
TRACK = """<script>
(function () {
  var w = window;
  if (w.opendisTrack) return;
  w.opendisTrack = true;
  w.opendisUsed = Date.now();
  function mark() { w.opendisUsed = Date.now(); }
  function hook(win) {
    try {
      if (!win || win.opendisHooked) return;
      win.opendisHooked = true;
      ["pointerdown", "keydown", "wheel", "touchstart"].forEach(function (name) {
        win.addEventListener(name, mark, {capture: true, passive: true});
      });
    } catch (e) {}
  }
  hook(w);
  setInterval(function () {
    var frames = w.document.querySelectorAll("iframe");
    for (var i = 0; i < frames.length; i++) hook(frames[i].contentWindow);
  }, 2000);
})();
</script>"""
# During the warning: Stay, by itself, if the page was used in the warning's time (RECENT_MS)
STAY = """<script>
(function () {
  var w = window;
  if (w.opendisStay) return;
  w.opendisStay = setInterval(function () {
    var button = w.document.querySelector(".st-key-opendis_stay button");
    var used = Date.now() - (w.opendisUsed || 0) < RECENT_MS;
    if (!button || used) {
      clearInterval(w.opendisStay);
      w.opendisStay = null;
      if (button) button.click();
    }
  }, 1000);
})();
</script>"""
# The sign-off: a link that the tab's own page follows. On the Cloud, the app runs in a frame that
# may not move the tab itself, but the page around it, of the same origin, may. Else the frame
# leaves the app, which ends the session all the same.
SIGN_OFF = """<script>
(function () {{
  var url = {url};
  try {{
    var page = window.top.document, link = page.createElement("a");
    link.href = url;
    page.body.appendChild(link);
    link.click();
  }} catch (e) {{}}
  setTimeout(function () {{ window.location.replace(url); }}, 3000);
}})();
</script>"""


def first(folder: str, paths: str) -> str:
    """A search path, folders joined by ``os.pathsep``, with ``folder`` first and only once."""
    return os.pathsep.join([folder] + [p for p in paths.split(os.pathsep) if p and p != folder])


def reset_due() -> bool:
    """Whether RESET holds a line that this process has not seen since its first page run."""
    line = RESET.read_text(encoding="utf-8").strip() if RESET.is_file() else ""
    seen = os.environ.get(SEEN)
    os.environ[SEEN] = line
    return seen is not None and seen != line


def number() -> str | None:
    """This engine's number, ``07`` for the entry file ``opendis_07.py``; None for another."""
    stem = Path(getattr(sys.modules.get("__main__"), "__file__", None) or "").stem
    return stem[len("opendis_"):] if stem.startswith("opendis_") else None


def setting(name: str, default: float) -> float:
    """Minutes from the environment, else from the deployment's secrets, else the default."""
    value = os.environ.get(name)
    if value is None and st.secrets.load_if_toml_exists():
        value = st.secrets.get(name)
    try:
        return default if value is None else float(value)
    except (TypeError, ValueError):
        return default


def idle_state(idle_s: float, idle_min: float, warning_min: float) -> str:
    """``off`` after idle_min minutes without a click, ``warn`` in the warning_min minutes before,
    else ``active``. No sign-off when idle_min is 0."""
    if idle_min <= 0:
        return "active"
    if idle_s >= idle_min * 60:
        return "off"
    return "warn" if idle_s >= (idle_min - warning_min) * 60 else "active"


def sessions() -> int | None:
    """The engine's sessions at this moment: the gauge ``active_sessions`` of the front page."""
    try:
        from streamlit.runtime import Runtime  # pylint: disable=import-outside-toplevel
        found = Runtime.instance().stats_mgr.get_stats(family_names=["active_sessions"])
        return int(sum(stat.value for stat in found.get("active_sessions", [])))
    except Exception:  # pylint: disable=broad-exception-caught
        return None  # no runtime, as in a test, or another Streamlit


def badge_html(nn: str, count: int | None) -> str:
    """The banner's right end: "engine 00 · users 1", the numbers in the code font. The users
    turn amber when another tab shares the engine. Without a count, the number alone."""
    users = "" if count is None else (
        '<span class="sep">·</span><span class="label">users</span>'
        f'<span class="num{" shared" if count > 1 else ""}" '
        f'title="Browser tabs open on this engine now">{count}</span>')
    return (f'<div class="opendis-users"><span class="label">engine</span>'
            f'<span class="num">{nn}</span>{users}</div>')


def _touch() -> None:
    st.session_state[TOUCHED] = time.time()


def _watch(name: str) -> None:
    """The banner's number and users, then the idle check: a warning, then the sign-off."""
    st.markdown(CSS + badge_html(name.rsplit("-", 1)[-1], sessions()), unsafe_allow_html=True)
    idle_min, warning_min = setting(*IDLE), setting(*WARNING)
    idle_s = time.time() - st.session_state.get(TOUCHED, time.time())
    state = idle_state(idle_s, idle_min, warning_min)
    if state == "warn":
        left = math.ceil((idle_min * 60 - idle_s) / 60)
        with st.container(key="opendis_idle"):
            st.warning(f"No activity for {idle_s // 60:.0f} minutes. This engine signs you off "
                       f"in {left} minute{'' if left == 1 else 's'}, so that another user can "
                       "have it.")
            st.button("Stay signed in", key="opendis_stay", type="primary", on_click=_touch)
        recent_ms = str(round(warning_min * 60 * 1000))
        st.html(STAY.replace("RECENT_MS", recent_ms), unsafe_allow_javascript=True)
    elif state == "off":
        url = f"{FRONT}/?signed_off={name}&idle={idle_min:g}"
        st.html(SIGN_OFF.format(url=json.dumps(url)), unsafe_allow_javascript=True)


def main() -> None:
    """Run the app's page as ``app.py`` does, after a reset if one is due, with the engine's
    number and users in the banner, and the idle check.

    The app's folder goes first on sys.path, for the app's imports, and on PYTHONPATH:
    Streamlit reloads a changed module only if it lies in the entry file's folder or on
    PYTHONPATH. So a push updates an engine in place, as it updates app.py.
    """
    if APP_DIR not in sys.path:
        sys.path.insert(0, APP_DIR)
    os.environ["PYTHONPATH"] = first(APP_DIR, os.environ.get("PYTHONPATH", ""))
    if reset_due():
        from common import job  # pylint: disable=import-outside-toplevel
        deleted, kept = job.clear_runs()
        st.cache_data.clear()
        st.cache_resource.clear()
        print(f"[OpenDiS engine] reset: deleted the results of {deleted} run(s), kept {kept}",
              file=sys.stderr, flush=True)
    nn = number()
    if nn is not None:
        _touch()  # every page run comes from a click, or from a push
        with st.container(key="opendis_engine"):  # out of the page's flow: no gap above the app
            st.html(TRACK, unsafe_allow_javascript=True)
            st.fragment(_watch, run_every=WATCH_S)(f"opendis-{nn}")
    import app  # pylint: disable=import-outside-toplevel
    app.main()
    if nn is not None:
        st.set_page_config(page_title=f"OpenDiS {nn}")  # the browser tab; the rest stays the app's
