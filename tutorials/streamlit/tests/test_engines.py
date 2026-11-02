"""tests.test_engines: the engines' entry files and helper, the front page, and manage.py.

The engines' gauges come from a stand-in server on this computer (README.md, Tests).
"""
import http.server
import json
import os
import re
import sys
import tempfile
import threading
import time
import types

os.environ["OPENDIS_APP_RUNS"] = tempfile.mkdtemp(prefix="opendis_app_runs_")  # before the app

# pylint: disable=wrong-import-position
from pathlib import Path

from streamlit.testing.v1 import AppTest

from common import front_page, style
from common.testing import check, summary
from engines import engine, manage

ENGINES = Path(engine.__file__).parent
APP = str(ENGINES.parent / "app.py")
ENTRY = ('"""{name}.streamlit.app: one engine of the OpenDiS tutorial (README.md)."""\n'
         "from engine import main\nmain()\n")
# A page run of the idle check alone, its last click idle_s seconds ago
WATCH = """
import sys, time
sys.path.insert(0, {engines!r})
import streamlit as st
import engine
st.session_state.setdefault(engine.TOUCHED, time.time() - {idle_s})
engine._watch("opendis-07")
"""
# The gauges as the Cloud's opendis.streamlit.app published them on 2026-11-01, shortened
METRICS = """# TYPE cache_memory_bytes gauge
cache_memory_bytes{cache_type="st_session_state",cache=""} 27
# TYPE session_events counter
session_events_total{type="connect"} 6
# TYPE active_sessions gauge
# HELP active_sessions Current number of active sessions.
active_sessions 2
# EOF
"""


class Stub(http.server.BaseHTTPRequestHandler):
    """Engine n's gauges at /opendis-NN: 3 sessions for 00, none for 01, 404 for 02, 500 for 03,
    a page without gauges for 04, a wait longer than the timeout for 05, 1 session for others."""

    def do_GET(self):  # pylint: disable=invalid-name
        """Answer as engine n would."""
        n = int(self.path.rsplit("-", 1)[-1])
        if n == 5:
            time.sleep(1.5)
        body = {0: METRICS.replace("sessions 2", "sessions 3"),
                1: METRICS.replace("sessions 2", "sessions 0"),
                4: "<html>Zzzz</html>"}.get(n, METRICS.replace("sessions 2", "sessions 1"))
        try:
            if n in (2, 3):
                self.send_error(404 if n == 2 else 500)
                return
            self.send_response(200)
            self.end_headers()
            self.wfile.write(body.encode())
        except OSError:  # engine 5's reader gave up
            pass

    def log_message(self, *args):  # pylint: disable=arguments-differ
        """Quiet."""


def test_entry_files() -> list[bool]:
    """50 entry files, opendis_00.py to opendis_49.py, each the same three lines."""
    names = sorted(p.name for p in ENGINES.glob("opendis_*.py"))
    wanted = [f"opendis_{n:02d}.py" for n in range(front_page.COUNT)]
    wrong = [name for name in names
             if (ENGINES / name).read_text() != ENTRY.format(name=name[:-3].replace("_", "-"))]
    return [check("entry files: opendis_00.py to opendis_49.py", names == wanted, names[:3]),
            check("entry files: each the same three lines", not wrong, wrong[:3])]


def test_paths() -> list[bool]:
    """The app's folder first on a search path, once."""
    sep, app = os.pathsep, engine.APP_DIR
    return [check("first: on an empty path", engine.first(app, "") == app),
            check("first: before the others",
                  engine.first(app, f"/a{sep}/b") == f"{app}{sep}/a{sep}/b"),
            check("first: only once",
                  engine.first(app, f"/a{sep}{app}{sep}") == f"{app}{sep}/a"),
            check("APP_DIR: the folder of app.py", (Path(app) / "app.py").is_file(), app)]


def test_reset() -> list[bool]:
    """A new line of RESET is due once; the first page run of a process sees no reset."""
    saved, seen = engine.RESET, os.environ.pop(engine.SEEN, None)
    with tempfile.TemporaryDirectory() as tmp:
        engine.RESET = Path(tmp) / "RESET"
        engine.RESET.write_text("one\n")
        due = [engine.reset_due(), engine.reset_due()]
        engine.RESET.write_text("two\n")
        due += [engine.reset_due(), engine.reset_due()]
        engine.RESET.unlink()
        due.append(engine.reset_due())
    engine.RESET = saved
    os.environ.pop(engine.SEEN, None)
    if seen is not None:
        os.environ[engine.SEEN] = seen
    return [check("reset: due once per new line, not at the first run",
                  due == [False, False, True, False, True], due)]


def test_code() -> list[bool]:
    """A changed code file is seen once, and a run's script not at all; the app's modules are
    dropped, the page's own and the others kept."""
    seen = os.environ.pop(engine.CODE, None)
    with tempfile.TemporaryDirectory() as tmp:
        app = Path(tmp)
        for name in ("app.py", "common/spec.py", "engines/engine.py", "runs/a/b/script.py"):
            (app / name).parent.mkdir(parents=True, exist_ok=True)
            (app / name).write_text("x = 1\n")
        changed = [engine.code_changed(tmp), engine.code_changed(tmp)]
        later = time.time() + 5
        os.utime(app / "common" / "spec.py", (later, later))
        changed += [engine.code_changed(tmp), engine.code_changed(tmp)]
        os.utime(app / "runs" / "a" / "b" / "script.py", (later + 5, later + 5))
        changed.append(engine.code_changed(tmp))
        modules = {"__main__": types.SimpleNamespace(__file__=f"{tmp}/engines/opendis_07.py"),
                   "app": types.SimpleNamespace(__file__=f"{tmp}/app.py"),
                   "common": types.SimpleNamespace(__file__=f"{tmp}/common/__init__.py",
                                                   __path__=[f"{tmp}/common"]),
                   "common.spec": types.SimpleNamespace(__file__=f"{tmp}/common/spec.py"),
                   "engines": types.SimpleNamespace(__file__=None, __path__=[f"{tmp}/engines"]),
                   "json": json, "gone": None}
        dropped = engine.forget_app_modules(modules, tmp)
    os.environ.pop(engine.CODE, None)
    if seen is not None:
        os.environ[engine.CODE] = seen
    return [check("code: a changed file seen once, a run's script not, nothing at the first run",
                  changed == [False, False, True, False, False], changed),
            check("code: the app's modules dropped, the page's own and the others kept",
                  sorted(dropped) == ["app", "common", "common.spec", "engines"]
                  and sorted(modules) == ["__main__", "gone", "json"], (dropped, list(modules)))]


def test_idle() -> list[bool]:
    """The idle check's thresholds and settings, the engine's number and its count of users;
    then, with AppTest, the warning, Stay, and the sign-off."""
    states = [engine.idle_state(minutes * 60, 60, 5) for minutes in (0, 54.9, 55, 59.9, 60, 600)]
    name, saved = engine.IDLE[0], os.environ.pop(engine.IDLE[0], None)
    minutes = []
    for value in ("2.5", "soon", None):
        if value is not None:
            os.environ[name] = value
        minutes.append(engine.setting(*engine.IDLE))
        os.environ.pop(name, None)
    if saved is not None:
        os.environ[name] = saved
    real, fake = sys.modules["__main__"], types.ModuleType("__main__")
    sys.modules["__main__"] = fake  # number() reads the entry file of the page run
    fake.__file__ = str(ENGINES / "opendis_07.py")
    numbers = [engine.number()]
    fake.__file__ = APP
    numbers.append(engine.number())
    sys.modules["__main__"] = real
    runs = {idle_s: AppTest.from_string(WATCH.format(engines=str(ENGINES), idle_s=idle_s),
                                        default_timeout=30).run()
            for idle_s in (60, 56 * 60, 61 * 60)}
    warn = runs[56 * 60]
    before = [w.value for w in warn.warning]
    warn_html = "".join(str(e.proto) for e in warn.main if getattr(e, "type", "") == "html")
    stay = warn.button(key="opendis_stay").label if warn.button else None
    if stay:
        warn.button(key="opendis_stay").click().run()
    scripts = {s: "".join(str(e.proto) for e in at.main if getattr(e, "type", "") == "html")
               for s, at in runs.items()}
    text = runs[56 * 60].warning
    return [check("idle: warn from 55 minutes, off from 60",
                  states == ["active", "active", "warn", "warn", "off", "off"], states),
            check("idle: never when its minutes are 0", engine.idle_state(1e6, 0, 5) == "active"),
            check("idle: minutes from the environment, else the default",
                  minutes == [2.5, 60.0, 60.0], minutes),
            check("number: 07 from opendis_07.py, none from app.py", numbers == ["07", None],
                  numbers),
            check("badge: engine 07 alone, with 1 user, with 3 (shared)",
                  'class="num">07<' in engine.badge_html("07", None)
                  and ">users<" not in engine.badge_html("07", None)
                  and '>users</span><span class="num" ' in engine.badge_html("07", 1)
                  and ">1</span>" in engine.badge_html("07", 1)
                  and '<span class="num shared"' in engine.badge_html("07", 3)
                  and ">3</span>" in engine.badge_html("07", 3)),
            check("idle: nothing in the first minutes", not runs[60].warning
                  and not scripts[60], scripts[60][:80]),
            check("idle: the warning at 56 minutes, and Stay", stay == "Stay signed in"
                  and len(before) == 1 and before[0].startswith("No activity for 56 minutes. "
                                                                "This engine signs you off in 4 "
                                                                "minutes"), before),
            check("idle: with the warning, Stay by itself while the page is in use, no sign-off",
                  not runs[56 * 60].exception and "opendis_stay" in warn_html
                  and str(round(engine.WARNING[1] * 60 * 1000)) in warn_html
                  and "signed_off" not in warn_html,
                  warn_html[:120]),
            check("idle: Stay ends the warning", not text, [w.value for w in text]),
            check("idle: the sign-off, to the front page", "?signed_off=opendis-07&idle=60"
                  in scripts[61 * 60] and "link.click()" in scripts[61 * 60],
                  scripts[61 * 60][:120])]


def test_parse() -> list[bool]:
    """The number of sessions in the gauges, and None without it."""
    return [check("parse: the Cloud's gauges", front_page.parse(METRICS) == 2),
            check("parse: none",
                  front_page.parse(METRICS.replace("sessions 2", "sessions 0")) == 0),
            check("parse: a float", front_page.parse("active_sessions 4.0\n") == 4),
            check("parse: no gauge", front_page.parse("<html>Zzzz</html>") is None),
            check("parse: another family", front_page.parse("active_sessions_total 3\n") is None)]


def test_read() -> list[bool]:
    """Each state from the stand-in server, all 50 engines at once, within the timeout."""
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Stub)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    metrics = f"http://127.0.0.1:{server.server_address[1]}/{{name}}"
    try:
        start = time.monotonic()
        engines = front_page.read_all(metrics, timeout=1.0)
        took = time.monotonic() - start
        unused = front_page.read(0, f"http://127.0.0.1:{_free_port()}/{{name}}", timeout=1.0)
    finally:
        server.shutdown()
    states = [(e.n, e.users, e.state, e.label) for e in engines[:7]]
    return [check("read: each state", states == [
        (0, 3, "busy", "3 users"), (1, 0, "free", "free"), (2, None, "missing", "not deployed"),
        (3, None, "asleep", "asleep"), (4, None, "asleep", "asleep"),
        (5, None, "asleep", "asleep"), (6, 1, "busy", "1 user")], states),
        check("read: 50 engines, in order", [e.n for e in engines] == list(range(50))),
        check("read: all at once, within the timeout", took < 3.0, f"{took:.1f} s"),
        check("read: no server", unused.state == "asleep", unused)]


def _free_port() -> int:
    """A port where nothing listens."""
    server = http.server.HTTPServer(("127.0.0.1", 0), Stub)
    port = server.server_address[1]
    server.server_close()
    return port


def test_table() -> list[bool]:
    """5 columns of 10: opendis-00 to opendis-09 in the first; each cell's link and state."""
    engines = [front_page.Engine(n, users=(0, 1, 2, None)[n % 4], missing=n == 49)
               for n in range(50)]
    text = front_page.table_html(engines, page="http://localhost:87{n:02d}")
    rows = text.split("<tr>")[1:]
    first = [re.sub(r"<[^>]+>", "", name) for name in re.findall(r"<a [^>]*>(.*?)</a>", rows[0])]
    return [check("table: 10 rows of 5 cells", len(rows) == 10
                  and all(row.count("<td") == 5 for row in rows), [r.count("<td") for r in rows]),
            check("table: the first row", first == ["opendis-00", "opendis-10", "opendis-20",
                                                    "opendis-30", "opendis-40"], first),
            check("table: the names' prefix apart, which a narrow box hides",
                  '<span class="prefix">opendis-</span>07</a>' in text),
            check("table: a link that opens a tab", 'href="http://localhost:8707" target="_blank"'
                  in text),
            check("table: the states", all(part in text for part in (
                '<td class="free">', '<span class="state">1 user</span>',
                '<span class="state">2 users</span>', '<td class="asleep">',
                '<td class="missing">', '<span class="state">not deployed</span>')))]


def test_news() -> list[bool]:
    """The announcements: the file's text, stripped; nothing when it is empty or missing."""
    saved, texts = front_page.NEWS, []
    with tempfile.TemporaryDirectory() as tmp:
        front_page.NEWS = Path(tmp) / "announcements.md"
        texts.append(front_page.news())
        front_page.NEWS.write_text("  \n\n")
        texts.append(front_page.news())
        front_page.NEWS.write_text("\n**Tutorial** on [MMM](https://example.org).\n")
        texts.append(front_page.news())
    front_page.NEWS = saved
    return [check("news: none, none, the text",
                  texts == ["", "", "**Tutorial** on [MMM](https://example.org)."], texts)]


def test_switch() -> list[bool]:
    """OPENDIS_FRONT_PAGE from the environment; no secrets file here."""
    saved, results = os.environ.pop(front_page.SWITCH, None), []
    for value, wanted in (("true", True), ("1", True), ("YES", True), ("false", False),
                          ("", False), (None, False)):
        if value is None:
            os.environ.pop(front_page.SWITCH, None)
        else:
            os.environ[front_page.SWITCH] = value
        results.append(front_page.wanted() == wanted)
    os.environ.pop(front_page.SWITCH, None)
    if saved is not None:
        os.environ[front_page.SWITCH] = saved
    return [check("switch: true, 1 or yes; else the app", all(results), results)]


def test_manage() -> list[bool]:
    """manage.py's table and links."""
    engines = [front_page.Engine(n, users=1 if n == 7 else 0) for n in range(50)]
    lines = manage.status_text(engines).splitlines()
    links = manage.links_text().splitlines()
    return [check("manage: status, 10 rows", lines[0].startswith("opendis-00  free")
                  and "opendis-17" in lines[7] and len(lines) == 12, lines[:2]),
            check("manage: status, the counts",
                  lines[-1].startswith("49 free, 1 in use, 0 asleep"), lines[-1]),
            check("manage: links", len(links) == 50 and links[7] == (
                "opendis-07  https://github.com/OpenDiS/OpenDiS/blob/main/tutorials/streamlit/"
                "engines/opendis_07.py"), links[7])]


def test_pages() -> list[bool]:
    """With Streamlit's AppTest: the front page when switched on, else the app; and the app
    from an engine's entry file. No engine answers here."""
    names = (front_page.SWITCH, "PYTHONPATH", engine.SEEN, engine.CODE)
    saved, metrics, path = {k: os.environ.get(k) for k in names}, front_page.METRICS, sys.path[:]
    title = style.BANNER_TITLE
    news, tmp = front_page.NEWS, tempfile.TemporaryDirectory()
    front_page.METRICS = f"http://127.0.0.1:{_free_port()}/{{name}}"
    front_page.NEWS = Path(tmp.name) / "announcements.md"
    front_page.NEWS.write_text("**A test announcement**, with [a link](https://example.org).\n")
    try:
        os.environ[front_page.SWITCH] = "true"
        front = AppTest.from_file(APP, default_timeout=60)
        front.query_params["signed_off"], front.query_params["idle"] = "opendis-07", "60"
        front.run()
        text = "".join(m.value for m in front.markdown)
        del os.environ[front_page.SWITCH]
        app = AppTest.from_file(APP, default_timeout=120)
        app.run()
        sys.path.insert(0, str(ENGINES))  # as Streamlit does for an entry file
        one = AppTest.from_file(str(ENGINES / "opendis_00.py"), default_timeout=120)
        one.run()
    finally:
        front_page.METRICS, front_page.NEWS, sys.path[:] = metrics, news, path
        banner = "".join(m.value for m in one.markdown) if one.markdown else ""
        tracked = any("opendisTrack" in str(e.proto) for e in one.main
                      if getattr(e, "type", "") == "html")
        style.BANNER_TITLE = title
        tmp.cleanup()
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    return [check("pages: the front page, without errors", not front.exception, front.exception),
            check("pages: the front page's table, every engine asleep",
                  '<table class="opendis-engines">' in text and text.count('<td class="asleep">')
                  == 50, text[:200]),
            check("pages: the front page's announcement",
                  any("A test announcement" in m.value for m in front.markdown)),
            check("pages: the front page's policy: shared, released, signed off",
                  "Every engine is shared with anyone who opens it" in text
                  and "clear all results" in text
                  and "Closing the browser tab releases the engine" in text
                  and "After 60 minutes of idling, the engine signs you off" in text),
            check("pages: the note after a sign-off", bool(front.info) and front.info[0].value
                  .startswith("opendis-07 signed you off after 60 minutes"),
                  [i.value for i in front.info]),
            check("pages: the front page has no example menu", not front.selectbox),
            check("pages: the app, when not switched on", not app.exception
                  and app.selectbox(key="picker_example").value is not None, app.exception),
            check("pages: the app from opendis_00.py", not one.exception
                  and one.selectbox(key="picker_example").value is not None, one.exception),
            check("pages: the engine's banner, OpenDiS, with engine 00 at its right end",
                  "<h1>OpenDiS</h1>" in banner and '>engine</span><span class="num">00<' in banner,
                  banner[-300:]),
            check("pages: the engine notes use that stays in the browser", tracked)]


def main() -> int:
    """Run every check; the exit code is 0 when all pass."""
    return summary(test_entry_files() + test_paths() + test_reset() + test_code() + test_idle()
                   + test_parse()
                   + test_read()
                   + test_table() + test_news() + test_switch() + test_manage() + test_pages())


if __name__ == "__main__":
    raise SystemExit(main())
