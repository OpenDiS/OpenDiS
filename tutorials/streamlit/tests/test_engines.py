"""tests.test_engines: the engines' entry files and helper, the front page, and manage.py.

The engines' gauges come from a stand-in server on this computer (README.md, Tests).
"""
import http.server
import os
import sys
import tempfile
import threading
import time

os.environ["OPENDIS_APP_RUNS"] = tempfile.mkdtemp(prefix="opendis_app_runs_")  # before the app

# pylint: disable=wrong-import-position
from pathlib import Path

from streamlit.testing.v1 import AppTest

from common import front_page
from common.testing import check, summary
from engines import engine, manage

ENGINES = Path(engine.__file__).parent
APP = str(ENGINES.parent / "app.py")
ENTRY = ('"""{name}.streamlit.app: one engine of the OpenDiS tutorial (README.md)."""\n'
         "from engine import main\nmain()\n")
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
    first = [cell.split(">", 2)[2].split("<", 1)[0] for cell in rows[0].split("<td")[1:]]
    return [check("table: 10 rows of 5 cells", len(rows) == 10
                  and all(row.count("<td") == 5 for row in rows), [r.count("<td") for r in rows]),
            check("table: the first row", first == ["opendis-00", "opendis-10", "opendis-20",
                                                    "opendis-30", "opendis-40"], first),
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
    names = (front_page.SWITCH, "PYTHONPATH", engine.SEEN)
    saved, metrics, path = {k: os.environ.get(k) for k in names}, front_page.METRICS, sys.path[:]
    news, tmp = front_page.NEWS, tempfile.TemporaryDirectory()
    front_page.METRICS = f"http://127.0.0.1:{_free_port()}/{{name}}"
    front_page.NEWS = Path(tmp.name) / "announcements.md"
    front_page.NEWS.write_text("**A test announcement**, with [a link](https://example.org).\n")
    try:
        os.environ[front_page.SWITCH] = "true"
        front = AppTest.from_file(APP, default_timeout=60)
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
            check("pages: the front page has no example menu", not front.selectbox),
            check("pages: the app, when not switched on", not app.exception
                  and app.selectbox(key="picker_example").value is not None, app.exception),
            check("pages: the app from opendis_00.py", not one.exception
                  and one.selectbox(key="picker_example").value is not None, one.exception)]


def main() -> int:
    """Run every check; the exit code is 0 when all pass."""
    return summary(test_entry_files() + test_paths() + test_reset() + test_parse() + test_read()
                   + test_table() + test_news() + test_switch() + test_manage() + test_pages())


if __name__ == "__main__":
    raise SystemExit(main())
