"""engines/manage.py: the 50 engines from a terminal, within the Cloud's terms (README.md).

Usage, with the app's Python, from this folder:

    python manage.py status        # each engine's users, as on the front page
    python manage.py open [--all]  # open in the browser each engine that does not answer (or all)
    python manage.py links         # each engine's name and the GitHub link of its entry file
    python manage.py reset         # a new line in RESET: a push of it resets every engine
"""
import argparse
import datetime
import sys
import time
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import front_page  # noqa: E402  pylint: disable=wrong-import-position

GITHUB = "https://github.com/OpenDiS/OpenDiS/blob/main/tutorials/streamlit/engines/{file}"
OPEN_GAP_S = 0.5  # between two tabs, so that the browser keeps up


def status_text(engines: list[front_page.Engine]) -> str:
    """The front page's table as text, then a count of the engines in each state."""
    lines = []
    for row in range(front_page.ROWS):
        lines.append("   ".join(f"{e.name}  {e.label:<12}" for e in engines[row::front_page.ROWS]))
    counts = {state: sum(e.state == state for e in engines)
              for state in ("free", "busy", "asleep", "missing")}
    lines.append(f"\n{counts['free']} free, {counts['busy']} in use, {counts['asleep']} asleep "
                 f"(or starting, or failing), {counts['missing']} not deployed")
    return "\n".join(lines)


def links_text() -> str:
    """For each engine: its name, which is its App URL, and the GitHub link of its entry file."""
    return "\n".join(f"{front_page.NAME.format(n=n)}  {GITHUB.format(file=f'opendis_{n:02d}.py')}"
                     for n in range(front_page.COUNT))


def main() -> int:
    """Run one command; see the module's docstring."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("command", choices=("status", "open", "links", "reset"))
    parser.add_argument("--all", action="store_true", help="open: every engine, not only those "
                        "that do not answer")
    args = parser.parse_args()
    if args.command == "links":
        print(links_text())
    elif args.command == "reset":
        line = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        (HERE / "RESET").write_text(line + "\n", encoding="utf-8")
        print(f"RESET now reads: {line}\nCommit tutorials/streamlit/engines/RESET and push it to "
              "GitHub's main. Each engine then clears its results at its next page run.")
    else:
        engines = front_page.read_all()
        if args.command == "status":
            print(status_text(engines))
        else:
            chosen = [e for e in engines if args.all or e.state == "asleep"]
            for engine in chosen:
                webbrowser.open_new_tab(engine.url())
                time.sleep(OPEN_GAP_S)
            print(f"Opened {len(chosen)} engine(s): {' '.join(e.name for e in chosen)}"
                  if chosen else "Every engine answers: none to open.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
