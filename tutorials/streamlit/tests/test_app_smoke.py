"""tests.test_app_smoke: the whole page, headless, with Streamlit's AppTest (README.md, Tests)."""
import json
import os
import sys
import tempfile
import threading
import time

os.environ["OPENDIS_APP_RUNS"] = tempfile.mkdtemp(prefix="opendis_app_runs_")  # before the app

# pylint: disable=wrong-import-position
from pathlib import Path

from streamlit.testing.v1 import AppTest

from common import job, network3d, paths
from common.picker import STEP_LIST
from common.testing import check, summary

APP = str(Path(__file__).resolve().parents[1] / "app.py")


def page() -> AppTest:
    """The page, run once."""
    at = AppTest.from_file(APP, default_timeout=120)
    at.run()
    return at


def click(at: AppTest, label: str) -> None:
    """Click the button with that label, and run."""
    next(b for b in at.button if b.label == label).click().run()


def choose(at: AppTest, label: str, code: str = "PyDiS", variant: str | None = None) -> None:
    """Choose an example and its code in the sidebar's two menus, and its sub-example in step
    1's menu, and run."""
    at.selectbox(key="picker_example").select(label).run()
    menu = at.selectbox(key=f"picker_core/{label}")
    if menu.value != code:  # a menu of one code is disabled
        menu.select(code).run()
    if variant is not None:
        at.selectbox(key=f"picker_variant_menu/{label}").select(variant).run()


def shown(at: AppTest) -> list[dict]:
    """What the 3-D window's component draws: its arguments."""
    return [json.loads(c.proto.json_args) for c in at.get("component_instance")
            if c.proto.component_name.endswith("network3d")]


def main() -> int:  # pylint: disable=too-many-locals,too-many-statements
    """Run every check."""
    at = page()
    picker = [s for s in at.sidebar.selectbox if s.label == "Example"]
    codes = [(s.options, s.value, s.disabled) for s in at.sidebar.selectbox if s.label == "Core"]
    text = "\n".join(m.value for m in at.markdown)
    outputs = [c.value for c in at.code if "No run yet" in c.value]
    initial = shown(at)
    results = [
        check("the page loads without an exception", not at.exception, at.exception),
        check("the banner says OpenDiS, with the logo", "<h1>OpenDiS</h1>" in text
              and "data:image/png;base64," in text),
        check("in the sidebar, the Example menu: Frank-Read src, Binary junction, Uniaxial "
              "tension, Initial configurations", len(picker) == 1 and picker[0].options ==
              ["Frank-Read src", "Binary junction", "Uniaxial tension", "Initial configurations"],
              [s.options for s in at.selectbox]),
        check("below it, the Core menu of the Frank-Read source: PyDiS, ExaDiS",
              codes == [(["PyDiS", "ExaDiS"], "PyDiS", False)], codes),
        check("no menu in the page's main part", not list(at.main.selectbox),
              [s.label for s in at.main.selectbox]),
        check("below 🗑 clear all results, ExaDiS's status",
              any(c.value.startswith("ExaDiS is") for c in at.caption),
              [c.value for c in at.caption][-3:]),
        check("the right column is the simulation output, not a chat",
              "Simulation output" in text and "Chat" not in text and outputs, outputs),
        check("seven steps", [o for r in at.radio if r.label == "Steps" for o in r.options]
              == [s.label for s in STEP_LIST]),
        check("the 3-D window draws the initial configuration, in the app's first view",
              len(initial) == 1 and [f["label"] for f in initial[0]["frames"]] ==
              ["initial configuration"] and initial[0]["layout"]["scene"]["camera"]["eye"]["y"] < 0,
              initial[:1]),
    ]

    # every step of each example, sub-example and code, by next › (the binary junction last)
    for example, code, variant in (
            ("Frank-Read src", "PyDiS", None), ("Frank-Read src", "ExaDiS", None),
            ("Uniaxial tension", "PyDiS", None), ("Uniaxial tension", "ExaDiS", None),
            ("Initial configurations", "PyDiS", None), ("Initial configurations", "ExaDiS", None),
            ("Binary junction", "PyDiS", "FCC"), ("Binary junction", "ExaDiS", "FCC"),
            ("Binary junction", "ExaDiS", None), ("Binary junction", "PyDiS", None)):
        at = page()
        choose(at, example, code, variant)
        example = f"{example} {variant}" if variant else example
        failures = []
        labels = next(r for r in at.radio if r.label == "Steps").options
        for label in labels:
            if at.exception or at.warning:
                failures.append((label, [e.value for e in at.exception],
                                 [w.value for w in at.warning]))
            if label != labels[-1]:
                click(at, "next ›")
        titles = "\n".join(m.value for m in at.markdown)
        results += [check(f"{example} with {code}: every step draws, cleanly", not failures,
                          failures),
                    check(f"{example} with {code}: next › reaches the Run step",
                          any(b.label == "▶ Run" for b in at.button)),
                    check(f"{example} with {code}: the code panel is {code}'s",
                          f"{code} code" in titles, titles[-500:])]
    junction = [(s.options, s.disabled) for s in at.selectbox if s.label == "Core"]
    results.append(check("the binary junction has two codes, PyDiS and ExaDiS",
                         junction == [(["PyDiS", "ExaDiS"], False)], junction))

    # the binary junction's crystal structures, in a menu of step 1: FCC is a sub-example of its
    # own, at its step 1, and the choice stays while step 1 is not shown
    at = page()
    choose(at, "Binary junction")
    menu = at.selectbox(key="picker_variant_menu/Binary junction")
    offered = (menu.label, menu.options, menu.value)
    menu.select("FCC").run()
    titles = "\n".join(m.value for m in at.markdown)
    step = next(r for r in at.radio if r.label == "Steps").value
    junction_types = [s.options for s in at.selectbox if s.label == "Junction type"]
    click(at, "next ›")
    click(at, "‹ prev")
    kept = at.selectbox(key="picker_variant_menu/Binary junction").value
    results += [
        check("in step 1 of the binary junction, the Crystal structure menu: BCC, FCC",
              offered == ("Crystal structure", ["BCC", "FCC"], "BCC"), offered),
        check("FCC is a sub-example of its own: its title, at step 1, with its junction types",
              "Binary junction · FCC · PyDiS" in titles and step == STEP_LIST[0].label
              and junction_types == [["lomer", "glissile", "hirth", "colinear"]],
              (step, junction_types)),
        check("the choice stays when step 1 is left and shown again", kept == "FCC", kept),
    ]

    # a changed selection marks its step; reset brings the example's values back
    at = page()
    next(n for n in at.number_input if n.label.startswith("Box length")).set_value(2000.0).run()
    marked = next(r for r in at.radio if r.label == "Steps").format_func(STEP_LIST[0].label)
    for _ in STEP_LIST[1:]:
        click(at, "next ›")
    click(at, "↺ reset to the example's values")
    for _ in STEP_LIST[1:]:
        click(at, "‹ prev")
    box = next(n for n in at.number_input if n.label.startswith("Box length")).value
    results += [check("a changed box length marks step 1 with •", marked.endswith("•")),
                check("reset brings the box length back to 1000", box == 1000.0, box)]

    # the binary junction, 20 steps; its replay plays in the 3-D window's component, in the
    # browser (its behaviour is components/network3d/network3d.js's, checked with Node.js)
    at = page()
    choose(at, "Binary junction")
    for _ in STEP_LIST[1:]:
        click(at, "next ›")
    next(n for n in at.number_input if n.label.startswith("max_step")).set_value(20).run()
    written, write_live = [], network3d.write_live  # the snapshots that the 3-D window fetches
    network3d.write_live = lambda run_id, label, net: (written.append(label),
                                                       write_live(run_id, label, net))
    try:
        click(at, "▶ Run")
    finally:
        network3d.write_live = write_live
    left = sorted(p.name for p in network3d.LIVE.glob("*")) if network3d.LIVE.exists() else []
    run_error = list(at.exception)
    output = "\n".join(c.value for c in at.code)
    successes = [s.value for s in at.success]
    downloads = [d.label for d in at.get("download_button")]
    view = next(r for r in at.radio if r.label == "3-D view").value
    replay = shown(at)
    next(r for r in at.radio if r.label == "3-D view").set_value("history").run()
    history = [json.loads(c.proto.spec)["data"] for c in at.get("plotly_chart")]
    results += [
        check("the 3-D window replays the run's 3 snapshots, from the final one",
              len(replay) == 1 and [f["label"] for f in replay[0]["frames"]] ==
              ["step 0", "step 10", "final"] and replay[0]["start"] == "last",
              [[f["label"] for f in r["frames"]] for r in replay]),
        check("during the run, its newest snapshots are written for the 3-D window to fetch",
              written and set(written) <= {"step 0", "step 10"}, written),
        check("the newest snapshot's file is removed when the run ends", not left, left),
        check("the run ends without an exception", not run_error, run_error),
        check("the output column shows PyDiS's lines",
              "init_two_disl_lines: z0 = 1.000000" in output and "step = 10" in output,
              output[-300:]),
        check("the history: the nodes and the density of the run's 2 snapshots",
              len(history) == 1 and [t["x"] for t in history[0]] == [[0, 10], [0, 10]]
              and all(y > 0 for t in history[0] for y in t["y"]),
              [(t.get("x"), t.get("y")) for data in history for t in data]),
        check("the output ends with the app's line on how the job ended",
              "[app] done after" in output, output[-200:]),
        check("the Run step says Done", any(s.startswith("Done after") for s in successes),
              successes),
        check("the view switched to the simulation", view == "simulation", view),
        check("the final network, the whole run and the log are offered for download",
              {"⬇ final network (JSON)", "⬇ the whole run (zip)", "⬇ log.txt"} <= set(downloads),
              downloads),
        check("no job is left running", not job.running()),
    ]

    # the Python console, on the run's objects; ↺ restart ends its session
    next(r for r in at.radio if r.label == "Code panel").set_value("Python console").run()
    at.text_area(key="console_source").input("G.num_nodes(), sim.max_step")  # sent on submit
    next(b for b in at.button if b.label == "▶ execute").click().run()
    final = next(Path(os.environ["OPENDIS_APP_RUNS"]).glob("binary_junction/*/output/*final.json"))
    nodes = len(json.loads(final.read_text())["nodes"]["tags"])
    transcript = "\n".join(c.value for c in at.code)
    session = at.session_state["console"]
    click(at, "↺ restart")
    restarted = [c.value for c in at.code if ">>>" in c.value]
    results += [
        check("the console runs a command on the run's objects",
              f">>> G.num_nodes(), sim.max_step\n({nodes}, 20)" in transcript, transcript[-300:]),
        check("↺ restart ends the session and clears its commands",
              not session.alive() and not restarted, restarted),
    ]

    # 🗑 clear all results, after the run: a confirmation in place, which Cancel closes, and
    # Delete too
    click(at, "🗑 clear all results")
    asked = [b.label for b in at.button if b.key in ("clear_results", "clear_cancel")]
    click(at, "Cancel")
    kept = [p for p in Path(os.environ["OPENDIS_APP_RUNS"]).glob("*/*") if p.is_dir()]
    cancelled = [b.label for b in at.button if b.key in ("clear_results", "clear_cancel")]
    click(at, "🗑 clear all results")
    click(at, "Delete")
    left = [p for p in Path(os.environ["OPENDIS_APP_RUNS"]).glob("*/*") if p.is_dir()]
    toasts = [t.value for t in at.toast]
    output = "\n".join(c.value for c in at.code)
    after = [b.label for b in at.button if b.key in ("clear_results", "clear_cancel", "clear_ask")]
    results += [
        check("🗑 clear all results asks first: Delete, or Cancel",
              asked == ["Delete", "Cancel"], asked),
        check("Cancel closes the question and deletes nothing", not cancelled and len(kept) == 1,
              (cancelled, kept)),
        check("Delete deletes the run's folder, and closes the question",
              not left and after == ["🗑 clear all results"], (left, after)),
        check("it says so, and the output column shows no run",
              toasts == ["Deleted the results of 1 run(s)."] and "No run yet" in output,
              (toasts, output[-200:])),
    ]

    # the Frank-Read source with ExaDiS, 20 steps: its snapshots, as PyDiS's, play
    at = page()
    choose(at, "Frank-Read src", "ExaDiS")
    for _ in STEP_LIST[1:]:
        click(at, "next ›")
    next(n for n in at.number_input if n.label.startswith("max_step")).set_value(20).run()
    click(at, "▶ Run")
    output = "\n".join(c.value for c in at.code)
    replay = shown(at)
    successes = [s.value for s in at.success]
    results += [
        check("ExaDiS: the run ends without an exception, and says Done", not at.exception
              and any(s.startswith("Done after") for s in successes),
              ([e.value for e in at.exception], successes)),
        check("ExaDiS: the output column shows ExaDiS's lines",
              "init_frank_read_src_loop: length = 125.000000" in output and "step = 20" in output
              and "RUN TIME" in output, output[-300:]),
        check("ExaDiS: the 3-D window replays its snapshots and the final network",
              len(replay) == 1 and [f["label"] for f in replay[0]["frames"]] ==
              ["step 0", "step 10", "step 20", "final"],
              [[f["label"] for f in r["frames"]] for r in replay]),
        check("no job is left running after ExaDiS's", not job.running()),
    ]

    # the initial configurations: each of the tutorial's networks draws in step 1, and Run
    # builds one with ExaDiS, which the run view shows
    at = page()
    choose(at, "Initial configurations", "ExaDiS")
    drawn = {}
    for label in ("simple dislocation source", "ParaDiS file", "infinite lines",
                  "prismatic loops"):
        next(s for s in at.selectbox if s.label == "Network").select(label).run()
        frames = shown(at)
        drawn[label] = (not at.exception and len(frames) == 1
                        and [f["label"] for f in frames[0]["frames"]] == ["initial configuration"])
    two = [o.removesuffix(" •") for r in at.radio if r.label == "Steps" for o in r.options]
    click(at, "next ›")
    click(at, "▶ Run")
    output = "\n".join(c.value for c in at.code)
    replay = shown(at)
    results += [
        check("initial configurations: two steps", two == ["1. Initial configuration", "2. Run"],
              two),
        check("initial configurations: each of the tutorial's networks draws in step 1",
              all(drawn.values()), drawn),
        check("initial configurations: Run builds the prismatic loops with ExaDiS",
              "Create initial configuration of prismatic loops" in output
              and any(s.value.startswith("Done after") for s in at.success), output[-300:]),
        check("initial configurations: the run view shows the network that the run built",
              len(replay) == 1 and [f["label"] for f in replay[0]["frames"]] == ["final"],
              [[f["label"] for f in r["frames"]] for r in replay]),
    ]

    # until ExaDiS is built for this Python, its examples cannot run; while the app builds it,
    # a bar shows how far it is (a stand-in for the build: a thread of its name and progress)
    built, paths.EXADIS_MODULE = paths.EXADIS_MODULE, "pyexadis.not-built.so"
    done = threading.Event()
    builder = threading.Thread(target=done.wait, name=paths.BUILD_THREAD, daemon=True)
    builder.exadis = {"phase": "compiling", "percent": 45, "file": "exadis_pybind.cpp",
                      "start": time.monotonic() - 90}
    try:
        at = page()
        choose(at, "Frank-Read src", "ExaDiS")
        for _ in STEP_LIST[1:]:
            click(at, "next ›")
        run = next(b for b in at.button if b.label == "▶ Run")
        warnings = [w.value for w in at.warning]
        bar = [c.value for c in at.caption if c.value.startswith("ExaDiS is")]
        builder.start()
        at.run()
        building = [(p.value, p.proto.text) for p in at.get("progress")]
        run_building = next(b for b in at.button if b.label == "▶ Run").disabled
    finally:
        done.set()
        paths.EXADIS_MODULE = built
    results += [
        check("without ExaDiS, its Run is off, and the Run step says why",
              run.disabled and any("ExaDiS is not built" in w for w in warnings),
              (run.disabled, warnings)),
        check("... and so does the status below 🗑 clear all results",
              any("not built" in c for c in bar), bar),
        check("while ExaDiS builds, a bar shows its phase, percent, file and time, and Run is off",
              building == [(45, "Building ExaDiS: compiling, 45 %, exadis_pybind.cpp · 1.5 min "
                                "so far")] and run_building, (building, run_building)),
    ]

    # a changed source file makes Streamlit load the app's modules again, with new classes:
    # the session's spec, of the old ones, must still work (last, as job.running() is the old
    # module's then)
    at = page()
    next(n for n in at.number_input if n.label == "Arm length / Lbox").set_value(0.2).run()
    for name in [m for m in sys.modules if m == "common" or m.startswith("common.")]:
        del sys.modules[name]
    at.run()
    arms = [n.value for n in at.number_input if n.label == "Arm length / Lbox"]
    results += [
        check("after the app's modules are loaded again, the page draws, with the changed value",
              not at.exception and arms == [0.2], ([e.value for e in at.exception], arms)),
    ]
    return summary(results)


if __name__ == "__main__":
    sys.exit(main())
