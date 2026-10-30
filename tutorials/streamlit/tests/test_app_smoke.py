"""tests.test_app_smoke: the whole page, headless, with Streamlit's AppTest (README.md, Tests)."""
import json
import os
import sys
import tempfile

os.environ["OPENDIS_APP_RUNS"] = tempfile.mkdtemp(prefix="opendis_app_runs_")  # before the app

# pylint: disable=wrong-import-position
from pathlib import Path

from streamlit.testing.v1 import AppTest

from common import job
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


def shown(at: AppTest) -> list[dict]:
    """What the 3-D window's component draws: its arguments."""
    return [json.loads(c.proto.json_args) for c in at.get("component_instance")
            if c.proto.component_name.endswith("network3d")]


def main() -> int:  # pylint: disable=too-many-locals,too-many-statements
    """Run every check."""
    at = page()
    picker = [s for s in at.selectbox if s.label == "Example"]
    text = "\n".join(m.value for m in at.markdown)
    outputs = [c.value for c in at.code if "No run yet" in c.value]
    initial = shown(at)
    results = [
        check("the page loads without an exception", not at.exception, at.exception),
        check("the banner says OpenDiS, with the logo", "<h1>OpenDiS</h1>" in text
              and "data:image/png;base64," in text),
        check("one drop-down menu: Frank-Read src, Binary junction",
              len(picker) == 1 and picker[0].options == ["Frank-Read src", "Binary junction"],
              [s.options for s in at.selectbox]),
        check("the right column is the simulation output, not a chat",
              "Simulation output" in text and "Chat" not in text and outputs, outputs),
        check("seven steps", [o for r in at.radio if r.label == "Steps" for o in r.options]
              == [s.label for s in STEP_LIST]),
        check("the 3-D window draws the initial configuration, in the app's first view",
              len(initial) == 1 and [f["label"] for f in initial[0]["frames"]] ==
              ["initial configuration"] and initial[0]["layout"]["scene"]["camera"]["eye"]["y"] < 0,
              initial[:1]),
    ]

    # every step of each example, by next ›
    for example in ("Frank-Read src", "Binary junction"):
        at = page()
        at.selectbox(key="picker_example").select(example).run()
        failures = []
        for step in STEP_LIST:
            if at.exception or at.warning:
                failures.append((step.label, [e.value for e in at.exception],
                                 [w.value for w in at.warning]))
            if step != STEP_LIST[-1]:
                click(at, "next ›")
        results += [check(f"{example}: every step draws, cleanly", not failures, failures),
                    check(f"{example}: next › reaches the Run step",
                          any(b.label == "▶ Run" for b in at.button))]

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
    at.selectbox(key="picker_example").select("Binary junction").run()
    for _ in STEP_LIST[1:]:
        click(at, "next ›")
    next(n for n in at.number_input if n.label.startswith("max_step")).set_value(20).run()
    click(at, "▶ Run")
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

    # 🗑 clear all results, after the run
    at.button(key="clear_results").click().run()
    left = [p for p in Path(os.environ["OPENDIS_APP_RUNS"]).glob("*/*") if p.is_dir()]
    toasts = [t.value for t in at.toast]
    output = "\n".join(c.value for c in at.code)
    results += [
        check("🗑 clear all results deletes the run's folder", not left, left),
        check("it says so, and the output column shows no run",
              toasts == ["Deleted the results of 1 run(s)."] and "No run yet" in output,
              (toasts, output[-200:])),
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
