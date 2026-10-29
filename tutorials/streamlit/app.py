"""OpenDiS: the Streamlit app (entry point).

Run it from this folder with ``./launch_app.sh``, or with
``python -m streamlit run app.py``. The page has the pinned banner, with the example picker in
it; the sidebar (the steps, prev / next); the chosen step's selections in the 2nd left column;
the 3-D window and the code panel in the centre; and the simulation's screen output on the
right. Run launches the PyDiS job (``common/job.py``) as the last thing in a script run, so that
Stop, or any other click, interrupts it and kills the job.
"""
import functools
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from common import (catalog, code_panel, codegen, job, network, network3d, output_panel, picker,
                    results, steps, style, views)
from common.catalog import Example
from common.paths import ASSETS, fetch_graph
from common.picker import RUN_STEP, Step
from common.spec import Spec

LOGO = ASSETS / "opendis_logo.png"
FAVICON = ASSETS / "opendis_favicon.png"  # the website's favicon
SIDEBAR_WIDTH_PX = 300  # initial width; Streamlit clamps it to 200-600
VIEWS = ("initial configuration", "simulation", "history")
FRAME_MS = 300  # the replay's time per snapshot
# The 3-D window's height, whatever it shows (the figure and the replay's controls, a run's newest
# snapshot, or a note), so that the page below it does not jump
VIEW_PX = 620
PLOTLY_CONFIG = {"displaylogo": False, "responsive": True}
SELECTIONS_OPEN = "selections_open"  # session key: whether the 2nd left column is shown
COLUMNS_OPEN, COLUMNS_SHUT = [0.95, 1.7, 0.95], [1.7, 0.95]


# ── Session state: one set of keys per example ──

def _key(example: Example, name: str) -> str:
    """A session-state key of the example, e.g. ``spec/frank_read_src``."""
    return f"{name}/{example.name}"


def _spec(example: Example) -> Spec:
    """The example's current spec; its file's at first."""
    key = _key(example, "spec")
    if key not in st.session_state:
        st.session_state[key] = example.default
    return st.session_state[key]


def _replace_spec(example: Example, spec: Spec) -> None:
    """Code replaces the spec (a reset): the widgets start again."""
    st.session_state[_key(example, "spec")] = spec
    st.session_state[_key(example, "gen")] = st.session_state.get(_key(example, "gen"), 0) + 1


def _widget_prefix(example: Example) -> str:
    """The prefix of the step widgets' keys; it changes with every spec that code sets."""
    return f"{example.name}/g{st.session_state.get(_key(example, 'gen'), 0)}"


def _status(example: Example, folder: job.RunFolder | None) -> job.Status | None:
    """This session's last job status, if it belongs to the folder of the current spec."""
    run_id, status = st.session_state.get(_key(example, "status"), (None, None))
    return status if folder is not None and run_id == folder.run_id else None


def _folder(example: Example, spec: Spec) -> job.RunFolder | None:
    """The run folder of the spec's script, or None while the spec is not valid."""
    if spec.validate():
        return None
    return job.RunFolder.for_script(example.name, example.script_name, example.final_json,
                                    codegen.full_script(spec, example.final_json))


def _follow_step(example: Example, step: Step) -> None:
    """When the step changes, show its selections and its 3-D view.

    This runs before the view radio and the 2nd left column exist in this run.
    """
    last = _key(example, "laststep")
    if st.session_state.get(last) != step.key:
        st.session_state[last] = step.key
        st.session_state[SELECTIONS_OPEN] = True  # a newly chosen step shows its selections
        st.session_state[f"view_{example.name}"] = (
            "simulation" if step.key == "run" else "initial configuration")


# ── Cached data: what the 3-D window draws ──

@st.cache_data(show_spinner=False, max_entries=32)
def _initial_payload(spec: Spec) -> network3d.Payload:
    """The initial configuration of the spec, as PyDiS builds it, for the 3-D window."""
    return network3d.frames_payload([("initial configuration", network.initial(spec))])


@st.cache_data(show_spinner=False, max_entries=8)
def _run_payload(files: tuple[tuple[str, str, float], ...]) -> network3d.Payload:
    """A run's snapshots, for its files (label, path, mtime), for the 3-D window to replay."""
    return network3d.frames_payload([(label, network.read_json(Path(path)))
                                     for label, path, _ in files])


# ── The panels ──

def _panel_title(text: str) -> None:
    """A centred serif panel title."""
    st.markdown(style.panel_title_html(text), unsafe_allow_html=True)


def _mode_radio(label: str, options: tuple[str, ...], key: str) -> str:
    """A centred horizontal radio that picks a panel's mode (label collapsed)."""
    _, middle, _ = st.columns([0.5, 5, 0.5])
    with middle:
        return st.radio(label, options, key=key, horizontal=True, label_visibility="collapsed")


def _chart(target: Any, fig: go.Figure, key: str | None = None) -> None:
    """Show a figure in the column (or slot), at the column's width (a run's newest snapshot,
    the history); the networks of the 3-D window are otherwise ``network3d``'s."""
    target.plotly_chart(fig, width="stretch", theme=None, config=PLOTLY_CONFIG, key=key)


def _simulation_view(example: Example, slot: DeltaGenerator,
                     folder: job.RunFolder | None) -> None:
    """The run's snapshots, to replay, or what to do to get them."""
    files = results.frames(folder) if folder is not None else []
    with slot.container(height=VIEW_PX, border=False):
        if not files:
            st.info(f"Press Run in step {RUN_STEP.label.split('.')[0]} to simulate these "
                    "selections." if folder is not None
                    else "Fix the selections first; the code panel lists the problems.")
            return
        if not folder.results.exists():
            st.warning("This run did not finish: these are its snapshots so far.")
        payload = _run_payload(tuple((label, str(path), path.stat().st_mtime)
                                     for label, path in files))
        network3d.show(payload, key=f"net3d/run/{example.name}", scope=example.name,
                       frame_ms=FRAME_MS)


def _history_view(example: Example, spec: Spec, folder: job.RunFolder | None) -> None:
    """The number of nodes, and the density, against the time step, from the run's snapshots."""
    steps, nodes, density = (results.history(folder, spec.material.burgmag) if folder
                             else ([], [], []))
    if not steps:
        st.info("The history appears with a run: the number of nodes and the density of each "
                "snapshot of the network that it writes (output/disnet_<step>.json).")
        return
    _chart(st, views.history_figure(steps, nodes, density), key=f"history/{example.name}")


def _view_panel(example: Example, spec: Spec,
                folder: job.RunFolder | None) -> DeltaGenerator | None:
    """The 3-D window; returns the slot that a run draws into, in the simulation view."""
    view = _mode_radio("3-D view", VIEWS, f"view_{example.name}")
    if view == "initial configuration":
        if spec.validate():
            st.info("Fix the selections first; the code panel lists the problems.")
            return None
        try:
            payload = _initial_payload(spec)
        except Exception as exc:  # pylint: disable=broad-exception-caught
            st.error(f"PyDiS cannot build this initial configuration: {exc}")
            return None
        with st.container(height=VIEW_PX, border=False):
            network3d.show(payload, key=f"net3d/initial/{example.name}", scope=example.name)
        return None
    if view == "history":
        _history_view(example, spec, folder)
        return None
    slot = st.empty()
    _simulation_view(example, slot, folder)
    return slot


# ── The 2nd left column: the chosen step's selections, or the run controls ──

def _selections(example: Example, step: Step, spec: Spec) -> Spec:
    """The step's widgets; the user's changes become the example's spec."""
    ctx = steps.StepContext(key=_widget_prefix(example), example=example)
    new = (steps.run_settings(spec, ctx) if step.key == "run"
           else steps.render(step.key, spec, ctx))
    if new != spec:
        st.session_state[_key(example, "spec")] = new
    errors = new.validate()
    if errors:
        st.warning("To fix before a run:\n\n" + "\n".join(f"- {e}" for e in errors))
    return new


def _request_run(example: Example) -> None:
    """The Run button: run after this script run's widgets, and show the run meanwhile."""
    st.session_state[_key(example, "run")] = True
    st.session_state[f"view_{example.name}"] = "simulation"


def _show_status(slot: DeltaGenerator, example: Example, folder: job.RunFolder | None) -> None:
    """The last run's outcome, or whether results are ready, in the status slot."""
    status = _status(example, folder)
    if folder is None:
        slot.warning("Fix the selections first; the code panel lists the problems.")
    elif status is not None and status.state in output_panel.FINAL:
        text = f"{status.state.capitalize()} after {status.elapsed_s:.1f} s: {status.message}."
        {"done": slot.success, "stopped": slot.warning}.get(status.state, slot.error)(text)
    elif folder.results.exists():
        slot.success("The results of these selections are ready: see the simulation view.")
    else:
        slot.info(f"Press Run to simulate these selections. {example.run_hint}")


def _downloads(example: Example, folder: job.RunFolder) -> None:
    """Download buttons for a run's results; the files are read only when clicked."""
    stem = f"{example.name}_{folder.run_id}"
    st.markdown("**Download the results**")
    st.download_button("⬇ final network (JSON)", data=folder.results.read_bytes,
                       file_name=f"{stem}_{example.final_json}", mime="application/json",
                       on_click="ignore", key=f"{example.name}/download_final", width="stretch")
    st.download_button("⬇ the whole run (zip)", data=functools.partial(results.output_zip, folder),
                       file_name=f"{stem}.zip", mime="application/zip", on_click="ignore",
                       key=f"{example.name}/download_zip", width="stretch")
    st.caption("The zip holds the script, its spec and log, and its output/ folder: a snapshot "
               "of the network every write_freq steps (disnet_<step>.json, as "
               "DisNetManager.write_json writes it), the final network, and any other file that "
               "the script writes there.")


def _run_controls(example: Example, folder: job.RunFolder | None) -> DeltaGenerator:
    """The Run step's buttons: Run, Stop, the status and reset; returns the slot for progress."""
    k = example.name
    requested = bool(st.session_state.get(_key(example, "run")))
    left, right = st.columns(2)
    left.button("▶ Run", key=f"{k}/run", type="primary", width="stretch",
                disabled=folder is None or requested, on_click=_request_run, args=(example,))
    right.button("■ Stop", key=f"{k}/stop", width="stretch", disabled=not requested)
    slot = st.empty()
    _show_status(slot, example, folder)
    if folder is not None and folder.results.exists() and not requested:
        _downloads(example, folder)
    st.button("↺ reset to the example's values", key=f"{k}/reset", width="stretch",
              on_click=_replace_spec, args=(example, example.default), disabled=requested)
    st.caption("Changing any selection during a run stops it.")
    return slot


def _show_selections(show: bool) -> None:
    """Hide or show the 2nd left column (a button callback)."""
    st.session_state[SELECTIONS_OPEN] = show


def _selections_column(example: Example, step: Step, spec: Spec
                       ) -> tuple[Spec, job.RunFolder | None, DeltaGenerator | None]:
    """The step's title, then its selections (and on the Run step the run controls)."""
    title, hide = st.columns([5, 1], vertical_alignment="center")
    with title:
        _panel_title(step.label)
    hide.button("«", key="selections_hide", help="Hide this column",
                on_click=_show_selections, args=(False,))
    with st.container(border=True, key="opendis_selections"):  # styled in common/style.py
        spec = _selections(example, step, spec)
        folder = _folder(example, spec)
        return spec, folder, _run_controls(example, folder) if step.key == "run" else None


# ── The run ──

def _progress(slot: DeltaGenerator, status: job.Status, max_step: int) -> None:
    """The progress bar: the time step reached, and the elapsed time."""
    step = status.step
    text = (f"step {step} of {max_step}" if step is not None else "starting PyDiS")
    slot.progress(min((step or 0) / max_step, 1.0), text=f"{text} · {status.elapsed_s:.1f} s")


def _run(example: Example, spec: Spec, folder: job.RunFolder | None,
         slots: tuple[DeltaGenerator | None, DeltaGenerator, DeltaGenerator | None]) -> None:
    """Run the job for the spec (or reuse its results), then rerun to show the outcome."""
    st.session_state.pop(_key(example, "run"), None)
    if folder is None:
        return
    progress_slot, output_slot, view_slot = slots
    progress_slot = progress_slot or st.sidebar.empty()
    if not folder.results.exists():
        job.kill_leftovers()
        folder.prepare(codegen.full_script(spec, example.final_json), spec.to_json())
        drawn = {"step": None}

        def update(status: job.Status) -> None:
            st.session_state[_key(example, "status")] = (folder.run_id, status)
            _progress(progress_slot, status, spec.run.max_step)
            output_panel.show_output(output_slot, status.lines, output_panel.no_run())
            found = network.snapshots(folder.output)
            if view_slot is not None and found and found[-1][0] != drawn["step"]:
                step, path = found[-1]
                try:
                    net = network.read_json(path)
                except ValueError:  # still being written
                    return
                drawn["step"] = step
                with view_slot.container(height=VIEW_PX, border=False):
                    _chart(st, views.network_figure(net, f"step {step} · {views.describe(net)}"))

        status = job.run(folder, update)
        st.session_state[_key(example, "status")] = (folder.run_id, status)
    st.rerun()


def _show_problems(examples: list[Example], broken: dict[str, list[str]]) -> bool:
    """List the input files that have problems; True if there is no example at all to show."""
    if broken:
        with st.expander(f"Input files with problems ({len(broken)})", expanded=not examples):
            for name, messages in broken.items():
                st.error(f"**{name}**\n\n" + "\n".join(f"- {message}" for message in messages))
    if not examples:
        st.warning("There is no example to show: add an input file to the app's `examples/` "
                   "folder, or fix the ones above.")
    return not examples


def main() -> None:
    """Build the page: config, CSS, banner, sidebar, the three columns, then any run."""
    st.set_page_config(page_title=style.BANNER_TITLE, page_icon=str(FAVICON), layout="wide",
                       initial_sidebar_state=SIDEBAR_WIDTH_PX)
    st.markdown(style.css(), unsafe_allow_html=True)
    st.markdown(style.banner_html(style.BANNER_TITLE, LOGO), unsafe_allow_html=True)
    fetch_graph()
    examples, broken = catalog.scan()
    if _show_problems(examples, broken):
        return
    with st.container(key="opendis_picker", horizontal=True, vertical_alignment="center"):
        example = picker.example_picker(examples)  # pinned into the banner by the CSS
    step = picker.current_step(example)
    _follow_step(example, step)
    with st.sidebar:
        steps_box = st.container()  # the steps radio, drawn after the selections
        picker.prev_next(example)
    spec, folder, progress_slot = _spec(example), None, None
    shown = st.session_state.get(SELECTIONS_OPEN, True)
    if shown:
        second, center, right = st.columns(COLUMNS_OPEN)
        with second:
            spec, folder, progress_slot = _selections_column(example, step, spec)
    else:
        center, right = st.columns(COLUMNS_SHUT)
        folder = _folder(example, spec)
    with steps_box:
        picker.steps_radio(example, spec.changed_parts(example.default))
    status = _status(example, folder)
    with center:
        if not shown:
            st.button(f"» show {step.label}", key="selections_show",
                      on_click=_show_selections, args=(True,))
        _panel_title(example.title)
        view_slot = _view_panel(example, spec, folder)
        _panel_title("PyDiS code")
        mode = _mode_radio("Code panel", code_panel.MODES, f"code_{example.name}")
        code_panel.show(mode, example, spec, step)
    with right:
        _panel_title("Simulation output")
        output_slot = output_panel.panel(example, folder, status)
    if st.session_state.get(_key(example, "run")):
        _run(example, spec, folder, (progress_slot, output_slot, view_slot))


if __name__ == "__main__":
    main()
