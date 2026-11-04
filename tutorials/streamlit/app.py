"""OpenDiS: the Streamlit app (entry point).

Run it from this folder with ``./launch_app.sh``, or with
``python -m streamlit run app.py``. The page has the pinned banner; the sidebar (the example and
its code, the steps, prev / next, two links); the chosen step's selections in the 2nd left
column; the 3-D window and the code panel in the centre; and on the right the simulation's
screen output, then the status of ExaDiS's build, or Build ExaDiS. Run launches the PyDiS or
ExaDiS job (``common/job.py``) as the last thing in a script run, so that Stop, or any other
click, interrupts it and kills the job.
"""
import functools
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from common import (catalog, code_panel, codegen, front_page, job, network, network3d,
                    output_panel, picker, results, steps, style, views)
from common.catalog import Example
from common.paths import (ASSETS, build_libraries_once, exadis_status, fetch_graph,
                          start_exadis_build)
from common.picker import Step, run_step
from common.spec import ENGINE_NAMES, LINE_SLIP, Precipitates, Spec, UserStressLine

LOGO = ASSETS / "opendis_logo.png"
FAVICON = ASSETS / "opendis_favicon.png"  # the website's favicon
SIDEBAR_WIDTH_PX = 300  # initial width; Streamlit clamps it to 200-600
# The sidebar's last lines: the app online, and on your own computer
LINKS = ("A [list of engines](https://opendis.streamlit.app/) to run the app online\n\n"
         "[Instructions](https://github.com/OpenDiS/OpenDiS/blob/main/tutorials/streamlit/"
         "README.md) on how to run the app locally")
VIEWS = ("initial configuration", "simulation", "history")
BUILD_VIEWS = ("initial configuration", "run")  # a script that builds a network: no history
FRAME_MS = 300  # the replay's time per snapshot
# The 3-D window's height, whatever it shows (the figure and the replay's controls, a run's newest
# snapshot, or a note), so that the page below it does not jump
VIEW_PX = 620
PLOTLY_CONFIG = {"displaylogo": False, "responsive": True}
SELECTIONS_OPEN = "selections_open"  # session key: whether the 2nd left column is shown
CLEARED = "results_cleared"  # session key: the note on the last clear of all results, shown once
CONFIRM_CLEAR = "confirm_clear"  # session key: whether that clear's confirmation is shown
EXADIS_STATE = "exadis_state"  # session key: the state of ExaDiS's build at the last page run
BUILD_POLL_S = 5  # the build's status bar refreshes itself this often while ExaDiS builds
# Why an ExaDiS example cannot run, by the state of ExaDiS's build
NOT_READY = {
    "unbuilt": "ExaDiS is not built on this server yet. Press Build ExaDiS at the bottom of the "
               "right column. Run turns on when ExaDiS is ready.",
    "building": "ExaDiS is being built on this server. Run turns on when it is ready: its status "
                "bar is in the right column.",
    "failed": "ExaDiS could not be built on this server: the right column says why. Build ExaDiS "
              "there tries again.",
    "missing": "ExaDiS is not built for this Python. Build it with CMake, as for the tutorials.",
}
BUILD_HELP = ("Builds ExaDiS on this server, for all its users, in some minutes. PyDiS runs "
              "meanwhile.")
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


def _renew_specs() -> None:
    """Rebuild the session's specs after Streamlit loaded the app's modules again.

    A changed source file (an edit, or an update of the app on the Cloud) makes Streamlit load
    every module of the app again, with new classes. A spec that the session kept is then of
    the old classes, and fails the code's ``isinstance`` tests; its JSON gives the new ones.
    """
    for key, value in list(st.session_state.items()):
        if str(key).startswith("spec/") and not isinstance(value, Spec):
            try:
                st.session_state[key] = Spec.from_json(value.to_json())
            except Exception:  # pylint: disable=broad-exception-caught
                del st.session_state[key]  # the example's own values, then


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


def _views(example: Example) -> tuple[str, ...]:
    """The 3-D window's views of the example: its second is the run's."""
    return VIEWS if example.default.config.simulates else BUILD_VIEWS


def _follow_step(example: Example, step: Step) -> None:
    """When the step changes, show its selections and its 3-D view.

    This runs before the view radio and the 2nd left column exist in this run.
    """
    last = _key(example, "laststep")
    if st.session_state.get(last) != step.key:
        st.session_state[last] = step.key
        st.session_state[SELECTIONS_OPEN] = True  # a newly chosen step shows its selections
        st.session_state[f"view_{example.name}"] = _views(example)[1 if step.key == "run" else 0]


# ── Cached data: what the 3-D window draws ──

@st.cache_resource(show_spinner="Drawing the 3-D window's background…", max_entries=8)
def _background(kind: str, *args) -> list:
    """What the 3-D window draws under the networks, by its key (``_background_key``): the same
    traces for every session, which nothing changes. The precipitates' spheres; or the user
    stress's field on its line's slip plane, at the offset (b): the frozen field is the link's,
    and user_stress.py's is drawn at t = 0, which each frame recolours at its time
    (``_plane_colors``)."""
    if kind == "spheres":
        spec, = args
        return views.spheres(*network.precipitates(spec), spec.config.box_length)
    url, box, offset = args
    if kind == "space_time":
        return views.slip_plane(network.space_time_field(box), box, LINE_SLIP, offset,
                                views.SPACE_TIME_CLIP, views.SPACE_TIME_TITLE,
                                views.SPACE_TIME_POINTS, 0)
    return views.slip_plane(network.user_field(url, box), box, LINE_SLIP, offset)


def _plane_colors(key: tuple | None, times: list[float]) -> list | None:
    """The colours of the slip plane's surfaces at each time (s), for a field that changes in
    time, user_stress.py's; None for a frozen field, the precipitates, or nothing: the
    background stays as drawn."""
    if key is None or key[0] != "space_time":
        return None
    _, _, box, offset = key
    return [views.plane_colors(network.space_time_field(box, t), box, LINE_SLIP, offset,
                               views.SPACE_TIME_POINTS, 0) for t in times]


def _on_plane(key: tuple | None) -> bool:
    """Whether the networks are drawn as over the user stress's slip plane: black, in the
    plane's first view (``views.network_figure``)."""
    return key is not None and key[0] != "spheres"


def _background_key(spec: Spec) -> tuple | None:
    """What the 3-D window draws under the networks, by what it depends on: the precipitates,
    by the spec; or the user stress's slip plane, by its field's kind and link, the box and the
    plane's offset, when its field is on and can be read; else None, with a warning if the
    field cannot be read."""
    if isinstance(spec.config, Precipitates):
        return ("spheres", spec)
    f = spec.forces
    if not isinstance(spec.config, UserStressLine) or not f.stress_field:
        return None
    box = spec.config.box_length
    key = (f.stress_kind, f.stress_url if f.stress_kind == "frozen" else None, box,
           views.plane_offset(box, LINE_SLIP, codegen.line_origin(box, spec.config.seed)))
    try:
        _background(*key)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        st.warning(f"The user stress's field cannot be drawn on its slip plane: {exc}")
        return None
    return key


@st.cache_data(show_spinner=False, max_entries=32)
def _initial_payload(spec: Spec, key: tuple | None = None) -> network3d.Payload:
    """The initial configuration of the spec, as PyDiS builds it, for the 3-D window, over its
    background, if given, as at t = 0. (Streamlit keys this cache on this code, not on
    _background's: a change there must change this too, or a server keeps old backgrounds.)"""
    return network3d.frames_payload([("initial configuration", network.initial(spec))],
                                    _background(*key) if key else None, _on_plane(key))


@st.cache_data(show_spinner=False, max_entries=8)
def _run_payload(files: tuple[tuple[str, str, float], ...], key: tuple | None = None,
                 times: tuple[float, ...] = ()) -> network3d.Payload:
    """A run's snapshots, for its files (label, path, mtime), for the 3-D window to replay, over
    its background, if given, which the frames' times (s) recolour when the user stress's field
    changes in time."""
    return network3d.frames_payload([(label, network.read_json(Path(path)))
                                     for label, path, _ in files],
                                    _background(*key) if key else None, _on_plane(key),
                                    colors=_plane_colors(key, list(times)) if times else None)


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
    """Show a figure in the column (or slot), at the column's width (the history); the networks
    of the 3-D window are ``network3d``'s."""
    target.plotly_chart(fig, width="stretch", theme=None, config=PLOTLY_CONFIG, key=key)


def _simulation_view(example: Example, spec: Spec, slot: DeltaGenerator,
                     folder: job.RunFolder | None) -> None:
    """The run's snapshots, to replay, or what to do to get them."""
    files = results.frames(folder) if folder is not None else []
    with slot.container(height=VIEW_PX, border=False):
        if not files:
            what = ("simulate these selections" if example.default.config.simulates
                    else f"build this network with {ENGINE_NAMES[example.engine]}")
            number = run_step(example).label.split(".")[0]
            st.info(f"Press Run in step {number} to {what}." if folder is not None
                    else "Fix the selections first; the code panel lists the problems.")
            return
        if not folder.results.exists():
            st.warning("This run did not finish: these are its snapshots so far.")
        key, times = _background_key(spec), ()
        if key is not None and key[0] == "space_time":  # the frames' times, from the log
            motion = results.displacement(folder)
            times = tuple(results.frame_times(
                folder.log.read_text(encoding="utf-8").splitlines() if folder.log.exists()
                else [], [label for label, _ in files], motion[0][-1] if motion else None))
        payload = _run_payload(tuple((label, str(path), path.stat().st_mtime)
                                     for label, path in files), key, times)
        network3d.show(payload, key=f"net3d/run/{example.name}", scope=example.name,
                       frame_ms=FRAME_MS)


def _motion(spec: Spec, folder: job.RunFolder) -> tuple | None:
    """The user stress's line: the time and its displacement, and its drag velocity M b τ, for
    FCC_0's mobility of its character (None with another mobility law)."""
    found = results.displacement(folder)
    if found is None:
        return None
    f = spec.forces
    mobility = (f.Mscrew if spec.config.theta == 0.0 else f.Medge) \
        if f.mobility_law == "FCC_0" else None
    return (*found, None if mobility is None else mobility * spec.material.burgmag
            * spec.loading.tau)


def _history_view(example: Example, spec: Spec, folder: job.RunFolder | None) -> None:
    """The number of nodes, and the density, against the time step, from the run's snapshots;
    the stress against the strain, or the user stress's displacement against the time."""
    steps, nodes, density = (results.history(folder, spec.material.burgmag) if folder
                             else ([], [], []))
    if not steps:
        st.info("The history appears with a run: the number of nodes and the density of each "
                "snapshot of the network that it writes (output/disnet_<step>.json)"
                + (", and the line's displacement against time." if spec.loading.tau is not None
                   else "."))
        return
    curve = (results.stress_strain(folder) if spec.loading.loading_mode == "strain_rate"
             else None)
    motion = _motion(spec, folder) if spec.loading.tau is not None else None
    _chart(st, views.history_figure(steps, nodes, density, curve, motion),
           key=f"history/{example.name}")


def _view_panel(example: Example, spec: Spec,
                folder: job.RunFolder | None) -> DeltaGenerator | None:
    """The 3-D window; returns the slot that a run draws into, in the simulation view."""
    view = _mode_radio("3-D view", _views(example), f"view_{example.name}")
    if view == "initial configuration":
        if spec.validate():
            st.info("Fix the selections first; the code panel lists the problems.")
            return None
        try:
            payload = _initial_payload(spec, _background_key(spec))
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
    _simulation_view(example, spec, slot, folder)
    return slot


# ── The 2nd left column: the chosen step's selections, or the run controls ──

def _selections(example: Example, step: Step, spec: Spec) -> Spec:
    """The step's widgets; the user's changes become the example's spec."""
    ctx = steps.StepContext(key=_widget_prefix(example), example=example)
    if step.key == "config":
        picker.variant_picker(example, catalog.examples())
    if step.key == "run" and spec.run is None:
        st.caption(f"This script builds a network and runs no simulation. Run runs it with "
                   f"{ENGINE_NAMES[example.engine]}.")
        new = spec
    else:
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
    st.session_state[f"view_{example.name}"] = _views(example)[1]


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
    st.caption("The zip holds the script, its spec and log, and its output/ folder: "
               + ("a snapshot of the network every write_freq steps (disnet_<step>.json, as "
                  "DisNetManager.write_json writes it), the final network, "
                  if example.default.config.simulates else "the network that it built, ")
               + "and any other file that the script writes there.")


def _not_ready(example: Example) -> str:
    """Why the example's code cannot run yet, or "" when it can: ExaDiS before its build."""
    return "" if example.engine != "exadis" else NOT_READY.get(exadis_status()["state"], "")


def _run_controls(example: Example, folder: job.RunFolder | None) -> DeltaGenerator:
    """The Run step's buttons: Run, Stop, the status and reset; returns the slot for progress."""
    k = example.name
    requested = bool(st.session_state.get(_key(example, "run")))
    not_ready = _not_ready(example)
    left, right = st.columns(2)
    left.button("▶ Run", key=f"{k}/run", type="primary", width="stretch",
                disabled=folder is None or requested or bool(not_ready), on_click=_request_run,
                args=(example,))
    right.button("■ Stop", key=f"{k}/stop", width="stretch", disabled=not requested)
    slot = st.empty()
    if not_ready and not (folder is not None and folder.results.exists()):
        slot.warning(not_ready)
    else:
        _show_status(slot, example, folder)
    if folder is not None and folder.results.exists() and not requested:
        _downloads(example, folder)
    st.button("↺ reset to the example's values", key=f"{k}/reset", width="stretch",
              on_click=_replace_spec, args=(example, example.default), disabled=requested)
    st.caption("Changing any selection during a run stops it.")
    return slot


def _ask_clear(ask: bool) -> None:
    """Show or hide the confirmation of 🗑 clear all results (a button callback)."""
    st.session_state[CONFIRM_CLEAR] = ask


def _clear_results() -> None:
    """Delete the results of all runs but those that run (a button callback)."""
    deleted, kept = job.clear_runs()
    for key in [k for k in st.session_state if str(k).startswith("status/")]:
        del st.session_state[key]  # this session's last runs are gone too
    st.session_state[CONFIRM_CLEAR] = False
    st.session_state[CLEARED] = (f"Deleted the results of {deleted} run(s)."
                                 + (f" {kept} run(s) in progress are kept." if kept else ""))


def _clear_button(example: Example) -> None:
    """The button that deletes the results of all past runs, after a confirmation in place."""
    note = st.session_state.pop(CLEARED, None)
    if note:
        st.toast(note)
    running = bool(st.session_state.get(_key(example, "run")))
    if not st.session_state.get(CONFIRM_CLEAR):
        st.button("🗑 clear all results", key="clear_ask", width="stretch", disabled=running,
                  on_click=_ask_clear, args=(True,))
        return
    st.warning("This deletes the results of all past runs: those of every example, for everyone "
               "who uses this app. A run in progress is kept.")
    delete, cancel = st.columns(2)
    delete.button("Delete", key="clear_results", type="primary", width="stretch",
                  disabled=running, on_click=_clear_results)
    cancel.button("Cancel", key="clear_cancel", width="stretch", on_click=_ask_clear,
                  args=(False,))


def _exadis_bar(running: bool) -> None:
    """The status of ExaDiS's build, below 🗑 clear all results, and Build ExaDiS while the app
    can build ExaDiS and has not: the app builds it only at that press. It is a fragment, which
    refreshes itself while ExaDiS builds; when the build's state changes, the whole page runs
    again, so that the Run buttons follow."""
    status = exadis_status()
    if status["state"] != st.session_state.get(EXADIS_STATE, status["state"]):
        st.rerun()
    if status["state"] != "unbuilt":
        with st.container(key="opendis_build_status"):  # its text is small (common/style.py)
            _exadis_lines(status)
    if status["state"] in ("unbuilt", "failed"):  # not while a run of this session goes on
        st.button("Build ExaDiS", key="build_exadis", width="stretch", disabled=running,
                  on_click=start_exadis_build, help=BUILD_HELP)


def _exadis_lines(status: dict) -> None:
    """The status bar's text, by the state of ExaDiS's build."""
    if status["state"] == "building":
        file = f" {status['file']} ·" if status.get("file") else ""
        st.progress(status["percent"] / 100,
                    text=f"Building ExaDiS: {status['phase']}, {status['percent']} %,{file} "
                         f"{status['elapsed_s'] / 60:.1f} min so far")
    elif status["state"] == "ready":
        st.caption(f"ExaDiS is ready ({status['module']}).")
    elif status["state"] == "failed":
        st.error(f"ExaDiS was not built: {status['error']}. The app's log says more.")
    else:
        st.caption("ExaDiS is not built for this Python. The app builds it on Linux, as on "
                   "Streamlit Community Cloud, at the press of Build ExaDiS. Elsewhere, build it "
                   "with CMake.")


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

def _progress(slot: DeltaGenerator, status: job.Status, max_step: int | None,
              code: str) -> None:
    """The progress bar: the time step reached, and the elapsed time (a build has no steps)."""
    step = status.step
    if max_step is None:
        slot.progress(0.0, text=f"{code} builds the network · {status.elapsed_s:.1f} s")
        return
    text = (f"step {step} of {max_step}" if step is not None else f"starting {code}")
    slot.progress(min((step or 0) / max_step, 1.0), text=f"{text} · {status.elapsed_s:.1f} s")


def _run(example: Example, spec: Spec, folder: job.RunFolder | None,
         slots: tuple[DeltaGenerator | None, DeltaGenerator, DeltaGenerator | None]) -> None:
    """Run the job for the spec (or reuse its results), then rerun to show the outcome."""
    st.session_state.pop(_key(example, "run"), None)
    if folder is None or _not_ready(example):
        return
    progress_slot, output_slot, view_slot = slots
    progress_slot = progress_slot or st.sidebar.empty()
    if not folder.results.exists():
        job.kill_leftovers()
        folder.prepare(codegen.full_script(spec, example.final_json), spec.to_json())
        drawn = {"step": None}
        key = _background_key(spec) if view_slot is not None else None
        background = _background(*key) if key else []

        def update(status: job.Status) -> None:
            st.session_state[_key(example, "status")] = (folder.run_id, status)
            _progress(progress_slot, status, spec.run.max_step if spec.run else None,
                      ENGINE_NAMES[example.engine])
            output_panel.show_output(output_slot, status.lines, output_panel.no_run(example))
            found = network.snapshots(folder.output)
            if view_slot is not None and found and found[-1][0] != drawn["step"]:
                step, path = found[-1]
                try:
                    net = network.read_json(path)
                except ValueError:  # still being written
                    return
                colors = _plane_colors(key, results.frame_times(status.lines, [f"step {step}"]))
                network3d.write_live(folder.run_id, f"step {step}", net, on_plane=_on_plane(key),
                                     colors=colors[0] if colors else None)
                if drawn["step"] is None:  # the 3-D window then fetches each new snapshot
                    with view_slot.container(height=VIEW_PX, border=False):
                        network3d.show_live(folder.run_id, f"net3d/live/{example.name}",
                                            example.name, background)
                drawn["step"] = step

        try:
            status = job.run(folder, update)
        finally:
            network3d.end_live(folder.run_id)
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
    if front_page.wanted():  # the front page's deployment (common/front_page.py)
        front_page.main()
        return
    st.set_page_config(page_title=style.BANNER_TITLE, page_icon=str(FAVICON), layout="wide",
                       initial_sidebar_state=SIDEBAR_WIDTH_PX)
    _renew_specs()
    st.markdown(style.css(), unsafe_allow_html=True)
    st.markdown(style.banner_html(style.BANNER_TITLE, LOGO), unsafe_allow_html=True)
    fetch_graph()
    build_libraries_once()
    examples, broken = catalog.scan()
    if _show_problems(examples, broken):
        return
    with st.sidebar:
        with st.container(key="opendis_picker"):
            example = picker.example_picker(examples)
        steps_box = st.container()  # the steps radio, drawn after the selections
        picker.prev_next(example)
        with st.container(key="opendis_links"):  # small and grey (common/style.py)
            st.caption(LINKS)
    step = picker.current_step(example)
    _follow_step(example, step)
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
        _panel_title(f"{ENGINE_NAMES[example.engine]} code")
        mode = _mode_radio("Code panel", code_panel.MODES, f"code_{example.name}")
        code_panel.show(mode, example, spec, step, folder)
    with right:
        _panel_title("Simulation output")
        output_slot = output_panel.panel(example, folder, status)
        _clear_button(example)
        state = exadis_status()["state"]
        st.session_state[EXADIS_STATE] = state
        st.fragment(_exadis_bar, run_every=BUILD_POLL_S if state == "building" else None)(
            bool(st.session_state.get(_key(example, "run"))))
    if st.session_state.get(_key(example, "run")):
        _run(example, spec, folder, (progress_slot, output_slot, view_slot))


if __name__ == "__main__":
    main()
