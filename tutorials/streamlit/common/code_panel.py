"""common.code_panel: the code panel's content, this step's code, the full script, or the console.

``app.py`` draws the panel's title and its mode radio, and this module fills the panel. The
Python console is this browser tab's own session on the run's objects (``common/console.py``).
"""
import streamlit as st

from common import codegen, console
from common.catalog import Example
from common.job import RunFolder
from common.picker import Step
from common.spec import Spec

CONSOLE = "Python console"
MODES = ("this step", "full script", CONSOLE)
HEIGHT = 460  # px; the code scrolls inside
SESSION, LOG, SOURCE = "console", "console_log", "console_source"  # session keys
SHOWN = 20  # the console shows its last commands
INTRO = ("Your own Python session, to look into this run's objects and build understanding. "
         "`net`, `sim`, `state` and `G` are as after `python -i` on the full script, and "
         "`snapshot(step)` loads one of the run's snapshots (`snapshot_steps` lists them). The "
         "session runs in a process of its own, on a copy of the run's results: what you do "
         "here changes neither the runs nor the other panels, and nobody else sees it. A "
         f"command stops after {console.COMMAND_S:g} s, and the session ends after "
         f"{console.IDLE_S / 60:g} min without one.")


def show(mode: str, example: Example, spec: Spec, step: Step, folder: RunFolder | None) -> None:
    """Fill the panel for its mode."""
    if mode == CONSOLE:
        _console(folder)
        return
    errors = spec.validate()
    if errors:
        st.error("Fix these selections first:\n\n" + "\n".join(f"- {e}" for e in errors))
        return
    if mode == "this step":
        st.code(codegen.step_code(spec, step.key, example.final_json), language="python",
                line_numbers=True, height=HEIGHT)
        return
    script = codegen.full_script(spec, example.final_json)
    st.download_button(f"⬇ download {example.script_name}", script,
                       file_name=example.script_name, mime="text/x-python", on_click="ignore",
                       key=f"download/{example.name}")
    st.code(script, language="python", line_numbers=True, height=HEIGHT)


# ── The Python console ──

def _end_session(forget: bool) -> None:
    """End this tab's session; ``forget`` also clears its commands (a button callback too)."""
    session = st.session_state.pop(SESSION, None)
    if session is not None:
        session.close()
    if forget:
        st.session_state[LOG] = []


def _execute(folder: RunFolder) -> None:
    """The ▶ execute button: run the command, in a session that starts if need be (a callback)."""
    source = st.session_state.get(SOURCE, "")
    if not source.strip():
        return
    log = st.session_state.setdefault(LOG, [])
    session = st.session_state.get(SESSION)
    if session is None or not session.alive():
        _end_session(forget=False)
        try:
            session = console.Console(folder)
        except console.Busy as exc:
            log.append((source, str(exc)))
            return
        log.append((None, "# a new session on this run's objects\n" + session.banner.output))
        if not session.banner.alive:
            return
        st.session_state[SESSION] = session
    reply = session.execute(source)
    log.append((source, reply.output))
    if not reply.alive:
        st.session_state.pop(SESSION, None)


def _transcript(log: list) -> str:
    """The commands and what they printed, as Python's prompt shows them."""
    lines = []
    for source, output in log[-SHOWN:]:
        if source is not None:
            first, *rest = source.rstrip().split("\n")
            lines += [">>> " + first] + ["... " + line for line in rest]
        if output.strip():
            lines += output.rstrip("\n").split("\n")
    return "\n".join(lines)


def _console(folder: RunFolder | None) -> None:
    """The console on the run of the current selections, once the run has results."""
    if folder is None or not folder.results.exists():
        st.info("The console opens on a run's results: press Run in step 7 first.")
        return
    session = st.session_state.get(SESSION)
    if session is not None and session.run_id != folder.run_id:
        _end_session(forget=True)  # it was on another run's objects
    st.caption(INTRO)
    log = st.session_state.get(LOG, [])
    if log:
        st.code(_transcript(log), language="python")
    with st.form("console_form", clear_on_submit=True, border=False):
        st.text_area("Python", key=SOURCE, height=100, label_visibility="collapsed",
                     placeholder="G.num_nodes()")
        st.form_submit_button("▶ execute", type="primary", on_click=_execute, args=(folder,))
    st.button("↺ restart", key="console_restart", on_click=_end_session, args=(True,),
              help="A new session, without the names that you made")
