"""common.output_panel: the right column, the simulation's screen output.

It shows what the job prints, as the script would print it in a terminal: during a run, the
lines as they come (the app writes them into the slot that ``panel`` returns); after a run, the
run folder's ``log.txt``, which ends with one ``[app]`` line on how the job ended.
"""
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from common import job
from common.catalog import Example
from common.picker import run_step

TAIL = 80  # the lines shown; log.txt has them all
FINAL = ("done", "failed", "timeout", "stopped")


def no_run(example: Example) -> str:
    """What the column says before a run of these selections."""
    number = run_step(example).label.split(".")[0]
    return f"No run yet for these selections: press Run in step {number}."


def output_lines(folder: job.RunFolder | None, status: job.Status | None) -> list[str]:
    """The job's lines: from this session's run of the folder, else from its log.txt."""
    if status is not None:
        lines = list(status.lines)
        return lines + [f"[app] {status.summary()}"] if status.state in FINAL else lines
    if folder is not None and folder.log.exists():
        return folder.log.read_text(encoding="utf-8", errors="replace").splitlines()
    return []


def show_output(slot: DeltaGenerator, lines: list[str], empty: str) -> None:
    """Write the job's last lines into the output slot, or ``empty`` if it has none."""
    if not lines:
        slot.code(empty, language="text", wrap_lines=True)
        return
    skipped = len(lines) - TAIL
    head = [f"... ({skipped} earlier lines are in log.txt)"] if skipped > 0 else []
    slot.code("\n".join(head + lines[-TAIL:]), language="text")


def panel(example: Example, folder: job.RunFolder | None,
          status: job.Status | None) -> DeltaGenerator:
    """The output, and the log to download; returns the slot that a run writes into."""
    with st.container(key="opendis_output"):  # styled in common/style.py
        slot = st.empty()
        show_output(slot, output_lines(folder, status), no_run(example))
        if (folder is not None and folder.log.exists()
                and (status is None or status.state in FINAL)):
            st.download_button("⬇ log.txt", data=folder.log.read_bytes, mime="text/plain",
                               file_name=f"{example.name}_{folder.run_id}_log.txt",
                               on_click="ignore", key=f"download_log/{example.name}",
                               width="stretch")
    return slot
