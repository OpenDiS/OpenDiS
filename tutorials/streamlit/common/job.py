"""common.job: run a generated script as a separate PyDiS job, with progress, time limits and a kill.

The job runs in its own process group, in its run folder, so that the script's ``output/`` goes
there. Matplotlib draws off screen (``MPLBACKEND=Agg``), so the script that runs is the one the
code panel shows, plots included, and no window opens; the one warning that this gives is
filtered out, so the job's output is what PyDiS prints. The kill is a safety net, used only for:

- a job that runs too long, or prints nothing for too long;
- a job the app stops (Stop, a rerun, a closed tab);
- the app shutting down;
- jobs left by an earlier app session.

This module imports neither Streamlit nor PyDiS.
"""
import atexit
import hashlib
import os
import queue
import re
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from common.paths import OPENDIS_DIR, RUNS

LOG, PID, SPEC = "log.txt", "job.pid", "spec.json"
OUTPUT = "output"  # the scripts' write_dir
POLL_S = 0.5       # the run loop wakes at least this often
_STEP = re.compile(r"\s*step\s*=\s*(\d+)\b")  # SimulateNetwork's print_freq line
_RUNNING: set[subprocess.Popen] = set()  # this app process's jobs, for the kill at shutdown
JOB_ENV = {"MPLBACKEND": "Agg",
           "PYTHONWARNINGS": "ignore:FigureCanvasAgg is non-interactive"}


@dataclass(frozen=True)
class Limits:
    """The time limits, in seconds."""

    first_line_s: float = 120.0  # the first line is the initial configuration's
    silence_s: float = 300.0     # no new line for this long
    total_s: float = 3600.0      # 1 h in total
    term_s: float = 2.0          # after SIGTERM, before SIGKILL


@dataclass(frozen=True)
class RunFolder:
    """A run's folder, ``runs/<example>/<run_id>/``; ``run_id`` hashes the script."""

    path: Path
    script_name: str
    final_json: str  # the final network, which the script writes into ``output/`` at its end

    @classmethod
    def for_script(cls, example: str, script_name: str, final_json: str, script_text: str,
                   runs: Path = RUNS) -> "RunFolder":
        """The folder of this script text: the same code is the same run."""
        run_id = hashlib.sha256(script_text.encode("utf-8")).hexdigest()[:12]
        return cls(runs / example / run_id, script_name, final_json)

    @property
    def run_id(self) -> str:
        """The short hash of the script text."""
        return self.path.name

    @property
    def script(self) -> Path:
        """The script the job runs."""
        return self.path / self.script_name

    @property
    def output(self) -> Path:
        """The script's ``output/`` folder: the snapshots, the results and the final network."""
        return self.path / OUTPUT

    @property
    def results(self) -> Path:
        """The final network, written last: the run is complete when it is there."""
        return self.output / self.final_json

    @property
    def log(self) -> Path:
        """Everything the job printed, then one ``[app]`` line on how it ended."""
        return self.path / LOG

    def prepare(self, script_text: str, spec_json: str | None = None) -> None:
        """A fresh folder with the script and, if there is one, ``spec.json`` in it."""
        self.path.mkdir(parents=True, exist_ok=True)
        for old in self.output.glob("*") if self.output.is_dir() else []:
            old.unlink()  # a stopped run's snapshots would mix with this run's
        self.script.write_text(script_text, encoding="utf-8")
        if spec_json is not None:
            (self.path / SPEC).write_text(spec_json, encoding="utf-8")


@dataclass
class Status:
    """A job's state while it runs and at its end."""

    state: str = "starting"  # starting, running; then done, failed, timeout, stopped
    lines: list[str] = field(default_factory=list)
    elapsed_s: float = 0.0
    exit_code: int | None = None
    killed: bool = False
    message: str = ""

    @property
    def ok(self) -> bool:
        """True when the job exited with code 0 and wrote the final network."""
        return self.state == "done"

    @property
    def step(self) -> int | None:
        """The time step of the last ``step = k ...`` line, if there is one: PyDiS prints
        ``step = k, nodes = n, ...`` or, before it printed the nodes, ``step = k dt = ...``."""
        for line in reversed(self.lines):
            match = _STEP.match(line)
            if match:
                return int(match.group(1))
        return None

    def summary(self) -> str:
        """One line on how the job ended, for ``log.txt``."""
        return f"{self.state} after {self.elapsed_s:.1f} s: {self.message}"


def _launch(folder: RunFolder, python: str) -> subprocess.Popen:
    """Start the job in its own process group, and note its pid and script in ``job.pid``."""
    env = {**os.environ, **JOB_ENV, "OPENDIS_DIR": str(OPENDIS_DIR)}
    # The process outlives this function, so it is not opened in a with block.
    proc = subprocess.Popen(  # pylint: disable=consider-using-with
        [python, "-u", str(folder.script)], cwd=folder.path, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
        start_new_session=True, env=env)
    _RUNNING.add(proc)
    (folder.path / PID).write_text(f"{proc.pid}\n{folder.script}\n", encoding="utf-8")
    return proc


def read_lines(proc: subprocess.Popen) -> queue.Queue:
    """A daemon thread that puts the job's lines on a queue, and ``None`` when they end."""
    lines: queue.Queue = queue.Queue()

    def read() -> None:
        for line in proc.stdout:
            lines.put(line.rstrip("\r\n"))
        lines.put(None)

    threading.Thread(target=read, daemon=True).start()
    return lines


def take_lines(lines: queue.Queue, timeout: float) -> tuple[list[str], bool]:
    """The lines that arrive within ``timeout``, and whether the output has ended."""
    taken: list[str] = []
    try:
        item = lines.get(timeout=timeout)
        while item is not None:
            taken.append(item)
            item = lines.get_nowait()
        return taken, True
    except queue.Empty:
        return taken, False


def kill(proc: subprocess.Popen, term_s: float = Limits.term_s) -> None:
    """SIGTERM to the job's process group, SIGKILL after ``term_s`` if it lives on, then wait."""
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=term_s)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def _over_limit(limits: Limits, start: float, last_line: float, now: float,
                started: bool) -> str:
    """Why the job must be stopped now, or ``""``."""
    if not started:
        if now - start > limits.first_line_s:
            return f"no output within {limits.first_line_s:g} s, so the app stopped the job"
        return ""
    if now - last_line > limits.silence_s:
        return f"no new output for {limits.silence_s:g} s, so the app stopped the job"
    if now - start > limits.total_s:
        return f"the job ran longer than {limits.total_s / 60:g} min, so the app stopped it"
    return ""


def _end(status: Status, folder: RunFolder, reason: str) -> None:
    """Set the final state: success is exit code 0 and the final network written."""
    if reason:
        status.state, status.message = "timeout", reason
    elif status.exit_code == 0 and folder.results.exists():
        status.state, status.message = "done", "the job exited by itself with code 0"
    elif status.exit_code == 0:
        status.state = "failed"
        status.message = f"the job exited with code 0 but wrote no {OUTPUT}/{folder.final_json}"
    else:
        status.state = "failed"
        status.message = f"the job exited with code {status.exit_code}"


def _record(status: Status, new: list[str], log: TextIO) -> None:
    """Add new lines to the status and the log."""
    for line in new:
        status.lines.append(line)
        log.write(line + "\n")
    if new:
        log.flush()


def _watch(proc: subprocess.Popen, folder: RunFolder, status: Status, log: TextIO,
           on_update: Callable[[Status], None] | None, limits: Limits) -> None:
    """The run loop: read lines, call ``on_update``, and stop the job at a limit."""
    lines = read_lines(proc)
    start = last_line = time.monotonic()
    while True:
        new, _ = take_lines(lines, POLL_S)
        _record(status, new, log)
        now = time.monotonic()
        last_line = now if new else last_line
        status.elapsed_s = now - start
        if proc.poll() is not None:  # it exited by itself: read what is left
            ended, deadline = False, now + 2.0
            while not ended and time.monotonic() < deadline:
                new, ended = take_lines(lines, 0.1)
                _record(status, new, log)
            status.exit_code = proc.returncode
            _end(status, folder, "")
            return
        reason = _over_limit(limits, start, last_line, now, bool(status.lines))
        if reason:
            kill(proc, limits.term_s)
            status.exit_code, status.killed = proc.returncode, True
            _end(status, folder, reason)
            return
        status.state = "running" if status.lines else "starting"
        if on_update is not None:
            on_update(status)


def run(folder: RunFolder, on_update: Callable[[Status], None] | None = None,
        limits: Limits = Limits(), python: str = sys.executable) -> Status:
    """Run the folder's script as a job and return how it ended.

    ``on_update`` is called at least every ``POLL_S`` seconds while the job
    runs. If it raises (Streamlit stops a script run that way), the job is
    killed and the exception goes on.
    """
    status = Status()
    with open(folder.log, "w", encoding="utf-8") as log:
        proc = _launch(folder, python)
        try:
            _watch(proc, folder, status, log, on_update, limits)
        finally:
            if proc.poll() is None:  # stopped from outside, by on_update
                kill(proc, limits.term_s)
                status.state, status.killed, status.exit_code = "stopped", True, proc.returncode
                status.message = "the app stopped the job (Stop, a rerun or a closed tab)"
            _RUNNING.discard(proc)
            (folder.path / PID).unlink(missing_ok=True)
            log.write(f"[app] {status.summary()}\n")
    return status


def command_line(pid: int) -> str:
    """The command line of a running process, or ``""`` if there is none."""
    result = subprocess.run(["ps", "-o", "command=", "-p", str(pid)], capture_output=True,
                            text=True, check=False)
    return result.stdout.strip()


def _kill_group(pid: int, script: str, term_s: float = Limits.term_s) -> None:
    """SIGTERM, then SIGKILL, to the process group of a job that is not this process's child."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(pid, sig)
        except ProcessLookupError:
            return
        deadline = time.monotonic() + term_s
        while time.monotonic() < deadline:
            if script not in command_line(pid):
                return
            time.sleep(0.1)


def kill_leftovers(runs: Path = RUNS) -> list[int]:
    """Kill the jobs that an earlier app session left running; returns their pids.

    Each run's ``job.pid`` names the job's pid and script. The process is killed
    only if it is alive and its command line contains that script. Jobs of this
    app process are left alone.
    """
    ours = {proc.pid for proc in _RUNNING}
    killed = []
    for pid_file in sorted(runs.glob(f"*/*/{PID}")):
        try:
            pid_text, script = pid_file.read_text(encoding="utf-8").split("\n")[:2]
            pid = int(pid_text)
        except (OSError, ValueError):
            pid_file.unlink(missing_ok=True)
            continue
        if pid in ours:
            continue
        if script and script in command_line(pid):
            _kill_group(pid, script)
            killed.append(pid)
        pid_file.unlink(missing_ok=True)
    return killed


def running() -> bool:
    """True while a job of this app process runs."""
    return any(proc.poll() is None for proc in _RUNNING)


@atexit.register
def _kill_running() -> None:
    """At app shutdown, kill every job that is still running."""
    for proc in list(_RUNNING):
        if proc.poll() is None:
            kill(proc)
