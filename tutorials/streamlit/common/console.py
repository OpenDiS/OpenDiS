"""common.console: a Python console on a run's objects, in a process of its own (README.md).

``Console`` copies a run's script and ``output/`` into a temporary folder and starts
``console_child.py`` there. The child runs the script with the run's final network loaded in
place of the run, so that ``net``, ``sim``, ``state`` and ``G`` are the run's, and then runs the
commands that ``execute`` sends. A command that runs too long ends the session. The session sees
only its own copy of the results, and an environment with no more than it needs; at most
``MAX_LIVE`` sessions run on this computer at once. This module imports neither Streamlit nor
PyDiS.
"""
import json
import os
import select
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from common.job import JOB_ENV, RunFolder, kill
from common.paths import OPENDIS_DIR

CHILD = Path(__file__).resolve().with_name("console_child.py")
PREFIX = "opendis_console_"  # the sessions' folders (console_child.py removes its own)
COMMAND_S = 30.0  # a command's time limit
IDLE_S = 900.0    # the session ends after this long without a command
START_S = 120.0   # the time to load the run's objects
MAX_LIVE = 4      # sessions at once, on this computer
ENV_KEPT = ("PATH", "LANG", "LC_ALL", "LC_CTYPE")


@dataclass(frozen=True)
class Reply:
    """What a command printed, and whether the session lives on."""

    output: str
    alive: bool


class Busy(Exception):
    """Too many sessions run already."""


def live_sessions() -> int:
    """How many console processes run on this computer, of any app process."""
    marker, proc = str(CHILD), Path("/proc")
    if proc.is_dir():  # Linux, as on the Cloud
        count = 0
        for cmdline in proc.glob("[0-9]*/cmdline"):
            try:
                count += marker in cmdline.read_bytes().replace(b"\0", b" ").decode(errors="replace")
            except OSError:
                pass
        return count
    result = subprocess.run(["ps", "-A", "-o", "command="], capture_output=True, text=True,
                            check=False)
    return sum(marker in line for line in result.stdout.splitlines())


class Console:  # pylint: disable=too-many-instance-attributes
    """One session on one run: the child process, its two pipes and its folder."""

    def __init__(self, folder: RunFolder, python: str = sys.executable) -> None:
        if live_sessions() >= MAX_LIVE:
            raise Busy(f"{MAX_LIVE} consoles are open already, the most at once. "
                       "Please try again in a few minutes.")
        self.run_id = folder.run_id
        self.workdir = Path(tempfile.mkdtemp(prefix=PREFIX))
        shutil.copytree(folder.output, self.workdir / "output")
        shutil.copy2(folder.script, self.workdir / folder.script_name)
        commands_r, commands_w = os.pipe()
        replies_r, replies_w = os.pipe()
        env = {name: os.environ[name] for name in ENV_KEPT if name in os.environ}
        env.update(JOB_ENV)
        env.update({"HOME": str(self.workdir), "TMPDIR": str(self.workdir),
                    "MPLCONFIGDIR": str(self.workdir), "OPENDIS_DIR": str(OPENDIS_DIR),
                    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
        # The process outlives this method, so it is not opened in a with block.
        self.proc = subprocess.Popen(  # pylint: disable=consider-using-with
            [python, "-I", "-B", "-u", str(CHILD), str(commands_r), str(replies_w), str(IDLE_S),
             folder.final_json, folder.script_name],
            cwd=self.workdir, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, pass_fds=(commands_r, replies_w), start_new_session=True)
        os.close(commands_r)
        os.close(replies_w)
        self._commands = os.fdopen(commands_w, "w", encoding="utf-8")
        self._replies = os.fdopen(replies_r, encoding="utf-8")
        self.banner = self._read(START_S)  # what loading printed

    def _read(self, timeout: float) -> Reply:
        """The child's next reply; without one within ``timeout``, the session ends."""
        try:
            ready = bool(select.select([self._replies], [], [], timeout)[0])
            line = self._replies.readline() if ready else ""
        except (OSError, ValueError):
            ready, line = True, ""
        if not line:
            self.close()
            why = (f"The command ran longer than {timeout:g} s, so the session was stopped."
                   if not ready else "The session has ended.")
            return Reply(why + " The next command starts a new one.", False)
        data = json.loads(line)
        if not data["ok"]:
            self.close()
        return Reply(data["output"], data["ok"])

    def execute(self, source: str, timeout: float = COMMAND_S) -> Reply:
        """Run one command in the session."""
        try:
            self._commands.write(json.dumps({"source": source}) + "\n")
            self._commands.flush()
        except (OSError, ValueError):
            self.close()
            return Reply("The session has ended. The next command starts a new one.", False)
        return self._read(timeout)

    def alive(self) -> bool:
        """True while the session's process runs."""
        return self.proc.poll() is None

    def close(self) -> None:
        """End the session: its process, its pipes and its folder."""
        if self.proc.poll() is None:
            kill(self.proc)
        for pipe in (self._commands, self._replies):
            try:
                pipe.close()
            except OSError:
                pass
        shutil.rmtree(self.workdir, ignore_errors=True)

    def __del__(self) -> None:
        try:  # a session that its browser tab left
            self.close()
        except Exception:  # pylint: disable=broad-exception-caught
            pass
