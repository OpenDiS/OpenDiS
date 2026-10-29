"""tests.test_job: the job runner, with small stand-in scripts (README.md, Tests)."""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import job
from common.testing import check, summary

FINAL = "final.json"
GOOD = """\
print("init: ready")
for istep in range(0, 30, 10):
    print("step = %d, nodes = 7, dt = 1.0e-08, time = 1.0e-08, elapsed = 0.0 sec" % istep)
import os
os.makedirs("output", exist_ok=True)
open("output/final.json", "w").write("{}")
"""
FAILS = 'print("step = 0, nodes = 7")\nraise SystemExit(3)\n'
NO_FINAL = 'print("step = 0, nodes = 7")\n'
SLOW = 'import time\nprint("step = 0, nodes = 7", flush=True)\ntime.sleep(60)\n'


class Stop(Exception):
    """What Streamlit raises in a script run that it stops."""


def prepared(runs: Path, name: str, text: str) -> job.RunFolder:
    """A prepared run folder for a stand-in script."""
    folder = job.RunFolder.for_script(name, "script.py", FINAL, text, runs=runs)
    folder.prepare(text)
    return folder


def main() -> int:  # pylint: disable=too-many-locals
    """Run every check."""
    with tempfile.TemporaryDirectory() as tmp:
        runs = Path(tmp)
        good = job.run(prepared(runs, "good", GOOD))
        log = next((runs / "good").glob("*/log.txt")).read_text().splitlines()
        fails = job.run(prepared(runs, "fails", FAILS))
        no_final = job.run(prepared(runs, "no_final", NO_FINAL))
        results = [
            check("a good run ends done", good.ok and good.exit_code == 0, good.summary()),
            check("its lines are kept, and the last step is read", good.lines[0] == "init: ready"
                  and good.step == 20, (good.lines[:2], good.step)),
            check("the step is read from PyDiS's line with or without the nodes",
                  job.Status(lines=["step = 10 dt = 1.000000e-08"]).step == 10
                  and job.Status(lines=["step = 20, nodes = 7, dt = 1.0e-08"]).step == 20),
            check("its log ends with the app's line", log[-1].startswith("[app] done after")),
            check("a script that exits with 3 fails", fails.state == "failed"
                  and "code 3" in fails.message, fails.summary()),
            check("a script that writes no final network fails", no_final.state == "failed"
                  and FINAL in no_final.message, no_final.summary()),
        ]

        # a stop from the caller (Streamlit's way), then a time limit: the job is killed
        calls = []

        def stop_soon(status: job.Status) -> None:
            calls.append(status.state)
            if status.lines:
                raise Stop

        stopped = prepared(runs, "stopped", SLOW)
        try:
            job.run(stopped, stop_soon)
            raised = False
        except Stop:
            raised = True
        log = stopped.log.read_text().splitlines()
        timed = job.run(prepared(runs, "timed", SLOW),
                        limits=job.Limits(first_line_s=10, silence_s=1.0, total_s=10))
        results += [
            check("a stop from the caller goes on, after the job is killed", raised and calls),
            check("the log says the app stopped it", log[-1].startswith("[app] stopped"),
                  log[-1:]),
            check("a silent job is stopped at the limit", timed.state == "timeout"
                  and timed.killed and "no new output" in timed.message, timed.summary()),
            check("no job is left running", not job.running()),
        ]

        # a job that an earlier session left, named in its job.pid
        leftover = prepared(runs, "leftover", SLOW)
        proc = subprocess.Popen([sys.executable, str(leftover.script)],  # pylint: disable=consider-using-with
                                cwd=leftover.path, stdout=subprocess.DEVNULL,
                                start_new_session=True)
        (leftover.path / job.PID).write_text(f"{proc.pid}\n{leftover.script}\n", encoding="utf-8")
        killed = job.kill_leftovers(runs)
        proc.wait(timeout=5)
        results += [
            check("a leftover job is killed", killed == [proc.pid] and proc.poll() is not None,
                  killed),
            check("its job.pid is removed", not (leftover.path / job.PID).exists()),
        ]

        # the run folders: the same script, the same folder; prepare() clears the output
        a = job.RunFolder.for_script("x", "s.py", FINAL, "print(1)\n", runs=runs)
        b = job.RunFolder.for_script("x", "s.py", FINAL, "print(1)\n", runs=runs)
        c = job.RunFolder.for_script("x", "s.py", FINAL, "print(2)\n", runs=runs)
        a.prepare("print(1)\n")
        a.output.mkdir()
        (a.output / "disnet_10.json").write_text("{}")
        a.prepare("print(1)\n")
        results += [
            check("the same script, the same run folder; another, another",
                  a.path == b.path and a.path != c.path),
            check("prepare() clears the output of an earlier run",
                  not list(a.output.glob("*"))),
        ]
    return summary(results)


if __name__ == "__main__":
    sys.exit(main())
