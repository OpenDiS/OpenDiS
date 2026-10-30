"""tests.test_console: the Python console, on 20-step runs of both examples (README.md, Tests)."""
import json
import os
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

from common import catalog, codegen, console, job
from common.testing import check, summary

SECRET = "OPENDIS_TEST_SECRET"  # in the app's environment, not in the console's


def short_run(example: catalog.Example, runs: Path) -> job.RunFolder:
    """A 20-step run of the example's own values."""
    spec = example.default.with_part("run", replace(example.default.run, max_step=20))
    script = codegen.full_script(spec, example.final_json)
    folder = job.RunFolder.for_script(example.name, example.script_name, example.final_json,
                                      script, runs=runs)
    folder.prepare(script, spec.to_json())
    job.run(folder)
    return folder


def main() -> int:  # pylint: disable=too-many-locals
    """Run every check."""
    os.environ[SECRET] = "not for the console"
    before = console.live_sessions()
    results, runs = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for example in catalog.examples():
            folder = short_run(example, Path(tmp))
            runs.append((example, folder))
            nodes = len(json.loads(folder.results.read_text())["nodes"]["tags"])
            session = console.Console(folder)
            got = session.execute("G.num_nodes(), sim.max_step, net.is_sane()")
            session.close()
            results.append(check(f"{example.label}: G, sim and net are the run's",
                                 session.banner.alive and got.output == f"({nodes}, 20, True)\n",
                                 (session.banner.output[-300:], got.output)))

        # the details, on the last example's run
        example, folder = runs[-1]
        session = console.Console(folder)
        replies = {name: session.execute(source) for name, source in (
            ("state", 'state["mu"]'),
            ("multi", "def twice(x):\n    return 2 * x\n\ntwice(21)"),
            ("print", "print('hello')"),
            ("error", "1/0"),
            ("kept", "twice(5)"),
            ("snapshot", "snapshot_steps, snapshot(10).get_disnet().num_nodes() > 0"),
            ("write", "open('output/mine.txt', 'w').write('x')"),
            ("env", f"import os; os.environ.get({SECRET!r})"),
            ("input", "input()"),
        )}
        in_copy = (session.workdir / "output" / "mine.txt").exists()
        in_run = (folder.output / "mine.txt").exists()
        limit, console.MAX_LIVE = console.MAX_LIVE, console.live_sessions()
        try:
            console.Console(folder).close()
            busy = False
        except console.Busy:
            busy = True
        console.MAX_LIVE = limit
        slow = session.execute("while True: pass", timeout=2.0)
        after = session.execute("1")
        mu = repr(float(example.default.material.mu))
        results += [
            check("state is the script's", replies["state"].output == mu + "\n",
                  replies["state"].output),
            check("several lines run, and the last expression shows its value",
                  replies["multi"].output == "42\n", replies["multi"].output),
            check("print shows", replies["print"].output == "hello\n", replies["print"].output),
            check("an error shows its traceback, from the user's code only, and the session "
                  "lives on", replies["error"].alive and "ZeroDivisionError" in
                  replies["error"].output and "console_child" not in replies["error"].output,
                  replies["error"].output),
            check("names stay from one command to the next", replies["kept"].output == "10\n",
                  replies["kept"].output),
            check("snapshot() loads a snapshot of the run",
                  replies["snapshot"].output == "([0, 10], True)\n", replies["snapshot"].output),
            check("what the session writes goes to its own copy, not to the run's folder",
                  in_copy and not in_run, (in_copy, in_run, replies["write"].output)),
            check("the app's environment is not the session's", replies["env"].output == "",
                  replies["env"].output),
            check("input() reads nothing", "EOFError" in replies["input"].output,
                  replies["input"].output),
            check("no more than MAX_LIVE sessions start", busy),
            check("a command that runs too long ends the session",
                  not slow.alive and "longer than 2 s" in slow.output and not session.alive(),
                  slow.output),
            check("the session's folder is removed", not session.workdir.exists()),
            check("then a command says that the session has ended",
                  not after.alive and "has ended" in after.output, after.output),
        ]
    results.append(check("no session is left running", console.live_sessions() == before))
    return summary(results)


if __name__ == "__main__":
    sys.exit(main())
