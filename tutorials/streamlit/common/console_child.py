"""common.console_child: the Python console's own process, which ``common/console.py`` starts.

It runs the run's script as ``python -i`` would, except that ``SimulateNetwork.run`` loads the
run's final network instead of simulating: then ``net``, ``sim``, ``state`` and ``G`` are the
run's. It then executes the commands that come on its command pipe, one JSON line each, and
answers each on its reply pipe with what the command printed, as Python's prompt shows it. It
ends after a time without a command, or when the app closes the pipe. It imports neither
Streamlit nor the app's modules: ``python -I`` runs it.

Usage: console_child.py <command fd> <reply fd> <idle seconds> <final network> <script>
"""
import ast
import atexit
import contextlib
import io
import json
import os
import re
import resource
import select
import shutil
import sys
import traceback

OUTPUT_CHARS = 20_000   # the most of one reply
LIMITS = {resource.RLIMIT_CPU: 600,           # seconds of CPU, in all
          resource.RLIMIT_FSIZE: 64 * 2**20,  # bytes in one file it writes
          resource.RLIMIT_NOFILE: 256}        # files open at once
PREFIX = "opendis_console_"  # the session's folder, which it removes at its end


def _limit() -> None:
    """Bound what one session can use: CPU time, the size of a file, the open files."""
    for kind, value in LIMITS.items():
        _, hard = resource.getrlimit(kind)
        top = value + 10 if kind == resource.RLIMIT_CPU else value  # SIGXCPU, then SIGKILL
        if hard != resource.RLIM_INFINITY:
            value, top = min(value, hard), min(top, hard)
        with contextlib.suppress(ValueError, OSError):
            resource.setrlimit(kind, (value, top))


def _snapshot(step: int):
    """A snapshot of the run, ``output/disnet_<step>.json``, as a DisNetManager."""
    from framework.disnet_manager import DisNetManager  # pylint: disable=import-outside-toplevel,import-error
    from pydis import DisNet  # pylint: disable=import-outside-toplevel,import-error
    manager = DisNetManager(DisNet())
    manager.read_json(os.path.join("output", f"disnet_{step}.json"))
    return manager


def load(script: str, final_json: str) -> dict:
    """The names of the run's script, once it ran with the final network in place of the run."""
    for path in ("python", "lib", "core/pydis/python"):  # as the scripts add them
        full = os.path.abspath(os.path.join(os.environ["OPENDIS_DIR"], path))
        if full not in sys.path:
            sys.path.append(full)
    import pydis  # pylint: disable=import-outside-toplevel,import-error

    def load_final(self, DM, state, *args, **kwargs):  # pylint: disable=invalid-name,unused-argument
        """In place of the run: the run's final network."""
        DM.read_json(os.path.join("output", final_json))

    pydis.SimulateNetwork.run = load_final
    namespace = {"__name__": "__main__", "__file__": script}
    with open(script, encoding="utf-8") as file:
        exec(compile(file.read(), script, "exec"), namespace)  # pylint: disable=exec-used
    namespace["snapshot"] = _snapshot
    namespace["snapshot_steps"] = sorted(
        int(match.group(1)) for match in (re.fullmatch(r"disnet_(\d+)\.json", name)
                                          for name in os.listdir("output")) if match)
    return namespace


def _run(source: str, namespace: dict) -> None:
    """Execute the source; if it ends with an expression, print its value, as the prompt does."""
    tree = ast.parse(source, "<console>")
    last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
    exec(compile(tree, "<console>", "exec"), namespace)  # pylint: disable=exec-used
    if last is not None:
        value = eval(compile(ast.Expression(last.value), "<console>", "eval"),  # pylint: disable=eval-used
                     namespace)
        if value is not None:
            namespace["_"] = value
            print(repr(value))


def execute(source: str, namespace: dict) -> str:
    """Run one command; returns what it printed, and its traceback if it failed."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        try:
            _run(source, namespace)
        except SystemExit:
            print("exit() does not end this session: ↺ restart starts a new one.")
        except BaseException as exc:  # pylint: disable=broad-exception-caught
            tb = None if isinstance(exc, SyntaxError) else exc.__traceback__
            while tb is not None and tb.tb_frame.f_code.co_filename == __file__:
                tb = tb.tb_next  # this module's frames are not the user's
            traceback.print_exception(type(exc), exc, tb, file=buffer)
    text = buffer.getvalue()
    if len(text) > OUTPUT_CHARS:
        text = text[:OUTPUT_CHARS] + f"\n... ({len(text) - OUTPUT_CHARS} more characters)"
    return text


def _reply(replies, output: str, ok: bool) -> None:
    replies.write(json.dumps({"output": output, "ok": ok}) + "\n")
    replies.flush()


def main() -> int:
    """Load the run's objects and say so, then answer each command, until the end."""
    commands_fd, replies_fd = int(sys.argv[1]), int(sys.argv[2])
    idle_s, final_json, script = float(sys.argv[3]), sys.argv[4], sys.argv[5]
    workdir = os.getcwd()
    if os.path.basename(workdir).startswith(PREFIX):
        atexit.register(shutil.rmtree, workdir, True)
    _limit()
    sys.stdin = io.StringIO("")  # input() reads nothing: the commands come on their own pipe
    with os.fdopen(commands_fd, encoding="utf-8") as commands, \
            os.fdopen(replies_fd, "w", encoding="utf-8") as replies:
        startup = io.StringIO()
        with contextlib.redirect_stdout(startup), contextlib.redirect_stderr(startup):
            try:
                namespace, ok = load(script, final_json), True
            except BaseException:  # pylint: disable=broad-exception-caught
                traceback.print_exc()
                namespace, ok = {}, False
        _reply(replies, startup.getvalue(), ok)
        while ok and select.select([commands], [], [], idle_s)[0]:
            line = commands.readline()
            if not line:  # the app closed the pipe
                break
            _reply(replies, execute(json.loads(line)["source"], namespace), True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
