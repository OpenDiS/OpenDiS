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
import contextlib
import io
import json
import os
import re
import resource
import select
import sys
import traceback

OUTPUT_CHARS = 20_000   # the most of one reply
LIMITS = {resource.RLIMIT_CPU: 600,           # seconds of CPU, in all
          resource.RLIMIT_FSIZE: 64 * 2**20,  # bytes in one file it writes
          resource.RLIMIT_NOFILE: 256}        # files open at once
MEMORY_ROOM = 512 * 2**20  # bytes for the commands, beyond what loading took (Linux)

# The guard rails, which an audit hook enforces once the run's objects are loaded. They refuse
# what looking into the objects does not need: other programs, threads and signals, the
# network, native libraries, the limits, a search of all objects (which would reach the hook)
BLOCKED = frozenset({
    "os.system", "os.exec", "os.spawn", "os.posix_spawn", "os.fork", "os.forkpty",
    "subprocess.Popen", "os.kill", "os.killpg", "signal.pthread_kill", "webbrowser.open",
    "_thread.start_new_thread",
    "socket.__new__", "socket.connect", "socket.bind", "socket.getaddrinfo",
    "socket.gethostbyname", "socket.gethostbyaddr", "socket.sendto", "socket.sendmsg",
    "ctypes.dlopen", "ctypes.dlsym", "ctypes.dlsym/handle", "ctypes.cdata",
    "ctypes.string_at", "ctypes.wstring_at",
    "resource.setrlimit", "resource.prlimit",
    "gc.get_objects", "gc.get_referrers", "gc.get_referents",
})
# ... changes to files outside the session's folder: each event's places of paths, and of
# dir_fd arguments, which make a path relative to a folder that the hook cannot see
CHANGES = {
    "os.remove": ((0,), (1,)), "os.rmdir": ((0,), (1,)), "os.mkdir": ((0,), (2,)),
    "os.rename": ((0, 1), (2, 3)), "os.replace": ((0, 1), (2, 3)), "os.link": ((0, 1), (2, 3)),
    "os.symlink": ((0, 1), (2,)), "os.chmod": ((0,), (2,)), "os.chown": ((0,), (3,)),
    "os.truncate": ((0,), ()), "os.utime": ((0,), (3,)), "os.mkfifo": ((0,), (2,)),
    "os.mknod": ((0,), (3,)), "os.chflags": ((0,), ()), "os.lchflags": ((0,), ()),
    "shutil.rmtree": ((0,), (1,)), "shutil.copyfile": ((1,), ()), "shutil.copymode": ((1,), ()),
    "shutil.copystat": ((1,), ()), "shutil.copytree": ((1,), ()), "shutil.chown": ((0,), ()),
    "shutil.make_archive": ((0,), ()), "shutil.unpack_archive": ((1,), ()),
}
SECRETS = (".streamlit", ".ssh", ".git", "secrets.toml")  # ... and reading these, or /proc
WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_APPEND | os.O_CREAT | os.O_TRUNC


def _limit() -> None:
    """Bound what one session can use: CPU time, the size of a file, the open files."""
    for kind, value in LIMITS.items():
        _, hard = resource.getrlimit(kind)
        top = value + 10 if kind == resource.RLIMIT_CPU else value  # SIGXCPU, then SIGKILL
        if hard != resource.RLIM_INFINITY:
            value, top = min(value, hard), min(top, hard)
        with contextlib.suppress(ValueError, OSError):
            resource.setrlimit(kind, (value, top))


def _limit_memory() -> None:
    """On Linux, bound the memory to what loading took, plus MEMORY_ROOM for the commands."""
    with contextlib.suppress(OSError, StopIteration, ValueError):
        with open("/proc/self/status", encoding="ascii") as status:
            size = next(int(line.split()[1]) * 1024 for line in status
                        if line.startswith("VmSize:"))
        resource.setrlimit(resource.RLIMIT_AS, (size + MEMORY_ROOM, size + MEMORY_ROOM))


def _guard(workdir: str) -> None:  # pylint: disable=too-many-locals
    """Install the guard rails: from now on the session refuses BLOCKED, CHANGES and SECRETS.

    The hook keeps its own references to what it calls, so that a command that replaces a
    function of a module changes nothing here. It is a guard against misuse, not a sandbox:
    Python's audit hooks cannot stop a determined expert.
    """
    import pathlib  # pylint: disable=import-outside-toplevel
    import posix  # pylint: disable=import-outside-toplevel,import-error
    import subprocess  # pylint: disable=import-outside-toplevel
    normpath = getattr(posix, "_path_normpath", None) or os.path.normpath
    getcwd, blocked, changes, secrets, write_flags = (posix.getcwd, BLOCKED, CHANGES, SECRETS,
                                                     WRITE_FLAGS)
    pure, is_a, size, refused = pathlib.PurePath, isinstance, len, PermissionError
    inside = normpath(workdir) + "/"

    def text(path) -> str | None:
        """The path as text, or None if it is not a plain path (which is refused)."""
        if is_a(path, str):
            return path
        if is_a(path, bytes):
            return path.decode("utf-8", "surrogateescape")
        if (is_a(path, pure) and type(path).__fspath__ is pure.__fspath__
                and type(path).__str__ is pure.__str__):
            return pure.__str__(path)
        return None

    def full(path: str) -> str:
        return normpath(path if path.startswith("/") else getcwd() + "/" + path)

    def refuse(what: str) -> None:
        raise refused(f"The console does not allow {what}: it is for looking into the "
                      "run's objects.")

    def check_change(event: str, args: tuple) -> None:
        paths, dir_fds = changes[event]
        if any(i < size(args) and args[i] is not None for i in dir_fds):
            refuse(f"{event} with dir_fd")
        for i in paths:
            if i < size(args) and not is_a(args[i], int):
                path = text(args[i])
                if path is None or not (full(path) + "/").startswith(inside):
                    refuse(f"{event} outside the session's folder")

    def hook(event: str, args: tuple) -> None:
        if event in blocked:
            refuse(event)
        elif event == "open" and not is_a(args[0], int):
            path = text(args[0])
            if path is None:
                refuse("opening a path of this kind")
            where = full(path)
            writing = (args[2] or 0) & write_flags or any(c in (args[1] or "") for c in "wax+")
            if writing and not (where + "/").startswith(inside):
                refuse(f"writing to {where}")
            if (any(part in secrets for part in where.split("/"))
                    or (where.startswith("/proc/") and not where.startswith("/proc/self/"))):
                refuse(f"reading {where}")
        elif event in changes:
            check_change(event, args)

    def no_fork_exec(*args, **kwargs):  # pylint: disable=unused-argument
        refuse("subprocess")

    # subprocess's own fork_exec reports no audit event, so it is replaced first
    subprocess._fork_exec = no_fork_exec  # pylint: disable=protected-access
    with contextlib.suppress(ImportError):
        import _posixsubprocess  # pylint: disable=import-outside-toplevel,import-error
        _posixsubprocess.fork_exec = no_fork_exec
    sys.addaudithook(hook)


def _snapshot(step: int):
    """A snapshot of the run, ``output/disnet_<step>.json``, as a DisNetManager."""
    from framework.disnet_manager import DisNetManager  # pylint: disable=import-outside-toplevel,import-error
    from pydis import DisNet  # pylint: disable=import-outside-toplevel,import-error
    manager = DisNetManager(DisNet())
    manager.read_json(os.path.join("output", f"disnet_{step}.json"))
    return manager


def load(script: str, final_json: str) -> dict:
    """The names of the run's script, once it ran with the final network in place of the run."""
    # as the scripts add them, ExaDiS's included
    for path in ("python", "lib", "core/pydis/python", "core/exadis/python"):
        full = os.path.abspath(os.path.join(os.environ["OPENDIS_DIR"], path))
        if full not in sys.path:
            sys.path.append(full)
    import pydis  # pylint: disable=import-outside-toplevel,import-error

    def load_final(self, DM, state, *args, **kwargs):  # pylint: disable=invalid-name,unused-argument
        """In place of the run: the run's final network."""
        DM.read_json(os.path.join("output", final_json))

    pydis.SimulateNetwork.run = load_final
    with open(script, encoding="utf-8") as file:
        source = file.read()
    if "import pyexadis" in source:  # an ExaDiS script: its run is replaced too, and its end
        import pyexadis  # pylint: disable=import-outside-toplevel,import-error
        import pyexadis_base  # pylint: disable=import-outside-toplevel,import-error
        pyexadis_base.SimulateNetwork.run = load_final
        pyexadis.finalize = lambda: None  # must not end ExaDiS: the session uses its objects
    namespace = {"__name__": "__main__", "__file__": script}
    exec(compile(source, script, "exec"), namespace)  # pylint: disable=exec-used
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
            shown = traceback.TracebackException(type(exc), exc, None if isinstance(
                exc, SyntaxError) else exc.__traceback__)
            shown.stack = traceback.StackSummary.from_list(  # this module's frames, the guard's
                [frame for frame in shown.stack if frame.filename != __file__])  # included
            buffer.write("".join(shown.format()))
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
        if ok:
            _limit_memory()
            _guard(os.getcwd())  # the session's folder, which the app removes at its end
        _reply(replies, startup.getvalue(), ok)
        while ok and select.select([commands], [], [], idle_s)[0]:
            line = commands.readline()
            if not line:  # the app closed the pipe
                break
            _reply(replies, execute(json.loads(line)["source"], namespace), True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
