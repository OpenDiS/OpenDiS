"""common.catalog: the examples, one input file each, ``examples/<name>.md``.

An input file has a front matter in TOML, between two lines of ``+++``, and a body in Markdown,
which step 1 shows. The front matter gives the example's label (the text of the drop-down menu),
its title, its script's file names, and the values of every step (``common/spec.py``). Adding an
example of a known kind of initial configuration is adding a file: the app's code does not
change. The folder is read on every call and a file again when it has changed, so an edited file
shows when the page is reloaded. A file with problems is not offered; ``problems()`` says why.
"""
import tomllib
from dataclasses import dataclass
from pathlib import Path

from common.paths import EXAMPLES, OPENDIS_DIR
from common.spec import PARTS, InputError, Spec

FENCE = "+++"
DEFAULT_ORDER = 100
META_KEYS = ("label", "title", "order", "script_name", "final_json", "tutorial_script",
             "tutorial_url", "run_hint")
REQUIRED = ("label", "title", "script_name", "final_json")


@dataclass(frozen=True)
class Example:  # pylint: disable=too-many-instance-attributes
    """An example: what its input file gives."""

    name: str                # the file's stem, e.g. ``frank_read_src``
    label: str               # the drop-down menu's text, e.g. ``Frank-Read src``
    title: str               # the centre panel's title
    order: int
    body: str                # the Markdown of the file: the problem, in words
    script_name: str         # the job script's file name
    final_json: str          # the final network, ``output/<final_json>``
    tutorial_script: str     # the tutorial's script, relative to the repository's root
    tutorial_url: str
    run_hint: str            # what the Run step says before a run
    default: Spec            # the selections that the file gives


def _split(text: str, name: str) -> tuple[str, str]:
    """The front matter and the body of an input file."""
    lines = text.splitlines()
    if not lines or lines[0].strip() != FENCE:
        raise InputError([f"{name}: the first line must be {FENCE}"])
    try:
        end = next(i for i in range(1, len(lines)) if lines[i].strip() == FENCE)
    except StopIteration:
        raise InputError([f"{name}: the front matter has no closing {FENCE}"]) from None
    return "\n".join(lines[1:end]), "\n".join(lines[end + 1:]).strip()


def load(path: Path) -> Example:
    """The example of an input file; raises ``InputError`` with every problem of the file."""
    front, body = _split(path.read_text(encoding="utf-8"), path.name)
    try:
        data = tomllib.loads(front)
    except tomllib.TOMLDecodeError as exc:
        raise InputError([f"front matter: {exc}"]) from None
    problems = [f"{key}: missing" for key in REQUIRED if key not in data]
    problems += [f"{key}: unknown key or table; the keys are {', '.join(META_KEYS + PARTS)}"
                 for key in data if key not in META_KEYS + PARTS]
    try:
        spec = Spec.from_tables(data)
    except InputError as exc:
        problems += exc.problems
    else:
        problems += [f"the selections: {message}" for message in spec.validate()]
    tutorial = data.get("tutorial_script", "")
    if tutorial and not (OPENDIS_DIR / tutorial).is_file():
        problems.append(f"tutorial_script: there is no file {tutorial} in {OPENDIS_DIR}")
    if problems:
        raise InputError(problems)
    return Example(name=path.stem, label=data["label"], title=data["title"],
                   order=int(data.get("order", DEFAULT_ORDER)), body=body,
                   script_name=data["script_name"], final_json=data["final_json"],
                   tutorial_script=tutorial, tutorial_url=data.get("tutorial_url", ""),
                   run_hint=data.get("run_hint", ""), default=spec)


_LOADED: dict[tuple[str, float], Example | InputError] = {}  # by the file's path and mtime


def _load_cached(path: Path) -> Example | InputError:
    """The example of a file, or what is wrong with it; kept while the file is unchanged."""
    key = (str(path), path.stat().st_mtime)
    if key not in _LOADED:
        try:
            _LOADED[key] = load(path)
        except InputError as exc:
            _LOADED[key] = exc
    return _LOADED[key]


def scan(root: Path = EXAMPLES) -> tuple[list[Example], dict[str, list[str]]]:
    """The examples in menu order (``order``, then name), and the problems of the other files."""
    examples, problems = [], {}
    for path in sorted(root.glob("*.md")):
        if path.stem.startswith("_") or path.stem == "README":
            continue
        loaded = _load_cached(path)
        if isinstance(loaded, InputError):
            problems[path.name] = loaded.problems
        else:
            examples.append(loaded)
    return sorted(examples, key=lambda e: (e.order, e.name)), problems


def examples(root: Path = EXAMPLES) -> list[Example]:
    """The examples, in menu order."""
    return scan(root)[0]


def problems(root: Path = EXAMPLES) -> dict[str, list[str]]:
    """What is wrong with each input file that has problems, by the file's name."""
    return scan(root)[1]
