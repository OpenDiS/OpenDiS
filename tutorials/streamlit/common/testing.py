"""common.testing: the PASS / FAIL helpers shared by the plain-script tests.

Each test is a script, ``python -m tests.test_<name>`` from ``tutorials/streamlit/``, that
prints one line per check and returns a non-zero exit code on failure.
"""

RED, GREEN, RESET = "\033[91m", "\033[92m", "\033[0m"


def check(name: str, ok: bool, detail: object = "") -> bool:
    """Print one result line and return ``ok``."""
    print(f"{name}: {GREEN}PASS{RESET}" if ok else f"{name}: {RED}FAIL{RESET} — {detail}")
    return ok


def summary(results: list[bool]) -> int:
    """Print the overall result; 0 when every check passed, else 1."""
    passed = all(results)
    print(f"\n{GREEN}ALL {len(results)} PASSED{RESET}" if passed
          else f"\n{RED}{results.count(False)} OF {len(results)} FAILED{RESET}")
    return 0 if passed else 1


def has(messages: list[str], *words: str) -> bool:
    """True if one message holds every word."""
    return any(all(word in message for word in words) for message in messages)
