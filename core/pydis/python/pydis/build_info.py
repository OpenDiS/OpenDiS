"""@package docstring
build_info: how the installed libpydis was built

Answers one question at present: whether calforce/SegSegForce.c was compiled
against the vendored pydis_log()/pydis_atan() rather than the platform libm,
which is what makes that kernel's output bitwise reproducible from one platform
to another. Every SYS whose name ends in _repro does so.

Reports the fact and stops there. What to do about it belongs to the caller;
the tests in tests/unit_tests/test1_node_force use it to choose between an
exact comparison and a tolerant one, and that choice is written where the
tolerances are, not here.

The library reports this itself, through SegSegForce_BitReproMath(), so the
answer is compiled from the same text under the same flags as the kernel it
describes and cannot be stale or disagree with it.
"""


def bitrepro_math() -> bool:
    """bitrepro_math: whether SegSegForce.c was built against the portable log/atan

    False when there is no compiled library to ask, and when the library
    predates the query function. Both mean "not known to be reproducible",
    which is the safe answer, and neither is worth raising over: a caller
    asking this wants to pick a tolerance, not to handle an exception.
    """
    try:
        pydis_lib = __import__('pydis_lib')
    except ImportError:
        return False
    query = getattr(pydis_lib, 'SegSegForce_BitReproMath', None)
    if query is None:
        return False
    return bool(query())


def build_description() -> str:
    """build_description: one line naming which of the two builds this is"""
    return ("libpydis: SegSegForce bitwise-reproducible math %s"
            % ("ON" if bitrepro_math() else "off"))
