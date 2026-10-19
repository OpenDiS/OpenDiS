"""FCC_slip_sys: slip-system tables, lookups, and Schmid factors.

Pins the port in pydis/util/fcc_slip_sys.py, including the five defects
fixed relative to the copy in this folder:

  1. get_slip_Snb returns the dominant system (the copy always gave -1)
  2. the junction placeholder SchmidFactor[12] = 1 never wins that search
  3. the (plane, Burgers) -> sysID table is built once, in __init__
  4. sysID_to_nb has no side effects on the object
  5. a table inconsistency raises ValueError instead of exiting

    python3 test_fcc_slip_sys.py
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict
from pydis.util.fcc_slip_sys import FCC_slip_sys, schmid_factor, JUNCTION_ID

A3_SCHMID_001 = 1.0 / np.sqrt(6.0)


def check_round_trip(s):
    """every system, with either sign on n and b, maps back to itself"""
    bad = []
    for i in range(12):
        n, b = (v.ravel() for v in s.sysID_to_nb(i))
        for sn in (1, -1):
            for sb in (1, -1):
                if s.nb_to_sysID(sn*n, sb*b)[0] != i:
                    bad.append((i, sn, sb))
    return report("round trip sysID -> (n, b) -> sysID, 12 systems x 4 signs",
                  not bad, suffix='' if not bad else ' %s' % bad)


def check_tables_agree(s):
    """slipSystem and sysID_to_plane_burgers describe the same 12 systems"""
    ok = True
    for i in range(12):
        n, b = (v.ravel() for v in s.sysID_to_nb(i))
        ok &= np.array_equal(n, np.asarray(s.slipSystem['slipPlane'][i], float))
        ok &= np.array_equal(b, s.slipSystem['burgVect'][i])
        ok &= s.SB_index_names[i] == s.slipSystem['SB_index'][i]
    return report("slipSystem and sysID_to_plane_burgers agree", ok)


def check_a3(s):
    """A3 is (-1,1,1)[101]/2 and has Schmid factor 1/sqrt(6) on [001]"""
    i = s.name_to_sysID('A3')
    n, b = (v.ravel() for v in s.sysID_to_nb(i))
    ok = report("name_to_sysID('A3') == 1", i == 1)
    ok &= report("A3 plane and Burgers vector",
                 np.array_equal(n, [-1, 1, 1]) and np.allclose(b, [0.5, 0, 0.5]))
    ok &= report_close("A3 Schmid factor on [001]",
                       schmid_factor([0, 0, 1], n, b), A3_SCHMID_001, atol=1e-15)
    return ok


def check_dominant_001(s):
    """[001] activates 8 systems equally; the lowest index, A2, is returned"""
    slip, dom = s.get_slip_Snb([0, 0, 1])
    f = np.array([slip['SchmidFactor'][i] for i in range(12)])
    ok = report("get_slip_Snb([001]) returns a system in 0..11, A2 on the tie",
                dom == 0, suffix=' (got %r)' % dom)
    ok &= report_close("schmid_factor on [001]", s.schmid_factor,
                       A3_SCHMID_001, atol=1e-15)
    ok &= report("8 systems at the maximum, 4 at zero",
                 np.sum(np.isclose(f, A3_SCHMID_001)) == 8
                 and np.sum(f == 0.0) == 4)
    return ok


def check_junction_excluded(s):
    """the placeholder SchmidFactor[12] = 1 must not be picked as dominant"""
    slip, dom = s.get_slip_Snb([0, 0, 1])
    return report("junction placeholder kept out of the dominant search",
                  slip['SchmidFactor'][JUNCTION_ID] == 1 and dom != JUNCTION_ID
                  and s.schmid_factor < 1)


def check_standard_triangle(s):
    """axes inside the [001]-[011]-[-111] triangle single-slip on B4"""
    dom = {tuple(e): s.get_slip_Snb(e)[1]
           for e in ([-1, 2, 3], [-1, 3, 5], [-2, 3, 8])}
    got = {e: s.SB_index_names.get(d, d) for e, d in dom.items()}
    return report("standard-triangle axes give primary system B4",
                  set(got.values()) == {'B4'}, suffix=' %s' % got)


def check_junction_lookup(s):
    """non-slip pairs map to 12, and 12 maps back to (None, None)"""
    ok = report("a non-slip pair maps to 12",
                s.nb_to_sysID([1, 0, 0], [1, 0, 0])[0] == JUNCTION_ID)
    ok &= report("sysID_to_nb(12) is (None, None)",
                 s.sysID_to_nb(JUNCTION_ID) == (None, None))
    try:
        s.sysID_to_nb(13)
        raised = False
    except ValueError:
        raised = True
    return ok & report("sysID_to_nb(13) raises ValueError", raised)


def check_no_side_effects():
    """sysID_to_nb leaves the object alone; the lookup table exists up front"""
    s = FCC_slip_sys()
    ok = report("nb2sysID_tab built in __init__", hasattr(s, 'nb2sysID_tab'))
    s.sysID_to_nb(1)
    return ok & report("sysID_to_nb sets no attributes",
                       not hasattr(s, 'slip_plane') and not hasattr(s, 'burg_vec'))


def check_mismatch_raises():
    """a corrupted table raises ValueError rather than exiting the process"""
    s = FCC_slip_sys()
    s.slipSystem['burgVect'][3] = s.slipSystem['burgVect'][4]
    try:
        s.get_slip_Snb([0, 0, 1])
        raised = False
    except ValueError:
        raised = True
    return report("table mismatch raises ValueError", raised)


def check_exported():
    """FCC_slip_sys is importable from the pydis package itself"""
    import pydis
    return report("pydis.FCC_slip_sys is pydis.util.fcc_slip_sys.FCC_slip_sys",
                  getattr(pydis, 'FCC_slip_sys', None) is FCC_slip_sys)


def main():
    print("FCC_slip_sys, pydis.util")
    s = FCC_slip_sys()
    checks = [check_exported(), check_round_trip(s), check_tables_agree(s), check_a3(s),
              check_dominant_001(s), check_junction_excluded(s),
              check_standard_triangle(s), check_junction_lookup(s),
              check_no_side_effects(), check_mismatch_raises()]
    return all(checks)


if __name__ == "__main__":
    passed = verdict("test_fcc_slip_sys", main())
    sys.exit(0 if passed else 1)
