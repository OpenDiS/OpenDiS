"""Checks that ExaDiS and PyDiS count the area swept in a merge of two connected nodes once."""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict, tier_tol
from framework.disnet_manager import DisNetManager
from pydis import DisNet, Cell
from pydis.util.plastic_strain import PlasticStrainAccum
import pyexadis
import pyexadis_base as ex

LBOX, B = 1000.0, np.array([1.0, 0.0, 0.0])
STATE = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
         "maxseg": 60.0, "minseg": 1.0, "rann": 3.0}
SWEPT_TOL = tier_tol(1e-11, 1e-9)            # relative to the area swept


def snapshot(N, name):
    """snapshot: (area functional F, name -> position) of the current network"""
    data = N.export_data()
    pos = np.asarray(data["nodes"]["positions"])
    ids = np.asarray(data["segs"]["nodeids"], dtype=int).reshape(-1, 2)
    b = np.asarray(data["segs"]["burgers"]).reshape(-1, 3)
    n = 0.5 * np.cross(pos[ids[:, 0]], pos[ids[:, 1]])
    F = (n[:, :, None] * b[:, None, :]).sum(axis=0) / N.cell.volume()
    return F, {name[tuple(t)]: p for t, p in zip(np.asarray(data["nodes"]["tags"]), pos)}


def sym(P):
    return 0.5 * (P + P.T)


def init_loop():
    """init_loop: (names, rn, links) of the initial configuration

    A closed loop, b along x, in a cell that is not periodic, with only s and m free;
    m has a 2 b arm to s in the xy-plane, and a 40 b arm to k3 in a plane tilted 45 degrees.
    """
    c, a, t = np.full(3, LBOX / 2), np.radians(30.0), np.radians(45.0)
    PIN, FREE = ex.NodeConstraints.PINNED_NODE, ex.NodeConstraints.UNCONSTRAINED
    loop = {"k1": (c + [-20.0, 40.0, 0.0], PIN),
            "s": (c + 2.0*np.array([np.cos(a), np.sin(a), 0.0]), FREE),
            "m": (c, FREE),
            "k3": (c + 40.0*np.array([np.cos(a), np.sin(a)*np.cos(t), np.sin(a)*np.sin(t)]), PIN),
            "k2": (c + [-40.0, -20.0, 60.0], PIN)}
    names = list(loop)
    rn = np.array([np.append(p, con) for p, con in loop.values()])
    links = []
    for k in range(len(rn)):
        k2 = (k + 1) % len(rn)
        pn = np.cross(B, rn[k2, :3] - rn[k, :3])
        links.append(np.concatenate(([k, k2], B, pn / np.linalg.norm(pn))))
    return names, rn, np.array(links)


def run_exadis(names, rn, links):
    """run_exadis: one call of ExaDiS' collision

    Returns (name -> position before, and after, the strain ExaDiS added, the area swept).
    """
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[False] * 3)
    N = DisNetManager(ex.ExaDisNet(cell, rn, links))
    name = {tuple(t): n for t, n in zip(np.array(N.get_disnet(ex.ExaDisNet).get_tags()), names)}
    F0, before = snapshot(N, name)
    ex.Collision(collision_mode='Retroactive', state=STATE).HandleCol(N, dict(STATE))
    F1, after = snapshot(N, name)
    added = np.array(N.get_disnet(ex.ExaDisNet).net.get_plastic_strain()[0]).reshape(3, 3)
    return before, after, added, sym(F0 - F1)


def run_pydis(names, rn, links, position):
    """run_pydis: the strain PyDiS accounts for merging s into m at position"""
    G = DisNet(cell=Cell(h=LBOX * np.eye(3), is_periodic=[False] * 3), rn=rn, links=links)
    tag = dict(zip(names, G.all_nodes_tags()))
    acc = PlasticStrainAccum()
    G.merge_node(tag["s"], tag["m"], position=position, plast_acc=acc)
    return acc.dEp


def main():
    print("the plastic strain of a merge of two connected nodes, ExaDiS against PyDiS")
    names, rn, links = init_loop()
    before, after, added, swept = run_exadis(names, rn, links)

    gone = sorted(n for n in before if n not in after)
    moved = sorted(n for n in after if np.abs(after[n] - before[n]).max() > 0.0)
    ok = report("ExaDiS' collision merged s and m, joined by a segment",
                {*gone, *moved} == {"m", "s"} and len(gone) == 1, suffix=' gone %s' % gone)
    if not ok:
        return ok
    merged = after["m"] if "m" in after else after["s"]
    line = (before["s"] - before["m"]) / np.linalg.norm(before["s"] - before["m"])
    off = np.linalg.norm(np.cross(merged - before["m"], line))
    ok &= report("the merge point lies off the line s-m, by %.3f b" % off, off > 0.1)

    dEp_pydis = run_pydis(names, rn, links, merged)
    scale = np.abs(swept).max()
    ok &= report_close("PyDiS' strain for the same merge against the area swept, relative",
                       dEp_pydis / scale, swept / scale, atol=SWEPT_TOL)
    # ExaDiS before 32082bb counted the triangle swept by the segment s-m twice
    ok &= report_close("ExaDiS' strain for the merge against the area swept, relative",
                       added / scale, swept / scale, atol=SWEPT_TOL)
    return ok


if __name__ == "__main__":
    pyexadis.initialize()
    passed = verdict("test_merge_plastic_strain_pydis_exadis", main())
    pyexadis.finalize()
    sys.exit(0 if passed else 1)
