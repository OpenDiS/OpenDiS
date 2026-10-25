"""Shows the bug in ExaDiS' node split, and fails until it is fixed (see bug_report/bug_report.md)."""

import sys, os
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python', 'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict, tier_tol
from pydis import DisNet, Cell
from pydis.util.plastic_strain import PlasticStrainAccum
import pyexadis
import pyexadis_base as ex

LBOX = 1000.0
B1, B2 = np.array([-1.0, 1.0, 1.0]), np.array([1.0, -1.0, 1.0])
STATE = {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
         "maxseg": 60.0, "minseg": 1.0, "rann": 3.0}
SWEPT_TOL = tier_tol(1e-11, 1e-9)            # relative to the area swept
OUTDIR = 'output'


def area_functional(nodes_array, segs_array):
    """area_functional: F = sum over segments of outer(0.5 r1 x r2, b) / V"""
    nodes_array, segs_array = np.asarray(nodes_array), np.asarray(segs_array)
    pos, ids, b = nodes_array[:, 2:5], segs_array[:, 0:2].astype(int), segs_array[:, 2:5]
    n = 0.5 * np.cross(pos[ids[:, 0]], pos[ids[:, 1]])
    return (n[:, :, None] * b[:, None, :]).sum(axis=0) / LBOX**3


def sym(P):
    return 0.5 * (P + P.T)


def init_junction():
    """init_junction: (names, rn, links) of the initial configuration

    Lines with b1 along [011] and b2 along [101] cross at the four-arm node c, in a cell
    that is not periodic, with their ends pinned 40 b away.
    """
    c = np.full(3, LBOX / 2)
    xi1, xi2 = np.array([0.0, 1.0, 1.0]) / np.sqrt(2), np.array([1.0, 0.0, 1.0]) / np.sqrt(2)
    pos = {"a1": c - 40.0*xi1, "a2": c + 40.0*xi1, "d1": c - 40.0*xi2, "d2": c + 40.0*xi2, "c": c}
    names = list(pos)
    PIN, FREE = ex.NodeConstraints.PINNED_NODE, ex.NodeConstraints.UNCONSTRAINED
    rn = np.array([np.append(pos[n], FREE if n == "c" else PIN) for n in names])
    links = []
    for n1, n2, b in [("a1", "c", B1), ("c", "a2", B1), ("d1", "c", B2), ("c", "d2", B2)]:
        pn = np.cross(b, pos[n2] - pos[n1])
        links.append(np.concatenate(([names.index(n1), names.index(n2)], b, pn / np.linalg.norm(pn))))
    return names, rn, np.array(links)


def run_exadis(names, rn, links, p0, p1):
    """run_exadis: ExaDiS' driver replays a record of the split of c, c to p0 and a new node to p1

    ExaDiS' topology never splits off one line, so the split is fed to execute_split this way.
    Returns (name -> position after, the Burgers vector from c to the new node, the strain
    ExaDiS added, the area swept).
    """
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[False] * 3)
    G = ex.ExaDisNet(cell, rn, links)
    serial = G.net._get_serial_network()
    assert len(serial._recycled_indices()) == 0     # else the new node's tag is not maxindex + 1
    tags = {n: tuple(int(x) for x in t) for n, t in zip(names, np.asarray(G.get_tags()))}
    tag_new = (tags["c"][0], serial._maxindex() + 1)
    F0 = area_functional(G.net.get_nodes_array(), G.net.get_segs_array())

    # initialize, zero the strain, split c with a2 and d2 going to the new node; the driver
    # clears its own folder, so the record is kept beside it
    os.makedirs(os.path.join(OUTDIR, 'split_record'), exist_ok=True)
    record = os.path.abspath(os.path.join(OUTDIR, 'split_record', 'oprec.split.exadis'))
    with open(record, 'w') as f:
        f.write("0 1.1 Serial 1 1 1\n")
        f.write("2" + " 0" * 19 + "\n")
        f.write("6 %d %d 2 %d %d %d %d %s %s %d %d\n" % (
            *tags["c"], *tags["a2"], *tags["d2"], " ".join("%.17g" % x for x in p0),
            " ".join("%.17g" % x for x in p1), *tag_new))

    system = pyexadis.System(G.net, ex.get_exadis_params(STATE))
    driver = pyexadis.Driver(system)
    driver.outputdir = os.path.abspath(os.path.join(OUTDIR, 'split_replay'))
    ctrl = pyexadis.Driver.Control()
    ctrl.loading = pyexadis.Driver.STRESS_CONTROL
    ctrl.appstress = np.zeros((3, 3))
    driver.oprec_replay(ctrl, record)

    nodes, segs = np.asarray(system.get_nodes_array()), np.asarray(system.get_segs_array())
    name = {t: n for n, t in tags.items()}
    name[tag_new] = "new"
    row = [name[(int(r[0]), int(r[1]))] for r in nodes]
    burgers = {(row[int(s[0])], row[int(s[1])]): s[2:5] for s in segs}
    b_new = burgers.get(("c", "new"), -burgers.get(("new", "c"), np.zeros(3)))
    added = np.array(system.get_plastic_strain()[0]).reshape(3, 3)
    return dict(zip(row, nodes[:, 2:5])), b_new, added, sym(F0 - area_functional(nodes, segs))


def run_pydis(names, rn, links, p0, p1):
    """run_pydis: the strain PyDiS accounts for the same split"""
    G = DisNet(cell=Cell(h=LBOX * np.eye(3), is_periodic=[False] * 3), rn=rn, links=links)
    tag = dict(zip(names, G.all_nodes_tags()))
    acc = PlasticStrainAccum()
    G.split_node(tag["c"], p0, p1, [tag["a2"], tag["d2"]], plast_acc=acc)
    return acc.dEp


def main():
    print("the plastic strain of a split that creates a segment, ExaDiS against PyDiS")
    names, rn, links = init_junction()
    pold = np.full(3, LBOX / 2)
    p0, p1 = pold + [-2.0, 1.0, 0.0], pold + [1.0, 2.0, 1.0]
    at, b_new, added, swept = run_exadis(names, rn, links, p0, p1)

    ok = report("ExaDiS' replay split c into c at p0 and a new node at p1",
                "new" in at and np.array_equal(at["c"], p0) and np.array_equal(at["new"], p1))
    ok &= report("they are joined by a segment with b1 + b2", np.allclose(b_new, B1 + B2))
    S = 0.5 * np.cross(p1 - pold, p0 - pold)
    ok &= report("p0, p1 and pold are not on one line: the new segment sweeps %.2f b^2"
                 % np.linalg.norm(S), np.linalg.norm(S) > 0.1)
    if not ok:
        return ok

    dEp_pydis = run_pydis(names, rn, links, p0, p1)
    scale = np.abs(swept).max()
    ok &= report_close("PyDiS' strain for the same split against the area swept, relative",
                       dEp_pydis / scale, swept / scale, atol=SWEPT_TOL)
    ok &= report_close("ExaDiS' strain for the split against the area swept, relative",
                       added / scale, swept / scale, atol=SWEPT_TOL)

    # the new segment, c -> new with b1 + b2, swept S; both of ExaDiS' calls see it at zero length
    missing = sym(np.outer(S, B1 + B2)) / LBOX**3
    print("    ExaDiS minus the area swept: %.3e in strain, %.1f%% of this split; the"
          % (np.abs(added - swept).max(), 100 * np.abs(added - swept).max() / scale))
    print("    new segment's triangle: %.3e; the two differ by %.1e"
          % (np.abs(missing).max(), np.abs(added - swept + missing).max()))
    return ok


if __name__ == "__main__":
    pyexadis.initialize()
    passed = verdict("test_split_plastic_strain_pydis_exadis", main())
    pyexadis.finalize()
    sys.exit(0 if passed else 1)
