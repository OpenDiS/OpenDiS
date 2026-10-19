"""Plastic strain swept by topological operations, not just time integration.

Every node move is accounted as segments sweeping area, P = outer(Sn, b)/V,
so plastic strain is not lost when collision, topology, or remesh move or
remove nodes. Prescribed moves, no force or mobility model:

0. move_node reproduces the time-integration sum bitwise, across a box face.
1. Inserting a node off a segment and removing it again sweeps nothing net.
2. Removing a two-arm node gives the same area whichever neighbour it
   collapses onto.
3. A merge of two connected nodes counts the segment between them once; a
   split counts the segment it creates between its products. ExaDiS gets
   the first twice and the second not at all; both are measured here.
4. A loop shrunk across a box face until remesh annihilates it sweeps its
   whole area, A b / V (Orowan), and F(before) - F(after) holds at every
   step, F the network's area functional. Without the topological
   accounting the same run leaks.

    python3 test_topological_plastic_strain.py
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[2]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]

import numpy as np
from framework.testing import report, report_close, verdict
from framework.disnet_manager import DisNetManager
from pydis import Cell, Remesh, DisNode
from pydis.util.plastic_strain import PlasticStrainAccum
import plast_strain_util as psu
import fcc_frank_read_util as fr

LBOX = 1000.0
B = np.array([1.0, 0.0, 0.0])


def cell():
    return Cell(h=LBOX*np.eye(3), is_periodic=[True, True, True])


def polygon(center, radius, n, z_wobble=0.0):
    """polygon: n vertices about center in the xy plane, optionally non-planar"""
    t = 2*np.pi*np.arange(n)/n
    pts = np.column_stack((radius*np.cos(t), radius*np.sin(t),
                           z_wobble*np.cos(2*t)))
    return cell().fold(center + pts)


def snapshot(G):
    return {tag: G.nodes(tag).R.copy() for tag in G.all_nodes_tags()}


def P_of(acc):
    return acc.dEp + acc.dWp


def tri_dyad(r1, r2, R, b):
    """tri_dyad: outer(0.5 (r2 - r1) x (R - r1), b) / V, unwrapped inputs"""
    return np.outer(0.5 * np.cross(r2 - r1, R - r1), b) / LBOX**3


def quad_dyad(r1, r2, r3, r4, b):
    """quad_dyad: segment (r3, r4) -> (r1, r2), unwrapped inputs, written out by hand"""
    return np.outer(0.5 * np.cross(r2 - r3, r1 - r4), b) / LBOX**3


def check_move_node_bitwise():
    """test 0: move_node equals the integration sum, near and across a box face"""
    ok = True
    for center, label in (([500.0, 500.0, 500.0], "interior"),
                          ([996.0, 500.0, 500.0], "across the x face")):
        G = fr.closed_loop(cell(), polygon(np.array(center), 80.0, 8), B)
        tag = list(G.all_nodes_tags())[0]
        new = G.cell.fold(G.nodes(tag).R + np.array([7.25, -3.5, 1.75]))
        G_ref = G.copy()
        old = snapshot(G_ref)
        G_ref.nodes(tag).R = new.copy()
        dEp, dWp = psu.plastic_strain_increment(G_ref, old)
        acc = PlasticStrainAccum()
        G.move_node(tag, new, plast_acc=acc)
        ok &= report_close("test 0, %s: move_node dEp == integration sum" % label,
                           acc.dEp, dEp, atol=0.0)
        ok &= report_close("test 0, %s: move_node dWp == integration sum" % label,
                           acc.dWp, dWp, atol=0.0)
    return ok


def check_insert_remove():
    """test 1: an off-segment insertion, then its removal, nets zero"""
    G = fr.closed_loop(cell(), polygon(np.array([500.0, 500.0, 500.0]), 80.0, 6), B)
    (t1, t2), e = next(iter(G.all_segments_dict().items()))
    r1, r2 = G.nodes(t1).R.copy(), G.nodes(t2).R.copy()
    R = 0.5*(r1 + r2) + np.array([0.0, 0.0, 9.0])
    acc = PlasticStrainAccum()
    new_tag = G.get_new_tag()
    G.insert_node_between(t1, t2, new_tag, R, plast_acc=acc)
    P_ins = P_of(acc).copy()
    want = tri_dyad(r1, r2, R, e.burg_vec_from(t1))
    scale = np.abs(want).max()
    ok = report_close("test 1: insertion sweeps triangle (r1, r2, R), relative",
                      P_ins / scale, want / scale, atol=1e-14)
    G.remove_two_arm_node(new_tag, plast_acc=acc)
    return ok & report_close("test 1: insert then remove nets zero, relative",
                             P_of(acc) / scale, np.zeros((3, 3)), atol=1e-14)


def bent_node_network():
    """bent_node_network: a loop with one node pushed off the line of its neighbours"""
    pos = polygon(np.array([500.0, 500.0, 500.0]), 80.0, 8)
    pos[3] += np.array([0.0, 0.0, 11.0])
    return fr.closed_loop(cell(), pos, B)


def check_route_independence():
    """test 2: collapsing a two-arm node onto either neighbour sweeps the same"""
    results = []
    for which in (0, 1):
        G = bent_node_network()
        tag = list(G.all_nodes_tags())[3]
        survivor = list(G.neighbors_tags(tag))[which]
        acc = PlasticStrainAccum()
        G.remove_two_arm_node(tag, survivor=survivor, plast_acc=acc)
        results.append(P_of(acc))
    scale = np.abs(results[0]).max()
    return report_close("test 2: removal onto either neighbour, relative",
                        results[1] / scale, results[0] / scale, atol=1e-14)


def merge_hand_values(G, dead, target, pos):
    """merge_hand_values: (correct, exadis) P for merging dead and target at pos"""
    def r(t):
        return G.nodes(t).R.copy()
    correct, exadis = np.zeros((3, 3)), np.zeros((3, 3))
    for (t1, t2), e in G.all_segments_dict().items():
        if not {t1, t2} & {dead, target}:
            continue
        b = e.burg_vec_from(t1)
        new = {t: (pos if t in (dead, target) else r(t)) for t in (t1, t2)}
        correct += quad_dyad(new[t1], new[t2], r(t1), r(t2), b)
        # exadis: each moved end against the other end's position before the merge
        for mover in {t1, t2} & {dead, target}:
            other = t2 if mover == t1 else t1
            moved = {mover: pos, other: r(other)}
            exadis += quad_dyad(moved[t1], moved[t2], r(t1), r(t2), b)
    return correct, exadis


def check_merge():
    """test 3a: merging two connected nodes at a point off their segment"""
    G = fr.closed_loop(cell(), polygon(np.array([500.0, 500.0, 500.0]), 80.0, 8), B)
    tags = list(G.all_nodes_tags())
    dead, target = tags[2], tags[3]
    pos = 0.5*(G.nodes(dead).R + G.nodes(target).R) + np.array([0.0, 0.0, 6.0])
    correct, exadis = merge_hand_values(G, dead, target, pos)
    acc = PlasticStrainAccum()
    G.merge_node(dead, target, position=pos, plast_acc=acc)
    scale = np.abs(correct).max()
    ok = report_close("test 3a: merge matches the hand sum, relative",
                      P_of(acc) / scale, correct / scale, atol=1e-14)
    excess = np.abs(exadis - correct).max() / scale
    print("    ExaDiS-style accounting differs from it by %.3f of its magnitude"
          % excess)
    return ok & report("test 3a: ExaDiS-style accounting double-counts the collapsing segment",
                       excess > 1e-3)


def three_arm_network():
    """three_arm_network: a free node c joined to three pinned nodes, b1 + b2 + b3 = 0"""
    c = np.array([500.0, 500.0, 500.0])
    arms = np.array([[60.0, 0.0, 3.0], [-30.0, 50.0, 4.0], [-30.0, -50.0, -4.0]])
    burgers = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, -1.0, 0.0]])
    rn = np.vstack((np.append(c, DisNode.Constraints.UNCONSTRAINED),
                    np.column_stack((c + arms, np.full(3, DisNode.Constraints.PINNED_NODE)))))
    links = np.array([np.concatenate(([0, k+1], burgers[k],
                                      fr.unit(np.cross(burgers[k], arms[k]))))
                      for k in range(3)])
    from pydis import DisNet
    return DisNet(cell=cell(), rn=rn, links=links)


def check_split():
    """test 3b: a split that creates a connecting segment counts it"""
    G = three_arm_network()
    tags = list(G.all_nodes_tags())
    c, k1, k2, k3 = tags
    pold = G.nodes(c).R.copy()
    pos1, pos2 = pold + np.array([3.0, 1.0, 0.5]), pold + np.array([-2.0, 0.5, -1.5])
    acc = PlasticStrainAccum()
    n1, n2 = G.split_node(c, pos1, pos2, [k2, k3], plast_acc=acc)
    b_conn = G.segments((n1, n2)).burg_vec_from(n1)
    first = G.segment_as_stored(n1, n2)[0]
    conn = (quad_dyad(pos1, pos2, pold, pold, b_conn) if first == n1
            else quad_dyad(pos2, pos1, pold, pold, -b_conn))
    arms = sum(quad_dyad(pos_new, G.nodes(k).R, pold, G.nodes(k).R,
                         G.segments((n, k)).burg_vec_from(n))
               for n, pos_new, ks in ((n1, pos1, [k1]), (n2, pos2, [k2, k3]))
               for k in ks)
    scale = np.abs(arms + conn).max()
    ok = report_close("test 3b: split matches the hand sum with the new segment, relative",
                      P_of(acc) / scale, (arms + conn) / scale, atol=1e-14)
    print("    the new segment is %.3f of the total; ExaDiS' execute_split omits it"
          % (np.abs(conn).max() / scale))
    return ok & report("test 3b: the connecting segment's share is not negligible",
                       np.abs(conn).max() / scale > 1e-3)


def shrink_step(G, center, lam):
    """shrink_step: every node to center + lam (r - center), folded, as one prescribed step"""
    for tag in list(G.all_nodes_tags()):
        r = G.cell.closest_image(Rref=center, R=G.nodes(tag).R)
        G.nodes(tag).R = G.cell.fold(center + lam*(r - center))


def shrink_to_nothing(topological):
    """shrink_to_nothing: (P total, F initial, worst identity error, steps, leak share)"""
    center = np.array([996.0, 500.0, 500.0])
    G = fr.closed_loop(cell(), polygon(center, 60.0, 6, z_wobble=5.0), B)
    state = {"maxseg": 40.0, "minseg": 10.0}
    remesh = Remesh(remesh_rule='LengthBased', state=state)
    acc = PlasticStrainAccum()
    state["plast_acc"] = acc if topological else None
    F0 = psu.network_area_dyad(G)
    total, topo_total, worst, steps = np.zeros((3, 3)), np.zeros((3, 3)), 0.0, 0
    while G.num_nodes() > 0 and steps < 200:
        old = snapshot(G)
        shrink_step(G, center, 0.9)
        dEp, dWp = psu.plastic_strain_increment(G, old)
        acc.zero()
        remesh.Remesh(DisNetManager(G), state)
        total += dEp + dWp + P_of(acc)
        topo_total += P_of(acc)
        F = psu.network_area_dyad(G) if G.num_nodes() > 0 else np.zeros((3, 3))
        worst = max(worst, np.abs(total - (F0 - F)).max() / np.abs(F0).max())
        steps += 1
    return total, F0, worst, steps, np.abs(topo_total).max() / np.abs(F0).max()


def check_orowan():
    """test 4: shrink a loop across a box face until remesh annihilates it"""
    total, F0, worst, steps, share = shrink_to_nothing(topological=True)
    scale = np.abs(F0).max()
    ok = report("test 4: loop annihilated, after %d steps" % steps, steps < 200)
    ok &= report_close("test 4: total swept equals the loop's A b / V, relative",
                       total / scale, F0 / scale, atol=1e-12)
    ok &= report("test 4: F(before) - F(after) at every step, worst relative %.1e" % worst,
                 worst < 1e-12)
    print("    topological operations carried %.3f of the total" % share)
    leak, _, worst_off, _, _ = shrink_to_nothing(topological=False)
    missing = np.abs(leak - F0).max() / scale
    print("    without them the run misses %.3f of A b / V" % missing)
    return ok & report("test 4: without topological accounting the same run leaks",
                       missing > 1e-3)


def main():
    print("plastic strain swept by topological operations")
    checks = [check_move_node_bitwise(), check_insert_remove(),
              check_route_independence(), check_merge(), check_split(),
              check_orowan()]
    return all(checks)


if __name__ == "__main__":
    passed = verdict("test_topological_plastic_strain", main())
    sys.exit(0 if passed else 1)
