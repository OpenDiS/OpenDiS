"""Topological operations on the 10-loop configuration, in PyDiS.

Two operations are applied to the loaded network, in order:

  remove_two_arm_node   on each tag in REMOVE_TAGS, which deletes the
                        node and joins its two arms into one segment
  insert_node_between   on the segment SPLIT_SEG, at its mid-point

The network is checked for sanity before and after (is_sane also
verifies Burgers vector conservation at every node), and the result is
written to OUT_DIR as JSON and compared against the stored reference,
which is what pins the outcome of the operations node by node.

The JSON is the interchange format both codes implement through
export_data / import_data, so the same reference can serve an exadis
counterpart; this file exercises the pydis side only.

Regenerate the reference with

    make loop_topol_op_ref
    cp output/loop_topol_op.json ref_data/

or by calling the script with --write-ref. A deliberate change to the
operations is installed that way, after checking that the differences
the test reports are the intended ones.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in
               ['python', 'lib', 'core/pydis/python']]
[sys.path.append(p) for p in pydis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np

from pydis.disnet import DisNet
from framework.disnet_manager import DisNetManager
from framework.testing import report

# the configuration, as node positions and connectivity
RN_FILE    = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the operations under test
REMOVE_TAGS = [(0, 5), (0, 6)]
SPLIT_SEG   = ((0, 8), (0, 9))

# where the result is written, and the reference it is compared against
OUT_DIR     = 'output'
RESULT_NAME = 'loop_topol_op.json'
REF_NAME    = 'loop_topol_op.json'


def init_loop_from_file(rn_file=RN_FILE, links_file=LINKS_FILE):
    """init_loop_from_file: build the loop configuration

    The link file carries no glide planes, but the topological
    operations propagate them and export_data writes them, so they have
    to be something meaningful rather than absent: n = b x t, the plane
    containing the Burgers vector and the line. That is well defined
    here (the smallest |b x t| over the 450 segments is 0.12), which it
    would not be for a configuration running exactly screw anywhere.
    """
    print("init_loop_from_file: rn_file = '%s', links_file = '%s'"
          % (rn_file, links_file))
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)

    n1 = links[:, 0].astype(int)
    n2 = links[:, 1].astype(int)
    planes = np.cross(links[:, 2:5], rn[n2] - rn[n1])
    norms = np.linalg.norm(planes, axis=1)
    if norms.min() < 1.0e-6:
        raise ValueError("init_loop_from_file: b x t vanishes on %d "
                         "segment(s); the glide plane is undefined "
                         "there" % int(np.sum(norms < 1.0e-6)))
    planes /= norms[:, None]

    G = DisNet()
    G.add_nodes_segments_from_list(rn, np.hstack((links[:, :5], planes)))
    return G


def set_glide_plane(G, tag1, tag2):
    """set_glide_plane: recompute a segment's glide plane as b x t

    Independent of which way round the segment is stored: reversing it
    negates both b and t, and (-b) x (-t) = b x t.
    """
    edge = G.segments((tag1, tag2))
    b = edge.burg_vec_from(tag1)
    t = G.nodes(tag2).R - G.nodes(tag1).R
    n = np.cross(b, t)
    edge.plane_normal = n / np.linalg.norm(n)


def remove_nodes(G, tags=REMOVE_TAGS):
    """remove_nodes: delete two-arm nodes, joining the arms they held

    The glide plane of the merged segment is recomputed rather than
    left as remove_two_arm_node leaves it. Remesh would normally not
    delete a two-arm node whose arms lie on different planes, which is
    the case all through this loop; deleting it anyway is the point of
    the test, and the segment it produces spans two planes, so the one
    inherited from whichever arm happened to be visited is not the
    plane of the new segment. It also would not be reproducible:
    neighbors_tags iterates a set of edge objects, so which arm is
    inherited from follows their addresses and varies between runs.

    What each removal should produce, the node gone and the two nodes
    it stood between connected to each other, is not asserted here: it
    is in the reference the result is compared against, which pins the
    whole network rather than one operation at a time.
    """
    for tag in tags:
        neighbors = tuple(G.neighbors_tags(tag))
        n_nodes, n_segs = G.num_nodes(), G.num_segments()

        G.remove_two_arm_node(tag)
        set_glide_plane(G, neighbors[0], neighbors[1])
        print("remove_two_arm_node(%s): neighbors %s now joined; "
              "nodes %d -> %d, segments %d -> %d"
              % (str(tag), str(neighbors), n_nodes, G.num_nodes(),
                 n_segs, G.num_segments()))


def insert_node(G, seg=SPLIT_SEG):
    """insert_node: split a segment by adding a node at its mid-point"""
    tag1, tag2 = seg
    n_nodes, n_segs = G.num_nodes(), G.num_segments()
    R = 0.5*(G.nodes(tag1).R + G.nodes(tag2).R)

    # get_new_tag recycles the tags freed by remove_nodes, so the new
    # node reuses one of them rather than extending the numbering
    new_tag = G.get_new_tag()
    G.insert_node_between(tag1, tag2, new_tag, R)
    print("insert_node_between(%s, %s): new node %s at %s; "
          "nodes %d -> %d, segments %d -> %d"
          % (str(tag1), str(tag2), str(new_tag),
             np.array2string(R, precision=4), n_nodes, G.num_nodes(),
             n_segs, G.num_segments()))


def write_result(N, out_dir=OUT_DIR):
    """write_result: save the network for comparison and for promotion"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / RESULT_NAME
    N.write_json(str(out_file))
    print("write_result: wrote %s (%d nodes, %d segments)"
          % (out_file, N.num_nodes(), N.num_segments()))
    print("")
    print("to install it, copy it into the ref_data/ next to this "
          "script:")
    print("    mkdir -p ref_data")
    print("    cp %s ref_data/%s" % (out_file, REF_NAME))
    print("")
    print("this reference encodes the operations this test applies")
    print("(remove %s, then split the segment %s);"
          % (', '.join(str(t) for t in REMOVE_TAGS),
             '-'.join(str(t) for t in SPLIT_SEG)))
    print("regenerate it whenever those change")
    return out_file


def load_ref():
    """load_ref: read the stored reference back into a pydis network

    Through DisNetManager.read_json, so the reference goes through the
    same import_data path any other consumer of the file would use.
    """
    N = DisNetManager(DisNet())
    N.read_json(str(ref_dir / REF_NAME))
    return N.get_disnet(DisNet)


def same_network(G, G_ref):
    """same_network: whether two networks describe the same disnet

    DisNet.is_equivalent does the comparison, and does it the way this
    test needs. It matches segments by their end tags rather than by
    position in an array, it reads the other network's Burgers vector
    from its own source tag, so a segment stored in either direction
    compares equal, and it already treats n and -n as one glide plane.

    Two things it does not do. It iterates only the nodes and segments
    of the network it is called on, so a reference holding more than
    the result would pass unnoticed, hence the counts below. And a
    segment the other network does not have makes it raise rather than
    return False, hence the guard.
    """
    if (G.num_nodes() != G_ref.num_nodes()
            or G.num_segments() != G_ref.num_segments()):
        print("counts differ: %d nodes, %d segments; the reference has "
              "%d and %d" % (G.num_nodes(), G.num_segments(),
                             G_ref.num_nodes(), G_ref.num_segments()))
        return False
    try:
        return bool(G.is_equivalent(G_ref))
    except (KeyError, AttributeError, TypeError) as err:
        print("the reference has no counterpart for %r" % (err,))
        return False


def main(write_ref_file=False):
    G = init_loop_from_file()
    print("nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    ok = report("network is sane as loaded", G.is_sane())
    remove_nodes(G)
    insert_node(G)
    ok &= report("network is sane after the operations", G.is_sane())

    N = DisNetManager(G)

    if write_ref_file:
        # a reference taken from a network that failed its own sanity
        # check would be worse than none
        if ok:
            write_result(N)
        else:
            print("the network is not sane, no reference written")
        return bool(ok)

    ref_file = ref_dir / REF_NAME
    if not ref_file.is_file():
        print("reference not found: ref_data/%s" % REF_NAME)
        print("generate it with")
        print("    make loop_topol_op_ref")
        print("then copy %s/%s into ref_data/"
              % (OUT_DIR, RESULT_NAME))
        return False

    ok &= report("resulting network matches the reference",
                 same_network(G, load_ref()))
    return bool(ok)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='topological operations on the 10-loop '
                    'configuration')
    parser.add_argument('--write-ref', dest='write_ref',
                        action='store_true', default=False,
                        help='save the resulting network as a reference '
                             'under output/ instead of comparing '
                             'against the stored one')
    args = parser.parse_args()

    passed = main(write_ref_file=args.write_ref)

    if not args.write_ref:
        tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
        print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
