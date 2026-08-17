"""Topological operations on the 10-loop configuration, PyDiS and ExaDiS.

Both codes apply the same two operations to the same configuration and
both results are compared against one stored reference, so they are
pinned to one network rather than each agreeing with itself:

  remove_two_arm_node   deletes each tag in REMOVE_TAGS, joining the
                        two arms it held into one segment
  insert_node_between   splits SPLIT_SEG at its mid-point

The exadis counterparts are reached through
ExaDisNet._get_serial_network(): merge_nodes_position, merging the node
into one of its neighbours at that neighbour's position so the others
stay put, followed by purge_network; and split_seg. Three things about
that interface are worth knowing before reading the code:
merge_nodes_position returns an error flag, so False means it worked;
it takes a Mat33 for the plastic strain increment, which a 3x3 numpy
array satisfies; and everything is addressed by index, so indices have
to be looked up from tags again after each purge_network, which
renumbers.

The glide plane of the merged segment is recomputed as b x t in both
codes rather than inherited from whichever arm survived. Remesh would
normally not delete a two-arm node whose arms lie on different planes,
which is the case all through this loop; deleting it anyway is the
point of the test, and the segment it produces spans two planes, so
neither arm's plane is the plane of the new segment.

Comparison runs through DisNet.is_equivalent, which matches nodes and
segments by tag rather than by array order, so it does not care that a
segment may be stored in either direction. The exadis result is
exported and imported into a DisNet for that; the network compared is
still the one exadis produced. The two codes agree on the tag given to
the new node only because both recycle freed tags LIFO.

Regenerate the reference with

    make loop_topol_op_ref
    cp output/loop_topol_op.json ref_data/

or by calling the script with --write-ref, which writes the pydis
result. A deliberate change to the operations is installed that way,
after checking that the differences the test reports are intended.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
opendis_paths = [str(opendis_root / p) for p in
                 ['python', 'lib', 'core/pydis/python',
                  'core/exadis/python']]
[sys.path.append(p) for p in opendis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np

from pydis.disnet import DisNet
from framework.disnet_manager import DisNetManager
from framework.testing import report, same_network

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

# cell edge for the exadis network. It holds the configuration and is
# not compared: it is non-periodic, so nothing here depends on it
BOX = 1000.0


def load_config(rn_file=RN_FILE, links_file=LINKS_FILE):
    """load_config: node positions and segments, shared by both codes

    The link file carries no glide planes, but the operations propagate
    them and export_data writes them, so they have to be something
    meaningful rather than absent: n = b x t, the plane containing the
    Burgers vector and the line. That is well defined here (the
    smallest |b x t| over the 450 segments is 0.12), which it would not
    be for a configuration running exactly screw anywhere.
    """
    print("load_config: rn_file = '%s', links_file = '%s'"
          % (rn_file, links_file))
    rn = np.loadtxt(input_dir / rn_file)[:, 1:]
    links = np.loadtxt(input_dir / links_file)

    n1 = links[:, 0].astype(int)
    n2 = links[:, 1].astype(int)
    planes = np.cross(links[:, 2:5], rn[n2] - rn[n1])
    norms = np.linalg.norm(planes, axis=1)
    if norms.min() < 1.0e-6:
        raise ValueError("load_config: b x t vanishes on %d segment(s); "
                         "the glide plane is undefined there"
                         % int(np.sum(norms < 1.0e-6)))
    planes /= norms[:, None]

    return rn, np.hstack((links[:, :5], planes))


def glide_plane(b, t):
    """glide_plane: b x t, normalized

    Independent of which way round the segment is stored: reversing it
    negates both b and t, and (-b) x (-t) = b x t.
    """
    n = np.cross(b, t)
    return n / np.linalg.norm(n)


# ---------------------------------------------------------------- pydis

def result_pydis():
    """result_pydis: the network after both operations, in pydis"""
    rn, segs = load_config()
    G = DisNet()
    G.add_nodes_segments_from_list(rn, segs)
    print("pydis : nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))
    if not G.is_sane():
        raise RuntimeError("result_pydis: the loaded network is not sane")

    for tag in REMOVE_TAGS:
        neighbors = tuple(G.neighbors_tags(tag))
        n_nodes, n_segs = G.num_nodes(), G.num_segments()
        G.remove_two_arm_node(tag)

        edge = G.segments((neighbors[0], neighbors[1]))
        edge.plane_normal = glide_plane(
            edge.burg_vec_from(neighbors[0]),
            G.nodes(neighbors[1]).R - G.nodes(neighbors[0]).R)
        print("pydis : remove_two_arm_node(%s): neighbors %s now "
              "joined; nodes %d -> %d, segments %d -> %d"
              % (str(tag), str(neighbors), n_nodes, G.num_nodes(),
                 n_segs, G.num_segments()))

    tag1, tag2 = SPLIT_SEG
    n_nodes, n_segs = G.num_nodes(), G.num_segments()
    R = 0.5*(G.nodes(tag1).R + G.nodes(tag2).R)
    new_tag = G.get_new_tag()
    G.insert_node_between(tag1, tag2, new_tag, R)
    print("pydis : insert_node_between(%s, %s): new node %s at %s; "
          "nodes %d -> %d, segments %d -> %d"
          % (str(tag1), str(tag2), str(new_tag),
             np.array2string(R, precision=4), n_nodes, G.num_nodes(),
             n_segs, G.num_segments()))
    return G


# --------------------------------------------------------------- exadis

def index_of(G, tag):
    """index_of: node index carrying a tag, which purging invalidates"""
    tags = np.array(G.get_tags())
    hit = np.where((tags[:, 0] == tag[0]) & (tags[:, 1] == tag[1]))[0]
    if hit.size != 1:
        raise ValueError("index_of: node %s not found" % str(tag))
    return int(hit[0])


def segment_index(serial, i, j):
    """segment_index: index of the segment joining two node indices"""
    for s in range(serial.number_of_segs()):
        seg = serial.segs(s)
        if {seg.n1, seg.n2} == {i, j}:
            return s
    raise ValueError("segment_index: no segment between %d and %d"
                     % (i, j))


def result_exadis():
    """result_exadis: the network after both operations, in exadis

    Returned as a DisNet, so that the same comparison serves both
    codes. Only the comparison goes through pydis; the network is the
    one exadis produced.
    """
    import pyexadis
    from pyexadis_base import ExaDisNet

    rn, segs = load_config()
    cell = pyexadis.Cell(h=BOX*np.eye(3), origin=-0.5*BOX*np.ones(3),
                         is_periodic=[False, False, False])
    G = ExaDisNet(cell, rn, segs)
    print("exadis: nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    dEp = np.zeros((3, 3))
    for tag in REMOVE_TAGS:
        serial = G.net._get_serial_network()
        i_drop = index_of(G, tag)
        conn = serial.conn(i_drop)
        if conn.num != 2:
            raise ValueError("result_exadis: node %s has %d arms"
                             % (str(tag), conn.num))
        tags = np.array(G.get_tags())
        keep_tag = tuple(tags[conn.node(0)])
        other_tag = tuple(tags[conn.node(1)])

        i_keep = conn.node(0)
        pos = np.array(serial.nodes(i_keep).pos)
        n_nodes, n_segs = G.num_nodes(), G.num_segments()
        # the return value is an error flag, not success
        if serial.merge_nodes_position(i_keep, i_drop, pos, dEp):
            raise RuntimeError("merge_nodes_position failed on %s"
                               % str(tag))
        serial.purge_network()

        # indices have been renumbered by the purge
        serial = G.net._get_serial_network()
        iseg = segment_index(serial, index_of(G, keep_tag),
                             index_of(G, other_tag))
        seg = serial.segs(iseg)
        seg.plane = glide_plane(np.array(seg.burg),
                                np.array(serial.nodes(seg.n2).pos)
                                - np.array(serial.nodes(seg.n1).pos))
        print("exadis: merge_nodes_position(%s into %s): nodes %d -> "
              "%d, segments %d -> %d"
              % (str(tag), str(keep_tag), n_nodes, G.num_nodes(),
                 n_segs, G.num_segments()))

    tag1, tag2 = SPLIT_SEG
    serial = G.net._get_serial_network()
    i1, i2 = index_of(G, tag1), index_of(G, tag2)
    iseg = segment_index(serial, i1, i2)
    R = 0.5*(np.array(serial.nodes(i1).pos)
             + np.array(serial.nodes(i2).pos))
    n_nodes, n_segs = G.num_nodes(), G.num_segments()
    inew = serial.split_seg(iseg, R)
    new_tag = tuple(np.array(G.get_tags())[inew])
    print("exadis: split_seg(%s, %s): new node %s at %s; nodes %d -> "
          "%d, segments %d -> %d"
          % (str(tag1), str(tag2), str(new_tag),
             np.array2string(R, precision=4), n_nodes, G.num_nodes(),
             n_segs, G.num_segments()))

    G_pydis = DisNet()
    G_pydis.import_data(G.export_data())
    return G_pydis


# ----------------------------------------------------- reference and run

def write_result(G, out_dir=OUT_DIR):
    """write_result: save the network for comparison and for promotion"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / RESULT_NAME
    N = DisNetManager(G)
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


def test_code(label, G, ref):
    """test_code: one code's result against the stored reference"""
    ok = report("%s: network is sane after the operations" % label,
                G.is_sane())
    ok &= report("%s: resulting network matches the reference" % label,
                 same_network(G, ref))
    return bool(ok)


def main(write_ref_file=False):
    G_pydis = result_pydis()

    if write_ref_file:
        # a reference taken from a network that failed its own sanity
        # check would be worse than none
        if not G_pydis.is_sane():
            print("the network is not sane, no reference written")
            return False
        write_result(G_pydis)
        return True

    ref_file = ref_dir / REF_NAME
    if not ref_file.is_file():
        print("reference not found: ref_data/%s" % REF_NAME)
        print("generate it with")
        print("    make loop_topol_op_ref")
        print("then copy %s/%s into ref_data/"
              % (OUT_DIR, RESULT_NAME))
        return False
    ref = load_ref()

    ok = test_code("pydis ", G_pydis, ref)

    try:
        import pyexadis
    except ImportError:
        print("pyexadis not available; the exadis half of this test "
              "did not run")
        return False

    pyexadis.initialize()
    ok &= test_code("exadis", result_exadis(), ref)
    pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='topological operations on the 10-loop '
                    'configuration, pydis and exadis')
    parser.add_argument('--write-ref', dest='write_ref',
                        action='store_true', default=False,
                        help='save the pydis result as a reference '
                             'under output/ instead of comparing '
                             'against the stored one')
    args = parser.parse_args()

    passed = main(write_ref_file=args.write_ref)

    if not args.write_ref:
        tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
        print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
