"""Topological operations on the 10-loop configuration, in ExaDiS.

The counterpart of test_loop_topol_op_pydis.py: the same two
operations, applied to the same configuration, compared against the
reference that test writes. So the reference pins the outcome of the
operations in both codes, which is what makes this a cross-code check
rather than each code agreeing with itself.

The operations map across as:

  remove_two_arm_node   merge_nodes_position, merging the node into one
                        of its neighbours at that neighbour's position,
                        then purge_network to drop what it emptied
  insert_node_between   split_seg at the mid-point

They are reached through ExaDisNet._get_serial_network(). Three things
about that interface are worth knowing before reading the code:
merge_nodes_position returns an error flag, so False means it worked;
it takes a Mat33 for the plastic strain increment, which a 3x3 numpy
array satisfies; and everything is addressed by index, so the indices
have to be looked up from tags again after each purge_network, which
renumbers.

As in the pydis test, the glide plane of the merged segment is
recomputed as b x t rather than inherited from whichever arm survived,
since the two arms of these nodes lie on different planes.

The comparison itself runs through pydis: the exadis result is exported
and imported into a DisNet so that DisNet.is_equivalent can match the
two networks by tag rather than by array order. Nothing of the exadis
side is bypassed by that, the network compared is the one exadis
produced.
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

import pyexadis
from pyexadis_base import ExaDisNet
from pydis.disnet import DisNet
from framework.disnet_manager import DisNetManager
from framework.testing import report, same_network

# the configuration, as node positions and connectivity
RN_FILE    = 'loop_rn.dat'
LINKS_FILE = 'loop_links.dat'

# the operations under test, the same ones the pydis test applies
REMOVE_TAGS = [(0, 5), (0, 6)]
SPLIT_SEG   = ((0, 8), (0, 9))

# the reference, written by test_loop_topol_op_pydis.py --write-ref
REF_NAME = 'loop_topol_op.json'

# cell edge. The cell holds the configuration and is not compared: it
# is non-periodic, so nothing here depends on it
BOX = 1000.0


def init_loop_from_file(rn_file=RN_FILE, links_file=LINKS_FILE):
    """init_loop_from_file: build the loop configuration in exadis

    The glide planes are b x t, as in the pydis test, since the link
    file carries none and the operations propagate them.
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

    cell = pyexadis.Cell(h=BOX*np.eye(3), origin=-0.5*BOX*np.ones(3),
                         is_periodic=[False, False, False])
    return ExaDisNet(cell, rn, np.hstack((links[:, :5], planes)))


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


def set_glide_plane(serial, iseg):
    """set_glide_plane: recompute a segment's glide plane as b x t

    Independent of which way round the segment is stored: reversing it
    negates both b and t, and (-b) x (-t) = b x t.
    """
    seg = serial.segs(iseg)
    b = np.array(seg.burg)
    t = (np.array(serial.nodes(seg.n2).pos)
         - np.array(serial.nodes(seg.n1).pos))
    n = np.cross(b, t)
    seg.plane = n / np.linalg.norm(n)


def remove_nodes(G, tags=REMOVE_TAGS):
    """remove_nodes: delete two-arm nodes, joining the arms they held

    Each node is merged into its first neighbour, at that neighbour's
    position, so the surviving nodes stay where they are and only the
    node between them goes, which is what remove_two_arm_node does on
    the pydis side.
    """
    dEp = np.zeros((3, 3))
    for tag in tags:
        serial = G.net._get_serial_network()
        i_drop = index_of(G, tag)
        conn = serial.conn(i_drop)
        if conn.num != 2:
            raise ValueError("remove_nodes: node %s has %d arms"
                             % (str(tag), conn.num))
        neighbors = [conn.node(j) for j in range(conn.num)]
        keep_tag = tuple(np.array(G.get_tags())[neighbors[0]])
        other_tag = tuple(np.array(G.get_tags())[neighbors[1]])

        i_keep = neighbors[0]
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
        set_glide_plane(serial, iseg)

        print("merge_nodes_position(%s into %s): nodes %d -> %d, "
              "segments %d -> %d"
              % (str(tag), str(keep_tag), n_nodes, G.num_nodes(),
                 n_segs, G.num_segments()))


def insert_node(G, seg=SPLIT_SEG):
    """insert_node: split a segment by adding a node at its mid-point"""
    tag1, tag2 = seg
    serial = G.net._get_serial_network()
    i1, i2 = index_of(G, tag1), index_of(G, tag2)
    iseg = segment_index(serial, i1, i2)
    R = 0.5*(np.array(serial.nodes(i1).pos)
             + np.array(serial.nodes(i2).pos))

    n_nodes, n_segs = G.num_nodes(), G.num_segments()
    inew = serial.split_seg(iseg, R)
    new_tag = tuple(np.array(G.get_tags())[inew])
    print("split_seg(%s, %s): new node %s at %s; nodes %d -> %d, "
          "segments %d -> %d"
          % (str(tag1), str(tag2), str(new_tag),
             np.array2string(R, precision=4), n_nodes, G.num_nodes(),
             n_segs, G.num_segments()))


def as_disnet(G):
    """as_disnet: the exadis network as a pydis DisNet, for comparison"""
    G_pydis = DisNet()
    G_pydis.import_data(G.export_data())
    return G_pydis


def load_ref():
    """load_ref: read the stored reference into a pydis network"""
    N = DisNetManager(DisNet())
    N.read_json(str(ref_dir / REF_NAME))
    return N.get_disnet(DisNet)


def main():
    G = init_loop_from_file()
    print("nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    remove_nodes(G)
    insert_node(G)

    G_result = as_disnet(G)
    ok = report("network is sane after the operations",
                G_result.is_sane())

    ref_file = ref_dir / REF_NAME
    if not ref_file.is_file():
        print("reference not found: ref_data/%s" % REF_NAME)
        print("generate it with")
        print("    make loop_topol_op_ref")
        print("then copy output/%s into ref_data/" % REF_NAME)
        return False

    ok &= report("resulting network matches the reference",
                 same_network(G_result, load_ref()))
    return bool(ok)


if __name__ == "__main__":
    pyexadis.initialize()
    passed = main()
    pyexadis.finalize()

    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
