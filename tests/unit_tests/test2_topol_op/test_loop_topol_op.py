"""Topological operations on the 10-loop configuration.

Two operations are applied to the loaded network, in order:

  remove_two_arm_node   on each tag in REMOVE_TAGS, which deletes the
                        node and joins its two arms into one segment
  insert_node_between   on the segment SPLIT_SEG, at its mid-point

Each is checked for the counts it should produce, for the connectivity
it should leave behind, and for network sanity (is_sane also verifies
Burgers vector conservation at every node). The result is written to
OUT_DIR as JSON, which is the interchange format both codes implement
through export_data / import_data, and is then compared against the
stored reference by loading that reference with each code in turn. So
the reference pins the outcome of the operations, and loading it twice
pins the claim that the file is readable by both.

Regenerate the reference with

    python3 test_loop_topol_op.py
    cp output/loop_topol_op.json ref_data/

The test always writes its result, so a deliberate change to the
operations is installed by copying the file, after checking that the
reported differences are the intended ones.
"""

import sys
from pathlib import Path
opendis_root = Path(__file__).resolve().parents[3]
pydis_paths = [str(opendis_root / p) for p in
               ['python', 'lib', 'core/pydis/python',
                'core/exadis/python']]
[sys.path.append(p) for p in pydis_paths if not p in sys.path]
input_dir = Path(__file__).resolve().parent / 'input_data'
ref_dir   = Path(__file__).resolve().parent / 'ref_data'

import numpy as np

from pydis.disnet import DisNet
from framework.disnet_manager import DisNetManager
from framework.testing import report, report_close

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

# positions survive the JSON round trip exactly, so this only absorbs
# a reference regenerated on a different platform
atol = 1.0e-12


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


def remove_nodes(G, tags=REMOVE_TAGS):
    """remove_nodes: delete two-arm nodes, joining the arms they held

    Returns whether every removal did what it should: the node is gone,
    the two nodes it stood between are now connected to each other, and
    the network is still sane.
    """
    ok = True
    for tag in tags:
        neighbors = tuple(G.neighbors_tags(tag))
        n_nodes, n_segs = G.num_nodes(), G.num_segments()

        G.remove_two_arm_node(tag)
        print("remove_two_arm_node(%s): neighbors %s now joined; "
              "nodes %d -> %d, segments %d -> %d"
              % (str(tag), str(neighbors), n_nodes, G.num_nodes(),
                 n_segs, G.num_segments()))

        ok &= report("  node %s removed" % str(tag),
                     not G.has_node(tag))
        ok &= report("  its neighbors %s are now connected"
                     % str(neighbors),
                     G.has_segment(neighbors[0], neighbors[1]))
        ok &= report("  counts dropped by one node and one segment",
                     G.num_nodes() == n_nodes - 1
                     and G.num_segments() == n_segs - 1)
        ok &= report("  network is sane", G.is_sane())
    return bool(ok)


def insert_node(G, seg=SPLIT_SEG):
    """insert_node: split a segment by adding a node at its mid-point

    Returns whether the split did what it should: the original segment
    is gone, the two halves exist, the new node sits at the mid-point,
    and the network is still sane.
    """
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

    ok = report("  segment %s-%s is gone" % (str(tag1), str(tag2)),
                not G.has_segment(tag1, tag2))
    ok &= report("  both halves exist",
                 G.has_segment(tag1, new_tag)
                 and G.has_segment(new_tag, tag2))
    ok &= report("  counts grew by one node and one segment",
                 G.num_nodes() == n_nodes + 1
                 and G.num_segments() == n_segs + 1)
    ok &= report_close("  new node at the mid-point",
                       G.nodes(new_tag).R, R, atol=0.0)
    ok &= report("  network is sane", G.is_sane())
    return bool(ok)


def write_result(N, out_dir=OUT_DIR):
    """write_result: save the network for comparison and for promotion"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / RESULT_NAME
    N.write_json(str(out_file))
    print("write_result: wrote %s (%d nodes, %d segments)"
          % (out_file, N.num_nodes(), N.num_segments()))
    return out_file


def load_ref(disnet_type):
    """load_ref: read the reference into a network of the given type

    This is what makes the format claim testable: the same file is read
    back through DisNetManager into a pydis DisNet and into an exadis
    ExaDisNet, both of which go through import_data.
    """
    ref_file = ref_dir / REF_NAME
    N = DisNetManager(disnet_type())
    N.read_json(str(ref_file))
    return N


def compare_networks(data, ref, label):
    """compare_networks: check two exported networks describe the same

    Tags, node ids and counts have to match exactly; positions, Burgers
    vectors and glide planes are compared with atol.
    """
    ok = True
    for group, exact_fields, close_fields in (
            ('nodes', ['tags', 'constraints'], ['positions']),
            ('segs', ['nodeids'], ['burgers', 'planes'])):
        for field in exact_fields:
            a = np.asarray(data[group][field])
            b = np.asarray(ref[group][field])
            ok &= report("%s %s.%s match" % (label, group, field),
                         a.shape == b.shape and np.array_equal(a, b))
        for field in close_fields:
            ok &= report_close("%s %s.%s" % (label, group, field),
                               data[group][field], ref[group][field],
                               atol=atol)
    return bool(ok)


def main():
    G = init_loop_from_file()
    print("nodes = %d, segments = %d"
          % (G.num_nodes(), G.num_segments()))

    ok = report("network is sane as loaded", G.is_sane())
    ok &= remove_nodes(G)
    ok &= insert_node(G)

    N = DisNetManager(G)
    write_result(N)
    data = N.export_data()

    ref_file = ref_dir / REF_NAME
    if not ref_file.is_file():
        print("reference not found: %s" % ref_file)
        print("install the result just written with")
        print("    cp %s/%s %s" % (OUT_DIR, RESULT_NAME, ref_dir))
        return False

    ok &= compare_networks(data, load_ref(DisNet).export_data(),
                           "pydis reload:")

    try:
        import pyexadis
        from pyexadis_base import ExaDisNet
    except ImportError:
        print("pyexadis not available; skipping the exadis reload. "
              "The format claim is only half tested in this run.")
        return bool(ok)

    pyexadis.initialize()
    ok &= compare_networks(data, load_ref(ExaDisNet).export_data(),
                           "exadis reload:")
    pyexadis.finalize()
    return bool(ok)


if __name__ == "__main__":
    passed = main()
    tag = '\033[32m' + 'PASSED' if passed else '\033[31m' + 'FAILED'
    print("test " + tag + '\033[0m')
    sys.exit(0 if passed else 1)
