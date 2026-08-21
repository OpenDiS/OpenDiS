"""@package docstring
testing: shared helpers for the OpenDiS test suite

Reporting, running an example script as a subprocess, and comparing two dislocation
configurations. Lives in framework/ rather than under tests/ because it is used by tests in
more than one place and depends on neither PyDiS nor ExaDiS. Import with the repository
located from __file__, e.g.

    sys.path.append(str(Path(__file__).resolve().parents[3] / 'python'))
    from framework.testing import report, run_script, compare_configs

Deliberately not named test_*.py: test discovery in the makefiles and in CMake globs
test_*.py, and this module is support code, not a test.
"""

import contextlib
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

GREEN = '\033[32m'
RED = '\033[31m'
RESET = '\033[0m'


class _Captured:
    """_Captured: carries text out of quiet_native_output's closed file"""

    def __init__(self):
        self.text = ''


@contextlib.contextmanager
def quiet_native_output(enabled=True):
    """quiet_native_output: swallow what C++ libraries write to the terminal

    Kokkos prints a long banner from pyexadis.initialize(), ExaDiS warns
    about Burgers vector conservation for any network of open segment
    pairs, and Intel's OpenMP runtime prints "OMP: Info #277" from inside
    ExaDiS on every call. All three land in the middle of test output and
    none of them indicates a problem.

    They come from C++, so contextlib.redirect_stdout does not catch them:
    that rebinds sys.stdout while the C++ writes to file descriptor 1.
    Redirecting the descriptors is what works.

    Yields a small object whose `.text` holds whatever was captured, set
    once the block finishes. It cannot be a file: the temporary is closed
    on the way out, so a caller reading it afterwards would find it shut.
    """
    captured = _Captured()
    if not enabled:
        yield captured
        return
    with tempfile.TemporaryFile(mode='w+') as tmp:
        sys.stdout.flush()
        sys.stderr.flush()
        saved_out, saved_err = os.dup(1), os.dup(2)
        os.dup2(tmp.fileno(), 1)
        os.dup2(tmp.fileno(), 2)
        try:
            yield captured
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_out, 1)
            os.dup2(saved_err, 2)
            os.close(saved_out)
            os.close(saved_err)
            tmp.seek(0)
            captured.text = tmp.read()


def kokkos_summary(text):
    """kokkos_summary: the lines of the Kokkos banner worth keeping

    Which device it selected and how large a thread pool it took. Those
    two decide whether a comparison is being made on the hardware you
    think it is, so they survive; the rest of the banner does not.
    """
    if not text:
        return []
    keep = []
    for line in text.splitlines():
        if 'thread_pool_topology' in line or ': Selected' in line:
            keep.append(line.strip())
    return keep


def array_digest(arr):
    """array_digest: a short stable fingerprint of a float array

    For pairing a reference against the fixture it was computed from, when
    the two live in different files and are blessed at different times. A
    reference carrying the digest of its inputs can say "these describe
    different pairs" instead of silently comparing forces to the wrong
    geometry.

    Bytes of the float64 representation, so it is exact rather than
    tolerant: any change of geometry at all, however small, is a different
    fixture and a stale reference.
    """
    a = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
    return hashlib.sha256(a.tobytes()).hexdigest()[:16]


def report(name, passed):
    """report: print a coloured PASSED/FAILED line; returns `passed` unchanged"""
    tag = GREEN + 'PASSED' + RESET if passed else RED + 'FAILED' + RESET
    print("%s ... %s" % (name, tag))
    return passed


def run_script(script, args=(), cwd=None, env_extra=None):
    """run_script: run a python script as a subprocess, returning True on exit status 0

    Uses sys.executable so the child runs under the same interpreter as the test, which
    matters because the interpreter that can import pyexadis is not always `python3`.
    cwd defaults to the caller's working directory; the example scripts write their output
    relative to it.
    """
    env = dict(os.environ)
    env.setdefault('MPLBACKEND', 'Agg')          # keep the examples headless
    if env_extra:
        env.update(env_extra)
    cmd = [sys.executable, str(script)] + [str(a) for a in args]
    print("run: %s" % ' '.join(cmd))
    result = subprocess.run(cmd, cwd=cwd, env=env)
    return result.returncode == 0


def load_config(json_file):
    """load_config: read a configuration written by DisNetManager.write_json

    Returns (positions, cell_h). Positions are returned as stored, not folded.
    """
    with open(json_file) as f:
        data = json.load(f)
    positions = np.array(data['nodes']['positions'], dtype=float)
    cell_h = np.array(data['cell']['h'], dtype=float)
    return positions, cell_h


def max_nearest_distance(ra, rb, h):
    """max_nearest_distance: greatest distance from any node of one configuration to the
    nearest node of the other, under the minimum image convention

    Deliberately order-independent. PyDiS and ExaDiS emit nodes in different orders, so an
    element-wise comparison of two identical configurations can report a difference of order
    the box size. Matching is done both ways, so a node present in one configuration but not
    the other cannot hide.

    Returns np.inf if either configuration is empty, including when both
    are. Two empty configurations do describe the same thing, but a run that
    lost every node is not a result worth passing silently, and the node
    counts are compared separately and would agree at zero.

    The emptiness test comes before the column slice: a configuration with no
    nodes reads back from write_json one-dimensional, so slicing it first
    raises IndexError rather than reporting anything.
    """
    ra = np.asarray(ra, dtype=float)
    rb = np.asarray(rb, dtype=float)
    if ra.size == 0 or rb.size == 0:
        return np.inf
    ra, rb = ra[:, :3], rb[:, :3]
    hinv = np.linalg.inv(np.asarray(h, dtype=float))
    d = ra[:, None, :] - rb[None, :, :]
    s = d @ hinv.T
    s -= np.rint(s)
    d = s @ np.asarray(h, dtype=float).T
    dist = np.linalg.norm(d, axis=-1)
    return max(dist.min(axis=1).max(), dist.min(axis=0).max())


def compare_configs(file_a, file_b):
    """compare_configs: compare two configurations written by write_json

    Returns (n_a, n_b, max_distance); see max_nearest_distance.
    """
    ra, h = load_config(file_a)
    rb, _ = load_config(file_b)
    return ra.shape[0], rb.shape[0], max_nearest_distance(ra, rb, h)


def write_ref_npz(net, filename, source=''):
    """write_ref_npz: save a configuration as a reference .npz

    Stores the node positions and the cell, which is all a comparison needs, plus the
    connectivity and the name of the code that produced it so the file is self-describing.
    """
    data = net.export_data()
    np.savez(filename,
             positions=np.array(data['nodes']['positions'], dtype=float),
             constraints=np.array(data['nodes']['constraints'], dtype=float),
             cell_h=np.array(data['cell']['h'], dtype=float),
             nodeids=np.array(data['segs']['nodeids'], dtype=int),
             burgers=np.array(data['segs']['burgers'], dtype=float),
             planes=np.array(data['segs']['planes'], dtype=float),
             source=np.array(source))
    print("write_ref_npz: wrote %s (%d nodes, source '%s')"
          % (filename, data['nodes']['positions'].shape[0], source))


def load_ref_npz(filename):
    """load_ref_npz: read a reference written by write_ref_npz

    Returns (positions, cell_h, source).
    """
    with np.load(filename, allow_pickle=False) as f:
        return f['positions'], f['cell_h'], str(f['source'])


def report_close(name, values, ref_values, atol, rtol=0.0):
    """report_close: report whether two arrays agree, showing the error

    The reporting counterpart of np.allclose: prints one PASSED/FAILED
    line carrying the largest deviation and the tolerance it was judged
    against, so a passing run still says how much margin it had.
    """
    values, ref_values = np.asarray(values), np.asarray(ref_values)
    if values.shape != ref_values.shape:
        return report("%s (shape %s vs reference %s)"
                      % (name, values.shape, ref_values.shape), False)
    max_err = float(np.max(np.abs(values - ref_values)))
    passed = bool(np.allclose(values, ref_values, rtol=rtol, atol=atol))
    return report("%s: max error %.4e, atol %.1e" % (name, max_err, atol),
                  passed)


def same_network(G, G_ref, verbose=True):
    """same_network: whether two networks describe the same disnet

    DisNet.is_equivalent does the comparison, and does it the way a
    test needs. It matches segments by their end tags rather than by
    position in an array, it reads the other network's Burgers vector
    from its own source tag, so a segment stored in either direction
    compares equal, and it treats n and -n as one glide plane.

    Two things it does not do. It iterates only the nodes and segments
    of the network it is called on, so a reference holding more than
    the result would pass unnoticed, hence the counts below. And a
    segment the other network does not have makes it raise rather than
    return False, hence the guard.
    """
    if (G.num_nodes() != G_ref.num_nodes()
            or G.num_segments() != G_ref.num_segments()):
        if verbose:
            print("counts differ: %d nodes, %d segments; the reference "
                  "has %d and %d" % (G.num_nodes(), G.num_segments(),
                                     G_ref.num_nodes(),
                                     G_ref.num_segments()))
        return False
    try:
        return bool(G.is_equivalent(G_ref))
    except (KeyError, AttributeError, TypeError) as err:
        if verbose:
            print("the reference has no counterpart for %r" % (err,))
        return False


def load_force_ref(npz_file, expected, regen_hint=''):
    """load_force_ref: read a nodal-force reference, checking its settings

    A force reference is only meaningful together with the constants it
    was computed from, so those are stored in the .npz and checked here
    against what the caller is about to compare. `expected` maps field
    name to value; a missing or differing field is reported and None is
    returned, so a stale reference is named as such rather than showing
    up as a physics failure.

    Returns the opened NpzFile, or None.
    """
    npz_file = Path(npz_file)
    # shown relative to the working directory where that is possible, so
    # the message reads the way the reader would type it
    try:
        shown = npz_file.relative_to(Path.cwd())
    except ValueError:
        shown = Path(npz_file.parent.name) / npz_file.name

    if not npz_file.is_file():
        print("reference not found: %s" % shown)
        if regen_hint:
            print(regen_hint)
        return None

    ref = np.load(npz_file, allow_pickle=False)
    for key, val in expected.items():
        if key not in ref.files:
            print("reference %s has no '%s' field; regenerate it"
                  % (shown, key))
            if regen_hint:
                print(regen_hint)
            return None
        if float(ref[key]) != float(val):
            print("reference %s was generated with %s = %g, this test "
                  "uses %g; regenerate it"
                  % (shown, key, float(ref[key]), float(val)))
            if regen_hint:
                print(regen_hint)
            return None

    source = str(ref['source']) if 'source' in ref.files else 'unknown'
    print("load_force_ref: '%s', source '%s'" % (npz_file.name, source))
    return ref


def compare_to_ref(json_file, npz_file):
    """compare_to_ref: compare a run's output against a stored reference

    Returns (n_run, n_ref, max_distance); see max_nearest_distance.
    """
    r_run, _ = load_config(json_file)
    r_ref, h, _ = load_ref_npz(npz_file)
    return r_run.shape[0], r_ref.shape[0], max_nearest_distance(r_run, r_ref, h)
