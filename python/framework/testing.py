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

import json
import os
import subprocess
import sys

import numpy as np

GREEN = '\033[32m'
RED = '\033[31m'
RESET = '\033[0m'


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

    Returns np.inf if either configuration is empty.
    """
    ra = np.asarray(ra, dtype=float)[:, :3]
    rb = np.asarray(rb, dtype=float)[:, :3]
    if ra.size == 0 or rb.size == 0:
        return np.inf
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


def compare_to_ref(json_file, npz_file):
    """compare_to_ref: compare a run's output against a stored reference

    Returns (n_run, n_ref, max_distance); see max_nearest_distance.
    """
    r_run, _ = load_config(json_file)
    r_ref, h, _ = load_ref_npz(npz_file)
    return r_run.shape[0], r_ref.shape[0], max_nearest_distance(r_run, r_ref, h)
