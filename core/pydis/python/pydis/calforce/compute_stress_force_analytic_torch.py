"""Segment-segment interaction forces on a GPU, through torch.

A port of compute_stress_force_analytic_python.py, the isotropic
non-singular (SBA) segment-segment kernel, onto torch tensors so the whole
pair set is evaluated on one device in a handful of fused kernels rather
than one batch at a time on the CPU.

    from pydis.calforce.compute_stress_force_analytic_torch import (
        torch_segseg_force_vec)
    f1, f2, f3, f4 = torch_segseg_force_vec(p1, p2, p3, p4, b12, b34,
                                            mu, nu, a)

The signature matches python_segseg_force_vec exactly, and numpy in gives
numpy out, so this is a drop-in replacement. Pass torch tensors instead
and you get tensors back on their own device, which is what avoids a
round trip when the caller already holds its data on the GPU.

WHY IT IS A REWRITE RATHER THAN A TRANSCRIPTION

The numpy original spells every broadcast out as np.multiply(A,
np.array([x, x, x]).T). That materialises a temporary of the full array
shape for each of roughly 500 such expressions, which costs little in
numpy and a great deal on a GPU, where each temporary is a separate
kernel launch and a separate allocation. Here the same arithmetic is
written as A * x[:, None], which broadcasts without materialising
anything. The algebra is unchanged and the variable names are kept
identical to the original so the two files can be read side by side.

BOTH BRANCHES ARE PORTED, AND THAT IS NOT OPTIONAL

The kernel switches to a lower-dimensional formula when the two segments
come within 1-c^2 < 1e-4 of parallel. It is tempting to port only the
general branch and leave the near-parallel one on the CPU, since it looks
like an edge case. It is not. Measured over 80 force calls of
examples/02_frank_read_src/test_frank_read_src_pydis_elast.py, 2778 of
17938 selected pairs take it, 15.5%. Neighbouring segments along a
dislocation line are nearly parallel with one another, and being
neighbours they are also the closest pairs, so they dominate the
short-range interactions rather than being rare. A CPU fallback would
have synchronised the device on every call and handed it a sixth of the
work.

ACCURACY

In float64 this agrees with the compiled ParaDiS kernel to the same
degree the numpy implementation does, which is the relevant comparison:
both are independent orderings of the same arithmetic. Near parallel that
is ~1e-8 absolute rather than ~1e-16, which is the kernel's own
conditioning floor rather than anything to do with torch; see the block
comment in tests/unit_tests/test1_node_force/test_segseg_force_pydis_exadis.py.

float32 DOES NOT WORK for this kernel and is refused unless the caller
passes allow_float32=True. This is not the usual "single precision is a
bit less accurate" caveat. Measured against the numpy kernel over 5000
random pairs, float32 gives a median relative error of 5.4e-03 where the
geometry is benign, 1-c^2 >= 1e-4, and 6.3e-03 near parallel, with worst
cases of 1.8e+04 and 6.6e+03 respectively. Not a single non-finite value
is produced, so nothing announces itself: the numbers come back looking
ordinary while carrying no correct digits at all in the worst rows and
two or three in the median.

The reason is that the kernel cancels quantities several orders larger
than its result, and float32 has only about seven digits to lose. There
is nowhere for the cancellation to happen. Halving the memory traffic is
not worth it at that price, so the knob exists for experiments only and
has to be asked for by name.
"""

import os
import sys

import numpy as np

try:
    import torch
except ImportError as exc:                              # pragma: no cover
    raise ImportError(
        "compute_stress_force_analytic_torch requires torch. Install it "
        "with 'pip install torch', choosing the build that matches the "
        "CUDA version on the machine you intend to run on."
    ) from exc


# Threshold on 1-c^2 below which the two segments are treated as parallel
# and handed to the special formula, which also receives it as ecrit. Must
# match the numpy implementation, the compiled ParaDiS kernel
# (core/pydis/c/calforce/SegSegForce.c: eps = 1e-4, ecrit = 1e-4) and
# exadis (core/exadis/src/force_types/force_common.h:1039, :611).
EPS_PARALLEL = 1e-04

# The "is this displacement nonzero at working precision" test that gates
# the two correction branches inside the near-parallel formula. Unlike
# EPS_PARALLEL this one is a property of the arithmetic rather than of the
# physics, so it scales with the dtype: 1e-16 is right for float64 and
# would never fire in float32, where the correction is needed just as much.
_CORRECTION_EPS = {torch.float64: 1e-16, torch.float32: 1e-7}

# The near-parallel formula calls the general one for its corrections, and
# those calls can in principle be near-parallel themselves. In practice
# cotanthetac places them just outside the threshold and the recursion
# terminates immediately, but a guard turns a pathological input into a
# clear error rather than a stack overflow.
_MAX_DEPTH = 4

# torch.compile support. Off unless asked for, by PYDIS_TORCH_COMPILE in the
# environment or enable_compile() in code.
#
#   'dynamic'  compile once, for symbolic shapes. The default.
#   'static'   compile separately for each input shape.
#
# Measured at 100k pairs, 4 threads, on a realistic geometry mix. Four
# independent runs, each with an empty Inductor cache and each compiled on the
# shape it then times:
#
#   round      eager   dynamic    static
#     1        239.8     101.6     116.0
#     2        271.0     106.0     124.5
#     3        251.6     109.6     115.1
#     4        270.9     118.4     123.2
#
# Compilation is worth having: roughly 258 ms to 109. And 'dynamic' came out
# ahead of 'static' in every round, by 4 to 17%, which is why it is the
# default. It also compiles once and serves any pair count afterwards, where
# 'static' pays around 13 s for each new one, and a simulation hands the
# kernel a different count almost every step. There is no case here for
# arranging fixed shapes: padding batches up to power-of-two buckets was
# measured at 795 ms against 611 for plain dynamic, since the padding wastes
# 1.36x the work, and splitting into fixed chunks came out level. Fusion
# shifts the answer by ~1.5e-12 relative, orders below the kernel's own noise
# floor.
#
# Two things that make this easy to measure wrongly, both of which produced
# badly misleading numbers before being controlled for:
#
#   - a dynamic-shape kernel is generated once and reused, and it inherits its
#     tuning from whichever shape compiled it. Warmed at 1000 pairs it then
#     runs 100k in 209 ms; warmed at 100k, 114 ms. Warm on a representative
#     pair count.
#   - Inductor caches generated kernels on disk, under TORCHINDUCTOR_CACHE_DIR,
#     and that cache outlives the process. A kernel built at some other shape
#     will be reused silently.
#
# The spread across rounds is 8 to 16%, so read differences of a few percent
# as noise.
_COMPILE_MODES = ('dynamic', 'static')
_compile_mode = os.environ.get('PYDIS_TORCH_COMPILE', '').strip().lower()
if _compile_mode in ('1', 'true', 'yes', 'on'):
    _compile_mode = 'dynamic'
elif _compile_mode not in _COMPILE_MODES:
    _compile_mode = None
_compiled_general = None
_compile_checked = False


def enable_compile(mode='dynamic'):
    """enable_compile: route the general branch through torch.compile

    mode is 'dynamic', 'static', or None/False to go back to eager. Takes
    effect on the next call; the compilation itself happens then, not here.

    Only the general branch is compiled. The classification above it and
    the near-parallel formula both branch on data, `int(spindex.sum())`
    and `nonzero()`, which Dynamo cannot trace without breaking the graph
    and which on a GPU would force a synchronisation. The general branch
    is pure elementwise arithmetic of fixed shape and is where nearly all
    the time goes on a realistic mix, so it is the part worth fusing.
    """
    global _compile_mode, _compiled_general, _compile_checked
    _compile_checked = False
    if mode in (None, False):
        _compile_mode = None
    elif mode in _COMPILE_MODES:
        _compile_mode = mode
    else:
        raise ValueError("enable_compile: mode must be one of %s, or None"
                         % (_COMPILE_MODES,))
    _compiled_general = None


def _general_dispatch(*args):
    """_general_dispatch: the general branch, compiled if that was asked for

    Compilation can fail in two places and neither should take a
    simulation down over what is only an optimisation.

    torch.compile itself raises when torch is too old for the running
    python, before tracing anything. And the first call raises when the
    backend cannot build: on a GPU, Inductor goes through Triton, which
    compiles a small CUDA helper using whatever $CC names. On a cluster
    that is often a vendor compiler such as nvc, which rejects the flags
    Triton passes, and the failure appears at the first force evaluation
    rather than at start-up. Setting CC to a plain gcc fixes it.

    Both paths fall back to eager, say so once, and carry on.
    """
    global _compiled_general, _compile_mode, _compile_checked
    if _compile_mode is None:
        return _general_branch(*args)

    if _compiled_general is None:
        try:
            _compiled_general = torch.compile(
                _general_branch, dynamic=(_compile_mode == 'dynamic'))
        except Exception as exc:
            _disable_compile("torch.compile is not available (%s). torch "
                             ">= 2.4 is needed on python 3.12 and >= 2.6 "
                             "on 3.13." % exc)
            return _general_branch(*args)

    if not _compile_checked:
        # guard the first call only: a later failure is a real error and
        # must not be swallowed by a blanket except on every evaluation
        try:
            out = _compiled_general(*args)
            _compile_checked = True
            return out
        except Exception as exc:
            _disable_compile(
                "the compiled kernel failed to build (%s).\n  If this "
                "mentions Triton and a compiler such as nvc, the usual "
                "cause is $CC naming a vendor compiler that rejects the "
                "flags Triton passes. Try CC=gcc."
                % str(exc).splitlines()[0][:160])
            return _general_branch(*args)

    return _compiled_general(*args)


def _disable_compile(reason):
    """_disable_compile: fall back to eager, once, with an explanation"""
    global _compile_mode, _compiled_general, _compile_checked
    sys.stderr.write("compute_stress_force_analytic_torch: %s\n"
                     "  Continuing without torch.compile.\n" % reason)
    _compile_mode = None
    _compiled_general = None
    _compile_checked = False


def resolve_device(device=None):
    """resolve_device: the device to compute on

    An explicit argument wins, then PYDIS_TORCH_DEVICE from the
    environment, then the best available: cuda, then mps, then cpu. The
    environment variable is there so a batch script can steer a run
    without editing the example it submits.
    """
    if device is not None:
        return torch.device(device)
    from_env = os.environ.get('PYDIS_TORCH_DEVICE')
    if from_env:
        return torch.device(from_env)
    if torch.cuda.is_available():
        return torch.device('cuda')
    mps = getattr(torch.backends, 'mps', None)
    if mps is not None and mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


def device_report(device=None):
    """device_report: one line naming the device and what it is"""
    dev = resolve_device(device)
    if dev.type == 'cuda':
        i = dev.index if dev.index is not None else torch.cuda.current_device()
        name = torch.cuda.get_device_name(i)
        gb = torch.cuda.get_device_properties(i).total_memory / 2**30
        return ("torch %s on %s (%s, %.1f GiB)"
                % (torch.__version__, dev, name, gb))
    return "torch %s on %s" % (torch.__version__, dev)


def _mf(f):
    """_mf: contract the four limits, the `@ [1,-1,-1,1]` of the original"""
    return f[:, 0] - f[:, 1] - f[:, 2] + f[:, 3]


def _remote_node_force(x1, x2, x3, x4, bp, b, a, mu, nu, depth=0):
    """_remote_node_force: forces on the four endpoints

    Splits the pairs by branch and sends each subset to the formula it
    needs, rather than evaluating the general formula on everything and
    overwriting the near-parallel rows afterwards.

    The numpy original does the latter, and so did this until it was
    measured. Numpy keeps it because subsetting consistently would mean
    touching every one of several hundred expressions; here the arithmetic
    already broadcasts, so the split is one index at the top and one
    scatter at the end.

    Do not expect much from it. The saving is bounded by the near-parallel
    fraction, and the gather of six inputs plus the scatter of four outputs
    takes most of that back. At 100k pairs, 4 threads: no measurable change
    on the 8% mix a real run produces, and about 13% on the two-thirds mix
    of the accuracy sweep. It is kept because it is strictly better on its
    own terms rather than because it is fast: nothing computed is
    discarded, 1-c^2 no longer needs clamping to keep discarded rows
    finite, and a batch that is all one branch pays no gather at all.

    What actually dominates, measured by disabling the correction branches:
    on near-parallel-heavy input, 74% of the time is the up to four
    recursive general-branch calls each near-parallel pair triggers for its
    corrections, which is the ParaDiS decomposition rather than anything
    about torch. On a realistic mix the corrections cost 1% and the whole
    near-parallel path 6%, so what this kernel spends its time on there is
    the materialised temporaries, not the branching.

    The classification needs t, tp and c, which the general formula needs
    too, so they are computed once here and passed down rather than
    recomputed on the subset.
    """
    n = x1.shape[0]
    if n == 0:
        z = torch.zeros_like(x1)
        return z, z.clone(), z.clone(), z.clone()

    Diff = x4 - x3
    oneoverL = 1.0 / torch.sqrt((Diff * Diff).sum(1))
    t = Diff * oneoverL[:, None]

    Diff = x2 - x1
    oneoverLp = 1.0 / torch.sqrt((Diff * Diff).sum(1))
    tp = Diff * oneoverLp[:, None]

    c = (t * tp).sum(1)
    spindex = (1.0 - c * c) < EPS_PARALLEL
    nsp = int(spindex.sum())

    # the common cases first: a batch that is entirely one branch pays no
    # gather and no scatter at all
    if nsp == 0:
        return _general_dispatch(x1, x2, x3, x4, bp, b, a, mu, nu,
                                 t, tp, c, oneoverL, oneoverLp)
    if nsp == n:
        return _special_remote_node_force(x1, x2, x3, x4, bp, b, a, mu, nu,
                                          EPS_PARALLEL, depth + 1)

    gi = (~spindex).nonzero(as_tuple=True)[0]
    si = spindex.nonzero(as_tuple=True)[0]
    fg = _general_dispatch(x1[gi], x2[gi], x3[gi], x4[gi], bp[gi], b[gi],
                           a, mu, nu, t[gi], tp[gi], c[gi],
                           oneoverL[gi], oneoverLp[gi])
    fs = _special_remote_node_force(x1[si], x2[si], x3[si], x4[si],
                                    bp[si], b[si], a, mu, nu,
                                    EPS_PARALLEL, depth + 1)
    out = []
    for g, sp in zip(fg, fs):
        f = torch.empty((n, 3), dtype=x1.dtype, device=x1.device)
        out.append(f.index_copy(0, gi, g).index_copy(0, si, sp))
    return tuple(out)


def _general_branch(x1, x2, x3, x4, bp, b, a, mu, nu,
                    t, tp, c, oneoverL, oneoverLp):
    """_general_branch: the general formula, for pairs known not to be parallel

    Every pair reaching here has 1-c^2 >= EPS_PARALLEL, so onemc2 needs no
    clamping and nothing computed is thrown away.
    """
    c2 = c * c
    onemc2 = 1.0 - c2
    onemc2inv = 1.0 / onemc2
    txtp = torch.linalg.cross(t, tp, dim=-1)

    R1 = x3 - x1
    R2 = x4 - x2
    d = (R1 * txtp).sum(1) * onemc2inv

    temp1 = torch.stack([(R1 * t).sum(1), (R2 * t).sum(1)], dim=1)
    temp2 = torch.stack([(R1 * tp).sum(1), (R2 * tp).sum(1)], dim=1)
    y = (temp1 - c[:, None] * temp2) * onemc2inv[:, None]
    z = (temp2 - c[:, None] * temp1) * onemc2inv[:, None]

    yin = torch.stack([y[:, 0], y[:, 0], y[:, 1], y[:, 1]], dim=1)
    zin = torch.stack([z[:, 0], z[:, 1], z[:, 0], z[:, 1]], dim=1)

    # the integrals
    a2 = a * a
    a2_d2 = a2 + d * d * onemc2
    y2 = yin * yin
    z2 = zin * zin
    cv = c[:, None]
    c2v = c2[:, None]
    oinv = onemc2inv[:, None]
    ad2 = a2_d2[:, None]

    Ra = torch.sqrt(ad2 + y2 + z2 + 2.0 * yin * zin * cv)
    Rainv = 1.0 / Ra

    Ra_Rdot_tp = Ra + zin + yin * cv
    Ra_Rdot_t = Ra + yin + zin * cv

    log_Ra_Rdot_tp = torch.log(Ra_Rdot_tp)
    ylog_Ra_Rdot_tp = yin * log_Ra_Rdot_tp
    log_Ra_Rdot_t = torch.log(Ra_Rdot_t)
    zlog_Ra_Rdot_t = zin * log_Ra_Rdot_t

    Ra2_R_tpinv = Rainv / Ra_Rdot_tp
    yRa2_R_tpinv = yin * Ra2_R_tpinv
    y2Ra2_R_tpinv = yin * yRa2_R_tpinv

    Ra2_R_tinv = Rainv / Ra_Rdot_t
    zRa2_R_tinv = zin * Ra2_R_tinv
    z2Ra2_R_tinv = zin * zRa2_R_tinv

    denom = 1.0 / torch.sqrt(onemc2 * a2_d2)
    cdenom = (1.0 + c) * denom

    f_003 = (-2.0 * denom[:, None]
             * torch.atan((Ra + yin + zin) * cdenom[:, None]))
    adf_003 = ad2 * f_003
    commonf223 = (cv * Ra - adf_003) * oinv

    f_103 = (cv * log_Ra_Rdot_t - log_Ra_Rdot_tp) * oinv
    f_013 = (cv * log_Ra_Rdot_tp - log_Ra_Rdot_t) * oinv
    f_113 = (cv * adf_003 - Ra) * oinv
    f_203 = zlog_Ra_Rdot_t + commonf223
    f_023 = ylog_Ra_Rdot_tp + commonf223

    commonf225 = f_003 - cv * Rainv
    commonf025 = cv * yRa2_R_tpinv - Rainv
    ycommonf025 = yin * commonf025
    commonf205 = cv * zRa2_R_tinv - Rainv
    zcommonf205 = zin * commonf205
    commonf305 = log_Ra_Rdot_t - (yin - cv * zin) * Rainv - c2v * z2Ra2_R_tinv
    zcommonf305 = zin * commonf305
    commonf035 = (log_Ra_Rdot_tp - (zin - cv * yin) * Rainv
                  - c2v * y2Ra2_R_tpinv)
    tf_113 = 2.0 * f_113

    f_005 = (f_003 - yRa2_R_tpinv - zRa2_R_tinv) / ad2
    f_105 = (Ra2_R_tpinv - cv * Ra2_R_tinv) * oinv
    f_015 = (Ra2_R_tinv - cv * Ra2_R_tpinv) * oinv
    f_115 = (Rainv - cv * (yRa2_R_tpinv + zRa2_R_tinv + f_003)) * oinv
    f_205 = (yRa2_R_tpinv + c2v * zRa2_R_tinv + commonf225) * oinv
    f_025 = (zRa2_R_tinv + c2v * yRa2_R_tpinv + commonf225) * oinv
    f_215 = (f_013 - ycommonf025 + cv * (zcommonf205 - f_103)) * oinv
    f_125 = (f_103 - zcommonf205 + cv * (ycommonf025 - f_013)) * oinv
    f_225 = (f_203 - zcommonf305 + cv * (y2 * commonf025 - tf_113)) * oinv
    f_305 = (y2Ra2_R_tpinv + cv * commonf305 + 2.0 * f_103) * oinv
    f_035 = (z2Ra2_R_tinv + cv * commonf035 + 2.0 * f_013) * oinv
    f_315 = (tf_113 - y2 * commonf025 + cv * (zcommonf305 - f_203)) * oinv
    f_135 = ((tf_113 - z2 * commonf205
              + cv * (yin * commonf035 - f_023)) * oinv)

    # Fintegrals[:, k-1] below, so this list is in the 1-based order the
    # original indexes it by
    Fi = torch.stack([_mf(f_003), _mf(f_103), _mf(f_013), _mf(f_113),
                      _mf(f_203), _mf(f_023), _mf(f_005), _mf(f_105),
                      _mf(f_015), _mf(f_115), _mf(f_205), _mf(f_025),
                      _mf(f_215), _mf(f_125), _mf(f_225), _mf(f_305),
                      _mf(f_035), _mf(f_315), _mf(f_135)], dim=1)

    # dot and cross products for the coefficients
    m4p = 0.25 * mu / np.pi
    m4pd = m4p * d
    m8p = 0.5 * m4p
    m8pd = m8p * d
    m4pn = m4p / (1.0 - nu)
    m4pnd = m4pn * d
    m4pnd2 = m4pnd * d
    m4pnd3 = m4pnd2 * d
    a2m4pnd = a2 * m4pnd
    a2m8pd = a2 * m8pd
    a2m4pn = a2 * m4pn
    a2m8p = a2 * m8p

    # The original also forms txbp, tpxb, txtpxbp, tpxtxb, txtpxbpxt and
    # tpxtxbxtp here and never uses any of them. In numpy that is six wasted
    # array temporaries; on a GPU it is six kernel launches and six
    # allocations per call, so they are dropped.
    tpxt = -txtp
    bxt = torch.linalg.cross(b, t, dim=-1)
    bpxtp = torch.linalg.cross(bp, tp, dim=-1)

    tdb = (t * b).sum(1)
    tdbp = (t * bp).sum(1)
    tpdb = (tp * b).sum(1)
    tpdbp = (tp * bp).sum(1)
    txtpdb = (txtp * b).sum(1)
    tpxtdbp = (tpxt * bp).sum(1)
    txbpdtp = tpxtdbp
    tpxbdt = txtpdb

    bpxtpdb = (bpxtp * b).sum(1)
    bxtdbp = (bxt * bp).sum(1)
    txbpdb = bxtdbp
    tpxbdbp = bpxtpdb

    txtpxt = tp - cv * t
    tpxtxtp = t - cv * tp
    txbpxt = bp - tdbp[:, None] * t
    tpxbxtp = b - tpdb[:, None] * tp
    bpxtpxt = tdbp[:, None] * tp - cv * bp
    bxtxtp = tpdb[:, None] * t - cv * b
    txtpxbpdtp = tdbp - tpdbp * c
    tpxtxbdt = tpdb - tdb * c
    txtpxbpdb = tdbp * tpdb - tpdbp * tdb
    tpxtxbdbp = txtpxbpdb

    # coefficients for f3 and f4
    temp1 = tdbp * tpdb + txtpxbpdb
    I00a = temp1[:, None] * tpxt
    I00b = bxt * txtpxbpdtp[:, None]
    temp1 = m4pnd * txtpdb
    temp2 = m4pnd * bpxtpdb
    I_003 = (m4pd[:, None] * I00a - m4pnd[:, None] * I00b
             + temp1[:, None] * bpxtpxt + temp2[:, None] * txtpxt)
    temp1 = m4pnd3 * txtpxbpdtp * txtpdb
    I_005 = (a2m8pd[:, None] * I00a - a2m4pnd[:, None] * I00b
             - temp1[:, None] * txtpxt)
    I10a = txbpxt * tpdb[:, None] - txtp * txbpdb[:, None]
    I10b = bxt * txbpdtp[:, None]
    temp1 = m4pn * tdb
    I_103 = temp1[:, None] * bpxtpxt + m4p * I10a - m4pn * I10b
    temp1 = m4pnd2 * (txbpdtp * txtpdb + txtpxbpdtp * tdb)
    I_105 = a2m8p * I10a - a2m4pn * I10b - temp1[:, None] * txtpxt
    I01a = txtp * bpxtpdb[:, None] - bpxtpxt * tpdb[:, None]
    temp1 = m4pn * tpdb
    temp2 = m4pn * bpxtpdb
    I_013 = m4p * I01a + temp1[:, None] * bpxtpxt - temp2[:, None] * txtp
    temp1 = m4pnd2 * txtpxbpdtp * tpdb
    temp2 = m4pnd2 * txtpxbpdtp * txtpdb
    I_015 = (a2m8p * I01a - temp1[:, None] * txtpxt + temp2[:, None] * txtp)
    temp1 = m4pnd * txbpdtp * tdb
    I_205 = -temp1[:, None] * txtpxt
    temp1 = m4pnd * txtpxbpdtp * tpdb
    I_025 = temp1[:, None] * txtp
    temp1 = m4pnd * (txtpxbpdtp * tdb + txbpdtp * txtpdb)
    temp2 = m4pnd * txbpdtp * tpdb
    I_115 = temp1[:, None] * txtp - temp2[:, None] * txtpxt
    temp1 = m4pn * txbpdtp * tdb
    I_215 = temp1[:, None] * txtp
    temp1 = m4pn * txbpdtp * tpdb
    I_125 = temp1[:, None] * txtp

    def _assemble(Fint, I003, I103, I013, I005, I105, I015,
                  I115, I205, I025, I215, I125, scale):
        f = (I003 * Fint[0][:, None] + I103 * Fint[1][:, None]
             + I013 * Fint[2][:, None] + I005 * Fint[3][:, None]
             + I105 * Fint[4][:, None] + I015 * Fint[5][:, None]
             + I115 * Fint[6][:, None] + I205 * Fint[7][:, None]
             + I025 * Fint[8][:, None] + I215 * Fint[9][:, None]
             + I125 * Fint[10][:, None])
        return f * scale[:, None]

    y0, y1 = y[:, 0], y[:, 1]
    f4 = _assemble(
        [Fi[:, 1] - y0*Fi[:, 0], Fi[:, 4] - y0*Fi[:, 1],
         Fi[:, 3] - y0*Fi[:, 2], Fi[:, 7] - y0*Fi[:, 6],
         Fi[:, 10] - y0*Fi[:, 7], Fi[:, 9] - y0*Fi[:, 8],
         Fi[:, 12] - y0*Fi[:, 9], Fi[:, 15] - y0*Fi[:, 10],
         Fi[:, 13] - y0*Fi[:, 11], Fi[:, 17] - y0*Fi[:, 12],
         Fi[:, 14] - y0*Fi[:, 13]],
        I_003, I_103, I_013, I_005, I_105, I_015, I_115, I_205, I_025,
        I_215, I_125, oneoverL)
    f3 = _assemble(
        [y1*Fi[:, 0] - Fi[:, 1], y1*Fi[:, 1] - Fi[:, 4],
         y1*Fi[:, 2] - Fi[:, 3], y1*Fi[:, 6] - Fi[:, 7],
         y1*Fi[:, 7] - Fi[:, 10], y1*Fi[:, 8] - Fi[:, 9],
         y1*Fi[:, 9] - Fi[:, 12], y1*Fi[:, 10] - Fi[:, 15],
         y1*Fi[:, 11] - Fi[:, 13], y1*Fi[:, 12] - Fi[:, 17],
         y1*Fi[:, 13] - Fi[:, 14]],
        I_003, I_103, I_013, I_005, I_105, I_015, I_115, I_205, I_025,
        I_215, I_125, oneoverL)

    # coefficients for f1 and f2
    temp1 = tpdb * tdbp + tpxtxbdbp
    I00a = temp1[:, None] * txtp
    I00b = bpxtp * tpxtxbdt[:, None]
    temp1 = m4pnd * tpxtdbp
    temp2 = m4pnd * bxtdbp
    I_003 = (m4pd[:, None] * I00a - m4pnd[:, None] * I00b
             + temp1[:, None] * bxtxtp + temp2[:, None] * tpxtxtp)
    temp1 = m4pnd3 * tpxtxbdt * tpxtdbp
    I_005 = (a2m8pd[:, None] * I00a - a2m4pnd[:, None] * I00b
             - temp1[:, None] * tpxtxtp)
    I01a = tpxt * tpxbdbp[:, None] - tpxbxtp * tdbp[:, None]
    I01b = -bpxtp * tpxbdt[:, None]
    temp1 = m4pn * tpdbp
    I_013 = -temp1[:, None] * bxtxtp + m4p * I01a - m4pn * I01b
    temp1 = m4pnd2 * (tpxbdt * tpxtdbp + tpxtxbdt * tpdbp)
    I_015 = a2m8p * I01a - a2m4pn * I01b + temp1[:, None] * tpxtxtp
    I10a = bxtxtp * tdbp[:, None] - tpxt * bxtdbp[:, None]
    temp1 = m4pn * tdbp
    temp2 = m4pn * bxtdbp
    I_103 = m4p * I10a - temp1[:, None] * bxtxtp + temp2[:, None] * tpxt
    temp1 = m4pnd2 * tpxtxbdt * tdbp
    temp2 = m4pnd2 * tpxtxbdt * tpxtdbp
    I_105 = (a2m8p * I10a + temp1[:, None] * tpxtxtp - temp2[:, None] * tpxt)
    temp1 = m4pnd * tpxbdt * tpdbp
    I_025 = -temp1[:, None] * tpxtxtp
    temp1 = m4pnd * tpxtxbdt * tdbp
    I_205 = temp1[:, None] * tpxt
    temp1 = m4pnd * (tpxtxbdt * tpdbp + tpxbdt * tpxtdbp)
    temp2 = m4pnd * tpxbdt * tdbp
    I_115 = temp1[:, None] * tpxt - temp2[:, None] * tpxtxtp
    temp1 = m4pn * tpxbdt * tpdbp
    I_125 = -temp1[:, None] * tpxt
    temp1 = m4pn * tpxbdt * tdbp
    I_215 = -temp1[:, None] * tpxt

    z0, z1 = z[:, 0], z[:, 1]
    f1 = _assemble(
        [Fi[:, 2] - z1*Fi[:, 0], Fi[:, 3] - z1*Fi[:, 1],
         Fi[:, 5] - z1*Fi[:, 2], Fi[:, 8] - z1*Fi[:, 6],
         Fi[:, 9] - z1*Fi[:, 7], Fi[:, 11] - z1*Fi[:, 8],
         Fi[:, 13] - z1*Fi[:, 9], Fi[:, 12] - z1*Fi[:, 10],
         Fi[:, 16] - z1*Fi[:, 11], Fi[:, 14] - z1*Fi[:, 12],
         Fi[:, 18] - z1*Fi[:, 13]],
        I_003, I_103, I_013, I_005, I_105, I_015, I_115, I_205, I_025,
        I_215, I_125, oneoverLp)
    f2 = _assemble(
        [z0*Fi[:, 0] - Fi[:, 2], z0*Fi[:, 1] - Fi[:, 3],
         z0*Fi[:, 2] - Fi[:, 5], z0*Fi[:, 6] - Fi[:, 8],
         z0*Fi[:, 7] - Fi[:, 9], z0*Fi[:, 8] - Fi[:, 11],
         z0*Fi[:, 9] - Fi[:, 13], z0*Fi[:, 10] - Fi[:, 12],
         z0*Fi[:, 11] - Fi[:, 16], z0*Fi[:, 12] - Fi[:, 14],
         z0*Fi[:, 13] - Fi[:, 18]],
        I_003, I_103, I_013, I_005, I_105, I_015, I_115, I_205, I_025,
        I_215, I_125, oneoverLp)

    return f1, f2, f3, f4


def _special_remote_node_force(x1, x2, x3, x4, bp, b, a, mu, nu, ecrit,
                               depth=0):
    """_special_remote_node_force: the near-parallel formula

    Mirrors SpecialRemoteNodeForce. Segments too close to parallel for the
    general expression are replaced by a segment lying exactly along t,
    tilted by cotanthetac, and the difference is added back through two
    calls to the general formula on the leftover pieces.

    Two departures from the original, neither of them algebraic:

      - the original swaps x1/x2 and negates bp in place when c < 0, then
        swaps back at the end. Here that is done functionally with
        torch.where, so a caller's tensors are never written through.
      - the original selects the rows needing the correction with a python
        loop over every pair. Here it is a boolean mask, which is the
        whole point of running on a device.
    """
    if depth > _MAX_DEPTH:
        raise RuntimeError(
            "near-parallel correction recursed deeper than %d levels; the "
            "input probably contains a degenerate segment pair" % _MAX_DEPTH)

    n = x1.shape[0]
    if n == 0:
        zz = torch.zeros_like(x1)
        return zz, zz.clone(), zz.clone(), zz.clone()

    eps = _CORRECTION_EPS.get(x1.dtype, 1e-16)
    cotanthetac = float(np.sqrt((1.0 - ecrit * 1.01) / (ecrit * 1.01)))

    Diff = x4 - x3
    oneoverL = 1.0 / torch.sqrt((Diff * Diff).sum(1))
    t = Diff * oneoverL[:, None]

    Diff = x2 - x1
    oneoverLp = 1.0 / torch.sqrt((Diff * Diff).sum(1))
    tp = Diff * oneoverLp[:, None]
    c = (t * tp).sum(1)

    # orient the first segment with the second, functionally
    flip = (c < 0)[:, None]
    x1o, x2o = x1, x2
    x1 = torch.where(flip, x2o, x1o)
    x2 = torch.where(flip, x1o, x2o)
    tp = torch.where(flip, -tp, tp)
    bp = torch.where(flip, -bp, bp)

    a2 = a * a
    m4p = 0.25 * mu / np.pi
    m8p = 0.5 * m4p
    m4pn = m4p / (1.0 - nu)
    a2m4pn = a2 * m4pn
    a2m8p = a2 * m8p

    # ---- first half: modify segment (x1,x2), forces on x3 and x4 ----
    temp = ((x2 - x1) * t).sum(1)
    x2mod = x1 + temp[:, None] * t
    diff = x2 - x2mod
    magdiff = torch.sqrt((diff * diff).sum(1))
    temp = (0.5 * cotanthetac * magdiff)[:, None] * t
    x1mod = x1 + 0.5 * diff + temp
    x2mod = x2mod + 0.5 * diff - temp

    R = x3 - x1mod
    Rdt = (R * t).sum(1)
    nd = R - Rdt[:, None] * t
    d2 = (nd * nd).sum(1)

    r4 = (x4 * t).sum(1)
    r3 = (x3 * t).sum(1)
    s2 = (x2mod * t).sum(1)
    s1 = (x1mod * t).sum(1)

    y = torch.stack([r3, r3, r4, r4], dim=1)
    z = torch.stack([-s1, -s2, -s1, -s2], dim=1)
    Fi = _parallel_integrals(y, z, a2, d2)

    tdb = (t * b).sum(1)
    tdbp = (t * bp).sum(1)
    nddb = (nd * b).sum(1)
    bxt = torch.linalg.cross(b, t, dim=-1)
    bpxt = torch.linalg.cross(bp, t, dim=-1)
    ndxt = torch.linalg.cross(nd, t, dim=-1)
    bpxtdb = (bpxt * b).sum(1)
    bpxtdnd = (bpxt * nd).sum(1)
    bpxtxt = tdbp[:, None] * t - bp

    I_003 = (m4pn * (nddb[:, None] * bpxtxt + bpxtdb[:, None] * ndxt
                     - bpxtdnd[:, None] * bxt)
             - (m4p * tdb * tdbp)[:, None] * nd)
    I_113 = ((m4pn - m4p) * tdb)[:, None] * bpxtxt
    I_005 = (-(a2m8p * tdb * tdbp)[:, None] * nd
             - (a2m4pn * bpxtdnd)[:, None] * bxt
             - (m4pn * bpxtdnd * nddb)[:, None] * ndxt)
    I_115 = (-(a2m8p * tdb)[:, None] * bpxtxt
             - (m4pn * bpxtdnd * tdb)[:, None] * ndxt)

    y0, y2_ = y[:, 0], y[:, 2]
    f4 = ((I_003 * (Fi[:, 1] - y0*Fi[:, 0])[:, None]
           + I_113 * (Fi[:, 4] - y0*Fi[:, 3])[:, None]
           + I_005 * (Fi[:, 7] - y0*Fi[:, 6])[:, None]
           + I_115 * (Fi[:, 10] - y0*Fi[:, 9])[:, None]) * oneoverL[:, None])
    f3 = ((I_003 * (y2_*Fi[:, 0] - Fi[:, 1])[:, None]
           + I_113 * (y2_*Fi[:, 3] - Fi[:, 4])[:, None]
           + I_005 * (y2_*Fi[:, 6] - Fi[:, 7])[:, None]
           + I_115 * (y2_*Fi[:, 9] - Fi[:, 10])[:, None]) * oneoverL[:, None])

    cond = (diff * diff).sum(1) > eps * ((x2mod * x2mod).sum(1)
                                         + (x1mod * x1mod).sum(1))
    if bool(cond.any()):
        k = cond.nonzero(as_tuple=True)[0]
        _, _, f3a, f4a = _remote_node_force(
            x1[k], x1mod[k], x3[k], x4[k], bp[k], b[k], a, mu, nu, depth + 1)
        _, _, f3b, f4b = _remote_node_force(
            x2mod[k], x2[k], x3[k], x4[k], bp[k], b[k], a, mu, nu, depth + 1)
        f3 = f3.index_add(0, k, f3a + f3b)
        f4 = f4.index_add(0, k, f4a + f4b)

    # ---- second half: modify segment (x3,x4), forces on x1 and x2 ----
    temp = ((x4 - x3) * tp).sum(1)
    x4mod = x3 + temp[:, None] * tp
    diff = x4 - x4mod
    magdiff = torch.sqrt((diff * diff).sum(1))
    temp = (0.5 * cotanthetac * magdiff)[:, None] * tp
    x3mod = x3 + 0.5 * diff + temp
    x4mod = x4mod + 0.5 * diff - temp

    R = x3mod - x1
    Rdtp = (R * tp).sum(1)
    nd = R - Rdtp[:, None] * tp
    d2 = (nd * nd).sum(1)

    r4 = (x4mod * tp).sum(1)
    r3 = (x3mod * tp).sum(1)
    s2 = (x2 * tp).sum(1)
    s1 = (x1 * tp).sum(1)

    y = torch.stack([r3, r3, r4, r4], dim=1)
    z = torch.stack([-s1, -s2, -s1, -s2], dim=1)
    Fi = _parallel_integrals(y, z, a2, d2)

    tpdb = (tp * b).sum(1)
    tpdbp = (tp * bp).sum(1)
    nddbp = (nd * bp).sum(1)
    bxtp = torch.linalg.cross(b, tp, dim=-1)
    bpxtp = torch.linalg.cross(bp, tp, dim=-1)
    ndxtp = torch.linalg.cross(nd, tp, dim=-1)
    bxtpdbp = (bxtp * bp).sum(1)
    bxtpdnd = (bxtp * nd).sum(1)
    bxtpxtp = tpdb[:, None] * tp - b

    I_003 = (m4pn * (nddbp[:, None] * bxtpxtp + bxtpdbp[:, None] * ndxtp
                     - bxtpdnd[:, None] * bpxtp)
             - (m4p * tpdbp * tpdb)[:, None] * nd)
    I_113 = ((m4pn - m4p) * tpdbp)[:, None] * bxtpxtp
    I_005 = (-(a2m8p * tpdbp * tpdb)[:, None] * nd
             - (a2m4pn * bxtpdnd)[:, None] * bpxtp
             - (m4pn * bxtpdnd * nddbp)[:, None] * ndxtp)
    I_115 = (-(a2m8p * tpdbp)[:, None] * bxtpxtp
             - (m4pn * bxtpdnd * tpdbp)[:, None] * ndxtp)

    z0, z1 = z[:, 0], z[:, 1]
    f2 = ((I_003 * (Fi[:, 2] - z0*Fi[:, 0])[:, None]
           + I_113 * (Fi[:, 5] - z0*Fi[:, 3])[:, None]
           + I_005 * (Fi[:, 8] - z0*Fi[:, 6])[:, None]
           + I_115 * (Fi[:, 11] - z0*Fi[:, 9])[:, None]) * oneoverLp[:, None])
    f1 = ((I_003 * (z1*Fi[:, 0] - Fi[:, 2])[:, None]
           + I_113 * (z1*Fi[:, 3] - Fi[:, 5])[:, None]
           + I_005 * (z1*Fi[:, 6] - Fi[:, 8])[:, None]
           + I_115 * (z1*Fi[:, 9] - Fi[:, 11])[:, None]) * oneoverLp[:, None])

    cond = (diff * diff).sum(1) > eps * ((x4mod * x4mod).sum(1)
                                         + (x3mod * x3mod).sum(1))
    if bool(cond.any()):
        k = cond.nonzero(as_tuple=True)[0]
        _, _, f1a, f2a = _remote_node_force(
            x3[k], x3mod[k], x1[k], x2[k], b[k], bp[k], a, mu, nu, depth + 1)
        _, _, f1b, f2b = _remote_node_force(
            x4mod[k], x4[k], x1[k], x2[k], b[k], bp[k], a, mu, nu, depth + 1)
        f1 = f1.index_add(0, k, f1a + f1b)
        f2 = f2.index_add(0, k, f2a + f2b)

    # undo the orientation swap: the forces follow the endpoints
    f1o, f2o = f1, f2
    f1 = torch.where(flip, f2o, f1o)
    f2 = torch.where(flip, f1o, f2o)
    return f1, f2, f3, f4


def _parallel_integrals(y, z, a2, d2):
    """_parallel_integrals: the twelve f_* integrals of the parallel formula

    Shared by both halves of _special_remote_node_force, which evaluate the
    same expressions with the roles of the two segments exchanged. Returns
    them contracted over the four limits, as (N, 12) in the 1-based order
    the original indexes.
    """
    a2_d2 = a2 + d2
    a2d2inv = (1.0 / a2_d2)[:, None]
    ypz = y + z
    ymz = y - z
    Ra = torch.sqrt(a2_d2[:, None] + ypz * ypz)
    Rainv = 1.0 / Ra
    Log_Ra_ypz = torch.log(Ra + ypz)

    f_003 = Ra * a2d2inv
    f_103 = -0.5 * (Log_Ra_ypz - ymz * Ra * a2d2inv)
    f_013 = -0.5 * (Log_Ra_ypz + ymz * Ra * a2d2inv)
    f_113 = -Log_Ra_ypz
    f_213 = z * Log_Ra_ypz - Ra
    f_123 = y * Log_Ra_ypz - Ra

    f_005 = a2d2inv * (2.0 * a2d2inv * Ra - Rainv)
    f_105 = a2d2inv * (a2d2inv * ymz * Ra - y * Rainv)
    f_015 = -a2d2inv * (a2d2inv * ymz * Ra + z * Rainv)
    f_115 = -a2d2inv * ypz * Rainv
    f_215 = Rainv - z * f_115
    f_125 = Rainv - y * f_115

    return torch.stack([_mf(f_003), _mf(f_103), _mf(f_013), _mf(f_113),
                        _mf(f_213), _mf(f_123), _mf(f_005), _mf(f_105),
                        _mf(f_015), _mf(f_115), _mf(f_215), _mf(f_125)],
                       dim=1)


def torch_segseg_force_vec(p1, p2, p3, p4, b1, b2, mu, nu, a,
                           device=None, dtype=None, allow_float32=False):
    """torch_segseg_force_vec: segment-segment forces, evaluated on a device

    Same arguments and same return as python_segseg_force_vec, plus an
    optional device and dtype. numpy arrays in gives numpy arrays out;
    torch tensors in gives tensors back on their own device, so a caller
    already holding device data pays no transfer.

    Anything other than float64 is refused unless allow_float32 is set;
    see the module docstring for the measurements behind that.
    """
    was_numpy = not torch.is_tensor(p1)
    if was_numpy:
        dev = resolve_device(device)
        dt = dtype or torch.float64
        cvt = lambda v: torch.as_tensor(np.ascontiguousarray(v), dtype=dt,
                                        device=dev)
    else:
        dev = p1.device if device is None else resolve_device(device)
        dt = dtype or p1.dtype
        cvt = lambda v: v.to(device=dev, dtype=dt)

    if dt != torch.float64 and not allow_float32:
        raise ValueError(
            "torch_segseg_force_vec: refusing to run in %s. This kernel "
            "cancels quantities orders larger than its result, and in "
            "single precision the answer is wrong in the fourth digit "
            "in the third digit even where the geometry is benign "
            "(median relative error 5.4e-03, worst 1.8e+04), with no "
            "non-finite value to signal it. Pass allow_float32=True if "
            "you are measuring that rather than relying on it."
            % dt)

    x1, x2, x3, x4 = cvt(p1), cvt(p2), cvt(p3), cvt(p4)
    bp, b = cvt(b1), cvt(b2)

    f1, f2, f3, f4 = _remote_node_force(x1, x2, x3, x4, bp, b,
                                        float(a), float(mu), float(nu))
    if was_numpy:
        return (f1.cpu().numpy(), f2.cpu().numpy(),
                f3.cpu().numpy(), f4.cpu().numpy())
    return f1, f2, f3, f4
