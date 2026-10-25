"""Frank-Read source on one FCC slip system, and helpers to inspect it.

The source is the five-node pinned loop of init_frank_read_src_loop in
examples/02_frank_read_src, built in the frame (b_hat, xi_hat, n_hat) of a
chosen slip system instead of (x, y, z). Only the free arm moves: the other
three segments run between pinned nodes, which remesh never refines.
"""

import numpy as np

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell
from pydis.util.fcc_slip_sys import unit

PINNED = DisNode.Constraints.PINNED_NODE
FREE = DisNode.Constraints.UNCONSTRAINED


def slip_frame(slip_sys, sysID):
    """slip_frame: rotation with columns (b_hat, xi_hat, n_hat) of system sysID

    xi_hat = n_hat x b_hat is the edge line direction, which makes the frame
    right-handed since b_hat x xi_hat = n_hat.
    """
    n, b = (unit(v) for v in slip_sys.sysID_to_nb(sysID))
    return np.column_stack((b, np.cross(n, b), n))


def frank_read_loop_local(arm_length):
    """frank_read_loop_local: node positions and constraints in the slip frame

    The free arm runs along local y on the glide plane z = 0; the loop closes
    through z = -arm_length.
    """
    h = arm_length / 2.0
    pos = np.array([[0.0, -h,  0.0],
                    [0.0, 0.0, 0.0],
                    [0.0,  h,  0.0],
                    [0.0,  h, -arm_length],
                    [0.0, -h, -arm_length]])
    constraints = np.array([PINNED, FREE, PINNED, PINNED, PINNED])
    return pos, constraints


def glide_normal(burg, dr):
    """glide_normal: unit normal of the plane holding burg and dr

    A screw segment, dr parallel to burg, gets some plane that holds burg.
    """
    pn = np.cross(burg, dr)
    if np.linalg.norm(pn) < 1e-12 * np.linalg.norm(burg) * np.linalg.norm(dr):
        pn = np.cross(burg, np.eye(3)[np.argmin(np.abs(burg))])
    return unit(pn)


def loop_links(cell, pos, burg):
    """loop_links: links of a closed loop through pos, glide planes from b x dr"""
    N = pos.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        dr = cell.closest_image(Rref=pos[i], R=pos[(i+1) % N]) - pos[i]
        links[i, :] = np.concatenate(([i, (i+1) % N], burg, glide_normal(burg, dr)))
    return links


def init_frank_read_src_on_slip_system(slip_sys, sysID, arm_length=500.0,
                                       box_length=1000.0, pbc=True,
                                       shift=np.zeros(3), pin_all=False):
    """init_frank_read_src_on_slip_system: edge Frank-Read source on system sysID

    Lengths are in units of |b|, so the Burgers vector is b_hat. The loop is
    centred in the cell, translated by shift, and folded into it. pin_all
    pins the free node too, for a network that cannot move at all.
    """
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc, pbc, pbc])
    Q = slip_frame(slip_sys, sysID)
    pos_local, constraints = frank_read_loop_local(arm_length)
    if pin_all:
        constraints = np.full_like(constraints, PINNED)
    pos = cell.fold(pos_local @ Q.T + cell.center() + shift)
    rn = np.column_stack((pos, constraints))
    links = loop_links(cell, pos, Q[:, 0])
    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))


def default_state(box_length=1000.0):
    """default_state: material and discretization of examples/02_frank_read_src"""
    return {"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0,
            "maxseg": 0.04*box_length, "minseg": 0.01*box_length, "rann": 3.0}


def dd_components(net, state, dt=1.0e-8):
    """dd_components: force, mobility, integrator, and topology handlers of 02"""
    from pydis import (CalForce, MobilityLaw, TimeIntegration, Topology,
                       Collision, Remesh, CellList)
    calforce = CalForce(force_mode='LineTension', state=state)
    mobility = MobilityLaw(mobility_law='SimpleGlide', state=state)
    return dict(
        calforce=calforce, mobility=mobility,
        timeint=TimeIntegration(integrator='EulerForward', dt=dt, state=state),
        topology=Topology(split_mode='MaxDiss', state=state,
                          force=calforce, mobility=mobility),
        collision=Collision(collision_mode='Retroactive', state=state,
                            nbrlist=CellList(cell=net.cell, n_div=[8, 8, 8])),
        remesh=Remesh(remesh_rule='LengthBased', state=state))


def uniaxial_voigt(sigma, edir):
    """uniaxial_voigt: sigma * e (x) e as Voigt [xx, yy, zz, yz, xz, xy]"""
    e = unit(edir)
    A = np.outer(e, e)
    return sigma * np.array([A[0, 0], A[1, 1], A[2, 2], A[1, 2], A[0, 2], A[0, 1]])


def free_segments(G):
    """free_segments: (tag1, tag2, edge) of segments with at least one free end"""
    def pinned(t):
        return G.nodes(t).constraint == PINNED
    return [(t1, t2, e) for (t1, t2), e in G.all_segments_dict().items()
            if not (pinned(t1) and pinned(t2))]


def free_arm_sysIDs(G, slip_sys):
    """free_arm_sysIDs: set of slip system indices found on the free segments"""
    return {int(slip_sys.nb_to_sysID(e.plane_normal, e.burg_vec_from(t1))[0])
            for t1, t2, e in free_segments(G)}


def free_line_length(G):
    """free_line_length: total length of the segments that can move"""
    total = 0.0
    for t1, t2, _ in free_segments(G):
        r1 = G.nodes(t1).R
        total += np.linalg.norm(G.cell.closest_image(Rref=r1, R=G.nodes(t2).R) - r1)
    return total


def closed_loop(cell, pos, burg, constraints=None):
    """closed_loop: DisNet of one loop through pos with Burgers vector burg"""
    N = pos.shape[0]
    if constraints is None:
        constraints = np.full(N, FREE)
    rn = np.column_stack((pos, constraints))
    return DisNet(cell=cell, rn=rn, links=loop_links(cell, pos, burg))


def random_loop(rng, nnode=10, box_length=1000.0):
    """random_loop: a polygon around (900, 500, 500) that crosses the x face, random b

    Nodes are jittered about even spacing so every segment stays well under
    half the box; a longer one has no unique minimum image.
    """
    cell = Cell(h=box_length*np.eye(3), is_periodic=[True, True, True])
    t = (np.arange(nnode) + rng.uniform(-0.3, 0.3, nnode)) * (2*np.pi / nnode)
    u, v = np.linalg.qr(rng.normal(size=(3, 3)))[0][:, :2].T
    pos = cell.fold(np.array([900.0, 500.0, 500.0])
                    + 300.0*(np.outer(np.cos(t), u) + np.outer(np.sin(t), v)))
    b = rng.normal(size=3)
    return closed_loop(cell, pos, b / np.linalg.norm(b))
