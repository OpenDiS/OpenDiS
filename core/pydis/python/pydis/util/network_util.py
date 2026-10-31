"""@package docstring
network_util: functions that build a dislocation network from lists of nodes and segments,
ports of ExaDiS' pyexadis_utils
"""

import numpy as np
from ..disnet import DisNode


def insert_frank_read_src(cell, nodes, segs, burg, plane, length, center, theta=0.0, linedir=None, numnodes=10):
    """insert_frank_read_src: append a Frank-Read source, a straight line of numnodes nodes pinned
       at its two ends, to the lists of nodes [x,y,z,constraint] and segments [node1,node2,burg,plane],
       as pyexadis_utils.insert_frank_read_src; the line is along linedir, or at theta degrees from
       burg in the plane
    """
    plane = plane / np.linalg.norm(plane)
    if np.abs(np.dot(burg, plane)) >= 1e-5:
        print('Warning: Burgers vector and plane normal are not orthogonal')

    if linedir is not None:
        ldir = np.array(linedir)
        ldir = ldir / np.linalg.norm(ldir)
    else:
        b = burg / np.linalg.norm(burg)
        y = np.cross(plane, b)
        y = y / np.linalg.norm(y)
        ldir = np.cos(theta*np.pi/180.0)*b + np.sin(theta*np.pi/180.0)*y

    istart = len(nodes)
    for i in range(numnodes):
        p = center - 0.5*length*ldir + i*length/(numnodes-1)*ldir
        pinned = i == 0 or i == numnodes-1
        constraint = DisNode.Constraints.PINNED_NODE if pinned else DisNode.Constraints.UNCONSTRAINED
        nodes.append(np.concatenate((p, [constraint])))
    for i in range(numnodes-1):
        segs.append(np.concatenate(([istart+i, istart+i+1], burg, plane)))

    return nodes, segs
