"""Reading per-node quantities out of the state dictionary, by tag.

Forces, velocities and old positions travel in `state` as arrays with a
companion array of tags. Every consumer therefore has to pair the two,
and doing it by array position is wrong: the network is reindexed by
export/import round trips and by any operation that adds or removes a
node, so entry i does not reliably describe node i.

These helpers turn each pair of arrays into a dict keyed by tag once, so
callers look nodes up by identity and never by offset.
"""

import numpy as np


def _by_tag(tags, values):
    """_by_tag: dict from a tag array and a matching value array"""
    tags = np.asarray(tags)
    values = np.asarray(values, dtype=float)
    if tags.shape[0] != values.shape[0]:
        raise ValueError("_by_tag: %d tags but %d values"
                         % (tags.shape[0], values.shape[0]))
    return {(int(t[0]), int(t[1])): v for t, v in zip(tags, values)}


def velocities_by_tag(state):
    """velocities_by_tag: nodal velocities keyed by tag, or None

    None when the mobility law has not run, which callers that need
    velocities should treat as an error rather than as zero velocity.
    """
    if "nodevels" not in state or "nodeveltags" not in state:
        return None
    return _by_tag(state["nodeveltags"], state["nodevels"])


def old_positions_by_tag(state):
    """old_positions_by_tag: positions at the start of the step, or None

    Accepts either of the two shapes in use: `oldnodes_dict`, as the
    ExaDiS driver saves it, carrying "tags" and "positions" arrays; or
    `oldpos_dict`, already keyed by tag.

    Returns None when neither is present. A rule that needs old positions
    should raise on that rather than substitute current positions, which
    silently turns a retroactive test into no test at all.
    """
    oldnodes = state.get("oldnodes_dict")
    if oldnodes is not None and "tags" in oldnodes:
        return _by_tag(oldnodes["tags"], oldnodes["positions"])

    oldpos = state.get("oldpos_dict")
    if oldpos is not None:
        return {tuple(k): np.asarray(v, dtype=float)
                for k, v in oldpos.items()}
    return None
