import numpy as np

####################################################################################################


def GetMinDist2_python(p1, v1, p2, v2, p3, v3, p4, v4, eps0 = 1.0e-12, epsM = 1e-6):
    """ this calculates the minimum distance between two line segments
    input:
        p1  Coordinates of the first endpoint of segment 1
        v1  Velocity of the node at point p1
        p2  Coordinates of the second endpoint of segment 1
        v2  Velocity of the node at point p2
        p3  Coordinates of the first endpoint of segment 2
        v3  Velocity of the node at point p3
        p4  Coordinates of the second endpoint of segment 2
        v4  Velocity of the node at point p4
    return:
        dist2       pointer to location in which to return the square
                    of the minimum distance between the two points
        ddist2dt    pointter to the location in which to return the
                    time rate of change of the distance between
                    L1 and L2
        L1          pointer to location at which to return the
                    normalized position on seg 1 closest to segment 2
        L2          pointer to location at which to return the
                    normalized position on seg 2 closest to segment 1
    """
    r1mr3 = p1 - p3
    r2mr1 = p2 - p1
    r4mr3 = p4 - p3

    seg1v = v2 - v1
    seg2v = v4 - v3

    M = np.zeros((2,2))
    M[0,0] = np.dot(r2mr1, r2mr1)
    M[1,0] =-np.dot(r4mr3, r2mr1)
    M[1,1] = np.dot(r4mr3, r4mr3)
    M[0,1] = M[1,0]

    rhs = np.array([-np.dot(r2mr1, r1mr3), np.dot(r4mr3, r1mr3)])
    detM = 1.0 - M[1,0] * M[1,0] / (M[0,0] * M[1,1])

    A = M[0,0]
    B = -2.0 * rhs[0]
    C = -2.0 * M[1,0]
    D = -2.0 * rhs[1]
    E = M[1,1]

    didDist2 = False

    if A < eps0: # seg1 is a point
        L1 = 0.0
        if E < eps0: L2 = 0.0
        else: L2 = -0.5 * D / E

    elif E < eps0: # seg2 is a point
        L2 = 0.0
        if A < eps0: L1 = 0.0
        else: L1 = -0.5 * B / A

    elif detM < epsM: # segments are parallel
        r4mr1 = p4 - p1
        r3mr2 = p3 - p2
        r4mr2 = p4 - p2

        dist = np.array([np.dot(r1mr3,r1mr3), np.dot(r4mr1,r4mr1), np.dot(r3mr2,r3mr2), np.dot(r4mr2,r4mr2)])
        pos = np.argmin(dist)
        dist2 = dist[pos]

        L1 = np.floor(pos/2.1)
        L2 = float(1 - pos%2)
        didDist2 = True

    else: # solve the general case
        detM *= M[0,0] * M[1,1]
        sol = np.zeros(2)
        sol[0] = ( M[1,1]*rhs[0] - M[1,0]*rhs[1]) / detM
        sol[1] = (-M[1,0]*rhs[0] + M[0,0]*rhs[1]) / detM

        if sol[0] >= 0 and sol[0] <= 1 and sol[1] >= 0 and sol[1] <= 1:
            # we are done here
            L1, L2 = sol

        else:
            # enumerate four cases
            trial = np.zeros((4,2))
            # alpha = 0
            icase = 0
            trial[icase,0] = 0
            trial[icase,1] = np.clip((rhs[1] - M[1,0]*trial[icase,0]) / M[1][1], 0, 1)

            # alpha = 1
            icase = 1
            trial[icase,0] = 1
            trial[icase,1] = np.clip((rhs[1] - M[1,0]*trial[icase,0]) / M[1][1], 0, 1)

            # beta = 0
            icase = 2
            trial[icase,1] = 0
            trial[icase,0] = np.clip((rhs[0] - M[0,1]*trial[icase,1]) / M[0][0], 0, 1)

            # beta = 1
            icase = 3
            trial[icase,1] = 1
            trial[icase,0] = np.clip((rhs[0] - M[0,1]*trial[icase,1]) / M[0][0], 0, 1)

            # compute distances for each trial case
            dist2_trial = np.zeros(4)
            for icase in range(4):
                L1, L2 = trial[icase]
                dist_vec = p1 + r2mr1 * trial[icase,0] - p3 - r4mr3 * trial[icase,1]
                dist2_trial[icase] = np.dot(dist_vec, dist_vec)

            # find minimum out of four trials
            pos = np.argmin(dist2_trial)
            dist2 = dist2_trial[pos]
            L1, L2 = trial[pos]
            didDist2 = True

    L1 = np.clip(L1, 0, 1)
    L2 = np.clip(L2, 0, 1)

    if not didDist2:
        dist_vec = p1 + r2mr1 * L1 - p3 - r4mr3 * L2
        dist2 = np.dot(dist_vec, dist_vec)

    ddist_vecdt = v1 + seg1v * L1 - v3 - seg2v * L2
    dist_vec    = p1 + r2mr1 * L1 - p3 - r4mr3 * L2

    ddist2dt = 2.0 * np.dot(ddist_vecdt, dist_vec)
    
    return dist2, ddist2dt, L1, L2