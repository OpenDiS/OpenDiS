/*****************************************************************************
 *
 *      Module:         Remesh.c
 *      Description:    This module contains functions common to more
 *                      than 1 of the supported version of remesh, plus
 *                      a generic entry function that invokes the
 *                      proper remesh version.
 *
 *      Included functions:
 *              CutSurfaceSegments()
 *              EstCoarsenForces()
 *              EstRefinementForces()
 *              Remesh()
 *
 *****************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "Home.h"
#include "Comm.h"
#include "Topology.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*-------------------------------------------------------------------------
 *
 *      Function:    EstRefineMentForces
 *      Description: Estimate the forces on the resultant segments
 *                   when an existing segment is split during
 *                   mesh refinement.
 *
 *      Arguments:
 *          node1    pointer to first endpoint of the original segment
 *          node2    pointer to second endpoint of the original segment
 *          newPos   coordinates of the point at which the node1/node2
 *                   segment will be split
 *          vec      vector from <node1> to <node2>
 *          f0Seg1   3 element array in which to return estimated segment
 *                   force on <node1> from a segment from <node1> to the
 *                   location <newPos>
 *          f1Seg1   3 element array in which to return estimated segment
 *                   force at location <newPos> from a segment from <node1>
 *                   to the location <newPos>
 *          f0Seg2   3 element array in which to return estimated segment
 *                   force at location <newPos> from a segment from <node2>
 *                   to the location <newPos>
 *          f1Seg2   3 element array in which to return estimated segment
 *                   force on <node2> from a segment from <node2> to the
 *                   location <newPos>
 *
 *------------------------------------------------------------------------*/
void EstRefinementForces(Home_t *home, Node_t *node1, Node_t *node2,
                         real8 newPos[3], real8 vec[3],
                         real8 f0Seg1[3], real8 f1Seg1[3],
                         real8 f0Seg2[3], real8 f1Seg2[3])
{
        int     arm12, arm21;
        real8   p1[3], p2[3], oldfp1[3], oldfp2[3], burg[3];

        arm12 = GetArmID(home, node1, node2);
        arm21 = GetArmID(home, node2, node1);

        burg[X] = node1->burgX[arm12];
        burg[Y] = node1->burgY[arm12];
        burg[Z] = node1->burgZ[arm12];

        oldfp1[X] = node1->armfx[arm12];
        oldfp1[Y] = node1->armfy[arm12];
        oldfp1[Z] = node1->armfz[arm12];

        oldfp2[X] = node2->armfx[arm21];
        oldfp2[Y] = node2->armfy[arm21];
        oldfp2[Z] = node2->armfz[arm21];

        p1[X] = node1->x;
        p1[Y] = node1->y;
        p1[Z] = node1->z;

        p2[X] = p1[X] + vec[X];
        p2[Y] = p1[Y] + vec[Y];
        p2[Z] = p1[Z] + vec[Z];

        FindSubFSeg(home, p1, p2, burg, oldfp1, oldfp2, newPos, f0Seg1,
                    f1Seg1, f0Seg2, f1Seg2);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    EstCoarsenForces
 *      Description: Estimate the forces on the resultant segment
 *                   when an existing discretization node is removed
 *                   during mesh coarsening.
 *
 *      Arguments:
 *          node1    pointer to first neighbor of node to be coarsened
 *          node2    pointer to node to be coarsened out
 *          node3    pointer to second neighbor of node to be coarsened
 *          f0Seg1   3 element array in which to return estimated segment
 *                   force on <node1> from a segment from <node1> to <node2>
 *          f1Seg1   3 element array in which to return estimated segment
 *                   force on <node2> from a segment from <node1> to <node2>
 *
 *------------------------------------------------------------------------*/
void EstCoarsenForces(Home_t *home, Node_t *node1, Node_t *node2,
                      Node_t *node3, real8 f0Seg[3], real8 f1Seg[3])
{
        int     i, arm21, arm23, arm12, arm32;
        real8   burg1[3], burg2[3];
        real8   p1[3], p2[3], p3[3];
        real8   fp0[3], fp1[3], fp2[3], fp3[3];
        real8   len1[3], len2[3];
        Param_t *param;

        param = home->param;

        arm21 = GetArmID(home, node2, node1);
        arm23 = GetArmID(home, node2, node3);
        arm12 = GetArmID(home, node1, node2);
        arm32 = GetArmID(home, node3, node2);

        p1[X] = node1->x; p2[X] = node2->x; p3[X] = node3->x;
        p1[Y] = node1->y; p2[Y] = node2->y; p3[Y] = node3->y;
        p1[Z] = node1->z; p2[Z] = node2->z; p3[Z] = node3->z;

/*
 *      If periodic boundaries are enabled, the nodes may
 *      be on opposite side of the problem space, so adjust
 *      the positions accordingly.
 */
        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);
        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p3[X], &p3[Y], &p3[Z]);

        burg1[X] = node1->burgX[arm12];
        burg1[Y] = node1->burgY[arm12];
        burg1[Z] = node1->burgZ[arm12];

        burg2[X] = node2->burgX[arm23];
        burg2[Y] = node2->burgY[arm23];
        burg2[Z] = node2->burgZ[arm23];

        fp0[X] = node1->armfx[arm12];
        fp0[Y] = node1->armfy[arm12];
        fp0[Z] = node1->armfz[arm12];

        fp1[X] = node2->armfx[arm21];
        fp1[Y] = node2->armfy[arm21];
        fp1[Z] = node2->armfz[arm21];

        fp2[X] = node2->armfx[arm23];
        fp2[Y] = node2->armfy[arm23];
        fp2[Z] = node2->armfz[arm23];

        fp3[X] = node3->armfx[arm32];
        fp3[Y] = node3->armfy[arm32];
        fp3[Z] = node3->armfz[arm32];

/*
 *      It is possible to we'll have to deal with extremely short
 *      segments here (zero length segments created during collisions and
 *      short segments being coarsened out).  Since the FindFSegComb()
 *      function does not handle these cases well, we do special handling
 *      of them here.
 *
 *      If we find a zero-length segment, the force on a segment between
 *      nodes 1 and 3 would equal the force on the non-zero length segment
 *      attached to node 2.  Otherwise, use FindFSegComb() to get the force.
 */
        for (i = 0; i < 3; i++) {
            len1[i] = p1[i] - p2[i];
            len2[i] = p2[i] - p3[i];
        }

/*
 *      Segments less than .01b in length are treated as zero length
 */
        if ((len1[X]*len1[X] + len1[Y]*len1[Y] + len1[Z]*len1[Z]) < 1.0e-4) {
            f0Seg[X] = node2->armfx[arm23];
            f0Seg[Y] = node2->armfy[arm23];
            f0Seg[Z] = node2->armfz[arm23];
            f1Seg[X] = node3->armfx[arm32];
            f1Seg[Y] = node3->armfy[arm32];
            f1Seg[Z] = node3->armfz[arm32];
        } else if ((len2[X]*len2[X]+len2[Y]*len2[Y]+len2[Z]*len2[Z])<1.0e-4) {
            f0Seg[X] = node1->armfx[arm12];
            f0Seg[Y] = node1->armfy[arm12];
            f0Seg[Z] = node1->armfz[arm12];
            f1Seg[X] = node2->armfx[arm21];
            f1Seg[Y] = node2->armfy[arm21];
            f1Seg[Z] = node2->armfz[arm21];
        } else {
            FindFSegComb(home, p1, p2, p3, burg1, burg2,
                         fp0, fp1, fp2, fp3, f0Seg, f1Seg);
        }

        return;
}


#ifdef _FEM
/*-------------------------------------------------------------------------
 *
 *      Function:       CutSurfaceSegments
 *      Description:    Loop through all the local nodes looking for 
 *                      surface nodes with connections to other nodes
 *                      on the same surface.
 *
 *-----------------------------------------------------------------------*/
static void CutSurfaceSegments(Home_t *home)
{
        int    i, j, thisDomain, globalOp = 1;
        real8  eps, tmp3[3];
        Node_t *node, *nbr;


        eps = 1.0e-12;
        thisDomain = home->myDomain;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          If the node is not a surface node, skip it.
 */
            if (node->constraint != SURFACE_NODE) {
                continue;
            }

            j = 0;

            while (j < node->numNbrs) {

                nbr = GetNeighborNode(home, node, j);
                if (nbr == (Node_t *)NULL) {
                    j++;
                    continue;
                }

/*
 *              if the neighboring node is not on the surface, skip
 *              to the next neighbor
 */
                if (nbr->constraint != SURFACE_NODE) {
                    j++;
                    continue;
                }
/*
 *              If the surface-normal vectors for the two nodes do not
 *              match the nodes are on different surfaces so we can skip
 *              to the next neighbor.
 */
                cross(node->fem_Surface_Norm, nbr->fem_Surface_Norm, tmp3);

                if ((tmp3[0]*tmp3[0] +
                     tmp3[1]*tmp3[1] +
                     tmp3[2]*tmp3[2]) > eps) {
                    j++;
                    continue;
                }
/*
 *              We found a segment we need to cut.  If the neighbor is
 *              in a remote domain, cut the connection in the lower
 *              numbered domain only.
 */
                if (thisDomain > nbr->myTag.domainID) {
                    j++;
                    continue;
                } 

                ChangeArmBurg(home, node, &node->nbrTag[j], 0, 0, 0,
                              0, 0, 0, globalOp, DEL_SEG_NONE);
                ChangeArmBurg(home, nbr, &node->myTag, 0, 0, 0,
                              0, 0, 0, globalOp, DEL_SEG_NONE);

                if ((nbr->myTag.domainID == thisDomain) &&
                    (nbr->numNbrs == 0)) {
                    RemoveNode(home, nbr, globalOp);
                }
            }
        }

        return;
}
#endif

#if 0
/*-------------------------------------------------------------------------
 *
 *      Function:       Remesh
 *      Description:    This is just a generic function to set up
 *                      for a remesh operation then select and execute
 *                      the proper remesh function.
 *
 *-----------------------------------------------------------------------*/
void Remesh(Home_t *home)
{
        int     i;
        Node_t  *node;
        Param_t *param;

        param = home->param;

        if (param->remeshRule < 0) return;

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, REMESH_START_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, REMESH_START_BARRIER);
#endif
#endif

        TimerStart(home, REMESH);
        ClearOpList(home);
        InitTopologyExemptions(home);

        switch(param->remeshRule) {
        case 2:
            RemeshRule_2(home);
            break;
        case 3:
            RemeshRule_3(home);
            break;
        default:
            Fatal("Remesh: undefined remesh rule %d", param->remeshRule);
            break;
        }
        TimerStop(home, REMESH);

#ifdef _FEM
/*
 *      If free-surfaces are in use rather than PBC, it is possible that
 *      some of the remesh operations have resulted in formation of
 *      segments between two nodes on the same surface.  Any such 
 *      segments must be cut here... This may orphan nodes, but those
 *      should be taken care of by the RemoveOrphanedNodes() call further
 *      on.
 */
        CutSurfaceSegments(home);
#endif

/*
 *      Send to the neighboring domains, a list of all local
 *      operations (add/delete nodes, relinking of arms,
 *      etc) that may affect nodes in remote domains, and process
 *      the remesh operations from the neighboring domains
 */
        TimerStart(home, SEND_REMESH);
        CommSendRemesh(home);
        TimerStop(home, SEND_REMESH);

        TimerStart(home, FIX_REMESH);
        FixRemesh(home);
        TimerStop(home, FIX_REMESH);

/*
 *      Under certain circumstances, parallel topological changes
 *      can create double links between nodes; links which can not be
 *      detected until after FixRemesh() is called... so, a quick
 *      check has to be done to clean up these potential double-
 *      links here, or they will cause problems later on.  Should
 *      only have to check nodes local to this domain.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
#ifdef _OP_REC
            (void)RemoveDoubleLinks(home, node, 2);
#else
            (void)RemoveDoubleLinks(home, node, 0);
#endif
            node->flags &= ~NODE_CHK_DBL_LINK;
        }

/*
 *      It is possible that remesh/coarsen has left an orphaned node.
 *      We need to get rid of any such nodes before the exchange of
 *      node information in CommSendGhosts().
 */
        RemoveOrphanedNodes(home);

#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, REMESH_END_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, REMESH_END_BARRIER);
#endif
#endif

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif

        return;
}
#endif
