#include <stdlib.h>
#include "Home.h"

/* Functions to be implemented in the future */

/* 
   Initialize            -- Initialize.c
   FMInit                -- FMComm.c
   GetCellDomainList     -- Decomp.c
   EvaluateMobility      -- Topology.c
   MergeNode             -- Topology.c
   SplitNode             -- Topology.c
   FindSubFSeg           -- LocalSegForces.c
   FindFSegComb          -- LocalSegForces.c
   FindPreciseGlidePlane -- FindPreciseGlidePlane.c
   PickScrewGlidePlane   -- PickScrewGlidePlane.c
*/

/*---------------------------------------------------------------------------
 *
 *      Function:     Initialize
 *      Description:  This is the driver routine for initialization,
 *                    handling some of the general initializations and
 *                    calling all the more specific initialization routines.
 *
 *      Last Modified:  01/09/08: Gregg Hommes - added call to
 *                                VerifyBurgersVectors() as a sanity check.
 *
 *-------------------------------------------------------------------------*/
void Initialize(Home_t *home,int argc, char *argv[])
{
    Param_t *param;
    printf("Stub.c: Initialize being implemented\n");

    home->ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    home->dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    //home->subcyc = (Subcyc_t *)calloc(1, sizeof(Subcyc_t));

    home->param = (Param_t *)calloc(1, sizeof(Param_t));
    param = home->param;

    CtrlParamInit(param, home->ctrlParamList);
    DataParamInit(param, home->dataParamList);

    SetBoxSize(param);
    SetRemainingDefaults(home);

    if (home->myDomain == 0) {
        DisableUnneededParams(home);
    }

    InitCellNatives(home);
    InitCellNeighbors(home);
    InitCellDomains(home);

    InitOpList(home);
    //InitOpRecList(home);

    //ReadRijm(home);
    //ReadRijmPBC(home);

    //SortNativeNodes(home);
    //if (param->useLabFrame) {}

}

/*---------------------------------------------------------------------------
 *
 *      Function:     FMInit
 *      Description:  This function controls allocation and initialization
 *                    of all the FM layers.  It should be called every
 *                    time the dynamic load balancing is invoked.
 *                    Some of the FM layer data is static and will
 *                    be calculated only the first time into this function,
 *                    the dynamic data will be recalculated as necessary.
 *
 *--------------------------------------------------------------------------*/
void FMInit(Home_t *home)
{
    printf("Stub.c: FMInit not yet implemented\n");
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetCellDomainList
 *      Description: Generic function to search the domain decomposition
 *                   for the list of domains intersecting the specified
 *                   cell.  This function will invoke the search function
 *                   appropriate to the type of domain decomposition
 *                   being used.
 *
 *      Arguments:
 *          cellID    ID of cell as returned by EncodeCellIdx().
 *          domCount  Location in which to return to caller the number of
 *                    domains intersecting the specified cell.
 *          domList   Location in which to return to caller the array
 *                    containing the IDs of all domains intersecting the
 *                    specified cell.
 *
 *------------------------------------------------------------------------*/
void GetCellDomainList(Home_t *home, int cellID, int *domCount, int **domList)
{
    printf("Stub.c: GetCellDomainList not yet implemented\n");
}

/*---------------------------------------------------------------------------
 *
 *	Function:	EvaluateMobility
 *	Description:	This function provides a way to invoke the proper
 *			mobility functions for specific nodes rather than
 *			looping through all the nodes native to the current
 *			domain.  This function is used solely to evaluate
 *			temporary nodes created to determine if a many-armed
 *			node need be split.
 *	Arguments:
 *		nodeA	Pointer to first node to evaluate
 *
 *
 *      Returns:
 *              0 if there were no errors
 *              1 if the mobility function failed
 *
 *-------------------------------------------------------------------------*/
int EvaluateMobility(Home_t *home, Node_t *nodeA)
{
    printf("Stub.c: EvaluateMobility not yet implemented\n");
    return(0);
}

/*---------------------------------------------------------------------------
 *
 *	Function:	MergeNode
 *	Description:	This function 'merges' two nodes by moving all
 *			arms of <deadNode> to <targetNode>, and then
 *			completely removing <deadNode>.  If the merge
 *			resulted in self-links from <targetNode> back
 *			to itself, or multiple links from <targetNode>
 *			to any other single node, these redundant links
 *			will be dealt with before returning to the caller.
 *
 *	Arguments:
 *              opClass         Flag indicating the class of topological
 *                              operation invoking this function.
 *                              initiated.  Valid types are:
 *
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *
 *              node1           pointer to first node to be merged
 *              node2           pointer to second node to be merged
 *              position        coordinates (x,y,z) at which final merged
 *                              node is to be placed.
 *              mergedNode      pointer to location in which to return
 *                              pointer to the node resulting from the merge.
 *                              A NULL pointer will be returned if the
 *                              the merge fails.
 *              status          pointer to location in which to return
 *                              a completion status to the caller.  Valid
 *                              statuses are the following, where all statuses
 *                              indicating success may be logically OR'ed
 *
 *                                  MERGE_SUCCESS
 *                                  MERGE_NO_REPOSITION
 *                                  MERGE_NODE_ORPHANED
 *                                  MERGE_NOT_PERMITTED
 *                                  MERGE_DOUBLE_LINK
 *
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *
 *-------------------------------------------------------------------------*/
void MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
               real8 *position, Node_t **mergedNode, int *status, int globalOp)
{
    printf("Stub.c: MergeNode not yet implemented\n");
}

/*---------------------------------------------------------------------------
 *
 *	Function:	SplitNode
 *	Description:	Create a new node and transfer the specified set of
 *			connections from an existing to to the new node.  If
 *			necessary, a new link will also be created between the
 *			existing node and the new node.
 *
 *	Arguments:
 *		node            Pointer to the node to be split
 *              pos1            coordinates at which splitNode1 will be
 *                              left after the split
 *              pos2            coordinates at which splitNode2 will be
 *                              left after the split
 *              vel1            velocity assigned to splitNode1 after the split
 *              vel2            velocity assigned to splitNode2 after the split
 *		armCount	number of arms of the original node selected
 *                              to be split off.
 *		armList		pointer to array of integers indicating the
 *				arms of existing node that are to be split off
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *              splitNode1      ptr to ptr to node to which all unselected arms
 *                              of the original node will be attached after the
 *                              split.  Returned to caller.
 *              splitNode2      ptr to ptr to node to which all selected arms
 *                              of the original node will be attached after the
 *                              after the split.  Returned to caller.
 *              flags           Bit field with additional processing flags
 *
 *	Returns:		1 if the split was successful
 *                              0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int SplitNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
              real8 *pos2, real8 *vel1, real8 *vel2, int armCount,
              int *armList, int globalOp, Node_t **splitNode1,
              Node_t **splitNode2, int flags)
{
    printf("Stub.c: SplitNode not yet implemented\n");
    return 0;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       FindSubFseg
 *      Description:    Given a segment p1-->p2 and the force at each endpoint
 *                      of the segment, estimate the resulting forces on the
 *                      segment pair created by bisecting p1-->p2 at newpos.
 *
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          burg     burgers vector from p1 to p2
 *          oldfp1   force of segment p1-->p2 at point p1
 *          oldfp2   force of segment p1-->p2 at point p2
 *          newpos   coordinates of position along the p1-->p2
 *                   segment at which a new node is to be added.
 *          f0seg1   Resulting forces at p1 on the p1-->newpos segment
 *          f1seg1   Resulting forces at newpos on the p1-->newpos segment
 *          f0seg2   Resulting forces at newpos on the newpos-->p2 segment
 *          f1seg2   Resulting forces at p2 on the newpos-->p2 segment
 *
 *-------------------------------------------------------------------------*/
void FindSubFSeg(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
                 real8 oldfp1[3], real8 oldfp2[3], real8 newpos[3],
                 real8 f0seg1[3], real8 f1seg1[3], real8 f0seg2[3],
                 real8 f1seg2[3])
{
    printf("Stub.c: FindSubFSeg not yet implemented\n");
}

/*---------------------------------------------------------------------------
 *
 *      Function:       FindFSegComb
 *      Description:
 *
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          p3       Coordinates of point 3 (corresponding to the periodic
 *                   image of p3 closest to point p2)
 *          burg1    burgers vector from p1 to p2
 *          burg2    burgers vector from p2 to p3
 * WARNING the fp*seg* arrays are modified!
 *          fp1seg1  force of segment p1-->p2 at point p1
 *          fp2seg1  force of segment p1-->p2 at point p2
 *          fp2seg2  force of segment p2-->p3 at point p2
 *          fp3seg2  force of segment p2-->p3 at point p3
 *          fsegnew  resulting forces at p1 and p3 from the
 *                   segment p1-->p3
 *
 *-------------------------------------------------------------------------*/
void FindFSegComb(Home_t *home, real8 p0[3], real8 p1[3], real8 p2[3],
                  real8 burg1[3], real8 burg2[3], real8 fp0seg1[3],
                  real8 fp1seg1[3], real8 fp1seg2[3], real8 fp2seg2[3],
                  real8 f0new[3], real8 f1new[3])
{
    printf("Stub.c: FindFSegComb not yet implemented\n");
}

/*-------------------------------------------------------------------------
 *
 *      Function:     FindGlidePlane
 *      Description:  If we need to enforce strict glide planes, then
 *                    this dispatch function will call an appropriate
 *                    calculate a glide plane normal from the burgers
 *                    vector and line direction, then explicitly select
 *                    from the permitted glide planes for the burgers
 *                    vector, the glide plane normal closest (in angle)
 *                    to the calculated glide plane and return the
 *                    selecetd plane normal to the caller.
 *
 *                    NOTE: This function converts (if necessary) burgers
 *                    vector and line dir from the crystal frame to the 
 *                    lab frame, so make sure the provided burgers vector
 *                    and linedir have not been previously converted
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          dirIn        Input array containing the 3 components of the
 *                       line direction. (This should not yet be shifted
 *                       from the lab frame to the crystal frame)
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plabne will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
void FindPreciseGlidePlane(Home_t *home, real8 burgVecIn[3], real8 dirIn[3],
                           real8 glidePlane[3])
{
    printf("Stub.c: FindPreciseGlidePlane not yet implemented\n");
}

/*-------------------------------------------------------------------------
 *
 *      Function:     PickScrewGlidePlane
 *      Description:  This is a generic dispatch function which calls an
 *                    appropriate material-specific routine to select, at
 *                    random, an appropriate glide plane for screw
 *                    dislocations in the specified type of material
 *                    crystal structure based on the burgers vector.
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
void PickScrewGlidePlane(Home_t *home, real8 burgVecIn[3], real8 glidePlane[3])
{
    printf("Stub.c: PickScrewGlidePlane not yet implemented\n");
}
