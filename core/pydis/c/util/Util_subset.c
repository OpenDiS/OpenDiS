#include <stdlib.h>
#include "Home.h"
#include "QueueOps.h"

/* Functions from ParaDiS Util.c */

/*-------------------------------------------------------------------------
 *
 *      Function:     cross
 *      Description:  Calculates the cross product of two vector
 *                    supplied as arrays by the caller.
 *
 *      Arguments:
 *          a    3 element array containing components of the first vector
 *          b    3 element array containing components of the second vector
 *          c    3 element array in which to return to the caller the
 *               components of the cross product of vectors <a> and <b>.
 *
 *------------------------------------------------------------------------*/
void cross(real8 a[3], real8 b[3], real8 c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Normalize
 *      Description:  Normalize vector a
 *
 *      Arguments:
 *          ax, ay, az  Pointers to the three elements of vector a.
 *                      The normalized values will be returned to the
 *                      caller via these same pointers.
 *
 *------------------------------------------------------------------------*/
void Normalize(real8 *ax, real8 *ay, real8 *az)
{
        real8 a2, a;

        a2 = ((*ax)*(*ax) + (*ay)*(*ay) + (*az)*(*az));

        if (a2 > 0.0) {
            a = sqrt(a2);
            *ax /= a;
            *ay /= a;
            *az /= a;
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     NormalizeVec
 *      Description:  Normalize vector a
 *
 *      Arguments:
 *          vec   Three-element vector to be normalized.  Contents
 *                are updated before control returned to the caller.
 *
 *------------------------------------------------------------------------*/
void NormalizeVec(real8 vec[3])
{
        real8 a2, a;

        a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        if (a2 > 0.0) {
            a = sqrt(a2);
            vec[0] /= a;
            vec[1] /= a;
            vec[2] /= a;
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     EncodeCellIdx
 *      Description:  Given the indices of a cell in each of the X, Y
 *                    and Z dimensions, return its ID encoded as
 *                    a single number. 
 *
 *                    IMPORTANT!  The base set of cells is padded
 *                    by 1 cell on each side of each dimension to
 *                    allow for periodic boundaries.  The encoded
 *                    index returned to the caller will have been
 *                    adjusted to take this layer of ghost cells
 *                    into account
 *
 *                    NOTE: It is assumed in the cell index encoding/
 *                    decoding functions that cells are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          xIndex    Base index of the cell in the X dimension
 *          yIndex    Base index of the cell in the Y dimension
 *          zIndex    Base index of the cell in the Z dimension
 *
 *      Returns:  Cell ID as single index into a 1-dimensional array.
 *
 *------------------------------------------------------------------------*/
int EncodeCellIdx(Home_t *home, int xIndex, int yIndex, int zIndex)
{
        int cellID;

        cellID = zIndex +
                 yIndex * (home->param->nZcells+2) +
                 xIndex * (home->param->nZcells+2) * (home->param->nYcells+2);

        return(cellID);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Connected
 *      Description:  Determines whether or not two nodes are connected
 *                    by a dislocation segment.
 *
 *      Arguments:
 *          node1     Pointer to first node
 *          node2     Pointer to second node
 *          armID     Pointer to location in which to return to the
 *                    caller the index (for node1) of the arm connecting
 *                    the two nodes.  (If the nodes are not connected
 *                    the contents of this pointer will be set to -1)
 *
 *      Returns:   0 if the nodes are not connected.
 *                 1 if the nodes are connected.
 *
 *------------------------------------------------------------------------*/
int Connected(Node_t *node1, Node_t *node2, int *armID)
{
        int i, dom2, idx2;

        *armID = -1;

        if ((node1 == NULL) || (node2 == NULL)) {
            return(0);
        }
    
        dom2 = node2->myTag.domainID;
        idx2 = node2->myTag.index;    

        for (i = 0; i < node1->numNbrs; i++) {
            if ((dom2 == node1->nbrTag[i].domainID) &&
                (idx2==node1->nbrTag[i].index)) {
                *armID = i;
                return(1);
            }
        }

        return 0;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       DomainOwnsSeg
 *      Description:    Determines if the specified domain owns
 *                      the segment beginning in the current domain
 *                      and terminating at the domain indicated by
 *                      <endTag>.
 *
 *                      Note:  Ownership of segments crossing
 *                             domain boundaries alternates on
 *                             and even/odd cycles.  Additionally,
 *                             the ownership is dependent on the
 *                             class of topological operation being
 *                             considered (i.e. during remesh operations
 *                             the ownership is the exact reverse of
 *                             ownership during collision handling.)
 *
 *      Arguments:
 *              opClass         Class of topological operation during which
 *                              this function is being invoked.  Valid
 *                              values are:
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *              thisDomain      domain containing the first endpoint of
 *                              the segment in question.
 *              endTag          pointer to the tag of the second endpoint
 *                              of the segment.
 *
 *      Returns:  1 if <thisDomain> owns the segment
 *                0 in all other cases.
 *
 *-------------------------------------------------------------------------*/
int DomainOwnsSeg(Home_t *home, int opClass, int thisDomain, Tag_t *endTag)
{
        int ownsSeg = 1;

/*
 *      If both endpoints are in the same domain, the domain owns
 *      the segment.
 */
        if (thisDomain == endTag->domainID) {
            return(1);
        }

/*
 *      For collision handling and node separations, ownership
 *      of segments crossing domain boundaries is the lower
 *      numbered domain on even numbered cycles and the higher
 *      numbered domain for odd numbered cycles.
 *
 *      For remesh operations, ownership rules are the opposite
 *      of those used for collision handling.
 */
        switch (opClass) {
        case OPCLASS_SEPARATION:
        case OPCLASS_COLLISION:
            if (home->cycle & 0x01)
                ownsSeg = (thisDomain > endTag->domainID);
            else
                ownsSeg = (thisDomain < endTag->domainID);
            break;
        case OPCLASS_REMESH:
            if (home->cycle & 0x01)
                ownsSeg = (thisDomain < endTag->domainID);
            else
                ownsSeg = (thisDomain > endTag->domainID);
            break;
        default:
            Fatal("Invalid opClass %d in DomainOwnsSeg()", opClass);
            break;
        }

        return(ownsSeg);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Fatal
 *      Description:  Prints a user specified message, aborts all
 *                    other parallel tasks (if any) and self-terminates
 *
 *      Arguments:    This function accepts a variable number of arguments
 *                    in the same fashion as printf(), with the first
 *                    option being the format string which determines
 *                    how the remainder of the arguments are to be
 *                    interpreted.
 *
 *------------------------------------------------------------------------*/
void Fatal(char *format, ...) 
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("Fatal: %s\n", msg);

#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        exit(1);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     FoldBox
 *      Description:  Given a set of coordinates, adjusts the 
 *                    coordinates to the corresponding point within
 *                    the primary (non-periodic) image of the problem
 *                    space.  If the provided coordinates are already
 *                    within the primary image, no adjustments are
 *                    made.
 *
 *      Arguments:
 *          param       Pointer to structure containing control parameters
 *          x, y, z     Pointers to components of coordinate.  On return
 *                      to the caller, these coordinates will have been
 *                      adjusted (if necessary) to the coordinates of the
 *                      corresponding point within the primary image of
 *                      the problem space.
 *
 *------------------------------------------------------------------------*/
void FoldBox(Param_t *param, real8 *x, real8 *y, real8 *z)
{
        real8 xc, yc, zc;

        xc = (param->maxSideX + param->minSideX) * 0.5;
        yc = (param->maxSideY + param->minSideY) * 0.5;
        zc = (param->maxSideZ + param->minSideZ) * 0.5;
    
        if (param->xBoundType == Periodic) {
            *x -= rint((*x-xc)*param->invLx) * param->Lx;
        }

        if (param->yBoundType == Periodic) {
            *y -= rint((*y-yc)*param->invLy) * param->Ly;
        }

        if (param->zBoundType == Periodic) {
            *z -= rint((*z-zc)*param->invLz) * param->Lz;
        }
    
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetArmID
 *      Description: Given two node pointers, return arm ID for node1
 *                   if node2 is its neighbor.
 *
 *      Arguments:
 *          node1    Pointer to first node structure
 *          node2    Pointer to second node structure
 *
 *      Returns:  Non-negative index of the arm of node1 terminating
 *                at node2 if the nodes are connected.
 *                -1 if the two nodes are not connected.
 *                
 *------------------------------------------------------------------------*/
int GetArmID (Home_t *home, Node_t *node1, Node_t *node2)
{
        int i, domID, index;

        if ((node1 == NULL) || (node2 == NULL)) {
            return(-1);
        }
    
        domID = node2->myTag.domainID;
        index = node2->myTag.index;

        for(i = 0; i < node1->numNbrs; i++) {
            if ((domID == node1->nbrTag[i].domainID) &&
                (index == node1->nbrTag[i].index)) {
                return(i);
            }
        }

        return(-1);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     FreeNodeArms
 *      Description:  Release all the specified node's arm related
 *                    arrays, zero out the pointers, and set the arm.
 *                    count to zero.
 *      Arguments:
 *          node    Pointer to node structure in which to free
 *                  all arm related arrays.
 *
 *------------------------------------------------------------------------*/
void FreeNodeArms(Node_t *node)
{
        if (node->numNbrs == 0) {
            return;
        }

        free(node->nbrTag);      node->nbrTag = (Tag_t *)NULL;
        free(node->burgX);       node->burgX = (real8 *)NULL;
        free(node->burgY);       node->burgY = (real8 *)NULL;
        free(node->burgZ);       node->burgZ = (real8 *)NULL;
        free(node->armfx);       node->armfx = (real8 *)NULL;
        free(node->armfy);       node->armfy = (real8 *)NULL;
        free(node->armfz);       node->armfz = (real8 *)NULL;
        free(node->nx);          node->nx = (real8 *)NULL;
        free(node->ny);          node->ny = (real8 *)NULL;
        free(node->nz);          node->nz = (real8 *)NULL;
        free(node->sigbLoc);     node->sigbLoc = (real8 *)NULL;
        free(node->sigbRem);     node->sigbRem = (real8 *)NULL;

        node->numNbrs = 0;

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     GetNeighborNode
 *      Description:  Given a node pointer, return the pointer to
 *                    the n'th neighbor of the node.
 *
 *      Arguments:
 *          node    Pointer to a node.
 *          n       Index in the node's neighbor list of the 
 *                  neighbor to return to the caller.
 *
 *      Returns:    Pointer to the requested neighbor if found.
 *                  NULL in all other cases.
 *
 * FIX ME!  Check if this function will ever be called with sparse arm list 
 *
 *------------------------------------------------------------------------*/
Node_t *GetNeighborNode(Home_t *home, Node_t *node, int n)
{
        int    i, j;
        Node_t *neighbor;
#if 0
/*
 *      New version which assumes no invalid arms in arm list
 */
        if (n >= node->numNbrs) {
            printf("GetNeighborNode: Error finding neighbor %d\n", n);
            PrintNode(node);
            return((Node_t *)NULL);
        }

        neighbor = GetNodeFromTag(home, node->nbrTag[n]);

        return(neighbor);
#else
/*
 *      Old version which assumes the arm list may be sparsely
 *      populated and returns the n'th valid neighbor, which may
 *      not be at index n.  
 */
        j = -1;
 
        for (i = 0; i < node->numNbrs; i++) {

            if (node->nbrTag[i].domainID >= 0) j++;

            if (j == n) {
                neighbor = GetNodeFromTag(home, node->nbrTag[i]);
                return(neighbor);
            }
        }

        printf("GetNeighborNode returning NULL for node (%d,%d) nbr %d\n",
               node->myTag.domainID, node->myTag.index, n);
        PrintNode(node);

        return((Node_t *)NULL);
#endif
}

/*-------------------------------------------------------------------------
 *
 *      Function:     GetNodeFromTag
 *      Description:  Given a node tag, returns a pointer to the
 *                    corresponding node structure.
 *
 *                    NOTE: If the specified tag is for a local node
 *                    and the local node is not found, the function
 *                    will force a code abort with a fatal error.
 *
 *      Arguments:
 *          tag    Tag identifying the desired node
 *
 *      Returns:    Pointer to the requested node if found.
 *                  NULL in all other cases.
 *
 *------------------------------------------------------------------------*/
Node_t *GetNodeFromTag (Home_t *home, Tag_t tag)
{
        Node_t         *node;
        RemoteDomain_t *remDom;
   
        if (tag.domainID < 0 || tag.index < 0) {
            Fatal("GetNodeFromTag: invalid tag (%d,%d)",
                  tag.domainID, tag.index);
        }

/*
 *      If the tag is for a local domain, look up the node in
 *      the local <nodeKeys> array.
 */
        if (tag.domainID == home->myDomain) {

            if (tag.index >= home->newNodeKeyPtr) {
                return((Node_t *)NULL);
            }

            if ((node = home->nodeKeys[tag.index]) == (Node_t *)NULL) {
                return((Node_t *)NULL);
            }

            return(node);

        } else {
/*
 *          If the node is in a remote domain, there are valid situations
 *          in which the current domain does not have information on
 *          either the remote doamin or the remote node.  Hence, it
 *          is not an error to return a NULL pointer.
 */
            remDom = home->remoteDomainKeys[tag.domainID];

            if (remDom == NULL) {
                return((Node_t *)NULL);
            }

            if (tag.index >= remDom->maxTagIndex) {
                return((Node_t *)NULL);
            }
      
            node = remDom->nodeKeys[tag.index];
            return(node);
        }
}

/*-------------------------------------------------------------------------
 *
 *      Function:     InitNodeArm
 *      Description:  Initialize all components of the specified arm
 *                    (segment) of a node.
 *
 *      Arguments:
 *          node   Pointer to the node
 *          armID  Index (0 offset) of the arm/segment to be initialized.
 *
 *------------------------------------------------------------------------*/
static void InitNodeArm(Node_t *node, int armID)
{
        node->nbrTag[armID].domainID = -1;
        node->nbrTag[armID].index = -1;
        node->burgX[armID] = 0.0;
        node->burgY[armID] = 0.0;
        node->burgZ[armID] = 0.0;
        node->armfx[armID] = 0.0;
        node->armfy[armID] = 0.0;
        node->armfz[armID] = 0.0;
        node->nx[armID] = 0.0;
        node->ny[armID] = 0.0;
        node->nz[armID] = 0.0;
        node->sigbLoc[3*armID  ] = 0.0;
        node->sigbLoc[3*armID+1] = 0.0;
        node->sigbLoc[3*armID+2] = 0.0;
        node->sigbRem[3*armID  ] = 0.0;
        node->sigbRem[3*armID+1] = 0.0;
        node->sigbRem[3*armID+2] = 0.0;

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    ReallocNodeArms
 *      Description: Reallocate memory for 'n' arms in the specified 
 *                   node structure preserving any previously 
 *                   existing arm data.
 *
 *      Arguments:
 *          node       Pointer to node for which to allocate arms
 *
 *------------------------------------------------------------------------*/
void ReallocNodeArms(Node_t *node, int n)
{
        int i, origNbrCnt;

        origNbrCnt = node->numNbrs;

        node->numNbrs = n;

        node->nbrTag = (Tag_t *) realloc(node->nbrTag, sizeof(Tag_t)*n);
        node->burgX = (real8 *) realloc(node->burgX, sizeof(real8  )*n);
        node->burgY = (real8 *) realloc(node->burgY, sizeof(real8  )*n);
        node->burgZ = (real8 *) realloc(node->burgZ, sizeof(real8  )*n);
        node->armfx = (real8 *) realloc(node->armfx, sizeof(real8  )*n);
        node->armfy = (real8 *) realloc(node->armfy, sizeof(real8  )*n);
        node->armfz = (real8 *) realloc(node->armfz, sizeof(real8  )*n);
        node->nx = (real8 *) realloc(node->nx, sizeof(real8  )*n);
        node->ny = (real8 *) realloc(node->ny, sizeof(real8  )*n);
        node->nz = (real8 *) realloc(node->nz, sizeof(real8  )*n);
        node->sigbLoc  = (real8 *) realloc(node->sigbLoc, n*3*sizeof(real8));
        node->sigbRem  = (real8 *) realloc(node->sigbRem, n*3*sizeof(real8));

/*
 *      And just initialize the newly allocated arms only, leaving
 *      the previously existing arms as they were.
 */
        for (i = origNbrCnt; i < node->numNbrs; i++) {
            InitNodeArm(node, i);
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ResetSegForces
 *      Description:  Reset the segment forces on nodeA for the
 *                    segment terminating at nodeB, then re-sum the
 *                    total nodal forces for nodeA
 *
 *      Arguments:
 *          nodeA       pointer to node for which to reset seg forces
 *          nodeBtag    pointer to the tag of the node at which the
 *                      segment in question terminates
 *          fx,fy,fz    segment forces 
 *          globalOp    set to 1 if this operation is an operation
 *                      that will have to be passed on to remote
 *                      domains
 *
 *------------------------------------------------------------------------*/
void ResetSegForces(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
                    real8 fx, real8 fy, real8 fz, int globalOp)
{
        int i;
    
/*
 *      If other domains need to be notified of this operation, add the
 *      action to the operation list
 */
        if (globalOp) {
            AddOp(home, RESET_SEG_FORCES,
                nodeA->myTag.domainID,
                nodeA->myTag.index,
                nodeBtag->domainID,
                nodeBtag->index,
                -1,-1,
                0.0, 0.0, 0.0, /* bx, by, bz */
                fx, fy, fz,
                0.0,0.0,0.0);  /* nx, ny, nz */
        }
    
/*
 *      Locate the segment of nodeA terminating at nodeb and update
 *      forces for that segment.
 */
        for (i = 0; i < nodeA->numNbrs; i++) {
            if ((nodeA->nbrTag[i].domainID == nodeBtag->domainID) &&
                (nodeA->nbrTag[i].index == nodeBtag->index)) {
                nodeA->armfx[i] = fx;
                nodeA->armfy[i] = fy;
                nodeA->armfz[i] = fz;
                break;
            }
        }
    
/*
 *      Reset the total forces for nodeA based on all its segment forces
 */
        nodeA->fX = 0;
        nodeA->fY = 0;
        nodeA->fZ = 0;
            
        for (i = 0; i < nodeA->numNbrs; i++) {
            nodeA->fX += nodeA->armfx[i];
            nodeA->fY += nodeA->armfy[i];
            nodeA->fZ += nodeA->armfz[i];
        }
    
        nodeA->flags |= NODE_RESET_FORCES;
    
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     MarkNodeForceObsolete
 *      Description:  Sets a flag for the specified node indicating that
 *                    the nodal force and velocity values are obsolete and
 *                    must be recalculated.  If the node is not local to
 *                    the current domain, the operation list going to
 *                    the remote domains will be modified to include a
 *                    operation to force the owning domain to recalculate
 *                    the node's force and velocity.
 *
 *      Arguments:
 *          node    pointer to node for which the force and velocity
 *                  values must be recomputed
 *
 *------------------------------------------------------------------------*/
void MarkNodeForceObsolete(Home_t *home, Node_t *node)
{
        node->flags |= NODE_RESET_FORCES;

/*
 *      If the node is locally owned, we're done.  If not, we
 *      need to let the owning domain know it needs to recompute
 *      the force/velocity.
 */
        if (node->myTag.domainID == home->myDomain) return;

        AddOp(home, MARK_FORCES_OBSOLETE,
            node->myTag.domainID,
            node->myTag.index,
            -1, -1,
            -1, -1,
            0.0, 0.0, 0.0, /* bx, by, bz */
            0.0, 0.0, 0.0, /* vx, vy, vz */
            0.0,0.0,0.0);  /* nx, ny, nz */

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     PBCPOSITION
 *      Description:  Finds the position of the image of point (x,y,z)
 *                    that is closest to point (x0,y0,z0).
 *
 *                    NOTE: The returned position is not required
 *                    to be in the primary image but may actually be in
 *                    one of the periodic images
 *                  
 *      Arguments:
 *          param       Pointer to structure containing control parameters
 *          x0, y0, z0  Components of position used as a base.
 *          x, y, z     Pointers to components of secondary point.  These
 *                      values will be overwritten with coordinates of the
 *                      image of this point closest to (x0,y0,z0).
 *
 *------------------------------------------------------------------------*/
void PBCPOSITION(Param_t *param, real8 x0, real8 y0, real8 z0,
                 real8 *x, real8 *y, real8 *z)
{
/*
 *      If periodic boundaries are not in use, the provided position
 *      of (x,y,z) will not be adjusted since there are no other
 *      images available.
 */
        if (param->xBoundType == Periodic) { 
            *x -= rint((*x-x0)*param->invLx) * param->Lx;
        }

        if (param->yBoundType == Periodic) {
            *y -= rint((*y-y0)*param->invLy) * param->Ly;
        }

        if (param->zBoundType == Periodic) {
            *z -= rint((*z-z0)*param->invLz) * param->Lz;
        }
    
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ZImage
 *      Description:  Finds the minimum image of (x,y,z) in the
 *                    problem box.
 *
 *                    Typical use of this function specifies (x,y,z)
 *                    as vector from a source point to a secondary
 *                    point.  Upon return to the caller, the vector
 *                    will have been adjusted to be a vector from the
 *                    source point to the closest image of the secondary
 *                    point.
 *                  
 *      Arguments:
 *          param       Pointer to structure containing control parameters
 *          x, y, z     Pointers to components of point (or vector).
 *
 *------------------------------------------------------------------------*/
void ZImage(Param_t *param, real8 *x, real8 *y, real8 *z)
{
/*
 *      If periodic boundaries are not in use, the provided position
 *      of (x,y,z) will not be adjusted since there are no other
 *      images available.
 */
        if (param->xBoundType == Periodic) {
            *x -= rint(*x * param->invLx) * param->Lx;
        }

        if (param->yBoundType == Periodic) {
            *y -= rint(*y * param->invLy) * param->Ly;
        }

        if (param->zBoundType == Periodic) {
            *z -= rint(*z * param->invLz) * param->Lz;
        }
    
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     PrintNode
 *      Description:  For the specified node, print some interesting
 *                    items of data.
 *
 *------------------------------------------------------------------------*/
void PrintNode(Node_t *node)
{
        int i;

        if (node == (Node_t *)NULL) return;

#if 1
        printf("  node(%d,%d) arms %d, ",
               node->myTag.domainID, node->myTag.index, 
               node->numNbrs);
#if 1
        for (i = 0; i < node->numNbrs; i++) {
            printf("(%d,%d) ", node->nbrTag[i].domainID,
                   node->nbrTag[i].index);
        }
#endif
        printf("\n");
#endif

#if 1
/*
 *      Print the nodal position
 */
        printf("  node(%d,%d) position = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->x, node->y, node->z);
#endif

#if 1
/*
 *      Print the nodal velocity and total node force
 */
        printf("  node(%d,%d) v = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->vX, node->vY, node->vZ);
        printf("  node(%d,%d) f = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->fX, node->fY, node->fZ);
#endif

#if 1
/*
 *      Print the arm specific forces
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) f = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->armfx[i], node->armfy[i], node->armfz[i]);
        }
#endif

#if 1
/*
 *      Print the burger's vector for each arm of the node
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) b = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif

#if 1
/*
 *      Print the glide plane normal for each arm of the node
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) n = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID,node->myTag.index, i,       
                   node->nbrTag[i].domainID,node->nbrTag[i].index,
                   node->nx[i],node->ny[i],node->nz[i]);
        }
#endif

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     RecalcSegGlidePlane
 *      Description:  Calculate and reset (if necessary) the glide plane
 *                    normal for a segment.  This function will only update
 *                    the information locally.  It is assumed the information
 *                    will be passed to remote domains elsewhere as needed.
 *
 *      Arguments:
 *          node1, node2   pointers to nodal endpoints of the segment
 *          ignoreIfScrew  Toggle.  If set to 1, the function will
 *                         leave the segment glide plane as is if the
 *                         segment is screw.  Otherwise it will pick
 *                         an appropriate glide plane for the bugers vector
 *
 *------------------------------------------------------------------------*/
void RecalcSegGlidePlane(Home_t *home, Node_t *node1, Node_t *node2,
                         int ignoreIfScrew)
{
        int node1SegID, node2SegID;
        real8 burg[3], lineDir[3], newPlane[3];

/*
 *      Just a couple quick sanity checks.
 */
        if ((node1 == (Node_t *)NULL) ||
            (node2 == (Node_t *)NULL) ||
            (node1 == node2)) {
            return;
        }

/*
 *      It is possible that the two nodes are not really connected.
 *      This can happen in cases such as MeshCoarsen() where a node is
 *      removed leaving two nodes doubly linked, and when the double links
 *      get reconciled they annihilate each other, leaving the two nodes
 *      unconnected.  In situations like that, there's nothing for this
 *      function to do...
 */
        if (!Connected(node1, node2, &node1SegID)) {
            return;
        }

        node2SegID = GetArmID(home, node2, node1);

        burg[X] = node1->burgX[node1SegID];
        burg[Y] = node1->burgY[node1SegID];
        burg[Z] = node1->burgZ[node1SegID];

        lineDir[X] = node2->x - node1->x;
        lineDir[Y] = node2->y - node1->y;
        lineDir[Z] = node2->z - node1->z;

        ZImage(home->param, &lineDir[X], &lineDir[Y], &lineDir[Z]);
        NormalizeVec(lineDir);

        FindPreciseGlidePlane(home, burg, lineDir, newPlane);

        if (DotProduct(newPlane, newPlane) < 1.0e-03) {
            if (ignoreIfScrew) return;
            PickScrewGlidePlane(home, burg, newPlane);
        }

        Normalize(&newPlane[X], &newPlane[Y], &newPlane[Z]);

        node1->nx[node1SegID] = newPlane[X];
        node1->ny[node1SegID] = newPlane[Y];
        node1->nz[node1SegID] = newPlane[Z];

        node2->nx[node2SegID] = newPlane[X];
        node2->ny[node2SegID] = newPlane[Y];
        node2->nz[node2SegID] = newPlane[Z];

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:      AddOp
 *      Description:   Add a topological operation to the list to
 *                     be sent to remote domains for processing.
 *
 *      Arguments:
 *          type        Type of operation to be performed
 *          dom1, idx1  Tag information for the 1st node involved
 *                      in the operation.
 *          dom2, idx2  Tag information for the 2nd node involved
 *                      in the operation.  (If no second node is
 *                      involved in this operation, these values
 *                      should be set to -1)
 *          dom3, idx3  Tag information for the 3rd node involved
 *                      in the operation.  (If no third node is
 *                      involved in this operation, these values
 *                      should be set to -1)
 *          bx, by, bz  Components of the burgers vector for the
 *                      operation.  (If not applicable to the
 *                      operation, these should be zeroes)
 *          x, y, z     Components of node's position (If not
 *                      applicable to the operation, should be zero'ed)
 *          nx, ny, nz  Components of glide plane normal. (If not
 *                      applicable to the operation, should be zero'ed)
 *          
 *------------------------------------------------------------------------*/
void AddOp (Home_t *home,
            OpType_t type, int dom1, int idx1,
            int dom2, int idx2, int dom3, int idx3,
            real8 bx, real8 by, real8 bz,
            real8 x, real8 y, real8 z,
            real8 nx, real8 ny, real8 nz) 
{
        Operate_t *op;

/*
 *      Make sure the buffer allocated for the operation list is
 *      large enough to contain another operation.  If not, increase
 *      the size of the buffer.
 */
        if (home->OpCount >= home->OpListLen) {
            ExtendOpList (home);
        }
        
        op = &home->opList[home->OpCount];

        op->type = type;
        op->dom1 = dom1;
        op->idx1 = idx1;
        op->dom2 = dom2;
        op->idx2 = idx2;
        op->dom3 = dom3;
        op->idx3 = idx3;
        op->bx = bx;
        op->by = by;
        op->bz = bz;
        op->x = x;
        op->y = y;
        op->z = z;
        op->nx = nx;
        op->ny = ny;
        op->nz = nz;
        home->OpCount++;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:      ClearOpList
 *      Description:   Zero's both the list of operations which is
 *                     communicated to remote domains during various
 *                     stages of topological changes and the count
 *                     of operations on the list .
 *
 *------------------------------------------------------------------------*/
void ClearOpList (Home_t *home)
{
        memset(home->opList, 0, sizeof(Operate_t)*home->OpListLen);
        home->OpCount = 0;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:      ExtendOpList
 *      Description:   Increases the buffer size dedicated to the
 *                     list of operations sent to remote domains
 *                     after topological changes.
 *
 *------------------------------------------------------------------------*/
void ExtendOpList (Home_t *home)
{
        int newSize;

        newSize = (home->OpListLen+OpBlock_Count) * sizeof(Operate_t);
        home->opList = (Operate_t *) realloc(home->opList, newSize);
        home->OpListLen += OpBlock_Count;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:      FreeOpList
 *      Description:   Deallocates the memory dedicated to the operation
 *                     list sent to remote domains after topological
 *                     changes.
 *
 *                     NOTE:  The current code no longer frees the
 *                     buffer associated with the list, instead allowing
 *                     it to persist through the entire execution for
 *                     efficiency reasons.  Hence, this function is
 *                     currently obsolete.
 *
 *------------------------------------------------------------------------*/
void FreeOpList (Home_t *home)
{
        free(home->opList);
        home->OpCount = 0;
        home->OpListLen = 0;
}


/*-------------------------------------------------------------------------
 *
 *      Function:      InitOpList
 *      Description:   Allocates an initial block of memory for the
 *                     list of operations this domain will send to
 *                     remote domains for processing during the
 *                     various stages of topological changes.  
 *
 *                     This function should only need to be called 
 *                     one time during initialization of the application.
 *                     Adding operations to the list as well as altering
 *                     the list size will be handled dynamically as
 *                     necessary
 *
 *------------------------------------------------------------------------*/
void InitOpList (Home_t *home)
{
        home->opList = (Operate_t *) malloc(sizeof(Operate_t)*OpBlock_Count);
        home->OpCount = 0;
        home->OpListLen = OpBlock_Count;
        home->rcvOpCount = 0;
        home->rcvOpList = 0;

        return;
}


/* from OpRec.c */

/*-------------------------------------------------------------------------
 *
 *      Function:     RequestNodeTag
 *
 *------------------------------------------------------------------------*/
int RequestNodeTag(Home_t *home, Tag_t *tag)
{
        if (tag->index >= home->newNodeKeyMax) {
/*
 *          The node index is higher than the current size
 *          of the node list. Reallocate the node list with
 *          the new size to allow the node index to be requested.
 */
            int NUM_INC = ceil((tag->index-home->newNodeKeyMax-1)/NEW_NODEKEY_INC);

            home->newNodeKeyMax += NUM_INC*NEW_NODEKEY_INC;
            home->nodeKeys = (Node_t **) realloc(home->nodeKeys,
            home->newNodeKeyMax * sizeof(Node_t *));

            home->newNodeKeyPtr = tag->index+1;

            return tag->index;

        } else {
/*
 *          Check if the requested node index is already used. If not
 *          make sure to update the newNodeKeyPtr variable if needed.
 */
            if (home->nodeKeys[tag->index] == (Node_t *)NULL) {

                if (tag->index >= home->newNodeKeyPtr) {
                    home->newNodeKeyPtr = tag->index+1;
                }
                return tag->index;

            } else {
                return -1;
            }
        }
}

/*-------------------------------------------------------------------------
 *
 *      Function:     RequestNewNativeNodeTag
 *
 *------------------------------------------------------------------------*/
Node_t *RequestNewNativeNodeTag(Home_t *home, Tag_t *tag)
{
        int     newIdx;
        Node_t *newNode;

        newNode = PopFreeNodeQ(home);
        newIdx = RequestNodeTag(home, tag);

        if (newIdx < 0) {
            Fatal("Cannot request node index (%d,%d)",
            tag->domainID, tag->index);
        }

        home->nodeKeys[newIdx] = newNode;

        newNode->myTag.domainID = home->myDomain;
        newNode->myTag.index    = newIdx;
        newNode->cellIdx        = -1;
        newNode->cell2Idx       = -1;
        newNode->cell2QentIdx   = -1;

/*
 *      Explicitly zero out velocity so we don't end up
 *      with garbage when we are calculating the velocity
 *      delta between timesteps.
 */
        newNode->vX = 0.0;
        newNode->vY = 0.0;
        newNode->vZ = 0.0;

        newNode->oldvX = 0.0;
        newNode->oldvY = 0.0;
        newNode->oldvZ = 0.0;

        newNode->flags = 0;

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
        newNode->multiNodeLife = 0;
#endif

#ifdef _FEM
        newNode->fem_Surface[0] = 0;
        newNode->fem_Surface[1] = 0;

        newNode->fem_Surface_Norm[0] = 0.0;
        newNode->fem_Surface_Norm[1] = 0.0;
        newNode->fem_Surface_Norm[2] = 0.0;
#endif

        return(newNode);
}

/* from InitSendDomains.c */

/*---------------------------------------------------------------------------
 *
 *      Function:     ExtendNodeKeys
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
static void ExtendNodeKeys(Home_t *home, int newLength)
{
        int i, oldLength;

        oldLength = home->newNodeKeyMax;
        home->newNodeKeyMax = newLength;

        home->nodeKeys = (Node_t **)realloc(home->nodeKeys,
                                            newLength * sizeof(Node_t *));

        for (i = oldLength; i < newLength; i++) {
            home->nodeKeys[i] = (Node_t *)NULL;
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     AddTagMapping
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag)
{
        int      newSize;
        TagMap_t *mapping;

        if (home->tagMapEnts >= home->tagMapSize) {
              home->tagMapSize += NEW_NODEKEY_INC;
              newSize = home->tagMapSize * sizeof(TagMap_t);
              home->tagMap = (TagMap_t *)realloc(home->tagMap, newSize);
        }

        mapping = &home->tagMap[home->tagMapEnts];

        mapping->oldTag.domainID = oldTag->domainID;
        mapping->oldTag.index    = oldTag->index;

        mapping->newTag.domainID = newTag->domainID;
        mapping->newTag.index    = newTag->index;

        home->tagMapEnts++;

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     GetNewTag
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
void GetNewTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
                      int *nextAvailableTag)
{
        int     nextTag, thisDomain;

        nextTag    = *nextAvailableTag;
        thisDomain = home->myDomain;

/*
 *      If the old tag belonged to a different domain, we must
 *      give the node a new tag.
 */
        if (oldTag->domainID != thisDomain) {

            for ( ; ; nextTag++) {

/*
 *              Extend the nodekeys array if necessary.
 */
                if (nextTag >= home->newNodeKeyMax) {
                    ExtendNodeKeys(home, home->newNodeKeyMax + NEW_NODEKEY_INC);
                }

                if (home->nodeKeys[nextTag] == (Node_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }
        } else {
/*
 *          The old tag belonged to this domain...  check if that tag
 *          is available for this run.
 *
 *          Extend the nodekeys array if necessary.
 */
            if (oldTag->index >= home->newNodeKeyMax) {
                ExtendNodeKeys(home, oldTag->index + NEW_NODEKEY_INC);
            }

/*
 *          If the old tag is still available, use it and return
 *          to the caller.
 */
            if (home->nodeKeys[oldTag->index] == (Node_t *)NULL) {
                newTag->domainID = oldTag->domainID;
                newTag->index    = oldTag->index;
                if (newTag->index >= home->newNodeKeyPtr) {
                    home->newNodeKeyPtr = newTag->index + 1;
                }
                return;
            }

/*
 *          The old tag is no longer available, so just use
 *          the next available tag.
 */
            for ( ; ; nextTag++) {

/*
 *              Extend node keys array if necessary
 */
                if (nextTag >= home->newNodeKeyMax) {
                    ExtendNodeKeys(home, home->newNodeKeyMax + NEW_NODEKEY_INC);
                }

                if (home->nodeKeys[nextTag] == (Node_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }
        }

        if (newTag->index >= home->newNodeKeyPtr) {
            home->newNodeKeyPtr = newTag->index + 1;
        }

        *nextAvailableTag = nextTag;

        return;
}

/* from Util.c */

void AddNodesFromArray(Home_t *home, real8 *buf)
{
        int        i, armID, numNbrs, nodesInBuf, bufIndex;
        Tag_t      tag;
        Node_t     *node;

/*
 *      Pull the node count out of the buffer then loop through
 *      all the nodal data provided.
 */
        bufIndex   = 0;

        nodesInBuf = (int)buf[bufIndex++];
        printf("AddNodesFromArray: nodesInBuf = %d\n", nodesInBuf);

        printf("AddNodesFromArray: ExtendNodeKeys = %d, %d\n", home->newNodeKeyMax, NEW_NODEKEY_INC);
        ExtendNodeKeys(home, home->newNodeKeyMax + NEW_NODEKEY_INC);

        for (i = 0; i < nodesInBuf; i++) {

            tag.domainID = (int)buf[bufIndex++];
            tag.index    = (int)buf[bufIndex++];

            node = RequestNewNativeNodeTag(home, &tag);

            node->x = buf[bufIndex++];
            node->y = buf[bufIndex++];
            node->z = buf[bufIndex++];

            numNbrs = (int)buf[bufIndex++];
            node->constraint = (int)buf[bufIndex++];

            ReallocNodeArms(node, numNbrs);
            printf("AddNodesFromArray: node(%d,%d) x=%g, y=%g, z=%g numNbrs=%d\n", 
                    node->myTag.domainID, node->myTag.index, node->x, node->y, node->z, node->numNbrs);

            for (armID = 0; armID < numNbrs; armID++) {
                node->nbrTag[armID].domainID = (int)buf[bufIndex++];
                node->nbrTag[armID].index    = (int)buf[bufIndex++];
                node->burgX[armID] = buf[bufIndex++];
                node->burgY[armID] = buf[bufIndex++];
                node->burgZ[armID] = buf[bufIndex++];
                node->nx[armID] = buf[bufIndex++];
                node->ny[armID] = buf[bufIndex++];
                node->nz[armID] = buf[bufIndex++];
            }

            PushNativeNodeQ(home, node);
        }  /* for (i = 0; i < nodesInBuf; ... ) */
        SortNativeNodes(home);
        printf("AddNodesFromArray: done\n");
        return;
}


/* From ReadRestart.c */
void FreeAllNodes(Home_t *home)
{
    Node_t       *node;
    int i;
    for (i = 0; i < home->newNodeKeyPtr; i++) {
        if (home->nodeKeys[i] != (Node_t *)NULL) {
            node =  home->nodeKeys[i];
            FreeNodeArms(node);
            /* node memory managed by NodeQ, not sure if this is going to cause memory leak */
            home->nodeKeys[i] = NULL;
        }
    }
    home->newNodeKeyPtr = 0;
}


/* From Initialize.c */
void SetBoxSize(Param_t *param){
/*
 *      Some of the parameters used in creating the nodal data file
 *      used for this run may not match the values to be used for
 *      this run.  We've completed processing the data file at this
 *      point, so update the data file parameters to match the values
 *      desired for this particular run.
 */
        param->dataDecompGeometry[X] = param->nXdoms;
        param->dataDecompGeometry[Y] = param->nYdoms;
        param->dataDecompGeometry[Z] = param->nZdoms;

        param->dataDecompType = param->decompType;
        param->dataFileVersion = NODEDATA_FILE_VERSION;
        param->numFileSegments = param->numIOGroups;

/*
 *      Calculate the length of the problem space in each
 *      of the dimensions
 */
        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

        param->invLx = 1.0 / param->Lx;
        param->invLy = 1.0 / param->Ly;
        param->invLz = 1.0 / param->Lz;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SetRemainingDefaults
 *      Description:  The default values of certain global parameters
 *                    are special in that they depend on values of
 *                    other global parameters.  If the user did not
 *                    specify values for these special parameters,
 *                    this function will calculate the necessary
 *                    defaults (as well as do some additional sanity
 *                    checks on some of the values).
 *
 *-------------------------------------------------------------------------*/
void SetRemainingDefaults(Home_t *home)
{
        real8   tmp, eps;
        real8   xCellSize, yCellSize, zCellSize, minCellSize;
        Param_t *param;

        param = home->param;

        param->delSegLength = 0.0;

        xCellSize = param->Lx / param->nXcells;
        yCellSize = param->Ly / param->nYcells;
        zCellSize = param->Lz / param->nZcells;

        minCellSize = MIN(xCellSize, yCellSize);
        minCellSize = MIN(minCellSize, zCellSize);

        eps = 1.0e-02;

/*
 *      The core radius and maximum segment length are required
 *      inputs.  If the user did not provide both values, abort
 *      now.
 */
        if (home->myDomain == 0) {
            if (param->rc < 0.0) {
                Fatal("The <rc> parameter is required but was not \n"
                      "    provided in the control file");
            }

            if (param->maxSeg < 0.0) {
                Fatal("The <maxSeg> parameter is required but was not \n"
                      "    provided in the control file");
            }
        }

/*
 *      If not provided, set position error tolerance based on <rc>
 */
        if (param->rTol <= 0.0) {
            param->rTol = 0.25 * param->rc;
        }

/*
 *      The deltaTT is set in the timestep integrator, but some
 *      mobility functions now use the deltaTT value, so it must
 *      be initialized before ParadisStep() is called since there
 *      is an initial mobility calculation done *before* timestep
 *      integration the first time into the function.
 */
        param->deltaTT = MIN(param->maxDT, param->nextDT);

        if (param->deltaTT <= 0.0) {
            param->deltaTT = param->maxDT;
        }

/*
 *      Set annihilation distance based on <rc>
 */
#ifdef _RETROCOLLISIONS
		if (param->rann<0) {
        	param->rann = 2.0 * param->rTol;
		}
#else
		param->rann = 2.0 * param->rTol;
#endif

/*
 *      Minimum area criteria for remesh is dependent on maximum
 *      and minumum segment lengths and position error tolerance.
 */
        param->remeshAreaMin = 2.0 * param->rTol * param->maxSeg;

        if (param->minSeg > 0.0) {
            param->remeshAreaMin = MIN(param->remeshAreaMin,
                                       (param->minSeg * param->minSeg *
                                        sqrt(3.0) / 4));
        }

/*
 *      Maximum area criteria for remesh is dependent on minimum area,
 *      and maximum segment length.
 */
        param->remeshAreaMax = 0.5 * ((4.0 * param->remeshAreaMin) +
                                      (0.25 * sqrt(3)) *
                                      (param->maxSeg*param->maxSeg));

/*
 *      If the user did not provide a minSeg length, calculate one
 *      based on the remesh minimum area criteria.
 */
        if (param->minSeg <= 0.0) {
            param->minSeg = sqrt(param->remeshAreaMin * (4.0 / sqrt(3)));
        }


/*
 *      If the user did not provide an Ecore value, set the default
 *      based on the shear modulus and rc values
 */
        if (param->Ecore < 0.0) {
            param->Ecore = (param->shearModulus / (4*M_PI)) *
                           log(param->rc/0.1);
        }

/*
 *      Now do some additional sanity checks.
 */
        if (home->myDomain == 0) {

/*
 *          First check for some fatal errors...
 */
            if (param->maxSeg <= param->rTol * (32.0 / sqrt(3.0))) {
                Fatal("Maximum segment length must be > rTol * 32 / sqrt(3)\n"
                      "    Current maxSeg = %lf, rTol = %lf",
                      param->maxSeg, param->rTol);
            }

            if (param->minSeg > (0.5 * param->maxSeg)) {
                Fatal("Minimum segment length must be < (0.5 * maxSeg)\n"
                      "    Current minSeg = %lf, maxSeg = %lf",
                      param->minSeg, param->maxSeg);
            }

            if (param->maxSeg <= param->minSeg) {
                Fatal("Max segment length (%e) must be greater than the\n"
                      "    minimum segment length (%e)", param->maxSeg,
                      param->minSeg);
            }

            if (param->maxSeg > (minCellSize * 0.9)) {
                Fatal("The maxSeg length must be less than the "
                      "minimum cell size * 0.9.  Current values:\n"
                      "    maxSeg    = %.1f\n    cellWidth = %.1f",
                      param->maxSeg, minCellSize);
            }

            if (param->remeshAreaMin > (0.25 * param->remeshAreaMax)) {
                Fatal("remeshAreaMin must be less than 0.25*remeshAreaMax\n"
                      "    Current remeshAreaMin = %lf, remeshAreaMax = %lf",
                      param->remeshAreaMin, param->remeshAreaMax);
            }

/*
 *          Now check for conditions that although not fatal, may result
 *          in undesired behaviour, and warn the user.
 */
            if (param->rc < 0.1) {
                fprintf(stderr, "WARNING: specified rc value (%e) will "
                                "yield a \nnegative core energy\n", param->rc);
            }

            tmp = (param->maxSeg * param->maxSeg * param->maxSeg);

            if (param->remeshAreaMax > (0.25 * sqrt(3) * tmp)) {
                fprintf(stderr, "WARNING: Area criteria will be unused "
                                "in remesh operations!\n");
                fprintf(stderr, "         rmeshAreaMax = %lf, maxSeg = %lf\n",
                                param->remeshAreaMax, param->maxSeg);
            }

            if (param->rann > (0.5 * param->rc + eps)) {
                fprintf(stderr, "WARNING: Separation distance is larger "
                                "than the core radius!\n");
                fprintf(stderr, "         rann = %lf, rc = %lf\n",
                                param->rann, param->rc);
            }

            if (param->rann > (2.0 * param->rTol)) {
                fprintf(stderr, "WARNING: Collision distance is outside the "
                                "position error tolerance!\n");
                fprintf(stderr, "         rann = %lf, rTol = %lf\n",
                                param->rann, param->rTol);
            }

#if 0
            tmp = param->remeshAreaMin - (2.0 * param->rTol * param->maxSeg);

            if (fabs(tmp) > eps) {
                fprintf(stderr, "WARNING: remesh minimum area != "
                                "2.0 * rTol * maxSeg\n");
                fprintf(stderr, "         remeshAreaMin = %lf, rTol = %lf"
                                "maxSeg = %lf\n", param->remeshAreaMin,
                                param->rTol, param->maxSeg);
            }
#endif

/*
 *          If free suraces are used but the specified surfaces are
 *          not within the primary bounding box, it's a problem.
 */
            if (((param->xBoundType == Free) &&
                ((param->xBoundMin < param->minCoordinates[X]) ||
                 (param->xBoundMax > param->maxCoordinates[X])))||
                ((param->yBoundType == Free) &&
                ((param->yBoundMin < param->minCoordinates[Y]) ||
                 (param->yBoundMax > param->maxCoordinates[Y])))||
                ((param->zBoundType == Free) &&
                ((param->zBoundMin < param->minCoordinates[Z]) ||
                 (param->zBoundMax > param->maxCoordinates[Z]))) ) {
                Fatal("Free surfaces are not within main bounding box!\n"
                      "    Surface min coordinates (%lf %lf %lf)\n"
                      "    Surface max coordinates (%lf %lf %lf)\n",
                      param->xBoundMin, param->yBoundMin, param->zBoundMin,
                      param->xBoundMax, param->yBoundMax, param->zBoundMax);
            }

#if !defined _FEM & !defined _FEMIMGSTRESS
/*
 *          If free surfaces are enabled but the finite element code
 *          is not linked in, results will not be accurate, so print
 *          a warning.
 */
            if ((param->xBoundType == Free) ||
                (param->yBoundType == Free) ||
                (param->zBoundType == Free)) {
                printf("***\n*** WARNING!  Use of free surfaces in ParaDiS "
                       "without the\n*** FEM/ParaDiS coupling is not "
                       "fully supported!\n***\n");
            }
#endif

        }  /* if domain == 0 */


/*
 *      If there are a mix of free surfaces and periodic boundaries,
 *      the boundary min/max values must default to the simulation
 *      boundaries in the dimensions without free surfaces.
 */
        if ((param->xBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->zBoundType == Free)) {

            if (param->xBoundType == Periodic) {
                param->xBoundMin = param->minSideX;
                param->xBoundMax = param->maxSideX;
            }
            if (param->yBoundType == Periodic) {
                param->yBoundMin = param->minSideY;
                param->yBoundMax = param->maxSideY;
            }
            if (param->zBoundType == Periodic) {
                param->zBoundMin = param->minSideZ;
                param->zBoundMax = param->maxSideZ;
            }
        }

/*
 *      Based on the mobility law selected in the control
 *      parameter file, set:
 *        1) the material type (BCC, FCC, etc.)
 *        2) the specific mobility type
 *        3) a pointer to the proper mobility function
 *        4) number of burgers vector groups used in
 *           tracking dislocation density per burgers vector
 *
 *      *************************************************
 *      ***                                           ***
 *      ***                  IMPORTANT!               ***
 *      ***   If you change any numBurgGroups value   ***
 *      ***   specified below, you must change the    ***
 *      ***   DENSITY_FILE_VERSION number defined     ***
 *      ***   in WriteProp.c!                         ***
 *      ***                                           ***
 *      *************************************************
 */
#if 0
        if (strcmp(param->mobilityLaw, "BCC_0") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_0;
            param->mobilityFunc = Mobility_BCC_0;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "BCC_0b") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_0B;
            param->mobilityFunc = Mobility_BCC_0b;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "BCC_glide") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_GLIDE;
            param->mobilityFunc = Mobility_BCC_glide;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "BCC_glide_0") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_GLIDE_0;
            param->mobilityFunc = Mobility_BCC_glide_0;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "FCC_0") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_0;
            param->mobilityFunc = Mobility_FCC_0;
            param->numBurgGroups = 7;
        } else if (strcmp(param->mobilityLaw, "FCC_0b") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_0B;
            param->mobilityFunc = Mobility_FCC_0b;
            param->numBurgGroups = 7;
        } else if (strcmp(param->mobilityLaw, "FCC_climb") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_CLIMB;
            param->mobilityFunc = Mobility_FCC_climb;
            param->numBurgGroups = 7;
            param->allowFuzzyGlidePlanes = 1;
        } else if (strcmp(param->mobilityLaw, "RELAX") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_RELAX;
            param->mobilityFunc = Mobility_Relax;
            param->numBurgGroups = 7;
        } else {
            Fatal("Unknown mobility function %s", param->mobilityLaw);
        }
#endif

#ifdef _GPU_SUBCYCLE
		if (param->mobilityType == MOB_FCC_0) {
			param->mobilityMatrixGPU = Mobility_FCC_0_matrix_GPU;
		} else if (param->mobilityType == MOB_BCC_0B) {
			// Don't do anything here.
		} else {
			Fatal("GPU mobility function is not available for %s", param->mobilityLaw);
		}
#endif

        param->partialDisloDensity =
                (real8 *)malloc(param->numBurgGroups * sizeof(real8));

/*
 *      Some types of mobility require the enforceGlidePlanes flag
 *      to be set.  Handle that here.
 */
        switch (param->mobilityType) {
            case MOB_BCC_GLIDE:
            case MOB_BCC_GLIDE_0:
            case MOB_FCC_0:
            case MOB_FCC_0B:
            case MOB_FCC_CLIMB:
                if (param->enforceGlidePlanes == 0) {
                    param->enforceGlidePlanes = 1;
                    if (home->myDomain == 0) {
                        printf("The specified mobility (%s) requires the "
                               "enforceGlidePlanes\ncontrol parameter "
                               "toggle to be set.  Enabling toggle now.\n",
                               param->mobilityLaw);
                    }
                }
                break;
        }

/*
 *      If type 1 domain decompositionis enabled, the DLBfreq
 *      value MUST be a multiple of 3.
 */
        if ((param->DLBfreq > 0) && (param->decompType == 1)) {
            param->DLBfreq = (param->DLBfreq + 2) / 3 * 3;
        }

/*
 *      If the cross slip flag has not been explicitly defined, give
 *      it a default setting based on the mobility being used.
 */
        if (param->enableCrossSlip < 0) {
            switch(param->mobilityType) {
                case MOB_BCC_GLIDE:
                case MOB_BCC_GLIDE_0:
                case MOB_FCC_0:
                case MOB_FCC_0B:
                case MOB_FCC_CLIMB:
                    param->enableCrossSlip = 1;
                    break;
                default:
                    param->enableCrossSlip = 0;
                    break;
            }
        }

#ifdef _MOBILITY_FIELD
	if (param->mobilityField || param->frictionField) {
		ReadMobilityField(home);
	}
#endif

	if (param->FricStress > 0.0 && param->mobilityType != MOB_FCC_0) {
		Fatal("Option FricStress is only available for MOB_FCC_0!");
	}
/*
 *      Several portions of the code need to calculate the total
 *      volume of the simulation and a volume factor used for
 *      converting dislocation length to density, so set those
 *      factors now.
 *
 *      If the FEM code is linked in, this is going to depend
 *      on the actual shape used within the primary image.  Otherwise
 *      we use the volume based on the free surfaces (a rectagular
 *      prism) or the full dimensions of the primary image.
 */
#ifdef _FEM
        switch (param->mesh_type) {
            case 1:
/*
 *              Shape is a rectangular prism
 */
                param->simVol = (param->xBoundMax-param->xBoundMin) *
                                (param->yBoundMax-param->yBoundMin) *
                                (param->zBoundMax-param->zBoundMin);
                break;
            case 2:
/*
 *              Shape is a cylinder
 */
                param->simVol = M_PI * param->fem_radius *
                                param->fem_radius * param->fem_height;
                break;
            default:
/*
 *              Unknown shape, so treat it as a rectangular prism
 */
                param->simVol = (param->xBoundMax-param->xBoundMin) *
                                (param->yBoundMax-param->yBoundMin) *
                                (param->zBoundMax-param->zBoundMin);
                break;
        }
#else
        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {

            param->simVol = (param->xBoundMax-param->xBoundMin) *
                            (param->yBoundMax-param->yBoundMin) *
                            (param->zBoundMax-param->zBoundMin);
        } else {
            param->simVol = param->Lx * param->Ly * param->Lz;
        }
#endif

        param->burgVolFactor = 1.0 / (param->burgMag * param->burgMag *
                                      param->simVol);

#ifdef FULL_N2_FORCES
/*
 *      To do full n^2 force calculations without remote forces, we need
 *      to explicitly set some flags.
 */
        param->elasticinteraction = 1;
        param->fmEnabled = 0;
        param->numDLBCycles = 0;
        param->DLBfreq = 0;
#endif

/*
 *      Deactivate FMM if line tension model
 */
		if (!param->elasticinteraction) {
			param->fmEnabled = 0;
		}

        return;
}
