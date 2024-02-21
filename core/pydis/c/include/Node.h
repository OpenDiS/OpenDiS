/*--------------------------------------------------------------------------
 *
 *	Node.h	Define the struct that holds all relevant data for a single
 *		node, either native or ghost
 *
 *		Notice: make depend if new files are added, otherwise always
 *			make clean after changing this file  Wei Cai 04/09/2002
 *
 *------------------------------------------------------------------------*/

#ifndef _Node_h
#define _Node_h

#include "Typedefs.h"
#include "Tag.h"
#include "ParadisThread.h"

/*
 *      Define the various node 'constraints' available.  Note: these
 *      constraints are mutually exclusive.
 */
#define UNCONSTRAINED  0
#define SURFACE_NODE   1
#define THINFILM_SURFACE_NODE   6
#define CYLINDER_SURFACE_NODE   6
#define PINNED_NODE    7

/*
 *      Define the bit flags that can be set for the node.  Note: these
 *      flags may be OR'ed together.
 */
#define NODE_RESET_FORCES    0x01
#define NODE_OUTSIDE_SURF    0x02
#define NO_COLLISIONS        0x04
#define NO_MESH_COARSEN      0x08
#define NODE_CHK_DBL_LINK    0x10

/*
 *      Used as bit flags to indicate the type of nodal data items
 *      preserved by calls to PreserveNodalData().
 */
#define NODE_POSITION 0x01
#define NODE_CURR_VEL 0x02
#define NODE_OLD_VEL  0x04

struct _node {
	int	flags;
#ifdef _SUBCYCLING
	int     subgroup;           /*group ID used for subcycling, 0 or 1*/
	int     G0_to_G4;
	int     newNode;
	int     CommSend[5];        /*flag used for parallel communication with subcycling*/
	real8   fxLong, fyLong, fzLong;

	int     fricPinned;
	real8   fricForce[3];
#endif
#ifdef _GPU_SUBCYCLE
	int     subindex;
	int     numInt;
#endif

	real8	x, y, z;		/* nodal position */
	real8	fX, fY, fZ;		/* nodal force: units=Pa*b^2) */
	real8	vX, vY, vZ;		/* nodal velocity: units=burgMag/sec */


	real8	oldx, oldy, oldz;	/* for strain increment, Moono.Rhee */
	real8	oldvX, oldvY, oldvZ;	/* nodal velocity at previous step */
        real8   currvX, currvY, currvZ; /* nodal velocity at beginning of the */
                                        /* current step.  Only used during    */
                                        /* timestep integration while vX, vY  */
                                        /* and vZ are being determined.       */
#if defined _SUBCYCLING | defined _RETROCOLLISIONS
	real8	olderx, oldery, olderz; /* for retroactive collision detection and strain increment when using subcycling*/
    real8	oldervX, oldervY, oldervZ;
	real8   olderfX, olderfY, olderfZ;
	real8   RKFx[6], RKFy[6], RKFz[6];
#endif

	Tag_t	myTag;

/*
 *	nbrTag and burgID are dynamically allocated to size numNbrs, when the
 *	node is allocated.
 */
	int	numNbrs;
	Tag_t	*nbrTag;
#ifdef _GPU_SUBCYCLE
	int     *armid;
	int     *armInt;
#endif

/*
 *	Arm information
 */
	real8	*armfx, *armfy, *armfz;	/* arm specific force contribution */
	real8	*burgX, *burgY, *burgZ;	/* burgers vector */
	real8	*nx, *ny, *nz;		/* glide plane */

	real8	*sigbLoc;		/* sig.b on arms (numNbr*3) */
	real8	*sigbRem;		/* sig.b on arms (numNbr*3) */

	int	*armCoordIndex;		/* Array of indices (1 per arm) into */
					/* the mirror domain's arrays of     */
					/* coordinates (armX, armY, armZ)    */
					/* of nodes' neighbors.  This array  */
					/* is only present and useful while  */
					/* task zero is downloading the data */
					/* from the remote domain for        */
					/* generating output                 */

	int	constraint;     /* constraint =  1 : any surface node   */
				/* constraint =  7 : Frank-Read end	*/
				/*                   points, fixed	*/
				/* constraint =  9 : Jog		*/
				/* constraint = 10 : Junction		*/
				/* constraint = 20 : cross slip(default)*/
				/*                   2-nodes non-deletable*/

	int	cellIdx;	/* cell node is currently sorted into */
	int	cell2Idx;	/* cell2 node is currently sorted into */
	int	cell2QentIdx;	/* Index of this node's entry in the */
				/* home->cell2QentArray.             */

	int	native;		/* 1 = native node, 0 = ghost node */

	Node_t	*next;		/* pointer to the next node in the queue */
				/* (ghost or free)			 */

	Node_t	*nextInCell;	/* used to queue node onto the current */
				/* containing cell		       */

	int	sgnv;		/* +1: if contribute to strain rate,	*/
				/* -1: if moving in opposite direction	*/

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
        int     multiNodeLife;
#endif

#ifdef _FEM
        int     fem_Surface[2];
        real8   fem_Surface_Norm[3];
#endif

#ifdef _OPENMP
        omp_lock_t nodeLock;
#endif
};

struct _nodeblock {
	NodeBlock_t	*next;
	Node_t		*nodes;

};

#endif
