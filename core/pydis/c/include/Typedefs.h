/***************************************************************************
 *
 *  Typedefs.h - This include file defines all the typedefs used for
 *               data structures in the code. It allows the typedefs
 *               to be used in various header files before the header
 *               file where the typedef'd structure is defined. It should
 *               be included at the top of any header file that references
 *               structures defined in other header files
 *
 *************************************************************************/

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include "Constants.h"

typedef double real8;

typedef struct _binnode BinNode_t;
typedef struct _binseg BinSeg_t;
typedef struct _cell Cell_t;
typedef struct _dcompdomain DcompDomain_t;
typedef struct _home Home_t;
typedef struct _indata InData_t;
typedef struct _innode InNode_t;
typedef struct _mirrordomain MirrorDomain_t;
typedef struct _node Node_t;
typedef struct _nodeblock NodeBlock_t;
typedef struct _operate Operate_t;

#ifdef _OP_REC
typedef struct _operaterec OperateRec_t;
#endif

typedef struct _param Param_t;
typedef struct _remotedomain RemoteDomain_t;
typedef struct _sortnode SortNode_t;
typedef struct _tag Tag_t;
typedef struct _timer Timer_t;
typedef struct _unmappedarm_t UnMappedArm_t;
typedef struct _segmentpair SegmentPair_t;

typedef enum {
	Periodic=0,
	Free=1,
	Reflecting=2
} BoundType_t;

/*
 *      Define the types of base operations that take place during
 *      the various topological changes.
 */
typedef enum {
	CHANGE_CONNECTION,
	INSERT_ARM,
	REMOVE_NODE,
	CHANGE_ARM_BURG,
	SPLIT_NODE,
	RESET_COORD,
	RESET_SEG_FORCES,
	RESET_SEG_FORCES2,
	MARK_FORCES_OBSOLETE,
        RESET_GLIDE_PLANE,
#ifdef _OP_REC
	START_STEP,
	TIME_INTEGRATION,
	START_SPLIT_MULTI,
	START_CROSS_SLIP,
	START_COLLISION,
	START_REMESH
#endif
} OpType_t;


/*
 *      Define a structure used when creating cell2 queues for
 *      collision handling.  An array of these structures will
 *      be allocated and initialized in SortNodesForCollision.c
 */
typedef struct {
	Node_t	*node;
	int	next;
} C2Qent_t;


/*
 *      Define a structure of data containing arrays and
 *      miscellaneous data needed for writing binary restart
 *      files.
 */
typedef struct {
        int   firstInGroup;
        int   lastInGroup;
        int   nodeCount;
        int   segCount;
        int   *nodeIndex;
        int   *nodeConstraint;
        int   *nodeNumSegs;
        int   *segTags;
        real8 *nodePos;
        real8 *burgersVec;
        real8 *glidePlane;
} BinFileData_t;

/*
 *      Define a couple structures used for managing the control
 *      and data file parameter lists
 */
typedef struct {
        char varName[MAX_STRING_LEN];
        int  valType;
        int  valCnt;
        int  flags;
        void *valList;
} VarData_t;


typedef struct {
        int       paramCnt;
        VarData_t *varList;
} ParamList_t;


#endif /* TYPEDEFS_H */
