/***************************************************************************
 *
 *	Topology.h	Define the struct that holds all relevant data for a
 *			topological event and an event list, plus prototypes 
 *			for various functions used in identifying and
 *			topological events and effecting the appropriate
 *			topology changes.
 *
 **************************************************************************/

#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"


/*
 *      Define the set of status codes that may be returned from
 *      the MergeNode() function.
 */
#define MERGE_SUCCESS              0x01  /* Merge succeeded               */
#define MERGE_NODE_ORPHANED        0x02  /* Merge succeeded, but left an  */
                                         /* orphaned node in a remote     */
                                         /* domain                        */
#define MERGE_NO_REPOSITION        0x04  /* Merge succeeded, but the      */
                                         /* resulting node was not able   */
                                         /* to be repositioned.           */
#define MERGE_NOT_PERMITTED        0x08  /* Merge failed because the      */
                                         /* requested merge would violate */
                                         /* segment ownership rules       */
#define MERGE_DOUBLE_LINK          0x10  /* Merge failed because it would */
                                         /* have resulted in doubly       */
                                         /* connected nodes in another    */
                                         /* domain                        */
/*
 *      Define the set of status codes that may be returned from
 *      the SplitNode() function.
 */
#define SPLIT_FAILED   0
#define SPLIT_SUCCESS  1

/*
 *      Define any processing flags that can be provided to SplitNode() to
 *      affect its behavior.
 */
#define SPLIT_DUP_SURFACE_PROP 0x01


#define OPCLASS_SEPARATION  1
#define OPCLASS_COLLISION   2
#define OPCLASS_REMESH      3
/*
 *	Prototypes for functions involved in altering topology
 */
int    EvaluateMobility(Home_t *home, Node_t *nodeA);
int    InitTopologyExemptions(Home_t *home);
void   MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
           real8 *position, Node_t **mergedNode, int *status, int globalOp);
int    NodeTopologyExemptions(Home_t *home, Node_t *node);
int    RemoveDoubleLinks(Home_t *home, Node_t *node, int globalOp);
void   RemoveOrphanedNodes(Home_t *home);
void   SplitMultiNodes(Home_t *home);
int    SplitNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
           real8 *pos2, real8 *vel1, real8 *vel2, int armCount,
           int *armList, int globalOp, Node_t **splitNode1,
           Node_t **splitNode2, int flags);

/*
 *     Some functions for preserving and restoring node info
 *     while attempting multi-node or surface-node splits
 */
void   BackupNode(Home_t *home, Node_t *origNode, Node_t *bkupNode);
void   FreeNodeArrays(Node_t *node);
void   GetForcesFromBkup(Home_t *home, Node_t *node, Node_t *bkupNodeList);
void   RestoreNode(Home_t *home, Node_t *origNode, Node_t *bkupNode);

#endif /* _TOPOLOGY_H */
