/*****************************************************************************
 *
 *  QueueOps.h   Define the Queue operation prototypes
 *
 ****************************************************************************/

#ifndef _QueueOps_h
#define _QueueOps_h

Node_t *PopFreeNodeQ(Home_t *home);
void PushFreeNodeQ(Home_t *home, Node_t *node);
void PushNativeNodeQ(Home_t *home, Node_t *node);
void PushGhostNodeQ(Home_t *home, Node_t *node);
void RecycleGhostNodes(Home_t *home);
void RemoveNodeFromCellQ(Home_t *home, Node_t *node);
void RemoveNodeFromCell2Q(Home_t *home, Node_t *node);

#endif
