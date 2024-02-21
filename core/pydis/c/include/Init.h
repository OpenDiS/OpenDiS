/*****************************************************************************
 *
 *	Init.h		Prototypes for the initialization routines
 *
 ****************************************************************************/

#ifndef _Init_h
#define _Init_h

#include "Home.h"
#include "InData.h"

void   SetBoxSize(Param_t *param);
void   InitRecycleNodeHeap(Home_t *home);
void   InitCellDomains(Home_t *home);
void   InitCellNatives(Home_t *home);
void   InitCellNeighbors(Home_t *home);
Home_t *InitHome(void);
void   Initialize(Home_t *home,int argc, char *argv[]);
int    OpenDir(Home_t *home);
void   ParadisInit(int argc, char *argv[], Home_t **homeptr);
void   RecvInitialNodeData(Home_t *home);
void   SendInitialNodeData(Home_t *home, InData_t *inData, int *msgCount,
           int **nodeLists, int *listCounts, int *nextAvailableTag);
void   SetRemainingDefaults(Home_t *home);

#endif
