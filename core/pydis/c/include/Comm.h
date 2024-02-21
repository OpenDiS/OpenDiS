/*************************************************************************
 *
 *  Comm.h - define various constants used during MPI communications
 *
 *************************************************************************/

#ifndef _Comm_h
#define _Comm_h

#include "Home.h"

#ifdef _FEM
#define FLTS_PER_GHOST_NODE      21
#define FLTS_PER_GHOST2_NODE     23
#else
#define FLTS_PER_GHOST_NODE      16
#define FLTS_PER_GHOST2_NODE     18
#endif

#define FLTS_PER_GHOST_ARM       11
#define FLTS_PER_GHOST_CELL       3
#define EXTRA_GHOST_FLTS          4

#define INIT_PARAMS         1
#define INIT_VALS_PER_NODE  7
#define INIT_VALS_PER_ARM   8

#define INTS_PER_TAG        2
/*
 *	For downloading entire problem space to domain zero.
 *	Used for various plotting functions and writing text
 * 	restart files.
 */
#define INTS_PER_MIRROR_NODE 4
#define INTS_PER_MIRROR_ARM  2
#define FLTS_PER_MIRROR_NODE 6
#define FLTS_PER_MIRROR_ARM  9
#define EXTRA_MIRROR_INTS    6

/*
 *	message tags
 */
#define MSG_INIT_LENS     1000
#define MSG_INIT_NODES    1002
#define MSG_GHOST_LEN     1020
#define MSG_GHOST         1030
#define MSG_MIG_LEN       1040
#define MSG_MIG           1050
#define MSG_OLDNEW_LEN    1060
#define MSG_OLDNEW        1070
#define MSG_MIRRORS       1090
#define	MSG_SEND_MIRROR   1091
#define MSG_GHOST2_REQ_LEN      1100
#define MSG_GHOST2_REQ          1101
#define MSG_GHOST2_RESPONSE_LEN 1102
#define MSG_GHOST2_RESPONSE     1103
#define MSG_VELOCITY      2000
#define MSG_VELOCITY_LEN  2001
#define MSG_OPLIST_LEN    2010
#define MSG_OPLIST        2020
#define MSG_TOKEN_RING    2030
#define MSG_REMESH_LEN    2040
#define MSG_REMESH        2050
#define	MSG_TAGREMAP_LEN  2060
#define	MSG_TAGREMAP	  2070
#define MSG_SEGDATA_LEN   2080
#define MSG_SEGDATA       2090
#define MSG_VISIT_COUNTS  2100

/*
 *	Prototypes
 */
void CommSendGhosts(Home_t *home);
void CommSendGhostPlanes(Home_t *home);
void CommSendRemesh(Home_t *home);
void CommSendSecondaryGhosts(Home_t *home);
void CommSendSegments(Home_t *home, int numRecvBufs, int numSendBufs,
         int *sendDomList, Segment_t **cellSegLists, int *cellSegCnts);
void CommSendVelocity(Home_t *home);

#ifdef _SUBCYCLING
void CommSendCoord(Home_t *home, int subGroup);
void CommSendVelocitySub(Home_t *home, int subGroup);
#endif

#endif
