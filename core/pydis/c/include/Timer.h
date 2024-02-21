/*************************************************************************
 *
 *  Timer.h - define the timing structures
 *
 ************************************************************************/

#ifndef _Timer_h
#define _Timer_h

#include "Typedefs.h"

struct _timer {
	real8	startTime; /* time at which most recent TimerStart called */
	real8	incr;      /* time of most recent event for this event type */
	real8	accum;     /* accumulated time for this event type */
	real8	save;      /* to save full force update times till */
			   /* next WriteStep */
	int	started;   /* 1 if event is active, 0 otherwise */
	char	*name;     /* label used during TimerPrint */
};


enum {
    TOTAL_TIME = 0,
    INITIALIZE,
	TIME_INTEGRATE,
#ifdef _SUBCYCLING
	SEG_LIST_MAKER,
#endif
#ifdef _GPU_SUBCYCLE
	SUBCYCLING_GPU,
#endif
    SORT_NATIVE_NODES,
    COMM_SEND_GHOSTS,
    GHOST_COMM_BARRIER,
    CELL_CHARGE,
    CELL_CHARGE_BARRIER,
    CALC_FORCE,
    LOCAL_FORCE,
    REMOTE_FORCE,
    CALC_FORCE_BARRIER,
    CALC_VELOCITY,
    CALC_VELOCITY_BARRIER,
    COMM_SEND_VELOCITY,
#ifdef _SUBCYCLING
	COMM_SEND_VELOCITYSUB,
    COMM_SEND_COORD,
#endif
    SPLIT_MULTI_NODES,
    COLLISION_HANDLING,
    POST_COLLISION_BARRIER,
    COL_SEND_REMESH,
    COL_FIX_REMESH,
    COL_FORCE_UPDATE,
    GENERATE_IO,
    PLOT,
    IO_BARRIER,
    REMESH_START_BARRIER,
    REMESH,
    SEND_REMESH,
    FIX_REMESH,
    FORCE_UPDATE_REMESH,
    REMESH_END_BARRIER ,
    MIGRATION,
    MIGRATION_BARRIER,
    LOADCURVE,
    LOAD_BALANCE,
    SEGFORCE_COMM,
    TIMER_BLOCK_SIZE  /* MUST BE LAST IN THE LIST */
};

/*
 *      Prototype the timer functions
 */
void TimeAtRestart(Home_t *home, int stage);
void TimerClear(Home_t *home, int index);
void TimerClearAll(Home_t *home);
void TimerInit(Home_t *home);
void TimerInitDLBReset(Home_t *home);
void TimerPrint(Home_t *home);
void TimerReinitialize(Home_t *home);
void TimerStart(Home_t *home, int index);
void TimerStop(Home_t *home, int index);
void TimerSave(Home_t *home, int index);

#endif
