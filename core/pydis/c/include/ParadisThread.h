/*--------------------------------------------------------------------------
 *
 *	ParadisThread.h	
 *
 *------------------------------------------------------------------------*/

#ifndef _ParadisThread_h
#define _ParadisThread_h

/*
 *      If we using OpenMP, define the locking mechanisms to be used
 *      when synchronizing access to a shared data item.  If OpenMP
 *      is not being used, the code is single threaded there's no need
 *      for locking, so calls to the locking mechanism are completely
 *      eliminated.
 */
#ifdef _OPENMP

#include <omp.h>

#define INIT_LOCK(a)    omp_init_lock((a))
#define LOCK(a)         omp_set_lock((a))
#define UNLOCK(a)       omp_unset_lock((a))
#define DESTROY_LOCK(a) omp_destroy_lock((a))

#else /* _OPENMP not defined */

#define INIT_LOCK(a)
#define LOCK(a)
#define UNLOCK(a)
#define DESTROY_LOCK(a)

#endif /* end ifdef _OPENMP */

/*
 *      Prototype any functions that we need to handle threads
 */
void GetThreadIterationIndices(int loopCount, int *threadID, int *startIndex,
        int *endIndex);

#endif /* ifndef _ParadisThread_h */
