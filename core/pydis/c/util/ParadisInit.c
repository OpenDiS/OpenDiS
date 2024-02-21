/**************************************************************************
 *
 *      Mopdule:     ParadisInit.c
 *      Description: Contains functions for doing one-time initializations
 *                   before entering the main loop of the application.
 *
 *      Includes functions:
 *          ParadisInit()
 *
 *************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include "Home.h"
#include "Init.h"
#include "Util.h"

#ifdef PARALLEL
#include "mpi.h"
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:    ParadisInit
 *      Description: Create the 'home' structure, setup timing categories
 *                   and initialize MPI.
 *
 *------------------------------------------------------------------------*/
void ParadisInit(int argc, char *argv[], Home_t **homeptr)
{
        Home_t         *home;

        home = InitHome();
        *homeptr = home;
    
        TimerInit(home);
    
#ifdef PARALLEL
        MPI_Init(&argc, &argv); 
        MPI_Comm_rank(MPI_COMM_WORLD, &home->myDomain);
        MPI_Comm_size(MPI_COMM_WORLD, &home->numDomains);
#else
        home->myDomain = 0;
        home->numDomains = 1;
#endif
    
        TimerStart(home, TOTAL_TIME);

        TimerStart(home, INITIALIZE);
        Initialize(home,argc,argv);  
        TimerStop(home, INITIALIZE);
    
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (home->myDomain == 0) printf("ParadisInit finished\n");

        return;
}
