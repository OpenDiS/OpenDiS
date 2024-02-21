#include <stdlib.h>
#include "Home.h"

/* Functions modified from ParaDiS to be used by pydis */

/* 
   ParadisInit            -- ParadisInit.c
*/

/*-------------------------------------------------------------------------
 *
 *      Function:    ParadisInit - (modified to not use argc, argv)
 *      Description: Create the 'home' structure, setup timing categories
 *                   and initialize MPI.
 *
 *------------------------------------------------------------------------*/
void ParadisInit_lean(Home_t **homeptr)
{
       Home_t         *home;

       home = InitHome();
       *homeptr = home;
    
       TimerInit(home);
    
       home->myDomain = 0;
       home->numDomains = 1;
    
       TimerStart(home, TOTAL_TIME);

       TimerStart(home, INITIALIZE);
       Initialize_lean(home);  
       TimerStop(home, INITIALIZE);
    
       if (home->myDomain == 0) printf("ParadisInit finished\n");

       return;
}


/* 
   Initialize            -- Initialize.c
*/

/*---------------------------------------------------------------------------
 *
 *      Function:     Initialize - (modified to not use argc, argv)
 *      Description:  This is the driver routine for initialization,
 *                    handling some of the general initializations and
 *                    calling all the more specific initialization routines.
 *
 *-------------------------------------------------------------------------*/
void Initialize_lean(Home_t *home)
{
    Param_t *param;
    printf("Initialize_lean:\n");

    home->ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    home->dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    //home->subcyc = (Subcyc_t *)calloc(1, sizeof(Subcyc_t));

    home->param = (Param_t *)calloc(1, sizeof(Param_t));
    param = home->param;

    CtrlParamInit(param, home->ctrlParamList);
    DataParamInit(param, home->dataParamList);

    SetBoxSize(param);

    // call this in python
    //SetRemainingDefaults(home);

    if (home->myDomain == 0) {
        DisableUnneededParams(home);
    }

    InitCellNatives(home);
    InitCellNeighbors(home);

    // call this in python
    //InitCellDomains(home);

    InitOpList(home);
    //InitOpRecList(home);

    //ReadRijm(home);
    //ReadRijmPBC(home);

    //SortNativeNodes(home);
    //if (param->useLabFrame) {}

    printf("Initialize_lean: completed\n");
}
