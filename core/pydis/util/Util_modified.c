#include <stdlib.h>
#include "Home.h"

/* Functions modified from ParaDiS to be used by pydis */

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
void Initialize(Home_t *home, int argc, char *argv[])
{
    Param_t *param;
    printf("Stub.c: Initialize being implemented\n");

    home->ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    home->dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
    //home->subcyc = (Subcyc_t *)calloc(1, sizeof(Subcyc_t));

    home->param = (Param_t *)calloc(1, sizeof(Param_t));
    param = home->param;

    CtrlParamInit(param, home->ctrlParamList);
    DataParamInit(param, home->dataParamList);

    SetBoxSize(param);
    SetRemainingDefaults(home);

    if (home->myDomain == 0) {
        DisableUnneededParams(home);
    }

    InitCellNatives(home);
    InitCellNeighbors(home);
    InitCellDomains(home);

    InitOpList(home);
    //InitOpRecList(home);

    //ReadRijm(home);
    //ReadRijmPBC(home);

    //SortNativeNodes(home);
    //if (param->useLabFrame) {}

}
