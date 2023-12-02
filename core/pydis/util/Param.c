/**************************************************************************
 *
 *      Module:       Param.c
 *      Description:  Contains functions for binding control and data
 *                    file parameter names to code variables.
 *
 *      Public functions:
 *          CtrlParamInit()
 *          DataParamInit()
 *          MarkParamDisabled()
 *          MarkParamEnabled()
 *
 *************************************************************************/
#include <stdarg.h>
#include <ctype.h>
#include "Home.h"
#include "Parse.h"


/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamDisabled
 *      Description:  Explicitly mark a control file parameter 
 *                    to be disabled.  Primarily used to
 *                    prevent writing to the restart file any
 *                    control file parameters that are not
 *                    applicable or appropriate to the current
 *                    execution of the code.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly disabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamDisabled(ParamList_t *CPList, char *name)
{
        int  paramIndex;

        if ((paramIndex = LookupParam(CPList, name)) >= 0) {
            CPList->varList[paramIndex].flags |= VFLAG_DISABLED;
        } else {
            Fatal("MarkParamDisabled: unknown parameter %s\n", name);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamEnabled
 *      Description:  Explicitly mark a control file parameter 
 *                    to be enabled.  Primarily used when
 *                    explicitly controlling whether a control
 *                    parameter is to be written to a restart
 *                    file.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly enabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamEnabled(ParamList_t *CPList, char *name)
{
        int  paramIndex;

        if ((paramIndex = LookupParam(CPList, name)) >= 0) {
            CPList->varList[paramIndex].flags &= ~VFLAG_DISABLED;
        } else {
            Fatal("MarkParamEnabled: unknown parameter %s\n", name);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     CtrlParamInit
 *      Description:  Bind all valid control file parameters to
 *                    the associated code variables.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *
 *----------------------------------------------------------------------*/
void CtrlParamInit(Param_t *param, ParamList_t *CPList)
{
/*
 *      Note: Parameters need only be initialized if their
 *      default values are non-zero.
 */
        CPList->paramCnt = 0;
        CPList->varList = (VarData_t *)NULL;

/*
 *      Simulation cell and processor setup
 */
        BindVar(CPList, "Simulation cell and processor setup",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "numXdoms", &param->nXdoms, V_INT, 1, VFLAG_NULL);
        param->nXdoms = 1;

        BindVar(CPList, "numYdoms", &param->nYdoms, V_INT, 1, VFLAG_NULL);
        param->nYdoms = 1;

        BindVar(CPList, "numZdoms", &param->nZdoms, V_INT, 1, VFLAG_NULL);
        param->nZdoms = 1;

        BindVar(CPList, "numXcells", &param->nXcells, V_INT, 1, VFLAG_NULL);
        param->nXcells = 3;   /* Must be >= 3 */

        BindVar(CPList, "numYcells", &param->nYcells, V_INT, 1, VFLAG_NULL);
        param->nYcells = 3;   /* Must be >= 3 */

#ifdef _BGP
        BindVar(CPList, "taskMappingMode", &param->taskMappingMode, V_INT, 1,
                VFLAG_NULL);
        param->taskMappingMode = 1;  /* use user-supplied decomp by default */
#endif

        BindVar(CPList, "numZcells", &param->nZcells, V_INT, 1, VFLAG_NULL);
        param->nZcells = 3;   /* Must be >= 3 */

        BindVar(CPList, "xBoundType", &param->xBoundType, V_INT, 1, VFLAG_NULL);
        param->xBoundType = Periodic;

        BindVar(CPList, "yBoundType", &param->yBoundType, V_INT, 1, VFLAG_NULL);
        param->yBoundType = Periodic;

        BindVar(CPList, "zBoundType", &param->zBoundType, V_INT, 1, VFLAG_NULL);
        param->zBoundType = Periodic;

        BindVar(CPList, "decompType", &param->decompType, V_INT, 1, VFLAG_NULL);
        param->decompType = 1;

        BindVar(CPList, "DLBfreq", &param->DLBfreq, V_INT, 1, VFLAG_NULL);
        param->DLBfreq = 3;

        BindVar(CPList, "xBoundMin", &param->xBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "xBoundMax", &param->xBoundMax, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "yBoundMin", &param->yBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "yBoundMax", &param->yBoundMax, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "zBoundMin", &param->zBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "zBoundMax", &param->zBoundMax, V_DBL, 1, VFLAG_NULL);

/*
 *      Simulation time and timestepping controls
 */
        BindVar(CPList, "Simulation time and timestepping controls",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "cycleStart", &param->cycleStart, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "maxstep", &param->maxstep, V_INT, 1, VFLAG_NULL);
        param->maxstep = 100;

        BindVar(CPList, "timeNow", &param->timeNow, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "timeStart", &param->timeStart, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "timestepIntegrator", param->timestepIntegrator,
                V_STRING, 1, VFLAG_NULL);
        strcpy(param->timestepIntegrator, "trapezoid");
#ifdef _SUBCYCLING
        BindVar(CPList, "subInteg0Integ1", param->subInteg0Integ1,
                V_STRING, 1, VFLAG_NULL);
        strcpy(param->subInteg0Integ1, "None");
#endif		
		BindVar(CPList, "deltaTT", &param->deltaTT, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "maxDT", &param->maxDT, V_DBL, 1, VFLAG_NULL);
        param->maxDT = 1.0e-07;

        BindVar(CPList, "nextDT", &param->nextDT, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "dtIncrementFact", &param->dtIncrementFact,
                V_DBL, 1, VFLAG_NULL);
        param->dtIncrementFact = 1.2;

        BindVar(CPList, "dtDecrementFact", &param->dtDecrementFact,
                V_DBL, 1, VFLAG_NULL);
        param->dtDecrementFact = 0.5;

        BindVar(CPList, "dtExponent", &param->dtExponent, V_DBL, 1, VFLAG_NULL);
        param->dtExponent = 4.0;

        BindVar(CPList, "dtVariableAdjustment",
                &param->dtVariableAdjustment, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "rTol", &param->rTol, V_DBL, 1, VFLAG_NULL);
        param->rTol = -1.0;

        BindVar(CPList, "rmax", &param->rmax, V_DBL, 1, VFLAG_NULL);
        param->rmax = 100; 

#ifdef _SUBCYCLING
		BindVar(CPList, "rTolrel" , &param->rTolrel , V_DBL, 1, VFLAG_NULL);
		param->rTolrel = 0.01;  //  param->rTolrel = 1e6;

		BindVar(CPList, "rTolth"  , &param->rTolth  , V_DBL, 1, VFLAG_NULL);
		param->rTolth = 0.10;  //   param->rTolth = 1e-6;

		BindVar(CPList, "renh"    , &param->renh    , V_DBL, 1, VFLAG_NULL);
		param->renh = 0.0;

		BindVar(CPList, "rg1"     , &param->rg1     , V_DBL, 1, VFLAG_NULL);
		param->rg1 = 0.0;

		BindVar(CPList, "rg2"     , &param->rg2     , V_DBL, 1, VFLAG_NULL);
		param->rg2 = 0.0;

		BindVar(CPList, "rg3"     , &param->rg3     , V_DBL, 1, VFLAG_NULL);
		param->rg3 = 0.0;

		BindVar(CPList, "rg4"     , &param->rg4     , V_DBL, 1, VFLAG_NULL);
		param->rg4 = 0.0;

		BindVar(CPList, "nTry"    , &param->nTry    , V_INT, 1, VFLAG_NULL);
		param->nTry = 0;

		BindVar(CPList, "sendSubGroupForc", &param->sendSubGroupForc, V_INT, 1, VFLAG_NULL);
		param->sendSubGroupForc = 0;
#endif

/*
 *      Discretization controls and controls for topological changes
 */
        BindVar(CPList, "Discretization and topological change controls",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "maxSeg", &param->maxSeg, V_DBL, 1, VFLAG_NULL);
        param->maxSeg = -1.0;

        BindVar(CPList, "minSeg", &param->minSeg, V_DBL, 1, VFLAG_NULL);
        param->minSeg = -1.0;

        BindVar(CPList, "remeshRule", &param->remeshRule, V_INT, 1, VFLAG_NULL);
        param->remeshRule = 2;

        BindVar(CPList, "splitMultiNodeFreq", &param->splitMultiNodeFreq,
                V_INT, 1, VFLAG_NULL);
        param->splitMultiNodeFreq = 1;

        BindVar(CPList, "collisionMethod", &param->collisionMethod, V_INT, 1,
                VFLAG_NULL);
#ifdef _RETROCOLLISIONS
        param->collisionMethod = 4;
#else
        param->collisionMethod = 2;
#endif

#ifdef _RETROCOLLISIONS
        BindVar(CPList, "rann", &param->rann, V_DBL, 1, VFLAG_NULL);
        param->rann = -1;
#endif

/*
 *      FMM controls
 */
        BindVar(CPList, "Fast Multipole Method controls", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "fmEnabled", &param->fmEnabled, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fmMPOrder", &param->fmMPOrder, V_INT, 1, VFLAG_NULL);
        param->fmMPOrder = 2;

        BindVar(CPList, "fmTaylorOrder", &param->fmTaylorOrder, V_INT,
                1, VFLAG_NULL);
        param->fmTaylorOrder = 5;

        BindVar(CPList, "fmCorrectionTbl", param->fmCorrectionTbl, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->fmCorrectionTbl,"inputs/fm-ctab.Ta.600K.0GPa.m2.t5.dat");

/*
 *      Identify tables needed for remote force calculations if
 *      FMM is not enabled
 */
        BindVar(CPList, "Tables for non-FMM far-field force calcs",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "Rijmfile", param->Rijmfile, V_STRING, 1, VFLAG_NULL);
        strcpy(param->Rijmfile,"inputs/Rijm.cube.out");

        BindVar(CPList, "RijmPBCfile", param->RijmPBCfile, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->RijmPBCfile,"inputs/RijmPBC.cube.out");

/*
 *      Loading condition parameters
 */
        BindVar(CPList, "Loading conditions", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "TempK", &param->TempK, V_DBL, 1, VFLAG_NULL);
        param->TempK = 600;

        BindVar(CPList, "loadType", &param->loadType, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "appliedStress", param->appliedStress, V_DBL, 6,
                VFLAG_NULL);

        BindVar(CPList, "appliedStressRate", param->appliedStressRate, V_DBL, 6,
                VFLAG_NULL);

        BindVar(CPList, "eRate", &param->eRate, V_DBL, 1, VFLAG_NULL);
        param->eRate = 1.0;

        BindVar(CPList, "indxErate", &param->indxErate, V_INT, 1, VFLAG_NULL);
        param->indxErate = 1;

        BindVar(CPList, "edotdir",param->edotdir, V_DBL, 3, VFLAG_NULL);
        param->edotdir[0] = 1.0;
        param->edotdir[1] = 0.0;
        param->edotdir[2] = 0.0;

        BindVar(CPList, "sRate", &param->sRate, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "cTimeOld", &param->cTimeOld, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "dCyclicStrain", &param->dCyclicStrain, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "netCyclicStrain", &param->netCyclicStrain,
                V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "numLoadCycle", &param->numLoadCycle, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "eAmp", &param->eAmp, V_DBL, 1, VFLAG_NULL);

/*
 *      To specify input sample axes in laboratory frame
 */
        BindVar(CPList, "useLabFrame", &param->useLabFrame, V_INT, 1,
                VFLAG_NULL);
        param->useLabFrame = 0;

        BindVar(CPList, "labFrameXDir", param->labFrameXDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameXDir[0] = 1.0;
        param->labFrameXDir[1] = 0.0;
        param->labFrameXDir[2] = 0.0;

        BindVar(CPList, "labFrameYDir", param->labFrameYDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameYDir[0] = 0.0;
        param->labFrameYDir[1] = 1.0;
        param->labFrameYDir[2] = 0.0;

        BindVar(CPList, "labFrameZDir", param->labFrameZDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameZDir[0] = 0.0;
        param->labFrameZDir[1] = 0.0;
        param->labFrameZDir[2] = 1.0;

/*
 *      Parameters for material specific constants and mobility values
 */
        BindVar(CPList, "Material and mobility parameters",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "mobilityLaw", param->mobilityLaw, V_STRING, 1,
                VFLAG_NULL);
        strcpy(param->mobilityLaw, "BCC_0");

        BindVar(CPList, "vacancyConc", &param->vacancyConc, V_DBL, 1,
                VFLAG_NULL);
        param->vacancyConc = -1.0; /* signifies not used */

        BindVar(CPList, "vacancyConcEquilibrium",
                &param->vacancyConcEquilibrium, V_DBL, 1, VFLAG_NULL);
        param->vacancyConcEquilibrium = -1.0; /* signifies not used */

        BindVar(CPList, "shearModulus", &param->shearModulus, V_DBL, 1,
                VFLAG_NULL);
        param->shearModulus = 6.488424e+10; /* Ta: units=Pa @t=600K, p=0Pa*/

        BindVar(CPList, "pois", &param->pois, V_DBL, 1, VFLAG_NULL);
        param->pois = 3.327533e-01;    /* Ta: temp=600K, pressure=0Pa */

        BindVar(CPList, "burgMag", &param->burgMag, V_DBL, 1, VFLAG_NULL);
        param->burgMag = 2.875401e-10;  /* Ta: units=m @t=600K, pressure=0Pa */

        BindVar(CPList, "YoungModulus", &param->YoungsModulus, V_DBL, 1,
                VFLAG_NULL);
        param->YoungsModulus=172.95e+09; /* Ta: units=Pa @t=600K, p=0Pa */

        BindVar(CPList, "rc", &param->rc, V_DBL, 1, VFLAG_NULL);
        param->rc = -1.0;

        BindVar(CPList, "Ecore", &param->Ecore, V_DBL, 1, VFLAG_NULL);
        param->Ecore = -1.0;

        BindVar(CPList, "MobScrew", &param->MobScrew, V_DBL, 1, VFLAG_NULL);
        param->MobScrew = 10.0;

        BindVar(CPList, "MobEdge", &param->MobEdge, V_DBL, 1, VFLAG_NULL);
        param->MobEdge  = 10.0;

        BindVar(CPList, "MobClimb", &param->MobClimb, V_DBL, 1, VFLAG_NULL);
        param->MobClimb = 1.0e-02;

        BindVar(CPList, "FricStress", &param->FricStress, V_DBL, 1, VFLAG_NULL);
        param->FricStress = 0.0;

//#ifdef _THERMAL_ACTIVATED_CROSSSLIP
        BindVar(CPList, "thermalCrossSlip", &param->thermalCrossSlip, V_INT, 1, VFLAG_NULL);
        param->thermalCrossSlip = 0;
        
        BindVar(CPList, "CRS_A", &param->CRS_A, V_DBL, 1, VFLAG_NULL);
        param->CRS_A =  2.100;
        
        BindVar(CPList, "CRS_T0", &param->CRS_T0, V_DBL, 1, VFLAG_NULL);
        param->CRS_T0 = 5.5949e9;
        
        BindVar(CPList, "pre_Cge", &param->pre_Cge, V_DBL, 1, VFLAG_NULL);
        param->pre_Cge = 1.0 ;
//#endif
	
		BindVar(CPList, "CRScount_bothScrew", &param->CRScount_bothScrew, V_INT, 1, VFLAG_NULL);
        param->CRScount_bothScrew = 0;
        BindVar(CPList, "CRScount_zipper", &param->CRScount_zipper, V_INT, 1, VFLAG_NULL);
        param->CRScount_zipper = 0;

/*
 *      List of burgers vectors/line directions to be considered sessile
 */
        BindVar(CPList, "sessileburgspec", param->sessileburgspec,
                V_DBL, 30, VFLAG_NULL);

        BindVar(CPList, "sessilelinespec", param->sessilelinespec,
                V_DBL, 30, VFLAG_NULL);

/*
 *      A couple values for including (if available) inertial terms in
 *      the mobility function.
 */
        BindVar(CPList, "includeInertia", &param->includeInertia, V_INT,
                1, VFLAG_NULL);

        BindVar(CPList, "massDensity", &param->massDensity, V_DBL,
                1, VFLAG_NULL);
        param->massDensity = -1.0;

/*
 *      Flux decomposition
 */
        BindVar(CPList, "Flux decomposition",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "totstraintensor", param->totstraintensor, V_DBL, 6,
                VFLAG_NULL);

        BindVar(CPList, "totpStn", param->totpStn, V_DBL, 6, VFLAG_NULL);

        BindVar(CPList, "totpSpn", param->totpSpn, V_DBL, 6, VFLAG_NULL);

        BindVar(CPList, "Ltot", param->Ltot, V_DBL, 16, VFLAG_NULL);
        BindVar(CPList, "FCC_Ltot", param->FCC_Ltot, V_DBL, 24, VFLAG_NULL);

        BindVar(CPList, "fluxtot", param->fluxtot, V_DBL, 28, VFLAG_NULL);
        BindVar(CPList, "FCC_fluxtot", param->FCC_fluxtot, V_DBL, 42,
                VFLAG_NULL);

/*
 *      Total system density: this is informational only and recalculated
 *      each time a restart file is written.
 */
        BindVar(CPList, "Total density. Informational only; ignored on input",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "disloDensity", &param->disloDensity, V_DBL, 1,
                VFLAG_NULL);

/*
 *      Velocity statistics
 */
        BindVar(CPList, "Velocity statistics",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "vAverage", &param->vAverage, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "vStDev", &param->vStDev, V_DBL, 1, VFLAG_NULL);


/*
 *      I/O controls and options
 */
        BindVar(CPList, "I/O controls and parameters", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "dirname", param->dirname, V_STRING, 1, VFLAG_NULL);

        BindVar(CPList, "writeBinRestart", &param->writeBinRestart, V_INT,
                1, VFLAG_NULL);

        BindVar(CPList, "skipIO", &param->skipIO, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "numIOGroups", &param->numIOGroups, V_INT, 1,
                VFLAG_NULL);
        param->numIOGroups = 1;

/* 
 *      segment/arm files
 */
        BindVar(CPList, "armfile", &param->armfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "armfilefreq", &param->armfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->armfilefreq = 100;

        BindVar(CPList, "armfiledt", &param->armfiledt, V_DBL, 1, VFLAG_NULL);
        param->armfiledt = -1.0;

        BindVar(CPList, "armfiletime", &param->armfiletime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "armfilecounter", &param->armfilecounter, V_INT, 1,
                VFLAG_NULL);


/* 
 *      flux decomposition files
 */
        BindVar(CPList, "fluxfile", &param->fluxfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fluxfreq", &param->fluxfreq, V_INT, 1, VFLAG_NULL);
        param->fluxfreq = 100;

        BindVar(CPList, "fluxdt", &param->fluxdt, V_DBL, 1, VFLAG_NULL);
        param->fluxdt = -1.0;

        BindVar(CPList, "fluxtime", &param->fluxtime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "fluxcounter", &param->fluxcounter, V_INT, 1,
                VFLAG_NULL);


/* 
 *      chain fragment files -- to be post-processed into a data format
 *      suitable for use with the VisIt visualization software.
 */
        BindVar(CPList, "fragfile", &param->fragfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fragfreq", &param->fragfreq, V_INT, 1, VFLAG_NULL);
        param->fragfreq = 100;

        BindVar(CPList, "fragdt", &param->fragdt, V_DBL, 1, VFLAG_NULL);
        param->fragdt = -1.0;

        BindVar(CPList, "fragtime", &param->fragtime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "fragcounter", &param->fragcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      gnuplot files
 */
        BindVar(CPList, "gnuplot", &param->gnuplot, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "gnuplotfreq", &param->gnuplotfreq, V_INT, 1,
                VFLAG_NULL);
        param->gnuplotfreq = 100;

        BindVar(CPList, "gnuplotdt", &param->gnuplotdt, V_DBL, 1, VFLAG_NULL);
        param->gnuplotdt = -1.0;

        BindVar(CPList, "gnuplottime", &param->gnuplottime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "gnuplotcounter", &param->gnuplotcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      pole figure files
 */
        BindVar(CPList, "polefigfile", &param->polefigfile, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "polefigfreq", &param->polefigfreq, V_INT, 1,
                VFLAG_NULL);
        param->polefigfreq = 100;

        BindVar(CPList, "polefigdt", &param->polefigdt, V_DBL, 1, VFLAG_NULL);
        param->polefigdt = -1.0;

        BindVar(CPList, "polefigtime", &param->polefigtime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "polefilecounter", &param->polefigcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      povray files
 */
        BindVar(CPList, "povray", &param->povray, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "povrayfreq", &param->povrayfreq, V_INT, 1, VFLAG_NULL);
        param->povrayfreq = 100;

        BindVar(CPList, "povraydt", &param->povraydt, V_DBL, 1, VFLAG_NULL);
        param->povraydt = -1.0;

        BindVar(CPList, "povraytime", &param->povraytime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "povraycounter", &param->povraycounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      atomeye files
 */
        BindVar(CPList, "atomeye", &param->atomeye, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "atomeyefreq", &param->atomeyefreq, V_INT, 1, VFLAG_NULL);
        param->povrayfreq = 100;

        BindVar(CPList, "atomeyedt", &param->atomeyedt, V_DBL, 1, VFLAG_NULL);
        param->povraydt = -1.0;

        BindVar(CPList, "atomeyetime", &param->atomeyetime, V_DBL, 1, VFLAG_NULL);
        BindVar(CPList, "atomeyesegradius", &param->atomeyesegradius, V_DBL, 1, VFLAG_NULL);
        param->atomeyesegradius = 500;

        BindVar(CPList, "atomeyecounter", &param->atomeyecounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      postscript files
 */
        BindVar(CPList, "psfile", &param->psfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "psfilefreq", &param->psfilefreq, V_INT, 1, VFLAG_NULL);
        param->psfilefreq = 100;

        BindVar(CPList, "psfiledt", &param->psfiledt, V_DBL, 1, VFLAG_NULL);
        param->psfiledt = -1.0;

        BindVar(CPList, "psfiletime", &param->psfiletime, V_DBL, 1, VFLAG_NULL);


/*
 *      restart files
 */
        BindVar(CPList, "savecn", &param->savecn, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "savecnfreq", &param->savecnfreq, V_INT, 1, VFLAG_NULL);
        param->savecnfreq = 100;

        BindVar(CPList, "savecndt", &param->savecndt, V_DBL, 1, VFLAG_NULL);
        param->savecndt = -1.0;

        BindVar(CPList, "savecntime", &param->savecntime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "savecncounter", &param->savecncounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      properties files
 */
        BindVar(CPList, "saveprop", &param->saveprop, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "savepropfreq", &param->savepropfreq, V_INT, 1,
                VFLAG_NULL);
        param->savepropfreq = 100;

        BindVar(CPList, "savepropdt", &param->savepropdt, V_DBL, 1, VFLAG_NULL);
        param->savepropdt = -1.0;

        BindVar(CPList, "saveproptime", &param->saveproptime, V_DBL, 1,
                VFLAG_NULL);


/*
 *      timer files
 */
        BindVar(CPList, "savetimers", &param->savetimers, V_INT, 1, VFLAG_NULL);
        param->savetimers = 0;

        BindVar(CPList, "savetimersfreq", &param->savetimersfreq, V_INT, 1,
                VFLAG_NULL);
        param->savetimersfreq = 100;

        BindVar(CPList, "savetimersdt", &param->savetimersdt, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "savetimerstime", &param->savetimerstime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "savetimerscounter", &param->savetimerscounter,
                V_INT, 1, VFLAG_NULL);
        param->savetimerscounter = 0;


/*
 *      tecplot files
 */
        BindVar(CPList, "tecplot", &param->tecplot, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "tecplotfreq", &param->tecplotfreq, V_INT, 1,
                VFLAG_NULL);
        param->tecplotfreq = 100;

        BindVar(CPList, "tecplotdt", &param->tecplotdt, V_DBL, 1, VFLAG_NULL);
        param->tecplotdt = -1.0;

        BindVar(CPList, "tecplottime", &param->tecplottime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "tecplotcounter", &param->tecplotcounter, V_INT, 1,
                VFLAG_NULL);
                
                
/*
 *      paraview files
 */
        BindVar(CPList, "paraview", &param->paraview, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "paraviewfreq", &param->paraviewfreq, V_INT, 1,
                VFLAG_NULL);
        param->paraviewfreq = 100;

        BindVar(CPList, "paraviewdt", &param->paraviewdt, V_DBL, 1, VFLAG_NULL);
        param->paraviewdt = -1.0;

        BindVar(CPList, "paraviewtime", &param->paraviewtime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "paraviewcounter", &param->paraviewcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      nodal velocity data
 */
        BindVar(CPList, "velfile", &param->velfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "velfilefreq", &param->velfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->velfilefreq = 100;

        BindVar(CPList, "velfiledt", &param->velfiledt, V_DBL, 1, VFLAG_NULL);
        param->velfiledt = -1.0;

        BindVar(CPList, "velfiletime", &param->velfiletime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "velfilecounter", &param->velfilecounter, V_INT, 1,
                VFLAG_NULL);

/*
 *      nodal force data
 */
        BindVar(CPList, "writeForce", &param->writeForce, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeForceFreq", &param->writeForceFreq, V_INT, 1,
                VFLAG_NULL);
        param->writeForceFreq = 100;

        BindVar(CPList, "writeForceDT", &param->writeForceDT, V_DBL, 1,
                VFLAG_NULL);
        param->writeForceDT = -1.0;

        BindVar(CPList, "writeForceTime", &param->writeForceTime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "writeForceCounter", &param->writeForceCounter, V_INT,
                1, VFLAG_NULL);

/*
 *      links data
 */
        BindVar(CPList, "linkfile", &param->linkfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "linkfilefreq", &param->linkfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->linkfilefreq = 100;

        BindVar(CPList, "linkfiledt", &param->linkfiledt, V_DBL, 1, VFLAG_NULL);
        param->linkfiledt = -1.0;

        BindVar(CPList, "linkfiletime", &param->linkfiletime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "linkfilecounter", &param->linkfilecounter, V_INT, 1,
                VFLAG_NULL);

#ifdef _OP_REC
/*
 *      topological operations record
 */
        BindVar(CPList, "oprec", &param->oprec, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "oprecmove", &param->oprecmove, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "oprecwritefreq", &param->oprecwritefreq, V_INT, 1,
                VFLAG_NULL);
        param->oprecwritefreq = 100;
	
        BindVar(CPList, "oprecfilefreq", &param->oprecfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->oprecfilefreq = 100;
	
        BindVar(CPList, "opreccounter", &param->opreccounter, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "runfromoprec", &param->runfromoprec, V_INT, 1,
                VFLAG_NULL);
        param->runfromoprec = 0;

        BindVar(CPList, "oprecfile", param->oprecfile, V_STRING, 1, VFLAG_NULL);
        strcpy(param->oprecfile, "oprec.dat");
#endif

/*
 *      Parameters related to creation of VisIt output files
 */
        BindVar(CPList, "writeVisit", &param->writeVisit, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitFreq", &param->writeVisitFreq,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitCounter", &param->writeVisitCounter,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitSegments", &param->writeVisitSegments,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitSegmentsAsText",
                &param->writeVisitSegmentsAsText, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitNodes", &param->writeVisitNodes,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitNodesAsText",
                &param->writeVisitNodesAsText, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitDT", &param->writeVisitDT,
                V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitTime", &param->writeVisitTime,
                V_DBL, 1, VFLAG_NULL);

        param->writeVisitFreq = 100;
        param->writeVisitDT = -1.0;

/*
 *      X-window display
 */
        BindVar(CPList, "winDefaultsFile", param->winDefaultsFile, V_STRING, 1,
                VFLAG_NULL);
        strcpy(param->winDefaultsFile, "inputs/paradis.xdefaults");


/*
 *      3D-dislocation density data
 */
        BindVar(CPList, "savedensityspec", param->savedensityspec, V_INT, 3,
                VFLAG_NULL);


/*
 *      Miscellaneous parameters
 */
        BindVar(CPList, "Miscellaneous parameters", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "enforceGlidePlanes", &param->enforceGlidePlanes,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "enableCrossSlip", &param->enableCrossSlip,
                V_INT, 1, VFLAG_NULL);
        param->enableCrossSlip = -1;

        BindVar(CPList, "TensionFactor", &param->TensionFactor, V_DBL, 1,
                VFLAG_NULL);
        param->TensionFactor = 1;

        BindVar(CPList, "elasticinteraction", &param->elasticinteraction,
                V_INT, 1, VFLAG_NULL);
        param->elasticinteraction = 1;


#ifdef _FEM
/*
 *      Add in some FEM specific stuff for M. Tang
 */
        BindVar(CPList, "BC_type", &param->BC_type, V_INT, 1, VFLAG_NULL);
        param->BC_type = 5;

        BindVar(CPList, "mesh_type", &param->mesh_type, V_INT, 1, VFLAG_NULL);
        param->mesh_type = 1;

        BindVar(CPList, "dirmax", &param->dirmax, V_INT, 1, VFLAG_NULL);
        param->dirmax = 2000;

        BindVar(CPList, "fem_nx", &param->fem_nx, V_INT, 1, VFLAG_NULL);
        param->fem_nx = 10;

        BindVar(CPList, "fem_ny", &param->fem_ny, V_INT, 1, VFLAG_NULL);
        param->fem_ny = 10;

        BindVar(CPList, "fem_nz", &param->fem_nz, V_INT, 1, VFLAG_NULL);
        param->fem_nz = 10;

        BindVar(CPList, "fem_radius", &param->fem_radius, V_DBL, 1, VFLAG_NULL);
        param->fem_radius = 1750.0;

        BindVar(CPList, "fem_height", &param->fem_height, V_DBL, 1, VFLAG_NULL);
        param->fem_height = 17500.0;

        BindVar(CPList, "fem_nr", &param->fem_nr, V_INT, 1, VFLAG_NULL);
        param->fem_nr = 5;

        BindVar(CPList, "fem_nh", &param->fem_nh, V_INT, 1, VFLAG_NULL);
        param->fem_nh = 10;

        BindVar(CPList, "fem_void_radius", &param->fem_void_radius, V_DBL,
                1, VFLAG_NULL);
        param->fem_void_radius = 500.0;

        BindVar(CPList, "fem_cube_edge_length", &param->fem_cube_edge_length,
                V_DBL, 1, VFLAG_NULL);
        param->fem_cube_edge_length = 1500.0;

        BindVar(CPList, "fem_numelm_arc", &param->fem_numelm_arc, V_INT,
                1, VFLAG_NULL);
        param->fem_numelm_arc = 5;

        BindVar(CPList, "fem_base_diagonal_half",
                &param->fem_base_diagonal_half, V_DBL, 1, VFLAG_NULL);
        param->fem_base_diagonal_half = 500.0;

        BindVar(CPList, "fem_grid_base", &param->fem_grid_base, V_INT,
                1, VFLAG_NULL);
        param->fem_grid_base = 3;

        BindVar(CPList, "fem_grid_vertical", &param->fem_grid_vertical,
                V_INT, 1, VFLAG_NULL);
        param->fem_grid_vertical = 3;

        BindVar(CPList, "fem_ageom_x", param->fem_ageom_x, V_DBL, 3,
                VFLAG_NULL);
        param->fem_ageom_x[0] = 1.0;
        param->fem_ageom_x[1] = 0.0;
        param->fem_ageom_x[2] = 0.0;

        BindVar(CPList, "fem_ageom_z", param->fem_ageom_z, V_DBL, 3,
                VFLAG_NULL);
        param->fem_ageom_z[0] = 0.0;
        param->fem_ageom_z[1] = 0.0;
        param->fem_ageom_z[2] = 1.0;
#endif

#ifdef _MOBILITY_FIELD
	BindVar(CPList, "mobilityField", &param->mobilityField, V_INT, 1, VFLAG_NULL);
        param->mobilityField = 0;
        
        BindVar(CPList, "frictionField", &param->frictionField, V_INT, 1, VFLAG_NULL);
        param->frictionField = 0;
        
        BindVar(CPList, "mobilityFieldFile", param->mobilityFieldFile, V_STRING, 1, VFLAG_NULL);
        strcpy(param->mobilityFieldFile,"mobility_field.dat");
        
        BindVar(CPList, "mobilityFieldScale", &param->mobilityFieldScale, V_DBL, 1, VFLAG_NULL);
        param->mobilityFieldScale = 1.0;
#endif

	BindVar(CPList, "strainTimeData", param->strainTimeData, V_STRING, 1, VFLAG_NULL);
        strcpy(param->strainTimeData,"");

        return;
}


void DataParamInit(Param_t *param, ParamList_t *DPList)
{

/*
 *      Note: Parameters need only be initialized if their
 *      default values are non-zero.
 */
        DPList->paramCnt = 0;
        DPList->varList = (VarData_t *)NULL;

        BindVar(DPList, "dataFileVersion", &param->dataFileVersion, V_INT,
                1, VFLAG_NULL);
        param->dataFileVersion = NODEDATA_FILE_VERSION;

        BindVar(DPList, "numFileSegments", &param->numFileSegments, V_INT,
                1, VFLAG_NULL);
        param->numFileSegments = 1;

        BindVar(DPList, "minCoordinates", param->minCoordinates, V_DBL, 3,
                VFLAG_NULL);

        BindVar(DPList, "maxCoordinates", param->maxCoordinates, V_DBL, 3,
                VFLAG_NULL);

        BindVar(DPList, "nodeCount", &param->nodeCount, V_INT, 1, VFLAG_NULL);

        BindVar(DPList, "dataDecompType", &param->dataDecompType, V_INT, 1,
                VFLAG_NULL);
        param->dataDecompType = 2;

        BindVar(DPList, "dataDecompGeometry", param->dataDecompGeometry,
                V_INT, 3, VFLAG_NULL);
        param->dataDecompGeometry[X] = 1;
        param->dataDecompGeometry[Y] = 1;
        param->dataDecompGeometry[Z] = 1;

        return;
}
