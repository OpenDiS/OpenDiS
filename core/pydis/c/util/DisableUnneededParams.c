/****************************************************************************
 *
 *      Module:       DisableUnneededParams.c
 *
 *      Description:  This function is called by task 0 during
 *                    initialization to explicitly mark certain
 *                    control parameters in order to prevent them
 *                    from being written out when creating restart
 *                    files.  This is done because not all the known
 *                    control parameters are applicable to every
 *                    simulation.  For instance, the set of mobility
 *                    parameters used during a given simulation
 *                    is dependent on the mobility law in use, and
 *                    will only be a subset of all the available
 *                    mobility parameters.
 * 
 ****************************************************************************/
#include "Home.h"

/*
 *      Need only be called by task zero
 */
void DisableUnneededParams(Home_t *home)
{
        Param_t  *param;

        param = home->param;

/*
 *      Not all mobility parameters are needed by all mobility functions.
 *      The easiest thing to do is first disable all of them and then only
 *      re-enable the ones appropriate to the selected mobility.
 */
        MarkParamDisabled(home->ctrlParamList, "MobScrew");
        MarkParamDisabled(home->ctrlParamList, "MobEdge");
        MarkParamDisabled(home->ctrlParamList, "MobClimb");
        MarkParamDisabled(home->ctrlParamList, "TempK");

        switch(param->mobilityType) {
            case MOB_BCC_0:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_BCC_0B:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_BCC_GLIDE:
            case MOB_BCC_GLIDE_0:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                break;
            case MOB_FCC_0:
            case MOB_FCC_0B:
            case MOB_FCC_CLIMB:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                break;
        }  /* end switch(mobilityType) */

/*
 *      If no free surfaces are used, we can drop the free surface
 *      boundary values.
 */
        if ((param->xBoundType == Periodic) &&
            (param->xBoundType == Periodic) &&
            (param->xBoundType == Periodic)) {
            MarkParamDisabled(home->ctrlParamList, "xBoundMin");
            MarkParamDisabled(home->ctrlParamList, "yBoundMin");
            MarkParamDisabled(home->ctrlParamList, "zBoundMin");
            MarkParamDisabled(home->ctrlParamList, "xBoundMax");
            MarkParamDisabled(home->ctrlParamList, "yBoundMax");
            MarkParamDisabled(home->ctrlParamList, "zBoundMax");
        }

/*
 *      If the dynamic load balancing is disabled, turn off the decomposition
 *      type.
 */
        if (param->DLBfreq < 1) {
            MarkParamDisabled(home->ctrlParamList, "decompType");
        }

/*
 *      If osmotic forces are not enabled, disable/enable any params 
 *      as needed.
 *
 *      Note: TempK was disabled along with mobility params above but
 *            it is needed to do osmotic forces
 */
        if (param->vacancyConcEquilibrium <= 0) {
            MarkParamDisabled(home->ctrlParamList, "vacancyConc");
            MarkParamDisabled(home->ctrlParamList, "vacancyConcEquilibrium");
        } else {
            MarkParamEnabled(home->ctrlParamList, "TempK");
        }

/*
 *      If FMM is not enabled, disable the corresponding parameters
 *      otherwise diable the parameters related to the Rijm* files
 */
        if (param->fmEnabled == 0) {
            MarkParamDisabled(home->ctrlParamList, "fmMPOrder");
            MarkParamDisabled(home->ctrlParamList, "fmTaylorOrder");
            MarkParamDisabled(home->ctrlParamList, "fmCorrectionTbl");
        } else {
            MarkParamDisabled(home->ctrlParamList, "Rijmfile");
            MarkParamDisabled(home->ctrlParamList, "RijmPBCfile");
        }

/*
 *      Not all of the timestep integrator parameters are used in
 *      all the integrator functions...
 */
        if (strcmp(param->timestepIntegrator, "forward-euler") == 0) {
            MarkParamDisabled(home->ctrlParamList, "dtDecrementFact");
            MarkParamDisabled(home->ctrlParamList, "dtExponent");
            MarkParamDisabled(home->ctrlParamList, "dtIncrementFact");
            MarkParamDisabled(home->ctrlParamList, "dtVariableAdjustment");
        } else {

            MarkParamDisabled(home->ctrlParamList, "rmax");
        }

/*
 *      The selected <loadType> affects which parameters are used.
 *      do that setup now.  As with the mobility parameters above,
 *      for some of the parameters, it's easier to disable a group
 *      of them and only re-enable them as needed.
 */
        MarkParamDisabled(home->ctrlParamList, "cTimeOld");
        MarkParamDisabled(home->ctrlParamList, "dCyclicStrain");
        MarkParamDisabled(home->ctrlParamList, "netCyclicStrain");
        MarkParamDisabled(home->ctrlParamList, "numLoadCycle");
        MarkParamDisabled(home->ctrlParamList, "eAmp");

        switch(param->loadType) {
            case 0:
                MarkParamDisabled(home->ctrlParamList, "indxErate");
                break;
            case 1:
                break;
            case 2:
                break;
            case 3:
                break;
/*
 *          For loadType's 4 and 5, we use re-enabled the same set
 *          of parameters
 */
            case 4:
            case 5:
                MarkParamEnabled(home->ctrlParamList, "cTimeOld");
                MarkParamEnabled(home->ctrlParamList, "dCyclicStrain");
                MarkParamEnabled(home->ctrlParamList, "netCyclicStrain");
                MarkParamEnabled(home->ctrlParamList, "numLoadCycle");
                MarkParamEnabled(home->ctrlParamList, "eAmp");
                break;
            case 6:
                break;
            case 7:
                break;
        }

/*
 *      If inertial terms are not being included, don't write parameters
 *      specific to that capability.
 */
        if (param->includeInertia == 0) {
            MarkParamDisabled(home->ctrlParamList, "massDensity");
        }

/*
 *      If no sessile burgers vectors have been specified, don't write
 *      the related arrays out.
 */
        if (param->sessileburgspec[0] <= 0) {
            MarkParamDisabled(home->ctrlParamList, "sessileburgspec");
            MarkParamDisabled(home->ctrlParamList, "sessilelinespec");
        }

/*
 *      If the geometries are not specified in a user-defined laboratory
 *      frame but use the standard crystalographic frame, don't write
 *      out any axes for a laboratory frame.
 */
        if (param->useLabFrame == 0) {
            MarkParamDisabled(home->ctrlParamList, "labFrameXDir");
            MarkParamDisabled(home->ctrlParamList, "labFrameYDir");
            MarkParamDisabled(home->ctrlParamList, "labFrameZDir");
        }

/*
 *      The flux data arrays always get dumped to the control file, but
 *      which arrays are written depends on the type of material, so
 *      only enable the ones we need.
 */
        MarkParamDisabled(home->ctrlParamList, "Ltot");
        MarkParamDisabled(home->ctrlParamList, "fluxtot");
        MarkParamDisabled(home->ctrlParamList, "FCC_Ltot");
        MarkParamDisabled(home->ctrlParamList, "FCC_fluxtot");

        switch (param->materialType) {
            case MAT_TYPE_BCC:
                MarkParamEnabled(home->ctrlParamList, "Ltot");
                MarkParamEnabled(home->ctrlParamList, "fluxtot");
                break;
            case MAT_TYPE_FCC:
                MarkParamEnabled(home->ctrlParamList, "FCC_Ltot");
                MarkParamEnabled(home->ctrlParamList, "FCC_fluxtot");
                break;
        }

/*
 *      Disable writing of the output-related parameters that do not
 *      apply to the types of output actually selected.
 */
        if (param->armfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "armfilefreq");
            MarkParamDisabled(home->ctrlParamList, "armfiledt");
            MarkParamDisabled(home->ctrlParamList, "armfiletime");
            MarkParamDisabled(home->ctrlParamList, "armfilecounter");
        };

        if (param->armfiledt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "armfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "armfiledt");
            MarkParamDisabled(home->ctrlParamList, "armfiletime");
        }


        if (param->fluxfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "fluxfreq");
            MarkParamDisabled(home->ctrlParamList, "fluxdt");
            MarkParamDisabled(home->ctrlParamList, "fluxtime");
            MarkParamDisabled(home->ctrlParamList, "fluxcounter");
        }
        

        if (param->armfiledt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "fluxfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "fluxdt");
            MarkParamDisabled(home->ctrlParamList, "fluxtime");
        }


        if (param->fragfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "fragfreq");
            MarkParamDisabled(home->ctrlParamList, "fragdt");
            MarkParamDisabled(home->ctrlParamList, "fragtime");
            MarkParamDisabled(home->ctrlParamList, "fragcounter");
        }

        if (param->fragdt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "fragfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "fragdt");
            MarkParamDisabled(home->ctrlParamList, "fragtime");
        }


        if (param->gnuplot == 0) {
            MarkParamDisabled(home->ctrlParamList, "gnuplotfreq");
            MarkParamDisabled(home->ctrlParamList, "gnuplotdt");
            MarkParamDisabled(home->ctrlParamList, "gnuplottime");
            MarkParamDisabled(home->ctrlParamList, "gnuplotcounter");
        }

        if (param->gnuplotdt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "gnuplotfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "gnuplotdt");
            MarkParamDisabled(home->ctrlParamList, "gnuplottime");
        }


        if (param->polefigfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "polefigfreq");
            MarkParamDisabled(home->ctrlParamList, "polefigdt");
            MarkParamDisabled(home->ctrlParamList, "polefigtime");
            MarkParamDisabled(home->ctrlParamList, "polefilecounter");
        }

        if (param->polefigdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "polefigfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "polefigdt");
            MarkParamDisabled(home->ctrlParamList, "polefigtime");
        }


        if (param->povray == 0) {
            MarkParamDisabled(home->ctrlParamList, "povrayfreq");
            MarkParamDisabled(home->ctrlParamList, "povraydt");
            MarkParamDisabled(home->ctrlParamList, "povraytime");
            MarkParamDisabled(home->ctrlParamList, "povraycounter");
        }

        if (param->povraydt > 0) {
            MarkParamDisabled(home->ctrlParamList, "povrayfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "povraydt");
            MarkParamDisabled(home->ctrlParamList, "povraytime");
        }


        if (param->psfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "psfilefreq");
            MarkParamDisabled(home->ctrlParamList, "psfiledt");
            MarkParamDisabled(home->ctrlParamList, "psfiletime");
        }

        if (param->psfile > 0) {
            MarkParamDisabled(home->ctrlParamList, "psfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "psfiledt");
            MarkParamDisabled(home->ctrlParamList, "psfiletime");
        }


        if (param->savecn == 0) {
            MarkParamDisabled(home->ctrlParamList, "savecnfreq");
            MarkParamDisabled(home->ctrlParamList, "savecndt");
            MarkParamDisabled(home->ctrlParamList, "savecntime");
            MarkParamDisabled(home->ctrlParamList, "savecncounter");
        }

        if (param->savecndt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savecnfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savecndt");
            MarkParamDisabled(home->ctrlParamList, "savecntime");
        }


        if (param->saveprop == 0) {
            MarkParamDisabled(home->ctrlParamList, "savepropfreq");
            MarkParamDisabled(home->ctrlParamList, "savepropdt");
            MarkParamDisabled(home->ctrlParamList, "saveproptime");
        }

        if (param->savepropdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savepropfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savepropdt");
            MarkParamDisabled(home->ctrlParamList, "saveproptime");
        }


        if (param->savetimers == 0) {
            MarkParamDisabled(home->ctrlParamList, "savetimersfreq");
            MarkParamDisabled(home->ctrlParamList, "savetimersdt");
            MarkParamDisabled(home->ctrlParamList, "savetimerstime");
            MarkParamDisabled(home->ctrlParamList, "savetimerscounter");
        }

        if (param->savetimersdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savetimersfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savetimersdt");
            MarkParamDisabled(home->ctrlParamList, "savetimerstime");
        }


        if (param->tecplot == 0) {
            MarkParamDisabled(home->ctrlParamList, "tecplotfreq");
            MarkParamDisabled(home->ctrlParamList, "tecplotdt");
            MarkParamDisabled(home->ctrlParamList, "tecplottime");
            MarkParamDisabled(home->ctrlParamList, "tecplotcounter");
        }

        if (param->tecplotdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "tecplotfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "tecplotdt");
            MarkParamDisabled(home->ctrlParamList, "tecplottime");
        }


        if (param->velfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "velfilefreq");
            MarkParamDisabled(home->ctrlParamList, "velfiledt");
            MarkParamDisabled(home->ctrlParamList, "velfiletime");
            MarkParamDisabled(home->ctrlParamList, "velfilecounter");
        }

        if (param->velfiledt > 0) {
            MarkParamDisabled(home->ctrlParamList, "velfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "velfiledt");
            MarkParamDisabled(home->ctrlParamList, "velfiletime");
        }


        if (param->writeForce == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeForceFreq");
            MarkParamDisabled(home->ctrlParamList, "writeForceDT");
            MarkParamDisabled(home->ctrlParamList, "writeForceTime");
            MarkParamDisabled(home->ctrlParamList, "writeForceCounter");
        }

        if (param->writeForceDT > 0) {
            MarkParamDisabled(home->ctrlParamList, "writeForceFreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "writeForceDT");
            MarkParamDisabled(home->ctrlParamList, "writeForceTime");
        }


        if ((param->savedensityspec[0] == 0) ||
            (param->savedensityspec[1] == 0) ||
            (param->savedensityspec[2] == 0)) {
            MarkParamDisabled(home->ctrlParamList, "savedensityspec");
        }

        if (param->writeVisit == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitFreq");
            MarkParamDisabled(home->ctrlParamList, "writeVisitDT");
            MarkParamDisabled(home->ctrlParamList, "writeVisitTime");
            MarkParamDisabled(home->ctrlParamList, "writeVisitCounter");
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodes");
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodesAsText");
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegments");
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegmentsAsText");
        }

        if (param->writeVisitDT > 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitFreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "writeVisitDT");
            MarkParamDisabled(home->ctrlParamList, "writeVisitTime");
        }

        if (param->writeVisitNodes == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodesAsText");
        }

        if (param->writeVisitSegments == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegmentsAsText");
        }

        return;
}
