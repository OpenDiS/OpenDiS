/****************************************************************************
 *
 *  Param.h  Define the Param struct, which holds or points to all control
 *           parameters for this simulation
 *
 ***************************************************************************/
#ifndef _PARAM_H
#define _PARAM_H

#include "Parse.h"
#include "Home.h"
#include "Mobility.h"
#include <assert.h>


/*
 *      Define a couple strings related to the nodal data files
 *      written along with the control file containing the global
 *      parameter values.
 */
#define HDF_DATA_FILE_SUFFIX   ".hdf"
#define NODEDATA_FILE_SUFFIX   ".data"
#define NODEDATA_FILE_VERSION  4

struct _param {
/*
 *      Defines the number of domains in each dimension of the
 *      entire problem space.
 */
        int nXdoms;
        int nYdoms;
        int nZdoms;
#ifdef _BGP
        int taskMappingMode;  /* Only valid for BG/P systems on which a  */
                              /* 3D hardware partition is allocated to   */
                              /* the users job.                          */
                              /*                                         */
                              /* Determines action taken when the user-  */
                              /* supplied domain decomposition is not    */
                              /* consistent with the 3D geometry of the  */
                              /* underlying hardware partition.  Allowed */
                              /* values are:                             */
                              /*     0 -- domain decomposition will be   */
                              /*          reset to be consistent with    */
                              /*          the hardware partition (if     */
                              /*          possible)                      */
                              /*     1 -- the user-supplied domain       */
                              /*          decomposition will be used     */
                              /*     2 -- the code will abort if there   */
                              /*          is a mismatch.                 */

#endif

/*
 *      Defines number of cells in each of the dimensions of the
 *      entire problem space.
 */
        int nXcells;
        int nYcells;
        int nZcells;

/*
 *      "natural" min and max cell indices for this domain. (i.e. indices
 *      before being incremented by 1 to allow for ghost cells)
 */
        int iCellNatMin;
        int iCellNatMax;
        int jCellNatMin;
        int jCellNatMax;
        int kCellNatMin;
        int kCellNatMax;

/*
 *      Specifies type of boundary in each dimension:
 *          0 == periodic
 *          1 == free surface
 */
        BoundType_t xBoundType;
        BoundType_t yBoundType;
        BoundType_t zBoundType;

/*
 *      When free surfaces are used, these define the upper and
 *      lower limits on dislocation coordinates.  (i.e. the planes
 *      defining the free surfaces)
 */
        real8 xBoundMin, xBoundMax;
        real8 yBoundMin, yBoundMax;
        real8 zBoundMin, zBoundMax;



/*
 *      Define min/max coordinate limits for each dimension of
 *      the problem space.  NOTE: These are redundant but kept
 *      until all code references to them have been removed.
 *      Use <minCoordinates> and <maxCoordinates> arrays instead.
 */
        real8 minSideX, maxSideX;
        real8 minSideY, maxSideY;
        real8 minSideZ, maxSideZ;

/*
 *      Domain decomposition and rebalance values
 */
        int   decompType;       /* Selects decomposition type */
        int   DLBfreq;          /* how often to load balance */
        int   numDLBCycles;     /* Number of initial load-balance-only */
                                /* cycles to be executed before main   */
                                /* processing loop is entered.  This   */
                                /* is a command line option value, not */
                                /* a control file parameter            */

/*
 *      Simulation time and timestepping controls
 */
        int   cycleStart;    /* Starting cycle number for the simulation */
        int   maxstep;       /* Cycles to execute before terminating */
        real8 timeStart;     /* Initial simulation time */
        real8 timeNow;       /* current simulation time */

        char  timestepIntegrator[MAX_STRING_LEN];
#ifdef _SUBCYCLING
        char  subInteg0Integ1[MAX_STRING_LEN];
#endif

        real8 deltaTT;       /* duration of previous timestep */
        real8 realdt;        /* *almost* obsolete */
        real8 nextDT;        /* timestep to attempt the next cycle */
        real8 maxDT;         /* Maximum timestep duration permitted*/
        real8 dtIncrementFact;  /* Maximum delta time increment allowed */
        real8 dtDecrementFact;  /* Factor by which delta time is        */
                                /* multiplied when cutting timestep.    */
                                /* value should be between 0 and 1      */
        real8 dtExponent;       /* Used in calculating variable delta   */
                                /* time adjustments                     */
        int   dtVariableAdjustment;  /* Flag indicating if timestep     */
                                     /* adjustments are of variable     */
                                     /* size or set percentages of the  */
                                     /* current dt.                     */
        real8 rTol;      /* Maximum error allowed in timestep */
        real8 rmax;      /* maximum migration distance per timestep */
                         /* for any node   */
#ifdef _SUBCYCLING
		real8 rTolth;           /* Threshold error tolerance, R. B. Sills 3-2-12*/
		real8 rTolrel;          /* Relative error tolerance, R. B. Sills 3-2-12*/
        real8 deltaTTsub;       /* duration of previous subcycle timestep */
        real8 realdtsub;        /* *almost* obsolete */
        real8 nextDTsub;        /* timestep to attempt the next subcycle step */
		real8 deltaTTsub2;      /* duration of previous subcycle timestep */
        real8 realdtsub2;       /* *almost* obsolete */
        real8 nextDTsub2;       /* timestep to attempt the next subcycle step */
		real8 deltaTTsub3;      /* duration of previous subcycle timestep */
        real8 realdtsub3;       /* *almost* obsolete */
        real8 nextDTsub3;       /* timestep to attempt the next subcycle step */
		real8 deltaTTsub4;      /* duration of previous subcycle timestep */
        real8 realdtsub4;       /* *almost* obsolete */
        real8 nextDTsub4;       /* timestep to attempt the next subcycle step */
		real8 renh;             /* elastic interaction enhancement radius for Jacobian*/
		real8 rg1;              /* subcycling grouping radius (first  layer) */
		real8 rg2;              /* subcycling grouping radius (second layer) */
		real8 rg3;              /* subcycling grouping radius (second layer) */
		real8 rg4;              /* Enlarged grouping radius*/
		int   nTry;             /* No of times we force the integrator to converge
		                           without cutting the global time-step-size by removing
								   some interactions from the global group. This feature
								   only works with force-based subcycling */
		int   sendSubGroupForc; /* Sending forces during subcycling steps for parallel
		                           processing. If it is equal to one forces will be sent,
								   otherwise each domain will calculate all necessary
								   forces during subcycling steps. */
#endif

		int   inSubcycling;      /* Flag indicating whether we are within the subcycling integration */
/*
 *      Discretization parameters and controls for topological changes
 */
        real8 minSeg;    /* min allowable segment length, before */
                         /* removing a node */
        real8 maxSeg;    /* max allowable segment length, before*/
                         /* adding a node*/
        int remeshRule;
        int collisionMethod;
        real8 remeshAreaMax; /* This is calculated from remeshAreaMin and  */
                             /* remeshAreaRatio and hence not user-provided*/
        real8 remeshAreaMin; /* This values is based on the minSeg value  */
                             /* and hence not specified by the user.      */
        int splitMultiNodeFreq;  /* Code will attempt to split multi-nodes */
                                 /* every cycle that is a multiple of this */
                                 /* value. */

/*
 *      Fast Multipole Method parameters
 */
        int fmEnabled;    /* Set to 1 remote forces are caculated */
                          /* using the Fast Multipole code        */

        int fmNumLayers;  /* Number of layers of cells defined for*/
                          /* the Fast Multipole remote force calcs*/

        int fmMPOrder;    /* Order used for multipole expansions  */

        int fmTaylorOrder;/* Order used for taylor expansions     */

        int fmNumPoints;  /* Number of points along each segment    */
                          /* at which remote stress will be         */
                          /* evaluated from the taylor expansion    */
                          /* coefficients. Automatically calculated */
                          /* not supplied by user.                  */

        char fmCorrectionTbl[MAX_STRING_LEN];

#ifdef _GPU_SUBCYCLE
		int *nFMCellLayer;
		int *nMkTaylorLayer;
		int **MkTaylorCellID;
		int **MkTaylorNCellID;
		real8 **MkTaylorR;
#endif

/*
 *      Names of tables for non-FMM far-field forces
 */
        char Rijmfile[MAX_STRING_LEN];
        char RijmPBCfile[MAX_STRING_LEN];

/*
 *      Loading condition parameters
 */
        real8 TempK;         /* Temperature in deg K */
        int   loadType;      /* 0 Creep test */
                             /* 1 Constant strain rate test */
                             /* 2 Displacement-controlled test */
                             /* 3 Junction unzipping jump test */
                             /* 4 Total strain controlled cyclic load test*/
                             /* 5  Plastic strain controlled cyclic load test*/
                             /* 6  Load-time curve test*/
                             /* 7  Constant stress rate test */

        real8 appliedStress[6]; /* External stress in units of Pa  */
                                /* as [sigma11, sigma22, sigma33,  */
                                /* sigma23, sigma31, sigma12] when */
                                /*  <loadType> == 0.               */

        real8 appliedStressRate[6]; /* External stress rate in units of Pa/s      */
                                    /* as [sigmaRate11, sigmaRate22, sigmaRate33, */
                                    /* sigmaRate23, sigmaRate31, sigmaRate12].    */

        real8 eRate;         /* Strain rate. Used when loadType == 1 */
        int   indxErate;     /* For normal loading, indxErate ==1 */
                             /* For shear loading, indxErate == 4 */ 
        real8 edotdir[3];    /* Uniaxial loading direction accompanying */
                             /* eRate                                   */

        real8 sRate;         /* Stress rate. Used when loadType == 7. Unit: Pa/s */              

        real8 cTimeOld;      /* Timestep related to cyclic loading */
        real8 netCyclicStrain; /* Net accumulated strain under cyclic load */
        real8 dCyclicStrain; /* Incremental strain under cyclic load */
        int   numLoadCycle;  /* Number of cyclic cycles */
        real8 eAmp;          /* Strain amplitude used with cyclic loading */

/*
 *      Values for specifying axes for a user-defined laboratory frame
 */
        int   useLabFrame;   /* 0 if standard crystalographic frame is to */
                             /* be used, 1 if user-supplied laboratory    */
                             /* frame is used                             */

        real8 labFrameXDir[3];  /* The Z direction is informational only */
        real8 labFrameYDir[3];  /* and is recalculated explicitly from   */
        real8 labFrameZDir[3];  /* the X and Y directions                */

/*
 *      Material and mobility parameters
 */
#if 0
        real8 meltTemp;      /* Material melting temperature in degrees K */
#endif

        char  mobilityLaw[MAX_STRING_LEN];
        int   mobilityType;  /* Integer value corresponding to the */
                             /* specified mobility law.  Redundant */
                             /* info, but easier to use in the code*/

        int   materialType;  /* Type of crystal structure (i.e. BCC, FCC)   */
                             /* This value is set within the code base on   */
                             /* the selected mobility law and is not user-  */
                             /* supplied                                    */

        real8 vacancyConc;            /* Concentration of vacancies in the */
                                      /* crystal: for calculating osmotic  */
                                      /* force (units ?)                   */

        real8 vacancyConcEquilibrium; /* Thermal equilibrium vacacy concen-*/
                                      /* ration in defect free crystal: for*/
                                      /* calculating osmotic force         */

        real8 shearModulus;
        real8 pois;
        real8 burgMag;

        real8 YoungsModulus;

        real8 rc;     /* core radius in elastic interaction calculation */

        real8 Ecore;  /* core energy (wrt. the choice of rc) in unit of Pa */

        int   enforceGlidePlanes;
        int   allowFuzzyGlidePlanes;
        int   enableCrossSlip;

        int   (*mobilityFunc)(Home_t *home, Node_t *node);  /* Set during */
                                   /* initialization to point to the      */
                                   /* appropriate mobility function       */
#ifdef _GPU_SUBCYCLE
        int   (*mobilityMatrixGPU)(Home_t *home, Node_t *node, double mobMatrix[3][3]);
#endif

        real8 MobScrew;
        real8 MobEdge;
        real8 MobClimb;
        real8 MobGlide;  /* floor on mobility for glide dislocations */
        real8 MobLine;
        real8 FricStress; /* Lattice friction stress */

/*
 *      Allow for specifying that dislocations with certain types
 *      of burgers vectors or line directions are immobile.
 */
        real8 sessileburgspec[30];
        real8 sessilelinespec[30];

/*
 *      Add some variables needed to include inertial terms into
 *      mobility functions.  If the selected mobility does not
 *      have support for inertial effects, these values will be
 *      quietly ignored.
 */
        int   includeInertia;  /* Toggle: 0 == no, 1 == yes */
        real8 massDensity;     /* Units are kg/m^3 */

//#ifdef _THERMAL_ACTIVATED_CROSSSLIP
/*
 *      Enables thermally-activated cross-slip.
 */
        int thermalCrossSlip;
		real8 CRS_A;
		real8 CRS_T0;
		real8 pre_Cge;
//#endif
		int CRScount_zipper;
		int CRScount_bothScrew;

/*
 *      Velocity statistics and parameters
 */
        real8 vAverage;        /* average nodal velocity */
        real8 vStDev;          /* St.Dev of nodal velocity */

/*
 *      I/O parameters
 */
        char dirname[MAX_STRING_LEN];

        int  writeBinRestart; /* if set, will write data portion of */
                              /* restart file in a binary format    */

        int  doBinRead;  /* If set, will attempt to read binary format */
                         /* restart file.  This flag is set internally */
                         /* and not specified by the user.             */

        int  numIOGroups;  /* number of groups into which to split tasks */
                           /* when doing parallel I/O                    */

        int  skipIO;    /* if set, all major IO is skipped regardless */
                        /* of what types of output generation would   */
                        /* be enabled by other control parameters     */

        int   armfile, armfilefreq, armfilecounter;
        real8 armfiledt, armfiletime;

        int   fluxfile, fluxfreq, fluxcounter;
        real8 fluxdt, fluxtime;

        int   fragfile, fragfreq, fragcounter;
        real8 fragdt, fragtime;

        int   gnuplot, gnuplotfreq, gnuplotcounter;
        real8 gnuplotdt, gnuplottime;

        int   polefigfile, polefigfreq, polefigcounter;
        real8 polefigdt, polefigtime;

        int   povray,  povrayfreq,  povraycounter;
        real8 povraydt, povraytime;

        int   atomeye,  atomeyefreq,  atomeyecounter;
        real8 atomeyedt, atomeyetime, atomeyesegradius;

        int   psfile,  psfilefreq;
        real8 psfiledt, psfiletime;

        int   savecn,  savecnfreq,  savecncounter;
        real8 savecndt, savecntime;

        int   saveprop, savepropfreq;
        real8 savepropdt, saveproptime;

        int   savetimers, savetimersfreq, savetimerscounter;
        real8 savetimersdt, savetimerstime;

        int   savedensityspec[3];

        int   tecplot, tecplotfreq, tecplotcounter;
        real8 tecplotdt, tecplottime;

        int   paraview, paraviewfreq, paraviewcounter;
        real8 paraviewdt, paraviewtime;

        int   velfile, velfilefreq, velfilecounter;
        real8 velfiledt, velfiletime;

        int   writeForce, writeForceFreq, writeForceCounter;
        real8 writeForceDT, writeForceTime;
        
        int   linkfile, linkfilefreq, linkfilecounter;
        real8 linkfiledt, linkfiletime;

#ifdef _OP_REC
        int   oprec, oprecmove;
        int   oprecwritefreq, oprecfilefreq, opreccounter;
        int   runfromoprec;
        char  oprecfile[MAX_STRING_LEN];
#endif

        int   writeVisit, writeVisitFreq, writeVisitCounter;
        int   writeVisitSegments, writeVisitSegmentsAsText;
        int   writeVisitNodes, writeVisitNodesAsText;
        real8 writeVisitDT, writeVisitTime;

        char  winDefaultsFile[MAX_STRING_LEN];



/*
 *      Lengths (and reciprocals) of each side of the problem
 *      space box.
 */
        real8 Lx, Ly, Lz;
        real8 invLx, invLy, invLz;


/*
 *      General stuff
 */
        real8 springConst;
        real8 rann;   /* closest distance before dislocations are */
                      /* considered in contact    */


        int  numBurgGroups;     /* Number of groups into which different  */
                                /* burgers vectors are organized in order */
                                /* to track dislocation density by burgers*/
                                /* vector.  This number is dependent on   */
                                /* the type of mobility used.             */

        real8 *partialDisloDensity;  /* Dynamically sized array to hold   */
                                     /* dislocation density for distinct  */
                                     /* burgers vector groupings.  This   */
                                     /* array will be dynamically al-     */
                                     /* located to the appropriate size   */
                                     /* depending on the type of mobility */
                                     /* being used.                       */
        real8 disloDensity;

        real8 delSegLength;  /* accumulated length of deleted segments    */
                             /* since most recent write to result_density */
                             /* property file                             */

        real8 densityChange[14];  /* For tracking density change by burgers */
                                  /* vector; currently only used for BCC.   */
                                  /* Values accumulated since most recent   */
                                  /* write of density change data to the    */
                                  /* property file.                         */

        real8 TensionFactor;
        int elasticinteraction;

        real8 delpStrain[6],delSig[6],totpStn[6];
        real8 delpSpin[6],totpSpn[6];
#ifdef _FIX_PLASTIC_STRAIN
        real8 delpStrainCorr[6], delpSpinCorr[6];
#endif

/*
 *      Added for strain decomposition and density flux decomp.
 */
        real8 totstraintensor[6];
        real8 totedgepStrain[6], totscrewpStrain[6];
        real8 dedgepStrain[6],   dscrewpStrain[6];

        real8 Ltot[4][4],  fluxtot[4][7];
        real8 dLtot[4][4], dfluxtot[4][7];

/*
 *      For strain decomposition and density flux decomp in FCC
 */
        real8 FCC_Ltot[6][4],  FCC_fluxtot[6][7];
        real8 FCC_dLtot[6][4], FCC_dfluxtot[6][7];

        int imgstrgrid[6];

/*
 *      The following two variables are used only temporarily and as such
 *      should not be specifiable in the control file, and hence, should
 *      not be bound to identifying strings via bindvar() as the other
 *      elements of this structure are.
 */
        char node_data_file[MAX_STRING_LEN];
                                  /* name of the file containing the    */
                                  /* nodal data for the run.  This data */
                                  /* is in the newer nodal data format, */
                                  /* (i.e. not contained as part of the */
                                  /* control file itself)               */

/*
 *      Define the parameters used in the nodal data file
 */
        int   dataFileVersion;       /* Version number of the data file */
        int   numFileSegments;       /* Number of files the nodal data  */
                                     /* is spread over.                 */
        int   nodeCount;             /* Total number of nodes in the    */
                                     /* data file (all file segments)   */
        int   dataDecompType;        /* Type of decomposition used when */
                                     /* generating the current data file*/
        int   dataDecompGeometry[3]; /* domain geometry used when gen-  */
                                     /* erating the current data file   */
        real8 minCoordinates[3];     /* minimum XYZ coordinates of sim  */
        real8 maxCoordinates[3];     /* maximum XYZ coordinates of sim  */

/*
 *      Define a couple factors used when calculating dislocation density.
 *      These will be calculated during initialization and will not be
 *      specified in the control parameter file.
 */
        real8 simVol;         /* Total volume of simulation */
        real8 burgVolFactor;  /* Volume factor used to convert dislocation */
                              /* length to dislocation density             */

/*
 *      Define the number of threads to use per MPI task. Ignored if
 *      threading has not been enabled.
 */
        int   maxNumThreads;

#ifdef _FEM
/*
 *      Include some FEM specific stuff for Meijie.  This is not the
 *      best place to add this in, but is the easiest way to add it
 *      quickly.
 */
        int mesh_type;          /* 1: rectangular prism                    */
                                /* 2: cylinder within cubic simulation box */
                                /* 3: void within cubic simulation box     */
                                /* 3: octahedron within cubic simulation box */

        int BC_type;            /* 1: only one free surface at z=zBoundMax */
                                /* 4: top and bottom along z are free      */
                                /*    surfaces.                            */
                                /* 5: all 6 surfaces are free              */

        int dirmax;             /* Iterative solver used if number of fem */
                                /* elements is larger than this value,    */
                                /* direct solver used otherwise.          */

/*
 *      Specify number of fem elements in the rectangular prism defined
 *      within the simulation (i.e.  mesh_type == 1)
 */
        int fem_nx;             /* number of fem elements along x, y and */
        int fem_ny;             /* z respectively                        */
        int fem_nz;

/*
 *      Some parameters needed when defining a cylindrical shape within
 *      the cubic simulation box (i.e mesh_type == 2)
 */
        real8 fem_radius;       /* radius of the cylinder (units of b)  */
        real8 fem_height;       /* height of the cylinder (units of b)  */
        int   fem_nr;           /* # of fem elements along circumference*/
                                /*   of the cylinder                    */
        int   fem_nh;           /* # of fem elements along height of    */
                                /*   the cylinder                       */

        real8 fem_void_radius;     /* radius of void */
        real8 fem_cube_edge_length;/* void???*/
        int   fem_numelm_arc;      /* void??*/

        real8 fem_base_diagonal_half; /* half diagonal distance of the base */
        int   fem_grid_base;          /* number of fem mesh nodes along     */
                                      /*   base cubic sides                 */
        int   fem_grid_vertical;      /* number of fem mesh nodes along     */
                                      /*   vertical distance between two    */
                                      /*   opposite vertices                */

        real8 fem_ageom_x[3];   /* vector defining X axis of cylinder   */
        real8 fem_ageom_z[3];   /* vector defining Z axis of cylinder   */
/*
 *
 */
        real8 fem_delSegLength; /* Used for tracking dislocation length  */
                                /* lost as dislocations move outside the */
                                /* simulation bounds                     */
#endif

#ifdef _MOBILITY_FIELD
		int     mobilityField;
		int     frictionField;
		char    mobilityFieldFile[MAX_STRING_LEN];
		double  mobilityFieldScale;
		int     mobilityFieldNx;
		int     mobilityFieldNy;
		int     mobilityFieldNz;
		double *mobilityFieldPhi;
#endif

		char strainTimeData[MAX_STRING_LEN];
};


#ifdef __cplusplus
extern "C" void CtrlParamInit(Param_t *param, ParamList_t *CPList);
extern "C" void DataParamInit(Param_t *param, ParamList_t *DPList);
#else
void CtrlParamInit(Param_t *param, ParamList_t *CPList);
void DataParamInit(Param_t *param, ParamList_t *DPList);
#endif

#endif /* _PARAM_H */
