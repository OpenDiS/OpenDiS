#ifndef _CONSTANTS_H
#define _CONSTANTS_H
/****************************************************************************
 *
 *      Module:      Constants.h  
 *      Description: Contains definitions for some constants that are used
 *                   by other include files.  
 *
 ***************************************************************************/

/*
 *      Define the numbers of nodes allocated when obtaining new node blocks
 *      within portions of the code.
 */
#define NEW_NODEKEY_INC     1000
#define RECYC_NODESTACK_INC 100
#define NODE_BLOCK_COUNT    500

/*
 *      Some definitions used to select whether to do a full or
 *      partial set of force calculations
 */
#define PARTIAL   1
#define FULL      2
#ifdef _SUBCYCLING
#define GROUP0    10
#define GROUP1    11
#define GROUP2    12
#define GROUP3    13
#define GROUP4    14
#endif

/*
 *      Define the subdirectories under which various types of
 *      output files will be created.  These directories are
 *      relative to the user-specified output directory in the
 *      control file.
 */
#define DIR_ARMDATA    "armdata"
#define DIR_FLUXDATA   "fluxdata"
#define DIR_FORCE      "force"
#define DIR_FRAGDATA   "visit"
#define DIR_GNUPLOT    "gnuplot"
#define DIR_POLEFIG    "polefig"
#define DIR_POVRAY     "povray"
#define DIR_ATOMEYE    "atomeye"
#define DIR_PROPERTIES "properties"
#define DIR_RESTART    "restart"
#define DIR_LINKDATA   "links"
#define DIR_TECPLOT    "tecplot"
#define DIR_PARAVIEW   "paraview"
#define DIR_TIMERS     "timers"
#define DIR_VELOCITY   "velocity"
#define DIR_VISIT      "visit"

#ifdef _OP_REC
#define DIR_OPREC      "oprec"
#endif

/*
 *      Some miscellaneous definitions used throughout the code
 */
#define MAX_NBRS            10 /* maximum number of segments for a node */
#define BOLTZMANNS_CONST    1.38e-23  /* joules/K */
#define FFACTOR_ORTH        0.05
#define FFACTOR_NORMAL      1.0e-4
#define FFACTOR_LMIN        1e-8  /* minimum allowed segment length */
#define FFACTOR_LMIN2       1e-16 /* minimum allowed segment length squared */

#define MAX_STRING_LEN      256
#define PRECISEGLIDEPLANE_THERESHOLD   0.999	/*dot product between cross(b, ksi) and nPlane from the table must exceed this number in order to copy the nPlane */

/*
 *      Define some values that can be used by GenerateOutput()
 *      to indicate the current stage of the execution.  Needed
 *      because output generation differs depending on the point
 *      the program is at.
 */
#define STAGE_INIT      1
#define STAGE_CYCLE     2
#define STAGE_TERM      3

#define FIRST_BLOCK     1
#define LAST_BLOCK      2

/*
 *      Define a set of flags corresponding to the types of
 *      output that may be generated.  Note: these are intended
 *      as bit values that can be or'ed together.
 */
#define GEN_ARM_DATA            0x00001
#define GEN_DENSITY_DATA        0x00002
#define GEN_FLUX_DATA           0x00004
#define GEN_GNUPLOT_DATA        0x00008
#define GEN_POLEFIG_DATA        0x00010
#define GEN_POVRAY_DATA         0x00020
#define GEN_POSTSCRIPT_DATA     0x00040
#define GEN_PROPERTIES_DATA     0x00080
#define GEN_RESTART_DATA        0x00100
#define GEN_TECPLOT_DATA        0x00200
#define GEN_TIMER_DATA          0x00400
#define GEN_VISIT_DATA          0x00800
#define GEN_VELOCITY_DATA       0x01000
#define GEN_XWIN_DATA           0x02000
#define GEN_FRAG_DATA           0x04000
#define GEN_FORCE_DATA          0x08000
#define GEN_ATOMEYE_DATA 	    0x10000
#define GEN_PARAVIEW_DATA 	    0x20000
#define GEN_LINK_DATA 	        0x40000

#ifdef _OP_REC
#define GEN_OPREC_DATA          0x80000
#endif

/*
 *      The DEL_SEG_* definitions are to be used as the <del_seg_factor>
 *      argument in invocations of the ChangeArmBurg() function.  This
 *      value indicates what portion of length of any segment deleted by
 *      the function will be accumulated in the <delSegLength> variable.
 */
#define DEL_SEG_NONE    0.0
#define DEL_SEG_HALF    0.5

#endif  /* _CONSTANTS_H */
