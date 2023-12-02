/****************************************************************************
 *
 *  Home.h  Define the Home struct which holds or points to all relevant
 *          data for this domain, plus define miscellaneous things
 *          that are used throughout the code.
 *
 *  Internal Units:
 *
 *    length in b (burgMag, read in meters, e.g. 2.725e-10)
 *    stress in Pa
 *    time   in second
 *
 *    force  in Pa*b^2
 *    force per unit length in Pa*b
 *    velocity in b/second
 *
 *    mobility in 1/Pa/second
 *
 ***************************************************************************/

#ifndef _Home_h
#define _Home_h

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <limits.h>

#ifdef USE_HDF
#include "hdf5.h"
#endif

#include "Constants.h"
#include "ParadisThread.h"
#include "Typedefs.h"
#include "ParadisProto.h"
#include "FM.h"
#include "Node.h"
#include "Param.h"
#include "Cell.h"
#include "RemoteDomain.h"
#include "Tag.h"
#include "MirrorDomain.h"
#include "Topology.h"
#include "OpList.h"
#include "Timer.h"
#include "Util.h"
#include "Init.h"
#include "Parse.h"
#include "Matrix.h"
#include "DebugFunctions.h"
#include "Force.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _SUBCYCLING
#include "SubCyc.h"
#endif

#ifdef _GPU_SUBCYCLE
#include "SubcycleGPU.h"
#endif

/*
 *      Don't let the code compile with both PARALLEL and FULL_N2_FORCES
 *      defined...
 */
#ifdef PARALLEL
#ifdef FULL_N2_FORCES
#error "Cannot define FULL_N2_FORCES with PARALLEL (See makefile.setup)"
#endif
#endif

/*
 *      If memory debugging is enabled, intercept calls to the
 *      memory management functions by redefining them to be
 *      locally defined functions.  This allows us to do some
 *      tracking of allocated memory and some checks to help
 *      identify memory corruption and memory leaks.
 */
#ifdef DEBUG_MEM

#define free(a)        ParadisFree(__FILE__, __LINE__, (a))
#define malloc(a)      ParadisMalloc(__FILE__, __LINE__, (a))
#define calloc(a, b)   ParadisCalloc(__FILE__, __LINE__, (a), (b))
#define realloc(a, b)  ParadisRealloc(__FILE__, __LINE__, (a), (b))

#ifdef __cplusplus
extern "C" void ParadisMemCheck(void);
extern "C" void ParadisFree(char *fileName, int lineNum, void *ptr);
extern "C" void *ParadisMalloc(char *fileName, int lineNum, size_t size);
extern "C" void *ParadisCalloc(char *fileName, int lineNum, size_t numElem,
                    size_t size);
extern "C" void *ParadisRealloc(char *fileName, int lineNum, void *ptr,
                    size_t size);
#else
void ParadisMemCheck(void);
void ParadisFree(char *fileName, int lineNum, void *ptr);
void *ParadisMalloc(char *fileName, int lineNum, size_t size);
void *ParadisCalloc(char *fileName, int lineNum, size_t numElem, size_t size);
void *ParadisRealloc(char *fileName, int lineNum, void *ptr, size_t size);
#endif

#endif /* ifdef DEBUG_MEM */


/*
 *      On BGL, rint was *really* slow, so here's a kludge from one
 *      of the IBM guru's
 */
#ifdef _BGL
#define rint(x) (6.7553994410557440e+15 + (x) - 6.7553994410557440e+15)
#endif


typedef struct {
        Tag_t  oldTag;
        Tag_t  newTag;
} TagMap_t;

#ifndef _SUBCYCLING
typedef struct {
        int        forcesSet;
        real8      f1[3];
        real8      f2[3];
        Node_t     *node1;
        Node_t     *node2;
#ifdef _OPENMP
        omp_lock_t segLock;
#endif
} Segment_t;
#endif

struct _segmentpair {
        Segment_t *seg1;
        Segment_t *seg2;
        int       setSeg1Forces;
        int       setSeg2Forces;
};


/*
 *     Define a structure to hold the burgers vectors and associated
 *     glide planes for material for which we dynamically calculate
 *     those values.
 */
typedef struct {
        int   numBurgVectors;       /* Total number of burgers vectors */
                                    /* in <burgList> below             */

        int   numPlanes;            /* Total number of planes in */
                                    /* <planeList> below         */

        int   *numPlanesPerBurg;    /* Number of planes associated with */
                                    /* each burgers vector in the list  */

        int   *burgFirstPlaneIndex; /* Index (0 offset) of the first plane */
                                    /* in <planeList> associated with each */
                                    /* burgers vector                      */

        real8 (*burgList)[3];       /* Array of burgers vectors */

        real8 (*planeList)[3];      /* Arraay of burgers vector planes */
} BurgInfo_t;




struct _home {

        int       myDomain;    /* encoded domain index for this domain */
        int       numDomains;
        int       cycle;       /* current cycle */
        int       lastCycle;   /* cycle to quit on */

        Param_t   *param;
#ifdef _SUBCYCLING
        Subcyc_t  *subcyc;
#endif

/*
 *      The following two pointers are used for registering control
 *      and data file parameters and when reading/writing restart
 *      files.  They will only be set and used in task zero.
 */
        ParamList_t *ctrlParamList;
        ParamList_t *dataParamList;

/*
 *      The following three pointers are queue-heads of singly-linked lists of
 *      nodes. Nodes are pushed onto the head of a queue, but may be removed
 *      from anywhere in the queue. Blocks of nodes are allocated at one time
 *      for efficiency. These blocks are themselves kept in the nodeBlockQ,
 *      so they can be found if they ever need to be.
 */
        Node_t    *nativeNodeQ;
        Node_t    *ghostNodeQ;
        Node_t    *freeNodeQ;

        Node_t    *lastFreeNode;
        Node_t    *lastGhostNode;

        NodeBlock_t  *nodeBlockQ;

/*
 *      the nodeKeys array contains pointers to Node_t structs. For a node
 *      with a given tag, nodeKeys[tag.index] points to the node's struct.
 *      The recycle node heap contains the indices of tags that were in
 *      use, but have been freed. When a node is created or moves into the
 *      domain, its tag is assigned from the recycle node heap if possible.
 *      Note: node tags retrieved from the heap will always be such that
 *      the tag retrieved is the lowest available tag.
 */
        Node_t    **nodeKeys;
        int       newNodeKeyPtr;
        int       newNodeKeyMax;

        int       *recycledNodeHeap;
        int       recycledNodeHeapSize;
        int       recycledNodeHeapEnts;

/*
 *      cellList keeps a list of all base cells of interest to this domain,
 *      whether native, ghost, or the base cell of a periodic ghost cell.
 *      The first nativeCellCount items in the list are native cells.
 *      cellKeys is long enough to include the encoded indices of all cells
 *      in the problem, including any periodic cells. For each cell (base
 *      or periodic) allocated for this domain cellKeys contains, at that
 *      cell's encoded index, a pointer to the cell.
 */
        int       *cellList;
        int       cellCount;
        int       nativeCellCount;
        int       firstTime_FMInit;

        Cell_t    **cellKeys;

/*
 *      There are two classes of remote domains. The primary remote domains
 *      are those associated with any primary ghost nodes (i.e. any
 *      domain that intersects any cell that intersects (or is an immediate
 *      neighbor of a cell native to the current domain).  The second class
 *      of remote domains are those from which the current domain only
 *      requires secondary ghost nodes.
 */
        int       remoteDomainCount; /* Number of primary remote domains */

        int       secondaryRemoteDomainCount;  /* Number of secondary */
                                               /* remote domains only */

        int       *remoteDomains;    /* encoded indices of the neighbor */
                                     /* domains.  Includes primary and  */
                                     /* secondary remote domains        */

        RemoteDomain_t  **remoteDomainKeys; /* pointers to RemoteDomain_t    */
                                            /* structs of neighboring remote */
                                            /* domains                       */
/*
 *      To allow for multiple types of domain decomposition, we
 *      maintain a generic pointer to decomposition data here.
 *      The pointer will be cast to the proper type of pointer
 *      based on the type of decomposition specified in the control
 *      file.
 */
        void      *decomp;

/*
 *      Also keep the domain boundaries for this domain, in more
 *      accessible form... saves us from looking these values up
 *      in the domain decomposition multiple times during each step.
 */
        real8     domXmin;
        real8     domXmax;
        real8     domYmin;
        real8     domYmax;
        real8     domZmin;
        real8     domZmax;

/*
 *      For the RCB domain decomposition, the number of levels in the
 *      decomposition hierarchy is dependent on the number of domains
 *      in each dimension.  So, at startup determine the maximum levels
 *      and store the values for future reference.
 */
        int       xMaxLevel;
        int       yMaxLevel;
        int       zMaxLevel;

/*
 *      Some MPI stuff. Might be better if it was in a comm sub-structure
 */
#ifdef PARALLEL
        int       maxPackSiz;     /* byte length required to accomodate */
                                  /* the largest communication buffer  */
                                  /* being sent.  Nearly obsolete: is  */
                                  /* now only used in mirror node comm */

        MPI_Request  *inRequests; /* used for asynchronous I/O */
        MPI_Request  *outRequests;
        MPI_Status   *inStatus;
        MPI_Status   *outStatus;

#endif

/*
 *      array of Burgers vectors -- almost obsolete
 */
        real8     *burgX;
        real8     *burgY;
        real8     *burgZ;

        int       nburg;

/*
 *      Used for various types of output generation when passing nodal
 *      information from ALL domains to domain 0
 */
        MirrorDomain_t  **mirrorDomainKeys;     /* image of entire domains */
        int             currentMirrors;         /* 1 if mirrors already sent */
                                                /* this cycle; else 0        */
        char      *inBuf;
        char      *outBuf;

/*
 *      Operation list for topological changes across domain
 */
        Operate_t *opList;
        int       OpCount;
        int       OpListLen;

#ifdef _OP_REC
        OperateRec_t *opRecList;
        int           OpRecCount;
        int           OpRecListLen;
#endif

        Operate_t *rcvOpList;
        int       rcvOpCount;

        Timer_t   *timers;

/*
 *      Define values related to cell2 grid overlaid on standard
 *      cell structure during collision handling.
 */
        int       *cell2;
        C2Qent_t  *cell2QentArray;

        int       cell2nx;
        int       cell2ny;
        int       cell2nz;

#ifdef _RETROCOLLISIONS
        real8	  cell2size;
#endif

/*
 *      array to hold net charge tensor of each cell
 */
        real8     *cellCharge;

/*
 *      During initialization and node migration, nodes may move
 *      between domains.  Nodes moved between domains receive new
 *      ID tags.  The tagMap pointer will be used to hold an array
 *      of mappings between old and new tags for these moving
 *      nodes with tagMapSize and tagMapEnts indicating the number
 *      of elements currently allocated in the tagMap and the
 *      number of those elements that are actually in use.
 */
        TagMap_t  *tagMap;
        int       tagMapSize;
        int       tagMapEnts;

/*
 *      Remote force calculations use a couple arrays for doing
 *      Gauss-Legendre integration, so we'll keep them here where
 *      they'll be available whenever needed.
 */
        real8     *glPositions;
        real8     *glWeights;

/*
 *      Define things needed for the fast-multipole stuff
 */
        FMLayer_t *fmLayer;         /* Pointer to array of structures */
                                    /* (1 per layer of the FM hier-   */
                                    /* archy).                        */

        int       fmNumMPCoeff;     /* Number of coefficients required*/
                                    /* for multipole expansions of the*/
                                    /* order specified in the control */
                                    /* file.                          */

        int       fmNumTaylorCoeff; /* Number of coefficients required*/
                                    /* for talyor expansions of the   */
                                    /* order specified in the control */
                                    /* file.                          */

/*
 *      Used only when load-balancing is done based on the number of
 *      segment to segment force calculations rather than the actual
 *      wallclock time spent.
 */
        int       cycleForceCalcCount;

/*
 *      For some simulations, the geometric information (burgers vectors,
 *      normal planes) may be provided in a user-specified laboratory frame.
 *      (See the <useLabFrame> control file parameter.  If this is being
 *      done, we'll need to rotate between the lab frame and the standard
 *      crystal frame during various portions of the simulation, so we
 *      compute the rotation matrix and invert once during initialization
 *      and store the matrices here for subsequent use.
 */
        real8 rotMatrix[3][3];
        real8 rotMatrixInverse[3][3];

/*
 *      For some of the material types, we dynamically calculate the list
 *      of possible burgers vectors (and associated glide planes) during
 *      initialization.  For those material types, we'll store all that
 *      information in the BurgListInfo_t structure for use throughout the
 *      code as needed.
 */
        BurgInfo_t burgData;

/*
 *      For parallel IO, each process needs to know things like
 *      which IO group it's in, the first and last processes in
 *      the group and so on.  Since those values are static, we
 *      can just calculate them once and store them in the Home_t
 *      structure for later use.
 */
        int   ioGroupNum;    /* IO group containing this process */

        int   firstInIOGroup, lastInIOGroup; /* first and last processes in */
                                             /* this IO group. */

        int   prevInIOGroup, nextInIOGroup;  /* Previous and next processes */
                                             /* in this IO group */

        int   isFirstInIOGroup;   /* Set to 1 if this process is the  */
                                  /* first in its IO group */
        int   isLastInIOGroup;    /* Set to 1 if this process is the  */
                                  /* last in its IO group */
#ifdef PARALLEL
        MPI_Comm commLastInIOGroup;  /* Communicator encompassing only   */
                                     /* the last process in each IO group*/
#endif

        struct timeval clock_time_beg;

#ifdef _FEM
/*
 *      When calculating Yoffe stress, we need to build a global list
 *      of segments that intersect the free surfaces.  These two
 *      variables are used for dealing with that list.  See code
 *      for description of contents of the segment list array.
 */
        int   surfaceSegCount;
        real8 *surfaceSegList;
#endif

#ifdef _GPU_SUBCYCLE
        int           deviceID;
        DeviceProp_t *deviceProp;
#endif
};

/* prototypes */
/* FIX ME!  Most of these prototypes should be moved elsewhere */

void    AddNode(Home_t *home, Node_t *nodeA, Node_t *newNode, Node_t *nodeB);
void    CommSendMirrorNodes(Home_t *home, int stage);
int     Connected(Node_t *node1, Node_t *node2, int *armidx);
Node_t *GetNewNativeNode(Home_t *home);
Node_t *GetNewGhostNode(Home_t *home, int domain, int index);
void    Gnuplot(Home_t *home, char *baseFileName, int ioGroup,
            int firstInGroup, int writePrologue, int writeEpilogue);

#endif
