/***************************************************************************
 *
 *      Module:         Decomp.h
 *      Description:    This contains structures and prototypes needed
 *                      in modules that deal with the generic functions
 *                      for creating, accessing and manipulating the
 *                      domain decomposition.
 *
 ***************************************************************************/
#ifndef _Decomp_h
#define _Decomp_h

/*
 *      The dynamic load balancing is done based on per-process load
 *      data.  Below are definitions for the various types of load
 *      information that might be used.
 */
#define DLB_USE_WALLCLK_TIME    0
#define DLB_USE_FORCECALC_COUNT 1

/*
 *      When allocating an RB decomposition structure, in some cases
 *      we want to allocate a new decomposition and initialize it from
 *      scratch, in other cases, we want to duplicate an existing
 *      decomposition.  Define values to indicate to the allocate
 *      function which behaviour is desired.
 */
#define ALLOC_NEW_DECOMP       0
#define ALLOC_DUPLICATE_DECOMP 1

/*
 *      The RSDecomp_t structure defines arrays which contain the domain
 *      boundaries for the entire problem based on the original (recursive
 *      sectioning) version of the domain decomposition. There are only
 *      nXdoms partitions in the x direction. For each such partition
 *      there are nYdoms partitions, not necessarily at the same place
 *      for different x partitions. Thus a 2D array is needed.  Similarly,
 *      a 3D array  is needed for the final partitioning in the z direction.
 *
 *      See comments in RSDecomp.c for more details.
 */
typedef struct {
        real8 *domBoundX;
        real8 **domBoundY;
        real8 ***domBoundZ;
} RSDecomp_t;


/*
 *      The RBDecomp_t structure contains information related to
 *      the recursive bisection decomposition of the problem space.
 *      This decomposition lends itself nicely to a tree-like structure
 *      where each structure represents a portion of the decomposition.
 *      If the portion corresponds to a single domain, the domID value
 *      will be >= zero and the subDecomp array will contain
 *      all NULL pointers.  If the structure does not represent
 *      a single domain, the domID value will be < zero and two
 *      or more of the subDecomp pointers will point to 
 *      subDecomposition at the next level down in the decomposition
 *      tree.  
 * 
 *      The numDoms array represents the number of domains (per dimension)
 *      included underneath this portion of the decomposition.  Needed
 *      for handling domain decompositions which do not have domain
 *      counts that are powers of two.
 */
#define MAX_DECOMP_LVLS 12

typedef struct _rbdecomp RBDecomp_t;

struct _rbdecomp {
        int          numDoms[3];
        int          domID;
        double       totLoad;
        double       cMin[3];
        double       cMax[3];
        RBDecomp_t   *subDecomp[8];
        char         decompID[MAX_DECOMP_LVLS];
};


/*
 *      Include prototypes for the set of generic function calls
 *      used to obtain information related to the domain decomposition.
 */
void BroadcastDecomp(Home_t *home, void *decomp);
void DLBfreeOld(Home_t *home);
int  FindCoordDomain(Home_t *home, int updateCoords, real8 *x,
         real8 *y, real8 *z);
void FreeDecomp(Home_t *home, void *decomp);
void GetAllDecompBounds(Home_t *home, real8 **decompBounds, int *numValues);
void GetCellDomainList(Home_t *home, int cellID, int *domCount,
         int **domList);
void GetLocalDomainBounds(Home_t *home, void *decomp);
void XPlotDecomp(Home_t *home, real8 xMin, real8 yMin, real8 zMin,
         real8 lMax, int color, real8 lineWidth);
void ReadDecompBounds(Home_t *home, void **filePtr, int doBinRead,
         int newDecompType, void **oldDecomp);
void Rebalance(Home_t *home, int criteria);
void UniformDecomp(Home_t *home, void **decomp);
void WriteDecompBounds(Home_t *home, FILE *fp);

#endif
