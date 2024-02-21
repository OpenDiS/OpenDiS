#ifndef _FM_h
#define _FM_h

#define X   0
#define Y   1
#define Z   2

#define NMAX        20    /* WARNING! MUST be at least param->fmOrder+4 */
#define NTMAX       (((NMAX+1)*(NMAX+2))>>1)
#define MAXORDER    (2*NMAX)

#define CELL_HASH_TABLE_SIZE 97

#define MP_COEFF     1
#define TAYLOR_COEFF 2

/*
 *  These macros are used in the FM code when calculating indices
 *  for all nearby neighbors of cells (where nearby is essentially
 *  all cells that are children of the parent (or children
 *  of the immediate neighbors of the parent) of cell <x>).
 *  Returned indices may be negative or larger than the dimensions
 *  of a block of cells, but adjustments for PBC are made elsewhere.
 */
#define SHIFTPOS(x)       (((((x) >> 1) + 1) << 1) + 1)
#define SHIFTNEG(x)       ((((x) >> 1) - 1) << 1)
/*
 *  Returns the index of <a> in a dimension of length <b> after
 *  correction for periodic boundary conditions.  If the initial
 *  condition is  0 <= a < b the returned value will be <a>,
 *  otherwise a will be adjusted.
 */
#define GETPBCINDEX(a,b)  (((a) % (b)) + ((b) * ((a) < 0)))


typedef double matrix[3][3];
typedef double vector[3];

typedef struct _fmcell FMCell_t;
typedef struct _fmlayer FMLayer_t;

struct _fmcell {
        int   cellID;     /* ID of the cell in this FM layer */
        int   owningDom;  /* Index of domain owning this cell.  There is  */
                          /* a single unique owner for each cell at each  */
                          /* FM layer.                                    */
        int   domCnt;     /* number of domains intersecting this cell.    */
                          /* this cell.                                   */

        int   *domList;   /* list of all domains that intersect this cell.*/
                          /* This list will only be created for cells at  */
                          /* the most refined FM layer, and is recomputed */
                          /* each time domain boundaries shift.           */

        real8 cellCtr[3]; /* Coordinates of the center of the cell */

        real8 *mpCoeff;   /* Array of multipole expansion coefficients */
                          /* for the cell.  Number of coefficients is  */
                          /* (m+3)*(m+2)*(m+1)/6*9 where m is the order*/
                          /* of the multipole expansion                */

        real8 *taylorCoeff; /* Array of taylor expansion coefficients    */
                            /* for the cell.  Number of coefficients is  */
                            /* (m+3)*(m+2)*(m+1)/6*9 where m is the order*/
                            /* of the multipole expansion                */
        FMCell_t *next;
        FMCell_t *prev;
};


struct _fmlayer {
        int   lDim[3];     /* Number of cells in each dimension at this layer */

        real8 cellSize[3]; /* dimensions of each cell at this layer */
    
        int   ownedCnt;    /* Total number of cells owned by this domain */
                           /* at this layer                              */

        int   ownedMin[3]; /* Min and max indices of cells owned by this */
        int   ownedMax[3]; /* domain.  This information is statically    */
                           /* defined during initialization.             */

        int   intersectCnt;     /* Number of cells intersecting this domain   */
        int   intersectMin[3];  /* Min and max indices of cells intersecting  */
        int   intersectMax[3];  /* this domain.  These are meaningless except */
                                /* at the most refined FM layer.  At the that */
                                /* layer, the intersection list will be       */
                                /* recomputed each cycle to handle shifting   */
                                /* domain boundries.                          */
        int        *domBuf;

/*
 *      Some lists of domains with which the current domain will 
 *      communictae during the updward and downward FM passes
 */
        int   fmUpPassSendDomCnt;
        int   *fmUpPassSendDomList;

        int   fmUpPassRecvDomCnt;
        int   *fmUpPassRecvDomList;

        int   fmDownPassSendDomCnt;
        int   *fmDownPassSendDomList;

        int   fmDownPassRecvDomCnt;
        int   *fmDownPassRecvDomList;

        int   *cellList; /* List of cellIDs that have been added to */
                         /* the cell table at this FM layer         */

        int   numCells;  /* Number of cells in the cell table at  */
                         /* this FM layer                         */

        FMCell_t *cellTable[CELL_HASH_TABLE_SIZE]; /* Hash table of  */
                         /* FM cell pointers. Contains all cells known  */
                         /* to the current domain at this FM layer      */
};


/*
 *      Prototypes for general FM for remote seg/seg forces
 */
void  CorrectionTableInit(Home_t *home);
void  CreateCorrectionTable(Param_t *param, int numLevels, int pbc[3],
          int mpi_np, int mpi_pid);
void  DecodeFMCellIndex(int dim[3], int cellID, int *x, int *y, int *z);
void  DoTableCorrection(Home_t *home);
int   EncodeFMCellIndex(int *dim, int x, int y, int z);
void  EvalTaylor(int norder, real8 *r, real8 *alpha, real8 sigma[3][3]);
void  FMCommUpPass(Home_t *home, int layerID);
void  FMCommDownPass(Home_t *home, int layerID);
void  FMDistTaylorExp(Home_t *home);
void  FMFree(Home_t *home);
void  FMInit(Home_t *home);
FMCell_t *LookupFMCell(FMCell_t *cellTable[], int cellID);
void  GetCellDomainList(Home_t *home, int cellIndex, int *domCount,
          int **domList);
void  MkTaylor(real8 MU, real8 NU, int norder, int uorder, int maxorder,
          real8 *r, real8 *eta, real8 *alpha);
void  FMSetRemoteStress(Home_t *home);
void  FMSetTaylorExpansions(Home_t *home);
void  FMShift(int norder, real8 *r, real8 *eta, real8 *neta);
void  FMSigma2(real8 mu, real8 nu, int norder, real8 Eeta[],
          matrix sigmatot, real8 pcterms[], int pcpows[][3]);
void  FMFindNearNbrs(Home_t *home, int layerID, int *inMin, int *inMax,
          int *outMin, int *outMax, int trimOverlap);
void  GaussQuadCoeff(int intOrder, real8 *positions, real8 *weights);
real8 ipow(real8 x, int n);
void  makeeta(int norder, real8 ksi0[3], real8 ksi[3], real8 b[3], real8 *Eeta);
void  makeftabs(real8 *fact, real8 *ifact, real8 *dfact);
void  makeqtab(real8 qtab[NMAX+1][NMAX+1]);
void  MeanStressCorrection(Home_t *home);
void  RemoteForceOneSeg(Home_t *home, Node_t *node1, Node_t *node2,
          real8 f1[3], real8 f2[3]);
void  SegForceFromTaylorExp(Home_t *home, int cellID, real8 *positions,
          real8 *weights, real8 *p1, real8 *p2, real8 *burg,
          real8 p1f[3], real8 p2f[3]);
void  TaylorShift(int norder, real8 *r, real8 *alpha, real8 *beta);
 
void  FMPrintMPCoeff(Home_t *home, int layerID);
void  FMPrint(Home_t *home);

#endif  /* FM_h */
