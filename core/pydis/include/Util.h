/****************************************************************************
 *
 *      Util.h  Contains miscellaneous definitions, and prototypes for
 *              functions defined in Util.c and used throughout the
 *              rest of the code.
 *
 ***************************************************************************/
#ifndef _Util_h
#define _Util_h

#include "Home.h"


#define DotProduct(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))


/***************************************************************************
 *
 *                              Prototypes
 *
 **************************************************************************/

/*
 *      3D to linear index mapping functions
 */
void   DecodeCellIdx(Home_t *home, int idx, int *iX, int *jY, int *kZ);
void   DecodeCell2Idx(Home_t *home, int idx, int *iX, int *jY, int *kZ);
void   DecodeDomainIdx(Home_t *home, int idx, int *iXdom, int *jYdom,
          int *kZdom);
int    EncodeCellIdx(Home_t *home, int iX, int jY, int kZ);
int    EncodeCell2Idx(Home_t *home, int iX, int jY, int kZ);
int    EncodeDomainIdx(Home_t *home, int iXdom, int jYdom, int kZdom);


/*
 *      Various numerical operations
 */
void   cross(real8 a[3], real8 b[3], real8 c[3]);
void   CSpline(real8 *x, real8 *y, real8 *y2, int numPoints);
void   CSplint(real8 *xa, real8 *ya, real8 *y2, int numPoints,
          real8 x, real8 *y);
void   DecompVec(real8 inVec[3], real8 vec1[3], real8 vec2[3], real8 ansVec[2]);
void   FindAbsMax(real8 *array, int numElements, real8 *maxVal,
          int *indexOfMax);
void   FindAbsMin(real8 *array, int numElements, real8 *minVal,
          int *indexOfMin);
void   FindMax(real8 *array, int numElements, real8 *maxVal,
          int *indexOfMax);
void   FindMin(real8 *array, int numElements, real8 *minVal,
          int *indexOfMin);
void   GetPlaneNormFromPoints(real8 p1[3], real8 p2[3], real8 p3[3],
          real8 normalVec[3]);
void   GetUnitVector(int unitFlag, real8 vx0, real8 vy0, real8 vz0,
          real8 vx1, real8 vy1, real8 vz1, real8 *ux, real8 *uy, real8 *uz,
          real8 *disMag);
void   InterpolateL(real8 *xa, real8 *ya, int numPoints, real8 x, real8 *y);
int    InterpolateBL(real8 *x, real8 *y, real8 *val, int numElem,
          real8 targetX, real8 targetY, real8 *result);
real8  Normal(real8 a[3]);
void   Normalize(real8 *ax, real8 *ay, real8 *az);
void   NormalizeVec(real8 vec[3]);
void   NormalizedCrossVector(real8 a[3], real8 b[3], real8 c[3]);
void   Orthogonalize(real8 *ax, real8 *ay, real8 *az,
          real8 bx, real8 by, real8 bz);
void   xvector(real8 ax, real8 ay, real8 az,
          real8 bx, real8 by, real8 bz,
          real8 *cx, real8 *cy, real8 *cz);


/*
 *      Nodal manipulations and operations
 */
void   AllocNodeArms(Node_t *node, int n);
void   FreeNode(Home_t *home, int index);
void   FreeNodeArms(Node_t *node);
#if __cplusplus
extern "C" Node_t *GetNeighborNode(Home_t *home, Node_t *node, int n);
#else
Node_t *GetNeighborNode(Home_t *home, Node_t *node, int n);
#endif
Node_t *GetNodeFromIndex(Home_t *home, int domID, int index);
Node_t *GetNodeFromTag(Home_t *home, Tag_t tag);
void   InsertArm(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
          real8 bx, real8 by, real8 bz,
          real8 nx, real8 ny, real8 nz, int Log);
void   MarkNodeForceObsolete(Home_t *home, Node_t *node);
void   PrintNode(Node_t *node);
void   ReallocNodeArms(Node_t *node, int n);
void   RemoveNode(Home_t *home, Node_t *node, int Log);
void   RepositionNode(Home_t *home, real8 newPos[3], Tag_t *tag, int globalOp);
void   ResetNodeArmForce(Home_t *home, Node_t *node);
void   SubtractSegForce(Home_t *home, Node_t *node1, Node_t *node2);


/*
 *      Recycled node handling functions
 */
int    GetFreeNodeTag(Home_t *home);
int    GetRecycledNodeTag(Home_t *home);
void   RecycleNodeTag(Home_t *home, int tagIdx);


/*
 *      PBC Image reconciliation functions
 */
void   FoldBox(Param_t *param, real8 *x, real8 *y, real8 *z);
void   PBCPOSITION(Param_t *param, real8 x0, real8 y0, real8 z0,
          real8 *x, real8 *y, real8 *z);
void   ZImage(Param_t *param, real8* x, real8* y, real8* z);


/*
 *      Topological operation and modification support functions
 */
void   ChangeArmBurg(Home_t *home, Node_t *node1, Tag_t *tag2,
          real8 bx, real8 by, real8 bz, real8 nx,
          real8 ny, real8 nz, int Log, real8 del_seg_factor);
int    ChangeConnection(Home_t *home, Node_t *node1, Tag_t *tag2,
          Tag_t *tag3, int Log);
void   CompressArmLists(Node_t *node);
#if __cplusplus
extern "C" int GetArmID(Home_t *home, Node_t *node1, Node_t *node2);
#else
int    GetArmID(Home_t *home, Node_t *node1, Node_t *node2);
#endif
void   GetBurgersVectorNormal(Home_t *home, Node_t *node1, Node_t *node2,
          real8 *bx, real8 *by, real8 *bz,
          real8 *nx, real8 *ny, real8 *nz);
void   RecalcSegGlidePlane(Home_t *home, Node_t *node1, Node_t *node2,
          int ignoreIfScrew);
void   ResetGlidePlane(Home_t *home, real8 newPlane[3], Tag_t *tag1,
          Tag_t *tag2, int globalOp);
void   ResetSegForces(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
          real8 fx, real8 fy, real8 fz, int globalOp);
void   ResetSegForces2(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
          real8 f1x, real8 f1y, real8 f1z,
          real8 f2x, real8 f2y, real8 f2z, int globalOp);
void   ResetNodalVelocity(Home_t *home, Node_t *nodeA, real8 vx, real8 vy,
          real8 vz, int globalOp);


/*
 *      Node ordering and sorting functions
 */
int    CollisionNodeOrder(Home_t *home, Tag_t *tagA, Tag_t *tagB);
int    DomainOwnsSeg(Home_t *home, int opClass, int thisDomain, Tag_t *endTag);
int    NodeCmpByTag(const void *, const void *);
#if __cplusplus
extern "C" int    OrderNodes(const void *a, const void *b);
#else
int    OrderNodes(const void *a, const void *b);
#endif
int    OrderTags(const void *a, const void *b);
void   SortNativeNodes(Home_t *home);


/*
 *      Prototypes for output related functions
 */
void   *CreateFragmentList(Home_t *home, int *totalFragmentCount);

#ifdef _CYGWIN
extern "C" void Fatal(const char *format, ...);
#else
#if __cplusplus
extern "C" void Fatal(const char *format, ...);
#else
void   Fatal(char *format, ...);
#endif
#endif

void   Plot(Home_t *home, int domIndex, int blkFlag);
void   WriteArms(Home_t *home, char *baseFileName, int ioGroup,
          int firstInGroup, int writePrologue, int writeEpilogue);
void   WriteDensFlux(char *fname, Home_t *home);
void   WriteDensityField(Home_t *home, char *fileName);
void   WriteFragments(Home_t *home, char *baseFileName, int ioGroup,
          int firstInGroup, int writePrologue, int writeEpilogue,
          void **fragmentList, int totalFragCount);
void   WritePoleFig(Home_t *home, char *baseFileName, int ioGroup,
          int firstInGroup, int writePrologue, int writeEpilogue);
void   WritePovray(Home_t *home, char *baseFileName, int ioGroup,
          int firstInGroup, int writePrologue, int writeEpilogue);
void   WriteAtomEye(Home_t *home, char *baseFileName, int ioGroup,
          int firstInGroup, int writePrologue, int writeEpilogue);


#if __cplusplus
extern "C" void   FindCellCenter(Param_t *param, real8 x, real8 y, real8 z, int type,
                                 real8 *xCenter, real8 *yCenter, real8 *zCenter);
#else
void   FindCellCenter(Param_t *param, real8 x, real8 y, real8 z, int type,
          real8 *xCenter, real8 *yCenter, real8 *zCenter);
#endif
void   LocateCell(Home_t *home, int *cellID, real8 coord[3]);
void   Meminfo(int *wss);
int    NodeHasSessileBurg(Home_t *home, Node_t *node);
int    NodePinned(Home_t *home, Node_t *node, int planeIndex,
          real8 glidedir[3][3]);
real8  randm(int *seed);
void   ReadTabulatedData(char *fileName, int numCols, real8 ***colData,
          int *numRows);
void   ReadRijm(Home_t *home);
void   ReadRijmPBC(Home_t *home);
int    Sign(real8 x);
void   testdeWitStress2();
void   Uniq(int*, int*);

void   Write_Node_Force_Vel(Home_t *home,char *info,int gids);

#endif
