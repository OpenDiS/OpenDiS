/***************************************************************************
 *
 *	Module:		ParadisProto.h
 *	Description:	This header is mainly just a dumping ground for
 *			the miscellaneous funtion prototypes that have
 *			not been included elsewhere.  This helps eliminate
 *			some of the compiler whines...
 *
 ***************************************************************************/
#ifndef _ParadisProto_h
#define _ParadisProto_h

#include "stdio.h"
#include "Tag.h"
#include "Home.h"

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/*
 *      Define some inline operations for handling 3-component vectors
 *
 * VECTOR_ADD:  add the components of the second vector to the first
 * VECTOR_COPY: copy the components of the second vector to the first
 * VECTOR_ZERO: Zero out the components of the vector
 */
#define VECTOR_ADD(a,b)  {(a)[0] += (b)[0]; (a)[1] += (b)[1]; (a)[2] += (b)[2];}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}
#define VECTOR_ZERO(a)   {(a)[0] = 0; (a)[1]  = 0; (a)[2] = 0;}

#ifdef __cplusplus
extern "C" void Getline(char *string, int len, FILE *fp);
#else
void Getline(char *string, int len, FILE *fp);
#endif

/* From NodeForce.c - stress due to a segment in coordinate-independent form */
void StressDueToSeg(real8 px, real8 py, real8 pz,
                    real8 p1x, real8 p1y, real8 p1z,
                    real8 p2x, real8 p2y, real8 p2z,
                    real8 bx, real8 by, real8 bz,
                    real8 a, real8 MU, real8 NU, real8 *stress);

/* From SegmentStress.c - stress due to a segment in coordinate-dependent form */
void SegmentStress(real8 MU, real8 NU,
                   real8 bX, real8 bY, real8 bZ,
                   real8 xA, real8 yA, real8 zA,
                   real8 xB, real8 yB, real8 zB,
                   real8 x0, real8 y0, real8 z0,
                   real8 a,
                   real8 Sigma[3][3]);

void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void GetVelocityStatistics(Home_t *home);
void AssignNodeToCell(Home_t *home, Node_t *node);
void TrapezoidIntegrator(Home_t *home);
#ifdef _SUBCYCLING
void RKFIntegrator(Home_t *home, int reqType);
void SubcycleIntegratorForceB(Home_t *home);
void FlagSubcycleNodes(Home_t *home, int subGroup);
#endif
void BroadcastDecomp(Home_t *home, void *decomp);
int  CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll);
void CellCharge(Home_t *home);
void CrossSlip(Home_t *home);
void CrossSlipBCC(Home_t *home);
void CrossSlipFCC(Home_t *home);
int ThermalActivation(Param_t *param, real8 Eb, real8 attfreq);
void EnergyBarrierCrossSlipFCC(Home_t *home, Node_t *node, real8 burg[3], real8 linedir[3], real8 newplane[3], real8 *Eb);
void DeltaPlasticStrain(Home_t *home);
void DeltaPlasticStrain_BCC(Home_t *home);
void DeltaPlasticStrain_FCC(Home_t *home);
void DistributeTagMaps(Home_t *home);
void FindPreciseGlidePlane(Home_t *home, real8 burgVecIn[3], real8 dirIn[3],
        real8 glidePlane[3]);
void FixGlideViolations(Home_t *home, Tag_t *tag, real8 oldpos[3]);
void FixRemesh(Home_t *home);
void ForwardEulerIntegrator(Home_t *home);
void FreeCellCenters(void);
void FreeCorrectionTable(void);
void FreeInitArrays(Home_t *home, InData_t *inData);
void FreeInNodeArray(InData_t *inData, int numNodes);
void FreeRijm(void);
void FreeRijmPBC(void);
void GenerateOutput(Home_t *home, int stage);
void GetDensityDelta(Home_t *home);
void GetMinDist(
        real8 p1x, real8 p1y, real8 p1z, real8 v1x, real8 v1y, real8 v1z,
        real8 p2x, real8 p2y, real8 p2z, real8 v2x, real8 v2y, real8 v2z,
        real8 p3x, real8 p3y, real8 p3z, real8 v3x, real8 v3y, real8 v3z,
        real8 p4x, real8 p4y, real8 p4z, real8 v4x, real8 v4y, real8 v4z,
        real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2);
void GetMinDist2(
        real8 p1x, real8 p1y, real8 p1z, real8 v1x, real8 v1y, real8 v1z,
        real8 p2x, real8 p2y, real8 p2z, real8 v2x, real8 v2y, real8 v2z,
        real8 p3x, real8 p3y, real8 p3z, real8 v3x, real8 v3y, real8 v3z,
        real8 p4x, real8 p4y, real8 p4z, real8 v4x, real8 v4y, real8 v4z,
        real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2);
void GetNbrCoords(Home_t *home, Node_t *node, int arm, real8 *x, real8 *y,
        real8 *z);
void GetParallelIOGroup(Home_t *home);
void HandleCollisions(Home_t *home);
void PredictiveCollisions(Home_t *home);
void ProximityCollisions(Home_t *home);
#ifdef _RETROCOLLISIONS
void RetroactiveCollisions(Home_t *home);
void RetroactiveCollisions2(Home_t *home);
int CollisionCriterion(real8 *dist2, real8 *L1ratio, real8 *L2ratio, const real8  mindist,
                       const real8 *x1t, const real8 *x1tau, const real8 *x2t, const real8 *x2tau,
                       const real8 *x3t, const real8 *x3tau, const real8 *x4t, const real8 *x4tau);
int HingeCollisionCriterion (real8 *L1ratio, const real8 *vec1, const real8 *vec2);
#endif
void HeapAdd(int **heap, int *heapSize, int *heapCnt, int value);
int  HeapRemove(int *heap, int *heapCnt);
void InitRemoteDomains(Home_t *home);
void InputSanity(Home_t *home);
int isCollinear(double nX, double nY, double nZ, double pX, double pY, double pZ);
void LoadCurve(Home_t *home, real8 deltaStress[3][3]);
void Migrate(Home_t *home);
int  NodeOwnsSeg(Home_t *home, Node_t *node1, Node_t *node2);
void ParadisStep(Home_t *home);
void ParadisFinish(Home_t *home);
void PickScrewGlidePlane(Home_t *home, real8 burgVec[3],
        real8 glidePlane[3]);
void ReadNodeDataFile(Home_t *home, InData_t *inData, char *dataFile);
void FreeAllNodes(Home_t *home);
void RecycleAllNodes(Home_t *home);
void ReadDataFile(Home_t *home, char *dataFile);
void RemapArmTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void Remesh(Home_t *home);
void RemeshRule_2(Home_t *home);
void RemeshRule_3(Home_t *home);
void ResetGlidePlanes(Home_t *home);
void SetLatestRestart(char *fileName);
void SortNodesForCollision(Home_t *home);
#ifdef _RETROCOLLISIONS
void SortNodes(Home_t *home, real8 *maxsep);
#endif
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
        int writePrologue, int writeEpilogue, int numSegs);
void TestGlidePlanes(Home_t *home, Node_t *node, real8 p[3], int coplanar, int *passtest);
void UniformDecomp(Home_t *home, void **decomp);
void WriteVelocity(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteForce(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteVisit(Home_t *home, char *baseFileName, int writePrologue,
        int writeEpilogue, int *nodesWritten, int *segsWritten);
void WriteVisitMetaDataFile(Home_t *home, char *baseFileName,
        int *groupDataCounts);
void WriteParaview(Home_t *home, char *baseFileName, int numSegs, int numLocSegs);

/*
 *      Some support function used by the various cross-slip modules
 */
#ifdef DEBUG_CROSSSLIP_EVENTS
void DumpCrossSlipEvent(Node_t *node, real8 newPlane[3], char *eventDesc);
#endif
void SaveCrossSlipInfo(Node_t *node, Node_t *nbr1, Node_t *nbr2,
        int nbr1ArmID, int nbr2ArmID, real8 segForceOrig[4][3],
        real8 nodePosOrig[3], real8 nbr1PosOrig[3],
        real8 nbr2PosOrig[3]);
void ResetPosition(Param_t *param, Node_t *node, real8 pos[3]);
void RestoreCrossSlipForce(Node_t *node, Node_t *nbr1, Node_t *nbr2,
        int nbr1ArmID, int nbr2ArmID,
        real8 segForceOrig[4][3]);

Node_t *RequestNewNativeNodeTag(Home_t *home, Tag_t *tag);
void AddNodesFromArray(Home_t *home, real8 *buf);
void ReleaseMemory(Home_t *home);

void ParadisInit_lean(Home_t **homeptr);
void Initialize_lean(Home_t *home);
#ifdef _GPU_SUBCYCLE
#ifdef __cplusplus
extern "C" {
void InitializeParadisGPU(Home_t *home);
}
#else
void InitializeParadisGPU(Home_t *home);
#endif
#endif

#ifdef CHECK_MEM_USAGE
void _CheckMemUsage(Home_t *home, char *msg);
#define CheckMemUsage(a,b) _CheckMemUsage((a),(b))
#else
#define CheckMemUsage(a,b) {}
#endif

#endif
