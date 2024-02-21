#ifndef _FORCE_H
#define _FORCE_H
/***************************************************************************
 *
 *	Module:		Force.h
 *	Description:	This header is primarily for prototypes of various
 *                      functions used in calculating or updating
 *                      forces or stresses.
 *
 ***************************************************************************/
#include "Home.h"
#ifdef _SUBCYCLING
#include "SubCyc.h"
#endif

void AddtoArmForce(Node_t *node, int arm, real8 f[3]);
void AddtoNodeForce(Node_t *node, real8 f[3]);
void ComputeForces(Home_t *home, Node_t *node1, Node_t *node2,
        Node_t *node3, Node_t *node4, real8 *f1, real8 *f2,
        real8 *f3, real8 *f4);
void ComputeSegSigbRem(Home_t *home, int reqType);
void deWitInteraction(real8 MU, real8 NU, real8  Sigma[][3],
        real8  px, real8  py, real8  pz,
        real8  tp1, real8  tp2, real8  tp3,
        real8  burgX, real8  burgY, real8  burgZ,
        real8  *cntx, real8  *cnty, real8  *cntz);
void dSegImgStress(Home_t *home, real8 Sigma[][3],
        real8 px, real8 py, real8 pz,
        real8 dlx, real8 dly, real8 dlz,
        real8 burgX, real8 burgY, real8 burgZ,
        real8 rx, real8 ry, real8 rz, int pbc);
void EstRefinementForces(Home_t *home, Node_t *node1, Node_t *node2,
        real8 newPos[3], real8 vec[3], real8 f0Seg1[3], real8 f1Seg1[3],
        real8 f0Seg2[3], real8 f1Seg2[3]);
void EstCoarsenForces(Home_t *home, Node_t *node1, Node_t *node2,
        Node_t *node3, real8 f0Seg[3], real8 f1Seg[3]);
void ExtPKForce(real8 str[3][3],
        real8 bx, real8 by, real8 bz, real8 X1, real8 Y1, real8 Z1,
        real8 X2, real8 Y2, real8 Z2, real8 f1[3], real8 f2[3]);
void FindFSegComb(Home_t *home, real8 p0[3], real8 p1[3], real8 p2[3],
        real8 burg1[3], real8 burg2[3], real8 fp0seg1[3],
        real8 fp1seg1[3], real8 fp1seg2[3], real8 fp2seg2[3],
        real8 f0new[3], real8 f1new[3]);
void FindSubFSeg(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
        real8 oldfp1[3], real8 oldfp2[3], real8 newpos[3],
        real8 f0seg1[3], real8 f1seg1[3], real8 f0seg2[3],
        real8 f1seg2[3]);
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z,
        real8 totStress[3][3]);
void GetFieldPointStressRem(Home_t *home, real8 x, real8 y, real8 z,
                            int cellX, int cellY, int cellZ,
                            real8 totStress[3][3]);
void LineTensionForce(Home_t *home, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
        real8 f1[3], real8 f2[3]);
void LocalSegForces(Home_t *home, int reqType);
void NodeForce(Home_t *home, int reqType);
void OsmoticForce(Home_t *home, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
        real8 f1[3], real8 f2[3]);
void PKForce(real8 sigb[3],
        real8 X1, real8 Y1, real8 Z1, real8 X2, real8 Y2, real8 Z2,
        real8 f1[3], real8 f2[3]);
void ReevaluateForces(Home_t *home);

#ifdef _FEM
void SegmentStress(real8 MU, real8 NU,
        real8 burgX, real8 burgY, real8 burgZ, real8 xA, real8 yA, real8 zA,
        real8 xB, real8 yB, real8 zB, real8 x0, real8 y0, real8 z0,
        real8 a, real8 Sigma[3][3]);
#endif
void SegmentSigb(real8 MU, real8 NU,
                 real8 bX, real8 bY, real8 bZ, real8 xA, real8 yA, real8 zA,
                 real8 xB, real8 yB, real8 zB, real8 x0, real8 y0, real8 z0,
                 real8 a, real8 bpX, real8 bpY, real8 bpZ, real8 Sigb[3] );

void SegSegForceIsotropic(real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 MU, real8 NU,
        int seg12Local, int seg34Local,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp2x, real8 *fp2y, real8 *fp2z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SegSegForce_SBN1(real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 MU, real8 NU,
        int Nint, real8 *quad_points, real8 *weights,
        int seg12Local, int seg34Local,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp2x, real8 *fp2y, real8 *fp2z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SegSegForce_SBN1_SBA(real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 MU, real8 NU, int Nint,
        int seg12Local, int seg34Local,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp2x, real8 *fp2y, real8 *fp2z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SegSegForce(real8 p1x, real8 p1y, real8 p1z,
	real8 p2x, real8 p2y, real8 p2z,
	real8 p3x, real8 p3y, real8 p3z,
	real8 p4x, real8 p4y, real8 p4z,
	real8 bpx, real8 bpy, real8 bpz,
	real8 bx, real8 by, real8 bz,
	real8 a, real8 MU, real8 NU,
	int seg12Local, int seg34Local,
	real8 *fp1x, real8 *fp1y, real8 *fp1z,
	real8 *fp2x, real8 *fp2y, real8 *fp2z,
	real8 *fp3x, real8 *fp3y, real8 *fp3z,
	real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SelfForce(int coreOnly, real8 MU, real8 NU,
        real8 bx, real8 by, real8 bz, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 a,  real8 Ecore,
        real8 f1[3], real8 f2[3]);
void SemiInfiniteSegSegForce(real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 MU, real8 NU,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SetOneNodeForce(Home_t *home, Node_t *node1);
void ZeroNodeForces(Home_t *home, int reqType);

#ifdef _SUBCYCLING
void SegSegListMaker(Home_t *home, int reqType);
void NodeForceList(Home_t *home, SegSeg_t *SegSegList, int SegSegListSize,
                   Segm_t *SegList, int SegListSize, int subGroup);
int CellPriority(Home_t *home, int cellID1, int cellID2);
#endif

#endif  /* ifndef _Force_h */
