/***************************************************************************
 *
 *	OpList.h	Define the struct that holds all data needed
 *			to handle a cross-domain topology change.
 *
 **************************************************************************/


#ifndef _OPLIST_H
#define _OPLIST_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"

#define OpBlock_Count 500

struct _operate {
	OpType_t	type;
	int		dom1;
	int		idx1;
	int		dom2;
	int		idx2;
	int		dom3;
	int		idx3;
	real8		bx;
	real8		by;
	real8		bz;
	real8		x;
	real8		y;
	real8		z;
	real8		nx;
	real8		ny;
	real8		nz;
};

#ifdef _OP_REC
struct _operaterec {
	OpType_t	type;
	int		i1;
	int		i2;
	int		i3;
	int		i4;
	int		i5;
	int		i6;
	real8	d1;
	real8	d2;
	real8	d3;
	real8	d4;
	real8	d5;
	real8	d6;
	real8	d7;
	real8	d8;
	real8	d9;
};
#endif

/*
 *      Prototype functions related to managing the remote operation list
 */
void AddOp(Home_t *home, OpType_t type, int dom1, int idx1,
        int dom2, int idx2, int dom3, int idx3,
        real8 bx, real8 by, real8 bz, real8 x, real8 y, real8 z,
        real8 nx, real8 ny, real8 nz);
void ClearOpList(Home_t *home);
void ExtendOpList(Home_t *home);
void FreeOpList(Home_t *home);
void InitOpList(Home_t *home);
void PrintOpList(Home_t *home);

#ifdef _OP_REC
void AddOpRec(Home_t *home, OpType_t type, int dom1, int idx1,
        int dom2, int idx2, int dom3, int idx3,
        real8 bx, real8 by, real8 bz, real8 x, real8 y, real8 z,
        real8 nx, real8 ny, real8 nz);
void StartStepOpRecList(Home_t *home, int step);
void NodeMoveOpRecList(Home_t *home);
void ExtendOpRecList(Home_t *home);
void FreeOpRecList(Home_t *home);
void InitOpRecList(Home_t *home);
void WriteOpRecList(Home_t *home, char *baseFileName);
void ReadOpRecList(Home_t *home, char *oprecFile);
void ExecuteOpRec(Home_t *home, OperateRec_t *op);
void RunFromOpRecList(Home_t *home);
void RunFromOpRecFile(Home_t *home);
#endif

#endif /* _OPLIST_H */
