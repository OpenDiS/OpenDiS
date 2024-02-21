#ifndef _Matrix_h
#define _Matrix_h
/****************************************************************************
 *
 *      Matrix.h  Contains prototypes for various matrix operation functions
 *
 ***************************************************************************/

void  MatrixMult(real8 *a, int aRows, int aCols, int aLD,
        real8 *b, int bCols, int bLD, real8 *c, int cLD);
void  MatrixMultArb(real8 *a, int aCols, int aRowOffset,
        int aColOffset, int aMultRows, int aMultCols,
        real8 *b, int bCols, int bRowOffset, int bColOffset, int bMultCols,
        real8 *c, int cCols, int cRowOffset, int cColOffset);
int   MatrixInvert(real8 *mat, real8 *invMat, int order, int lda);
real8 Matrix22Det(real8 a[2][2]);
void  Matrix22Invert(real8 A[2][2], real8 B[2][2]);
void  Matrix22Vector2Mult(real8 a[2][2], real8 b[2], real8 c[2]);
void  Matrix31Vector3Mult(real8 mat[3], real8 vec[3], real8 result[3][3]);
real8 Matrix33Det(real8 a[3][3]);
int   Matrix33Invert(real8 A[3][3], real8 B[3][3]);
void  Matrix33Mult33(real8 a[3][3], real8 b[3][3], real8 c[3][3]);
void  Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3]);
void  Matrix33Transpose(real8 mat[3][3], real8 trans[3][3]);
void  Matrix43Vector3Multiply(real8 A[4][3], real8 x[3], real8 y[4]);
void  Vec3TransposeAndMult(real8 vec[3], real8 result[3][3]);
void  Vector3Matrix33Mult(real8 vec[3], real8 mat[3][3], real8 result[3]);

#endif
