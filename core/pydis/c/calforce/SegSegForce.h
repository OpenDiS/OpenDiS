#include <math.h>
#define real8 double

int SegSegForce_BitReproMath(void);

void SpecialSegSegForce(real8 p1x, real8 p1y, real8 p1z,
                        real8 p2x, real8 p2y, real8 p2z,
                        real8 p3x, real8 p3y, real8 p3z,
                        real8 p4x, real8 p4y, real8 p4z,
                        real8 bpx, real8 bpy, real8 bpz,
                        real8 bx, real8 by, real8 bz,
                        real8 a, real8 MU, real8 NU, real8 ecrit,
                        int seg12Local, int seg34Local,
                        real8 *fp1x, real8 *fp1y, real8 *fp1z,
                        real8 *fp2x, real8 *fp2y, real8 *fp2z,
                        real8 *fp3x, real8 *fp3y, real8 *fp3z,
                        real8 *fp4x, real8 *fp4y, real8 *fp4z);


void SpecialSegSegForceHalf(real8 p1x, real8 p1y, real8 p1z,
                            real8 p2x, real8 p2y, real8 p2z,
                            real8 p3x, real8 p3y, real8 p3z,
                            real8 p4x, real8 p4y, real8 p4z,
                            real8 bpx, real8 bpy, real8 bpz,
                            real8 bx, real8 by, real8 bz,
                            real8 a, real8 MU, real8 NU, real8 ecrit,
                            real8 *fp3x, real8 *fp3y, real8 *fp3z,
                            real8 *fp4x, real8 *fp4y, real8 *fp4z);


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

/*
 *  SegSegForceList
 *
 *  Batched form of SegSegForce: evaluates n segment pairs in one call, looping over the
 *  existing scalar kernel internally.
 *
 *  SegSegForce takes 22 scalar doubles and 12 output pointers, so calling it from python
 *  costs one ctypes marshalling round trip per pair. That overhead was measured at ~8 us,
 *  against ~6 us of actual computation, so more than half the cost of a pydis force
 *  evaluation was spent crossing the language boundary rather than doing physics. This
 *  entry point pays it once per batch instead.
 *
 *  p1..p4, b12, b34 are n*3 arrays in row-major (x,y,z) order; f1..f4 are n*3 outputs.
 */
void SegSegForceList(int n,
                     const real8 *p1, const real8 *p2,
                     const real8 *p3, const real8 *p4,
                     const real8 *b12, const real8 *b34,
                     real8 a, real8 MU, real8 NU,
                     int seg12Local, int seg34Local,
                     real8 *f1, real8 *f2, real8 *f3, real8 *f4);

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
