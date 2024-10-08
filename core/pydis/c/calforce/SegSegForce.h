#include <math.h>
#define real8 double

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
