#include <math.h>
#define real8 double

void StressDueToSeg(real8 px, real8 py, real8 pz,
                    real8 p1x, real8 p1y, real8 p1z,
                    real8 p2x, real8 p2y, real8 p2z,
                    real8 bx, real8 by, real8 bz,
                    real8 a, real8 MU, real8 NU, real8 *stress);