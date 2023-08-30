#include <math.h>
#define real8 double

void SegmentStress(real8 MU, real8 NU,
                   real8 bX, real8 bY, real8 bZ,
                   real8 xA, real8 yA, real8 zA,
                   real8 xB, real8 yB, real8 zB,
                   real8 x0, real8 y0, real8 z0,
                   real8 a,
                   real8 Sigma[3][3] );
