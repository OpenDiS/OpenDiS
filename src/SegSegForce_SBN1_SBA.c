#include "SegSegForce_SBN1.h"
#include "SegSegForce.h"
#include "SegmentStress.h"
#include <stdio.h>
#include <stdlib.h>

void SegSegForce_SBN1_SBA(real8 p1x, real8 p1y, real8 p1z,
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
                          real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
    real8 Rc, L12, L34;
    Rc = sqrt( (p1x+p2x-p3x-p4x)*(p1x+p2x-p3x-p4x) + (p1y+p2y-p3y-p4y)*(p1y+p2y-p3y-p4y) + (p1z+p2z-p3z-p4z)*(p1z+p2z-p3z-p4z) )/2.0;
    L12 = sqrt((p2x-p1x)*(p2x-p1x)+(p2y-p1y)*(p2y-p1y)+(p2z-p1z)*(p2z-p1z));
    L34 = sqrt((p4x-p3x)*(p4x-p3x)+(p4y-p3y)*(p4y-p3y)+(p4z-p3z)*(p4z-p3z));
    if (Rc < L12*3.0 || Rc < L34*3.0)
    {
        SegSegForce( p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z,
                     bpx, bpy, bpz, bx, by, bz, a, MU, NU, seg12Local, seg34Local,
                     fp1x, fp1y, fp1z, fp2x, fp2y, fp2z, fp3x, fp3y, fp3z, fp4x, fp4y, fp4z); 
    } 
    else
    {
        SegSegForce_SBN1( p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z,
                          bpx, bpy, bpz, bx, by, bz, a, MU, NU, Nint, quad_points, weights, seg12Local, seg34Local,
                          fp1x, fp1y, fp1z, fp2x, fp2y, fp2z, fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
    }
}