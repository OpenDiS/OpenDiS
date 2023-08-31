#include "SegSegForce_SBN1.h"
#include "SegmentStress.h"
#include <stdio.h>
#include <stdlib.h>

/* copy ParaDiS function for computing stress of a straight segment */

/* from DD_Force.git/python/utils/compute_stress_force_numerical.py */

#define MAX_QUAD_POINTS 7

real8 Shape_Func_N1(real8 x)
{
    return 0.5 * (1.0 - x);
}

real8 Shape_Func_N2(real8 x)
{
    return 0.5 * (1.0 + x);
}

void SegSegForceHalf_SBN1(real8 p1x, real8 p1y, real8 p1z,
                            real8 p2x, real8 p2y, real8 p2z,
                            real8 p3x, real8 p3y, real8 p3z,
                            real8 p4x, real8 p4y, real8 p4z,
                            real8 bpx, real8 bpy, real8 bpz,
                            real8 bx, real8 by, real8 bz,
                            real8 a, real8 MU, real8 NU,
                            int Nint, real8 *quad_points, real8 *weights,
                            real8 *fp3x, real8 *fp3y, real8 *fp3z,
                            real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
    int i; 
    real8 N1_quad_points[MAX_QUAD_POINTS], N2_quad_points[MAX_QUAD_POINTS];
    real8 p34_half_x, p34_half_y, p34_half_z, p34_mid_x, p34_mid_y, p34_mid_z;
    real8 px, py, pz, sigma[3][3], sigb[3];

    if (Nint > MAX_QUAD_POINTS) {
        fprintf(stderr, "Nint > %d\n", MAX_QUAD_POINTS);
        exit(1);
    }
    p34_half_x = 0.5*(p4x-p3x); p34_half_y = 0.5*(p4y-p3y); p34_half_z = 0.5*(p4z-p3z);
    p34_mid_x  = 0.5*(p4x+p3x); p34_mid_y  = 0.5*(p4y+p3y); p34_mid_z  = 0.5*(p4z+p3z);

    for(i = 0; i < Nint; i++)
    {
        N1_quad_points[i] = Shape_Func_N1(quad_points[i]);
        N2_quad_points[i] = Shape_Func_N2(quad_points[i]);
    }
    *fp3x = 0.0; *fp3y = 0.0; *fp3z = 0.0;
    for(i = 0; i < Nint; i++)
    {
        px = p34_mid_x + p34_half_x * quad_points[i];
        py = p34_mid_y + p34_half_y * quad_points[i];
        pz = p34_mid_z + p34_half_z * quad_points[i];

        SegmentStress(MU, NU, bpx, bpy, bpz, p1x, p1y, p1z, p2x, p2y, p2z,
                      px, py, pz, a, sigma);

        sigb[0] = sigma[0][0] * bx + sigma[0][1] * by + sigma[0][2] * bz;
        sigb[1] = sigma[1][0] * bx + sigma[1][1] * by + sigma[1][2] * bz;
        sigb[2] = sigma[2][0] * bx + sigma[2][1] * by + sigma[2][2] * bz;

        *fp3x += (sigb[1]*p34_half_z - sigb[2]*p34_half_y) * N1_quad_points[i] * weights[i];
        *fp3y += (sigb[2]*p34_half_x - sigb[0]*p34_half_z) * N1_quad_points[i] * weights[i];
        *fp3z += (sigb[0]*p34_half_y - sigb[1]*p34_half_x) * N1_quad_points[i] * weights[i];

        *fp4x += (sigb[1]*p34_half_z - sigb[2]*p34_half_y) * N2_quad_points[i] * weights[i];
        *fp4y += (sigb[2]*p34_half_x - sigb[0]*p34_half_z) * N2_quad_points[i] * weights[i];
        *fp4z += (sigb[0]*p34_half_y - sigb[1]*p34_half_x) * N2_quad_points[i] * weights[i];
    }
}

/*-------------------------------------------------------------------------
 *
 *      Function:       SegSegForce_SBN1
 *      Description:    Used to calculate the interaction forces between
 *                      dislocation segments by numerically integrating the
 *                      stress field of one segment over the other segment.
 *
 *      Arguments:
 *              p1*,p2*      endpoints for first dislocation segment starting
 *                           at p1x,p1y,p1z and ending at p2x,p2y,p2z
 *              p3*,p4*      endpoints for seond dislocation segment starting
 *                           at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              seg12Local   1 if either node of segment p1->p2 is local to
 *                           the current domain, zero otherwise.
 *              seg34Local   1 if either node of segment p3->p4 is local to
 *                           the current domain, zero otherwise.
 *              fp1*,fp2*,   pointers to locations in which to return
 *              fp3*,fp4*    forces on nodes located at p1, p2, p3 and
 *                           p4 respectively
 *
 *-----------------------------------------------------------------------*/
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
                          real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
/*
 *          Only calculate the forces for segment p3->p4 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            if (seg34Local) {
                SegSegForceHalf_SBN1(p1x, p1y, p1z, p2x, p2y, p2z,
                                     p3x, p3y, p3z, p4x, p4y, p4z,
                                     bpx, bpy, bpz, bx, by, bz,
                                     a, MU, NU, Nint, quad_points, weights,
                                     fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);

            } /* if segment p3->p4 is "local" */

/*
 *          Only calculate the forces for segment p1->p2 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            if (seg12Local) {
                SegSegForceHalf_SBN1(p3x, p3y, p3z, p4x, p4y, p4z,
                                     p1x, p1y, p1z, p2x, p2y, p2z,
                                     bx, by, bz, bpx, bpy, bpz,
                                     a, MU, NU, Nint, quad_points, weights,
                                     fp1x, fp1y, fp1z, fp2x, fp2y, fp2z);
            } /* if segment p1->p2 is "local" */

       return;
}