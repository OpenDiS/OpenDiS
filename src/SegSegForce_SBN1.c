#include "SegSegForce_SBN1.h"

/* copy ParaDiS function for computing stress of a straight segment */


void SegSegForceHalf_SBN1(real8 p1x, real8 p1y, real8 p1z,
                            real8 p2x, real8 p2y, real8 p2z,
                            real8 p3x, real8 p3y, real8 p3z,
                            real8 p4x, real8 p4y, real8 p4z,
                            real8 bpx, real8 bpy, real8 bpz,
                            real8 bx, real8 by, real8 bz,
                            real8 a, real8 MU, real8 NU, int Nint,
                            real8 *fp3x, real8 *fp3y, real8 *fp3z,
                            real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
    /* hard code integration points and weights for Nint = 1,2,3 */

    /* calculate the stress field of segment 1-2 */

    /* numerically integrate Peach-Koehler force on segment 3-4 */
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
                          real8 a, real8 MU, real8 NU, int Nint,
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
                                     a, MU, NU, Nint,
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
                                     a, MU, NU, Nint,
                                     fp1x, fp1y, fp1z, fp2x, fp2y, fp2z);
            } /* if segment p1->p2 is "local" */

       return;
}