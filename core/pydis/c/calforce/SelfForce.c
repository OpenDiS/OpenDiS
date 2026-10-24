/* Translated from ParaDiS SelfForceIsotropic (NodeForce.c:2628) and SelfForce
 * (:2688). Every expression and every grouping is as the C has it, including
 * be*(S+ft) + fL*t, which adds S and ft before multiplying rather than scaling
 * be twice. That grouping is not what ExaDiS does, and the two differ by up to
 * 3.7e-09 on realistic segments; see tests/unit_tests/test1_node_force.
 *
 * The _NONSINGULAR_SELF_FORCE branch of the original is not reproduced: it is
 * commented out in ParaDiS' own makefile.setup, so the default build takes the
 * GetUnitVector path below, dividing by L rather than by La. ExaDiS takes the
 * same path.
 *
 * Signature differs from ParaDiS' only in returning six pointers rather than
 * two real8[3], which is what ctypesgen and numpy work with comfortably.
 */
#include "SelfForce.h"

void SelfForceIsotropic(int coreOnly, real8 MU, real8 NU,
                        real8 bx, real8 by, real8 bz,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 a, real8 Ecore,
                        real8 *f1x, real8 *f1y, real8 *f1z,
                        real8 *f2x, real8 *f2y, real8 *f2z)
{
        real8 tx, ty, tz, L, La, S;
        real8 bs, bs2, bex, bey, bez, be2, fL, ft;
        real8 ux, uy, uz;

        /* GetUnitVector(1, ...) with unitFlag 1: the length, then the
         * direction divided componentwise by it */
        ux = x2 - x1;
        uy = y2 - y1;
        uz = z2 - z1;
        L = sqrt(ux*ux + uy*uy + uz*uz);
        tx = ux / L;
        ty = uy / L;
        tz = uz / L;
        La = sqrt(L*L + a*a);

        bs = bx*tx + by*ty + bz*tz;
        bex = bx - bs*tx;  bey = by - bs*ty;  bez = bz - bs*tz;
        be2 = (bex*bex + bey*bey + bez*bez);
        bs2 = bs*bs;

        if (coreOnly) {
                S = 0.0;
        } else {
                S = (-(2*NU*La+(1-NU)*a*a/La-(1+NU)*a)/L +
                     (NU*log((La+L)/a)-(1-NU)*0.5*L/La))*MU/4/M_PI/(1-NU)*bs;
        }

        /* Ecore = MU/(4*pi) log(a/a0); the self force from the core energy
         * changing as the core radius goes from a0 to a */
        fL = -Ecore*(bs2 + be2/(1-NU));
        ft =  Ecore*2*bs*NU/(1-NU);

        /* CHANGED, and deliberately not what ParaDiS writes.
         *
         * ParaDiS assembles this as
         *
         *     *f2x = bex*(S+ft) + fL*tx;
         *     *f2y = bey*(S+ft) + fL*ty;
         *     *f2z = bez*(S+ft) + fL*tz;
         *
         * adding S and ft before scaling be. ExaDiS reaches the same algebra by
         * a different route: core_force returns -(ft*be + fL*t), self_force
         * returns -(S*be), and force_lt.h adds them, which is
         * (ft*be + fL*t) + S*be with the parenthesization below. The two differ
         * by up to 3.7e-09 on realistic segments with a non-zero Ecore.
         *
         * The grouping below is used because it is what makes pydis and ExaDiS
         * agree bitwise, which is the point of the _repro builds. It is a
         * departure from the reference C, so it is spelled out rather than left
         * to look like a transcription slip. Restore the three lines above to go
         * back to ParaDiS' own rounding.
         */
        *f2x = (ft*bex + fL*tx) + S*bex;
        *f2y = (ft*bey + fL*ty) + S*bey;
        *f2z = (ft*bez + fL*tz) + S*bez;

        *f1x = -*f2x;
        *f1y = -*f2y;
        *f1z = -*f2z;
}

void SelfForce(int coreOnly, real8 MU, real8 NU,
               real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1,
               real8 x2, real8 y2, real8 z2,
               real8 a, real8 Ecore,
               real8 *f1x, real8 *f1y, real8 *f1z,
               real8 *f2x, real8 *f2y, real8 *f2z)
{
        real8 dx, dy, dz, len2;

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;
        len2 = dx*dx + dy*dy + dz*dz;

        if (len2 < 1.0e-20) {
                *f1x = 0.0;  *f1y = 0.0;  *f1z = 0.0;
                *f2x = 0.0;  *f2y = 0.0;  *f2z = 0.0;
                return;
        }

        SelfForceIsotropic(coreOnly, MU, NU, bx, by, bz, x1, y1, z1,
                           x2, y2, z2, a, Ecore,
                           f1x, f1y, f1z, f2x, f2y, f2z);
}

/* One call per segment list, so python pays one crossing rather than nseg of
 * them. fseg is nseg rows of six: f1 then f2. */
void SelfForceList(int coreOnly, real8 MU, real8 NU, real8 a, real8 Ecore,
                   int nseg, real8 *burg, real8 *p1, real8 *p2, real8 *fseg)
{
        int i;
        for (i = 0; i < nseg; i++) {
                SelfForce(coreOnly, MU, NU,
                          burg[3*i], burg[3*i+1], burg[3*i+2],
                          p1[3*i], p1[3*i+1], p1[3*i+2],
                          p2[3*i], p2[3*i+1], p2[3*i+2],
                          a, Ecore,
                          &fseg[6*i],   &fseg[6*i+1], &fseg[6*i+2],
                          &fseg[6*i+3], &fseg[6*i+4], &fseg[6*i+5]);
        }
}
