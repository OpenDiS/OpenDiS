/* SelfForce: the self force on one dislocation segment
 *
 * ParaDiS SelfForceIsotropic and its SelfForce wrapper, NodeForce.c:2628 and
 * :2688. One function covers what ExaDiS splits between self_force and
 * core_force in force_types/force_common.h and force_core.h: the S term is
 * ExaDiS' self_force, the fL and ft terms are its core_force.
 *
 * coreOnly = 1 drops the S term and leaves only the core-energy part, which is
 * what ParaDiS uses where the elastic self force is accounted for elsewhere.
 */
#ifndef PYDIS_SELF_FORCE_H
#define PYDIS_SELF_FORCE_H

#include <math.h>
#define real8 double

void SelfForceIsotropic(int coreOnly, real8 MU, real8 NU,
                        real8 bx, real8 by, real8 bz,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 a, real8 Ecore,
                        real8 *f1x, real8 *f1y, real8 *f1z,
                        real8 *f2x, real8 *f2y, real8 *f2z);

void SelfForce(int coreOnly, real8 MU, real8 NU,
               real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1,
               real8 x2, real8 y2, real8 z2,
               real8 a, real8 Ecore,
               real8 *f1x, real8 *f1y, real8 *f1z,
               real8 *f2x, real8 *f2y, real8 *f2z);

void SelfForceList(int coreOnly, real8 MU, real8 NU, real8 a, real8 Ecore,
                   int nseg, real8 *burg, real8 *p1, real8 *p2, real8 *fseg);

#endif
