/* pydis_log / pydis_atan: platform-independent log() and atan()
 *
 * The system libm's log() and atan() are not required by the C standard
 * or IEEE 754 to be correctly rounded (unlike +,-,*,/,sqrt), so glibc and
 * Apple's libm can legitimately return answers a few ULP apart for the
 * same input. SegSegForce.c leans on log/atan inside a near-parallel
 * decomposition that cancels terms 1-2 orders of magnitude larger than
 * the final force, which turns that ULP-level libm disagreement into a
 * measurable (though still tiny) difference in the force computed on
 * different machines. These two functions replace the OS-supplied log
 * and atan with a fixed, portable implementation so the same binary
 * arithmetic runs on every platform.
 */
#ifndef PYDIS_PORTABLE_MATH_H
#define PYDIS_PORTABLE_MATH_H

double pydis_log(double x);
double pydis_atan(double x);

#endif
