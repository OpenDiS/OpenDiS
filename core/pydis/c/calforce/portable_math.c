/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft/SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 *
 * pydis_log() and pydis_atan() below are fdlibm's e_log.c and s_atan.c
 * (as carried by FreeBSD's lib/msun/src/), transcribed with the word-
 * extraction macros (GET_HIGH_WORD etc., from fdlibm's math_private.h)
 * replaced by the memcpy-based helpers just below, and u_int32_t
 * replaced by the standard uint32_t. The algorithms, comments and
 * every hexadecimal/decimal constant are unchanged from the original;
 * see portable_math.h for why this file exists at all.
 */

#include "portable_math.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "pydis_log/pydis_atan assume a little-endian double representation"
#endif

static inline void extract_words(int32_t *hi, uint32_t *lo, double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    *lo = (uint32_t)bits;
    *hi = (int32_t)(bits >> 32);
}

static inline int32_t get_high_word(double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (int32_t)(bits >> 32);
}

static inline uint32_t get_low_word(double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (uint32_t)bits;
}

static inline void set_high_word(double *x, uint32_t hi)
{
    uint64_t bits;
    memcpy(&bits, x, sizeof(bits));
    bits = ((uint64_t)hi << 32) | (bits & 0xffffffffu);
    memcpy(x, &bits, sizeof(bits));
}

/* ---- fdlibm e_log.c: log(x) ---- */

static const double
log_ln2_hi =  6.93147180369123816490e-01, /* 3fe62e42 fee00000 */
log_ln2_lo =  1.90821492927058770002e-10, /* 3dea39ef 35793c76 */
log_two54  =  1.80143985094819840000e+16, /* 43500000 00000000 */
Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

static const double log_zero = 0.0;
static volatile double log_vzero = 0.0;

double pydis_log(double x)
{
    double hfsq, f, s, z, R, w, t1, t2, dk;
    int32_t k, hx, i, j;
    uint32_t lx;

    extract_words(&hx, &lx, x);

    k = 0;
    if (hx < 0x00100000) {              /* x < 2**-1022  */
        if (((hx & 0x7fffffff) | lx) == 0)
            return -log_two54 / log_vzero;     /* log(+-0)=-inf */
        if (hx < 0) return (x - x) / log_zero; /* log(-#) = NaN */
        k -= 54; x *= log_two54;    /* subnormal number, scale up x */
        hx = get_high_word(x);
    }
    if (hx >= 0x7ff00000) return x + x;
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    i = (hx + 0x95f64) & 0x100000;
    set_high_word(&x, hx | (i ^ 0x3ff00000));  /* normalize x or x/2 */
    k += (i >> 20);
    f = x - 1.0;
    if ((0x000fffff & (2 + hx)) < 3) {  /* -2**-20 <= f < 2**-20 */
        if (f == log_zero) {
            if (k == 0) {
                return log_zero;
            } else {
                dk = (double)k;
                return dk * log_ln2_hi + dk * log_ln2_lo;
            }
        }
        R = f * f * (0.5 - 0.33333333333333333 * f);
        if (k == 0) return f - R; else {
            dk = (double)k;
            return dk * log_ln2_hi - ((R - dk * log_ln2_lo) - f);
        }
    }
    s = f / (2.0 + f);
    dk = (double)k;
    z = s * s;
    i = hx - 0x6147a;
    w = z * z;
    j = 0x6b851 - hx;
    t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
    t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
    i |= j;
    R = t2 + t1;
    if (i > 0) {
        hfsq = 0.5 * f * f;
        if (k == 0) return f - (hfsq - s * (hfsq + R)); else
            return dk * log_ln2_hi - ((hfsq - (s * (hfsq + R) + dk * log_ln2_lo)) - f);
    } else {
        if (k == 0) return f - s * (f - R); else
            return dk * log_ln2_hi - ((s * (f - R) - dk * log_ln2_lo) - f);
    }
}

/* ---- fdlibm s_atan.c: atan(x) ---- */

static const double atanhi[] = {
  4.63647609000806093515e-01, /* atan(0.5)hi 0x3FDDAC67, 0x0561BB4F */
  7.85398163397448278999e-01, /* atan(1.0)hi 0x3FE921FB, 0x54442D18 */
  9.82793723247329054082e-01, /* atan(1.5)hi 0x3FEF730B, 0xD281F69B */
  1.57079632679489655800e+00, /* atan(inf)hi 0x3FF921FB, 0x54442D18 */
};

static const double atanlo[] = {
  2.26987774529616870924e-17, /* atan(0.5)lo 0x3C7A2B7F, 0x222F65E2 */
  3.06161699786838301793e-17, /* atan(1.0)lo 0x3C81A626, 0x33145C07 */
  1.39033110312309984516e-17, /* atan(1.5)lo 0x3C700788, 0x7AF0CBBD */
  6.12323399573676603587e-17, /* atan(inf)lo 0x3C91A626, 0x33145C07 */
};

static const double aT[] = {
  3.33333333333329318027e-01, /* 0x3FD55555, 0x5555550D */
 -1.99999999998764832476e-01, /* 0xBFC99999, 0x9998EBC4 */
  1.42857142725034663711e-01, /* 0x3FC24924, 0x920083FF */
 -1.11111104054623557880e-01, /* 0xBFBC71C6, 0xFE231671 */
  9.09088713343650656196e-02, /* 0x3FB745CD, 0xC54C206E */
 -7.69187620504482999495e-02, /* 0xBFB3B0F2, 0xAF749A6D */
  6.66107313738753120669e-02, /* 0x3FB10D66, 0xA0D03D51 */
 -5.83357013379057348645e-02, /* 0xBFADDE2D, 0x52DEFD9A */
  4.97687799461593236017e-02, /* 0x3FA97B4B, 0x24760DEB */
 -3.65315727442169155270e-02, /* 0xBFA2B444, 0x2C6A6C2F */
  1.62858201153657823623e-02, /* 0x3F90AD3A, 0xE322DA11 */
};

static const double atan_one  = 1.0,
                     atan_huge = 1.0e300;

double pydis_atan(double x)
{
    double w, s1, s2, z;
    int32_t ix, hx, id;

    hx = get_high_word(x);
    ix = hx & 0x7fffffff;
    if (ix >= 0x44100000) {    /* if |x| >= 2^66 */
        uint32_t low = get_low_word(x);
        if (ix > 0x7ff00000 || (ix == 0x7ff00000 && (low != 0)))
            return x + x;      /* NaN */
        if (hx > 0) return  atanhi[3] + *(volatile double *)&atanlo[3];
        else        return -atanhi[3] - *(volatile double *)&atanlo[3];
    } if (ix < 0x3fdc0000) {   /* |x| < 0.4375 */
        if (ix < 0x3e400000) { /* |x| < 2^-27 */
            if (atan_huge + x > atan_one) return x;    /* raise inexact */
        }
        id = -1;
    } else {
    x = fabs(x);
    if (ix < 0x3ff30000) {     /* |x| < 1.1875 */
        if (ix < 0x3fe60000) { /* 7/16 <=|x|<11/16 */
            id = 0; x = (2.0 * x - atan_one) / (2.0 + x);
        } else {                /* 11/16<=|x|< 19/16 */
            id = 1; x  = (x - atan_one) / (x + atan_one);
        }
    } else {
        if (ix < 0x40038000) { /* |x| < 2.4375 */
            id = 2; x  = (x - 1.5) / (atan_one + 1.5 * x);
        } else {                /* 2.4375 <= |x| < 2^66 */
            id = 3; x  = -1.0 / x;
        }
    }}
    /* end of argument reduction */
    z = x * x;
    w = z * z;
    /* break sum from i=0 to 10 aT[i]z**(i+1) into odd and even poly */
    s1 = z * (aT[0] + w * (aT[2] + w * (aT[4] + w * (aT[6] + w * (aT[8] + w * aT[10])))));
    s2 = w * (aT[1] + w * (aT[3] + w * (aT[5] + w * (aT[7] + w * aT[9]))));
    if (id < 0) return x - x * (s1 + s2);
    else {
        z = atanhi[id] - ((x * (s1 + s2) - atanlo[id]) - x);
        return (hx < 0) ? -z : z;
    }
}
