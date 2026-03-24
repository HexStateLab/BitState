/*
 * arithmetic.h — IEEE-754 Constants & Fast Floating-Point Primitives
 *
 * Reverse-engineered constants from the physical substrate's FPU.
 * Direct port from HexState — these constants are architecture-level,
 * independent of qudit dimension.
 */

#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include <stdint.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * IEEE-754 DOUBLE LAYOUT
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define IEEE754_SIGN_MASK     0x8000000000000000ULL
#define IEEE754_EXP_MASK      0x7FF0000000000000ULL
#define IEEE754_MANT_MASK     0x000FFFFFFFFFFFFFULL
#define IEEE754_EXP_BIAS      1023
#define IEEE754_MANT_BITS     52

/* ═══════════════════════════════════════════════════════════════════════════════
 * MAGIC NUMBERS — Fast approximate operations via bit manipulation
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Quake III fast inverse sqrt magic constant */
#define MAGIC_ISQRT_F32       0x5F3759DF
#define MAGIC_ISQRT_F64       0x5FE6EB50C7B537A9ULL

/* Fast reciprocal */
#define MAGIC_RECIP_F64       0x7FDE623822FC16E6ULL

/* Fast sqrt */
#define MAGIC_SQRT_F64        0x1FF7A7EF9DB22D0EULL

/* ═══════════════════════════════════════════════════════════════════════════════
 * SUBSTRATE CONSTANTS
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define ARITH_PHI             1.6180339887498948482   /* Golden ratio         */
#define ARITH_PHI_INV         0.6180339887498948482   /* 1/φ                  */
#define ARITH_SQRT2           1.4142135623730950488   /* √2                   */
#define ARITH_SQRT2_INV       0.7071067811865475244   /* 1/√2 = √2/2         */
#define ARITH_DOTTIE          0.7390851332151606416   /* Fixed point of cos   */

/* ═══════════════════════════════════════════════════════════════════════════════
 * FAST OPERATIONS — Bit-level IEEE-754 hacks
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Fast approximate inverse sqrt: 1/√x
 * One Newton-Raphson iteration for ~23 bits of precision. */
static inline double arith_fast_isqrt(double x)
{
    union { double d; uint64_t u; } v;
    v.d = x;
    v.u = MAGIC_ISQRT_F64 - (v.u >> 1);
    /* Newton-Raphson: y = y * (1.5 - 0.5 * x * y * y) */
    double y = v.d;
    y = y * (1.5 - 0.5 * x * y * y);
    y = y * (1.5 - 0.5 * x * y * y);  /* Second iteration: ~46 bits */
    return y;
}

/* Fast approximate reciprocal: 1/x */
static inline double arith_fast_recip(double x)
{
    union { double d; uint64_t u; } v;
    v.d = x;
    v.u = MAGIC_RECIP_F64 - v.u;
    double y = v.d;
    y = y * (2.0 - x * y);  /* Newton */
    y = y * (2.0 - x * y);
    return y;
}

#endif /* ARITHMETIC_H */
