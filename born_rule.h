/*
 * born_rule.h — Born Rule Probability & Sampling for D=2
 *
 * Born rule: P(k) = |⟨k|ψ⟩|² = re[k]² + im[k]²
 *
 * For D=2, sampling is trivial: compute P(0), if rand < P(0) → outcome 0,
 * else → outcome 1. No CDF scan needed.
 *
 * Includes fast inverse sqrt for renormalization (from HexState arithmetic.h).
 */

#ifndef BORN_RULE_H
#define BORN_RULE_H

#include <math.h>
#include <stdint.h>
#include "arithmetic.h"

#define BORN_D  2

/* ═══════════════════════════════════════════════════════════════════════════════
 * BORN PROBABILITY — |amplitude|²
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Exact probability */
static inline double born_prob(double re, double im)
{
    return re * re + im * im;
}

/* Total probability of a D=2 state */
static inline double born_total_prob(const double *re, const double *im, int D)
{
    double total = 0;
    for (int k = 0; k < D; k++)
        total += re[k] * re[k] + im[k] * im[k];
    return total;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BORN SAMPLING — For D=2, just compare P(0) against random
 *
 * Returns outcome 0 or 1.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline int born_sample(const double *re, const double *im,
                              int D, double rand_01)
{
    (void)D;  /* Always 2 */
    double p0 = re[0] * re[0] + im[0] * im[0];
    return (rand_01 < p0) ? 0 : 1;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BORN COLLAPSE — Post-measurement state update
 *
 * After measuring outcome k:
 *   ψ[k] → 1 (with original phase)
 *   ψ[¬k] → 0
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void born_collapse(double *re, double *im, int D, int outcome)
{
    double mag = sqrt(re[outcome] * re[outcome] + im[outcome] * im[outcome]);
    if (mag < 1e-30) mag = 1.0;  /* fallback */

    double inv_mag = 1.0 / mag;
    double phase_re = re[outcome] * inv_mag;
    double phase_im = im[outcome] * inv_mag;

    for (int k = 0; k < D; k++) {
        re[k] = 0.0;
        im[k] = 0.0;
    }
    re[outcome] = phase_re;
    im[outcome] = phase_im;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * FAST INVERSE SQRT — For renormalization
 *
 * Uses Quake III magic constant with two Newton iterations.
 * ~46 bits of precision, ~2ns.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double born_fast_isqrt(double x)
{
    return arith_fast_isqrt(x);
}

#endif /* BORN_RULE_H */
