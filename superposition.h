/*
 * superposition.h — DFT₂ (Hadamard) Transform
 *
 * The D=2 Fourier transform is the Hadamard gate:
 *   H = (1/√2) [[1,  1],
 *                [1, -1]]
 *
 * H|0⟩ = |+⟩ = (|0⟩+|1⟩)/√2
 * H|1⟩ = |−⟩ = (|0⟩-|1⟩)/√2
 *
 * H² = I  (self-inverse — the key difference from DFT₆ where F⁴=I)
 *
 * Precomputed twiddle table for generic DFT₂ application.
 */

#ifndef SUPERPOSITION_H
#define SUPERPOSITION_H

#include <math.h>

#define SUP_D  2

/* ═══════════════════════════════════════════════════════════════════════════════
 * DFT₂ TWIDDLE TABLE
 *
 * W[j][k] = ω^(j·k) where ω = e^(2πi/2) = -1
 *   W[0][0] = 1   W[0][1] = 1
 *   W[1][0] = 1   W[1][1] = -1
 *
 * Normalized by 1/√2.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static const double DFT2_RE[2][2] = {
    { 0.7071067811865475244,  0.7071067811865475244 },  /* 1/√2 × {1, 1}  */
    { 0.7071067811865475244, -0.7071067811865475244 }   /* 1/√2 × {1, -1} */
};

static const double DFT2_IM[2][2] = {
    { 0.0, 0.0 },
    { 0.0, 0.0 }
};

/* ═══════════════════════════════════════════════════════════════════════════════
 * APPLY HADAMARD — In-place DFT₂ of 2 complex amplitudes
 *
 *   out[0] = (in[0] + in[1]) / √2
 *   out[1] = (in[0] - in[1]) / √2
 *
 * Cost: 2 adds + 2 multiplies. FMA-friendly.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void sup_apply_hadamard(double *re, double *im)
{
    double r0 = re[0], r1 = re[1];
    double i0 = im[0], i1 = im[1];

    static const double INV_SQRT2 = 0.7071067811865475244;  /* 1/√2 */

    re[0] = (r0 + r1) * INV_SQRT2;
    im[0] = (i0 + i1) * INV_SQRT2;
    re[1] = (r0 - r1) * INV_SQRT2;
    im[1] = (i0 - i1) * INV_SQRT2;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * INVERSE HADAMARD — H† = H (Hadamard is self-adjoint)
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void sup_apply_hadamard_inv(double *re, double *im)
{
    sup_apply_hadamard(re, im);  /* H† = H */
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * RENORMALIZE — Force total probability to 1.0
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void sup_renormalize(double *re, double *im, int D)
{
    double norm2 = 0;
    for (int k = 0; k < D; k++)
        norm2 += re[k] * re[k] + im[k] * im[k];

    if (norm2 > 1e-30 && fabs(norm2 - 1.0) > 1e-15) {
        double scale = 1.0 / sqrt(norm2);
        for (int k = 0; k < D; k++) {
            re[k] *= scale;
            im[k] *= scale;
        }
    }
}

#endif /* SUPERPOSITION_H */
