/*
 * entanglement.h — D=2 Joint State Storage & Entanglement Analysis
 *
 * A JointState stores 4 complex amplitudes for two entangled qubits:
 *   ψ(a,b) where a,b ∈ {0,1}
 *   idx = a * 2 + b
 *   |00⟩=0, |01⟩=1, |10⟩=2, |11⟩=3
 *
 * Memory: 64 bytes per pair (vs 576 bytes for D=6 HexState).
 */

#ifndef ENTANGLEMENT_H
#define ENTANGLEMENT_H

#include <math.h>
#include <string.h>
#include "statevector.h"
#include "superposition.h"

#define ENT_D   2
#define ENT_D2  4

/* ═══════════════════════════════════════════════════════════════════════════════
 * BELL STATES — (1/√2)(|00⟩ + |11⟩) = maximally entangled
 *
 * The standard Bell pair for qubits. 64 bytes.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void ent_bell(JointState *j)
{
    memset(j, 0, sizeof(*j));
    static const double INV_SQRT2 = 0.7071067811865475244;
    j->re[0] = INV_SQRT2;   /* |00⟩ */
    j->re[3] = INV_SQRT2;   /* |11⟩ */
}

/* Bell Φ⁻: (|00⟩ - |11⟩)/√2 */
static inline void ent_bell_minus(JointState *j)
{
    memset(j, 0, sizeof(*j));
    static const double INV_SQRT2 = 0.7071067811865475244;
    j->re[0] =  INV_SQRT2;  /* |00⟩ */
    j->re[3] = -INV_SQRT2;  /* |11⟩ */
}

/* Bell Ψ⁺: (|01⟩ + |10⟩)/√2 */
static inline void ent_bell_psi_plus(JointState *j)
{
    memset(j, 0, sizeof(*j));
    static const double INV_SQRT2 = 0.7071067811865475244;
    j->re[1] = INV_SQRT2;   /* |01⟩ */
    j->re[2] = INV_SQRT2;   /* |10⟩ */
}

/* Bell Ψ⁻: (|01⟩ - |10⟩)/√2 — the singlet */
static inline void ent_bell_psi_minus(JointState *j)
{
    memset(j, 0, sizeof(*j));
    static const double INV_SQRT2 = 0.7071067811865475244;
    j->re[1] =  INV_SQRT2;  /* |01⟩ */
    j->re[2] = -INV_SQRT2;  /* |10⟩ */
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRODUCT STATE — |ψ_a⟩ ⊗ |ψ_b⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void ent_product(JointState *j,
                               const QubitState *sa, const QubitState *sb)
{
    for (int a = 0; a < ENT_D; a++)
    for (int b = 0; b < ENT_D; b++) {
        int idx = a * ENT_D + b;
        j->re[idx] = sa->re[a] * sb->re[b] - sa->im[a] * sb->im[b];
        j->im[idx] = sa->re[a] * sb->im[b] + sa->im[a] * sb->re[b];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PARTIAL TRACE — Tr_B(ρ) → marginal probabilities for qubit A
 *
 * P_A(a) = Σ_b |ψ(a,b)|²
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void ent_marginal_a(const JointState *j, double *probs)
{
    for (int a = 0; a < ENT_D; a++) {
        double sum = 0;
        for (int b = 0; b < ENT_D; b++) {
            int idx = a * ENT_D + b;
            sum += j->re[idx] * j->re[idx] + j->im[idx] * j->im[idx];
        }
        probs[a] = sum;
    }
}

/* P_B(b) = Σ_a |ψ(a,b)|² */
static inline void ent_marginal_b(const JointState *j, double *probs)
{
    for (int b = 0; b < ENT_D; b++) {
        double sum = 0;
        for (int a = 0; a < ENT_D; a++) {
            int idx = a * ENT_D + b;
            sum += j->re[idx] * j->re[idx] + j->im[idx] * j->im[idx];
        }
        probs[b] = sum;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SCHMIDT RANK — Number of nonzero singular values
 *
 * For D=2: rank is 1 (separable) or 2 (entangled).
 * Uses the 2×2 determinant trick: if |det(ψ-matrix)| > ε, rank = 2.
 *
 * ψ as 2×2 matrix:  [[ψ00, ψ01],
 *                     [ψ10, ψ11]]
 * det = ψ00·ψ11 - ψ01·ψ10
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline int ent_schmidt_rank(const JointState *j)
{
    /* Complex determinant: det = ψ00·ψ11 - ψ01·ψ10 */
    double det_re = (j->re[0] * j->re[3] - j->im[0] * j->im[3])
                  - (j->re[1] * j->re[2] - j->im[1] * j->im[2]);
    double det_im = (j->re[0] * j->im[3] + j->im[0] * j->re[3])
                  - (j->re[1] * j->im[2] + j->im[1] * j->re[2]);
    double det_mag2 = det_re * det_re + det_im * det_im;

    return (det_mag2 > 1e-20) ? 2 : 1;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENTANGLEMENT ENTROPY — S = -Σ λ_k log₂(λ_k)
 *
 * For D=2, the reduced density matrix ρ_A has at most 2 eigenvalues.
 * Uses the 2×2 determinant to compute concurrence, then entropy.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double ent_entropy(const JointState *j)
{
    double probs[ENT_D];
    ent_marginal_a(j, probs);

    double S = 0;
    for (int k = 0; k < ENT_D; k++) {
        if (probs[k] > 1e-14)
            S -= probs[k] * log2(probs[k]);
    }
    return S;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TOTAL PROBABILITY — Sanity check (must equal 1.0)
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double ent_total_prob(const JointState *j)
{
    double sum = 0;
    for (int i = 0; i < ENT_D2; i++)
        sum += j->re[i] * j->re[i] + j->im[i] * j->im[i];
    return sum;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * RENORMALIZE — Force total probability to 1.0
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void ent_renormalize(JointState *j)
{
    double norm2 = ent_total_prob(j);
    if (norm2 > 1e-30 && fabs(norm2 - 1.0) > 1e-15) {
        double scale = 1.0 / sqrt(norm2);
        for (int i = 0; i < ENT_D2; i++) {
            j->re[i] *= scale;
            j->im[i] *= scale;
        }
    }
}

#endif /* ENTANGLEMENT_H */
