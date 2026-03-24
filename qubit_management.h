/*
 * qubit_management.h — Per-Qubit State Management Primitives
 *
 * Low-level operations on individual QubitState and JointState objects.
 * Used by the engine's entangle/disentangle routines.
 */

#ifndef QUBIT_MANAGEMENT_H
#define QUBIT_MANAGEMENT_H

#include <math.h>
#include <string.h>
#include "statevector.h"

#define QM_D  2

/* ═══════════════════════════════════════════════════════════════════════════════
 * INITIALIZATION
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Init to |0⟩ */
static inline void qm_init_zero(QubitState *s)
{
    memset(s, 0, sizeof(*s));
    s->re[0] = 1.0;
}

/* Init to |+⟩ = (|0⟩+|1⟩)/√2 */
static inline void qm_init_plus(QubitState *s)
{
    memset(s, 0, sizeof(*s));
    static const double INV_SQRT2 = 0.7071067811865475244;
    s->re[0] = INV_SQRT2;
    s->re[1] = INV_SQRT2;
}

/* Init to |k⟩ */
static inline void qm_init_basis(QubitState *s, int k)
{
    memset(s, 0, sizeof(*s));
    if (k >= 0 && k < QM_D)
        s->re[k] = 1.0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENTANGLEMENT PRIMITIVES
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Create Bell state (|00⟩+|11⟩)/√2 in a JointState */
static inline void qm_entangle_bell(JointState *j)
{
    memset(j, 0, sizeof(*j));
    static const double INV_SQRT2 = 0.7071067811865475244;
    j->re[0] = INV_SQRT2;   /* |00⟩ */
    j->re[3] = INV_SQRT2;   /* |11⟩ */
}

/* Create product state |ψ_a⟩⊗|ψ_b⟩ in a JointState */
static inline void qm_entangle_product(JointState *j,
                                       const QubitState *sa,
                                       const QubitState *sb)
{
    for (int a = 0; a < QM_D; a++)
    for (int b = 0; b < QM_D; b++) {
        int idx = a * QM_D + b;
        j->re[idx] = sa->re[a] * sb->re[b] - sa->im[a] * sb->im[b];
        j->im[idx] = sa->re[a] * sb->im[b] + sa->im[a] * sb->re[b];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TOTAL PROBABILITY
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double qm_total_prob(const QubitState *s)
{
    double sum = 0;
    for (int k = 0; k < QM_D; k++)
        sum += s->re[k] * s->re[k] + s->im[k] * s->im[k];
    return sum;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * RENORMALIZE
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void qm_renormalize(QubitState *s)
{
    double norm2 = qm_total_prob(s);
    if (norm2 > 1e-30 && fabs(norm2 - 1.0) > 1e-15) {
        double scale = 1.0 / sqrt(norm2);
        for (int k = 0; k < QM_D; k++) {
            s->re[k] *= scale;
            s->im[k] *= scale;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * COPY
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void qm_copy(QubitState *dst, const QubitState *src)
{
    memcpy(dst, src, sizeof(*dst));
}

#endif /* QUBIT_MANAGEMENT_H */
