/*
 * flat_qubit.h — Flat Representation Optimization for Qubits
 *
 * Two-tier representation that breathes between complexity levels:
 *
 *   FLAT_BASIS:   Single |k⟩ with phase — all Pauli gates O(1)
 *   QUANTUM_FULL: Full TrialityQubit — general gate application
 *
 * For D=2, there is no intermediate "subspace" tier (unlike D=6 with its
 * FLAT_SUBSPACE mode), because D=2 means subspace IS either basis or full.
 *
 * Auto-promotes when gates create superposition.
 * Auto-demotes after collapse.
 */

#ifndef FLAT_QUBIT_H
#define FLAT_QUBIT_H

#include <math.h>
#include <string.h>
#include <stdint.h>

#define FLAT_D  2

/* ═══════════════════════════════════════════════════════════════════════════════
 * REPRESENTATION TIERS
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef enum {
    FLAT_BASIS   = 0,   /* Single |k⟩ with phase — cheapest operations    */
    QUANTUM_FULL = 1    /* Full 2-amplitude state — general gates          */
} FlatMode;

/* ═══════════════════════════════════════════════════════════════════════════════
 * FLAT QUBIT
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    FlatMode mode;

    /* FLAT_BASIS mode: only these are used */
    uint8_t  basis_k;           /* Which basis state: 0 or 1               */
    double   phase_re;          /* Phase factor: re + i*im                 */
    double   phase_im;          /* (always unit magnitude)                 */

    /* QUANTUM_FULL mode: full amplitudes */
    double   re[FLAT_D];
    double   im[FLAT_D];

    /* Statistics */
    uint32_t promotions;        /* Times promoted to QUANTUM_FULL          */
    uint32_t demotions;         /* Times demoted to FLAT_BASIS             */
} FlatQubit;

/* ═══════════════════════════════════════════════════════════════════════════════
 * INITIALIZATION
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Init to |0⟩ in FLAT_BASIS mode */
static inline void flat_init(FlatQubit *fq)
{
    memset(fq, 0, sizeof(*fq));
    fq->mode = FLAT_BASIS;
    fq->basis_k = 0;
    fq->phase_re = 1.0;
    fq->phase_im = 0.0;
}

/* Init to |k⟩ with phase */
static inline void flat_init_basis(FlatQubit *fq, int k, double ph_re, double ph_im)
{
    memset(fq, 0, sizeof(*fq));
    fq->mode = FLAT_BASIS;
    fq->basis_k = (uint8_t)k;
    fq->phase_re = ph_re;
    fq->phase_im = ph_im;
}

/* Auto-detect mode from amplitudes */
static inline void flat_from_amplitudes(FlatQubit *fq,
                                        const double *re, const double *im)
{
    memset(fq, 0, sizeof(*fq));

    /* Check if it's a basis state */
    double p0 = re[0]*re[0] + im[0]*im[0];
    double p1 = re[1]*re[1] + im[1]*im[1];

    if (p0 > 1e-14 && p1 < 1e-14) {
        fq->mode = FLAT_BASIS;
        fq->basis_k = 0;
        double mag = sqrt(p0);
        fq->phase_re = re[0] / mag;
        fq->phase_im = im[0] / mag;
    } else if (p1 > 1e-14 && p0 < 1e-14) {
        fq->mode = FLAT_BASIS;
        fq->basis_k = 1;
        double mag = sqrt(p1);
        fq->phase_re = re[1] / mag;
        fq->phase_im = im[1] / mag;
    } else {
        fq->mode = QUANTUM_FULL;
        fq->re[0] = re[0]; fq->re[1] = re[1];
        fq->im[0] = im[0]; fq->im[1] = im[1];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PROMOTION — FLAT_BASIS → QUANTUM_FULL
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void flat_promote(FlatQubit *fq)
{
    if (fq->mode == QUANTUM_FULL) return;

    fq->re[0] = 0; fq->im[0] = 0;
    fq->re[1] = 0; fq->im[1] = 0;
    fq->re[fq->basis_k] = fq->phase_re;
    fq->im[fq->basis_k] = fq->phase_im;

    fq->mode = QUANTUM_FULL;
    fq->promotions++;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DEMOTION — QUANTUM_FULL → FLAT_BASIS (if possible)
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline int flat_try_demote(FlatQubit *fq)
{
    if (fq->mode == FLAT_BASIS) return 1;

    double p0 = fq->re[0]*fq->re[0] + fq->im[0]*fq->im[0];
    double p1 = fq->re[1]*fq->re[1] + fq->im[1]*fq->im[1];

    if (p0 > 1e-14 && p1 < 1e-14) {
        fq->mode = FLAT_BASIS;
        fq->basis_k = 0;
        double mag = sqrt(p0);
        fq->phase_re = fq->re[0] / mag;
        fq->phase_im = fq->im[0] / mag;
        fq->demotions++;
        return 1;
    }
    if (p1 > 1e-14 && p0 < 1e-14) {
        fq->mode = FLAT_BASIS;
        fq->basis_k = 1;
        double mag = sqrt(p1);
        fq->phase_re = fq->re[1] / mag;
        fq->phase_im = fq->im[1] / mag;
        fq->demotions++;
        return 1;
    }

    return 0;  /* Can't demote — genuine superposition */
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GATES — O(1) for FLAT_BASIS, O(D) for QUANTUM_FULL
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* X gate: |0⟩↔|1⟩ */
static inline void flat_apply_x(FlatQubit *fq)
{
    if (fq->mode == FLAT_BASIS) {
        fq->basis_k = 1 - fq->basis_k;  /* O(1) */
    } else {
        double tr = fq->re[0], ti = fq->im[0];
        fq->re[0] = fq->re[1]; fq->im[0] = fq->im[1];
        fq->re[1] = tr;         fq->im[1] = ti;
    }
}

/* Z gate: |0⟩→|0⟩, |1⟩→-|1⟩ */
static inline void flat_apply_z(FlatQubit *fq)
{
    if (fq->mode == FLAT_BASIS) {
        if (fq->basis_k == 1) {
            fq->phase_re = -fq->phase_re;  /* O(1) */
            fq->phase_im = -fq->phase_im;
        }
    } else {
        fq->re[1] = -fq->re[1];
        fq->im[1] = -fq->im[1];
    }
}

/* Phase gate: |1⟩→e^{iθ}|1⟩ */
static inline void flat_apply_phase(FlatQubit *fq, double theta)
{
    if (fq->mode == FLAT_BASIS) {
        if (fq->basis_k == 1) {
            double cs = cos(theta), sn = sin(theta);
            double pr = fq->phase_re, pi = fq->phase_im;
            fq->phase_re = pr * cs - pi * sn;  /* O(1) */
            fq->phase_im = pr * sn + pi * cs;
        }
    } else {
        double cs = cos(theta), sn = sin(theta);
        double r = fq->re[1], m = fq->im[1];
        fq->re[1] = r * cs - m * sn;
        fq->im[1] = r * sn + m * cs;
    }
}

/* T gate: |1⟩→e^{iπ/4}|1⟩ */
static inline void flat_apply_t(FlatQubit *fq)
{
    flat_apply_phase(fq, atan(1.0));
}

/* T† gate: |1⟩→e^{-iπ/4}|1⟩ */
static inline void flat_apply_td(FlatQubit *fq)
{
    flat_apply_phase(fq, -atan(1.0));
}

/* Hadamard: always promotes (creates superposition from basis state) */
static inline void flat_apply_hadamard(FlatQubit *fq)
{
    flat_promote(fq);

    static const double INV_SQRT2 = 0.7071067811865475244;
    double r0 = fq->re[0], i0 = fq->im[0];
    double r1 = fq->re[1], i1 = fq->im[1];

    fq->re[0] = (r0 + r1) * INV_SQRT2;
    fq->im[0] = (i0 + i1) * INV_SQRT2;
    fq->re[1] = (r0 - r1) * INV_SQRT2;
    fq->im[1] = (i0 - i1) * INV_SQRT2;

    /* Try to demote if it collapsed back to a basis state */
    flat_try_demote(fq);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * READ OUT — Get current amplitudes regardless of mode
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void flat_get_amplitudes(const FlatQubit *fq,
                                       double *re, double *im)
{
    if (fq->mode == FLAT_BASIS) {
        re[0] = 0; im[0] = 0;
        re[1] = 0; im[1] = 0;
        re[fq->basis_k] = fq->phase_re;
        im[fq->basis_k] = fq->phase_im;
    } else {
        re[0] = fq->re[0]; im[0] = fq->im[0];
        re[1] = fq->re[1]; im[1] = fq->im[1];
    }
}

#endif /* FLAT_QUBIT_H */
