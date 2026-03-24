/*
 * qubit_triality.h — Square-Geometry Triality for Qubits
 *
 * Three mutually-defining views of every qubit state, mapped to the
 * three Pauli axes on the Bloch sphere:
 *
 *   EDGE     (Z-basis): |0⟩, |1⟩        — computational basis
 *   VERTEX   (X-basis): |+⟩, |−⟩        — Hadamard basis
 *   DIAGONAL (Y-basis): |+i⟩, |−i⟩      — S†·H basis
 *
 * Square Geometry:
 *   A square has edges, vertices, and diagonals.
 *   Edge of square A = vertex of square B = diagonal of square C.
 *   This is the D=2 analog of HexState's hexagonal triality.
 *
 * Conversion cost: O(2) for all view changes (2 complex multiplies).
 * Key property: H² = I (unlike DFT₆ where F⁴=I).
 *
 * Gate affinity:
 *   Z, S, T, Phase, CZ → O(1) in EDGE view (only touch |1⟩)
 *   X                  → O(1) in VERTEX view (only touch |−⟩)
 *   Y                  → O(1) in DIAGONAL view (only touch |−i⟩)
 */

#ifndef QUBIT_TRIALITY_H
#define QUBIT_TRIALITY_H

#include <math.h>
#include <string.h>
#include <stdint.h>

#define TRI_D    2

/* ═══════════════════════════════════════════════════════════════════════════════
 * VIEW ENUMERATION
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef enum {
    VIEW_EDGE     = 0,   /* Z-basis: computational   */
    VIEW_VERTEX   = 1,   /* X-basis: Hadamard         */
    VIEW_DIAGONAL = 2,   /* Y-basis: S†·H             */
    VIEW_COUNT    = 3
} TrialityView;

static const char *TRIALITY_VIEW_NAMES[VIEW_COUNT] = {
    "Edge (Z)", "Vertex (X)", "Diagonal (Y)"
};

/* ═══════════════════════════════════════════════════════════════════════════════
 * TRIALITY QUBIT — State held in three simultaneous views
 *
 * Only the current view is guaranteed valid. Other views have dirty flags.
 * Lazy conversion: recompute dirty views only on access.
 *
 * 32 bytes × 3 views + metadata = ~112 bytes per TrialityQubit.
 * (HexState's TrialityQuhit is ~580 bytes.)
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double re[VIEW_COUNT][TRI_D];     /* Real parts per view               */
    double im[VIEW_COUNT][TRI_D];     /* Imaginary parts per view          */
    uint8_t dirty[VIEW_COUNT];        /* 1 = view needs recomputation      */
    TrialityView current_view;        /* Which view was last written        */

    /* Optimization flags */
    uint8_t is_eigenstate;            /* 1 = Pauli eigenstate (O(1) convert)*/
    uint8_t is_real;                  /* 1 = all imag parts zero           */
    uint8_t active_mask;              /* 2-bit mask: which basis states are
                                         nonzero (bit 0 = |0⟩, bit 1 = |1⟩) */

    /* Statistics */
    uint32_t conversions;             /* Total view conversions performed   */
    uint32_t gates_applied;           /* Total gates applied                */
    uint32_t conversion_savings;      /* Gates that avoided conversion      */
} TrialityQubit;

/* ═══════════════════════════════════════════════════════════════════════════════
 * CONVERSION MATRICES
 *
 * Edge→Vertex:   H = (1/√2) [[1,1],[1,-1]]
 * Edge→Diagonal: S†·H = (1/√2) [[1,1],[-i,i]]
 * Vertex→Diagonal: S† = [[1,0],[0,-i]]
 *
 * All inverses:
 * Vertex→Edge:    H (self-adjoint)
 * Diagonal→Edge:  H·S = (1/√2) [[1,i],[1,-i]]
 * Diagonal→Vertex: S = [[1,0],[0,i]]
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* ═══════════════════════════════════════════════════════════════════════════════
 * LAZY TRIALITY QUBIT — Heisenberg-picture deferred evaluation
 *
 * Gates accumulate as diagonal phase vectors.
 * DFTs (Hadamards) accumulate as a counter mod 2 (H²=I).
 *
 * Phase segment chain:
 *   Apply Z(θ): multiply |1⟩ by e^{iθ} → accumulate θ as phase offset
 *   Apply H:    increment dft_count mod 2
 *   Apply S:    accumulate π/2 phase
 *   Apply T:    accumulate π/4 phase
 *
 * Materializes only at measurement time.
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define LAZY_MAX_SEGMENTS 64

typedef struct {
    double phase;     /* Phase angle accumulated on |1⟩ */
} LazySegment;

typedef struct {
    LazySegment segments[LAZY_MAX_SEGMENTS];
    uint32_t    num_segments;
    uint32_t    dft_count;         /* Hadamard count mod 2 (H²=I)           */
    int         materialized;      /* 1 = amplitudes are live               */
    double      re[TRI_D];        /* Materialized amplitudes (if live)     */
    double      im[TRI_D];
} LazyTrialityQubit;

/* ═══════════════════════════════════════════════════════════════════════════════
 * API — qubit_triality.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Initialization */
void tri_init(TrialityQubit *tq);
void tri_init_state(TrialityQubit *tq, TrialityView view,
                    const double *re, const double *im);

/* View conversion */
void tri_ensure_view(TrialityQubit *tq, TrialityView target);
void tri_convert(TrialityQubit *tq, TrialityView from, TrialityView to);

/* Get current amplitudes in a given view (lazy recomputes if dirty) */
void tri_get_amplitudes(TrialityQubit *tq, TrialityView view,
                        double *re, double *im);

/* Gate application — routes to optimal view automatically */
void tri_apply_z(TrialityQubit *tq, double theta);
void tri_apply_x(TrialityQubit *tq);
void tri_apply_y(TrialityQubit *tq);
void tri_apply_hadamard(TrialityQubit *tq);
void tri_apply_s(TrialityQubit *tq);
void tri_apply_t(TrialityQubit *tq);

/* Triality rotation: Edge→Vertex→Diagonal→Edge (O(1) relabeling) */
void tri_rotate(TrialityQubit *tq);
void tri_rotate_back(TrialityQubit *tq);

/* Measurement in any view */
int tri_measure(TrialityQubit *tq, TrialityView view, double rand_01);

/* Statistics */
void tri_print_stats(const TrialityQubit *tq);

/* Lazy evaluation API */
void lazy_tri_init(LazyTrialityQubit *lq);
void lazy_tri_apply_phase(LazyTrialityQubit *lq, double theta);
void lazy_tri_apply_hadamard(LazyTrialityQubit *lq);
void lazy_tri_materialize(LazyTrialityQubit *lq);
int  lazy_tri_measure(LazyTrialityQubit *lq, double rand_01);

#endif /* QUBIT_TRIALITY_H */
