/*
 * clifford_exotic.h — Clifford Group "Exotic" Automorphism for D=2
 *
 * The D=2 analog of HexState's S₆ outer automorphism.
 *
 * S₆ is the ONLY symmetric group with a non-trivial outer automorphism.
 * S₂ has no such structure. But our TRIALITY creates 3 views (Z, X, Y),
 * and the single-qubit Clifford group acts on these views by permuting
 * and reflecting the 3 Pauli axes.
 *
 * The Clifford group C₁ has 24 elements (isomorphic to S₄ ≅ octahedral group):
 *   - 6 axis permutations of {Z, X, Y} (≅ S₃)
 *   - 4 sign choices per permutation
 *   - Total: 24 elements
 *
 * Each Clifford element maps (Z, X, Y) → (±P₁, ±P₂, ±P₃) where
 * {P₁, P₂, P₃} is a permutation of {Z, X, Y} with signs constrained
 * by the requirement that the product ZXY → P₁P₂P₃ preserves chirality.
 *
 * The "exotic" operations:
 *   - Apply a gate in a Clifford-conjugated basis
 *   - Dual measurement: standard vs Clifford-conjugate
 *   - Clifford invariant Δ_C: how much structure depends on Pauli frame
 */

#ifndef CLIFFORD_EXOTIC_H
#define CLIFFORD_EXOTIC_H

#include <stdint.h>
#include <math.h>

#define CLIFF_ORDER    24   /* |C₁| = 24 single-qubit Clifford elements   */
#define CLIFF_N_AXES   3    /* Z, X, Y                                     */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ═══════════════════════════════════════════════════════════════════════════════
 * CLIFFORD ELEMENT — Axis permutation + sign flip
 *
 * Maps (Z, X, Y) → (sign[0]·axes[0], sign[1]·axes[1], sign[2]·axes[2])
 * where axes[] is a permutation of {0=Z, 1=X, 2=Y}
 * and sign[] ∈ {+1, -1}.
 *
 * Constraint: det(permutation) × product(signs) = +1
 * (preserves the right-hand rule ZXY = iI)
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    int axes[3];     /* Permutation: where Z,X,Y map to (0=Z, 1=X, 2=Y) */
    int sign[3];     /* Sign: +1 or -1 for each axis                     */
} CliffordElement;

/* ═══════════════════════════════════════════════════════════════════════════════
 * CLIFFORD GATE — 2×2 unitary matrix representation
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double re[2][2];
    double im[2][2];
} CliffordGate;

/* ═══════════════════════════════════════════════════════════════════════════════
 * API — clifford_exotic.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Initialization (generates the 24 Clifford elements) */
void cliff_init(void);

/* Access the 24 elements and their gate matrices */
const CliffordElement *cliff_get_element(int idx);
const CliffordGate    *cliff_get_gate(int idx);

/* Group operations */
int cliff_compose(int a, int b);    /* Returns index of a∘b          */
int cliff_inverse(int idx);          /* Returns index of element⁻¹    */
int cliff_identity(void);            /* Returns index of identity     */

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC GATE APPLICATION
 *
 * Instead of applying gate G, apply C·G·C† for Clifford element C.
 * This transforms the gate into the Clifford-conjugated basis.
 *
 * Physical meaning: apply the "same" operation but in a different
 * Pauli frame. The gate "looks the same" from the rotated perspective.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void cliff_apply_exotic_gate(const double *in_re, const double *in_im,
                              double *out_re, double *out_im,
                              int clifford_idx);

/* ═══════════════════════════════════════════════════════════════════════════════
 * DUAL MEASUREMENT
 *
 * Returns probabilities in BOTH standard and Clifford-conjugate bases.
 *   Standard: probs_std[k] = |ψ[k]|²
 *   Exotic:   probs_exo[k] = |⟨k|C|ψ⟩|² for a probe Clifford C
 *
 * Cost: O(4) per Clifford element.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void cliff_dual_probabilities(const double *re, const double *im,
                               double *probs_std, double *probs_exo,
                               int clifford_idx);

/* ═══════════════════════════════════════════════════════════════════════════════
 * CLIFFORD INVARIANT Δ_C
 *
 * Δ_C(ψ) = (1/24) Σ_{C ∈ C₁} || P_std(ψ) - P_C(ψ) ||²
 *
 * where P_std(ψ) = (|⟨0|ψ⟩|², |⟨1|ψ⟩|²) = standard Born probabilities
 * and   P_C(ψ)   = (|⟨0|C|ψ⟩|², |⟨1|C|ψ⟩|²) = Clifford-rotated probs
 *
 * Δ_C = 0: state is symmetric under all Clifford conjugations
 *           (maximally "generic" — no preferred Pauli frame)
 * Δ_C > 0: state has a preferred Pauli frame
 *           (exploits triality structure)
 *
 * Analog of S₆'s exotic invariant Δ.
 *
 * Cost: O(24 × 8) = O(192).
 * ═══════════════════════════════════════════════════════════════════════════════ */

double cliff_exotic_invariant(const double *re, const double *im);

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC ENTROPY ΔS
 *
 * ΔS = S_std - S_exotic
 *   S_std = Shannon entropy of standard probabilities
 *   S_exotic = Shannon entropy of Clifford-conjugated probabilities
 *
 * ΔS > 0: more ordered in exotic view
 * ΔS < 0: more ordered in standard view
 * ═══════════════════════════════════════════════════════════════════════════════ */

double cliff_exotic_entropy(const double *re, const double *im,
                              int clifford_idx);

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC FINGERPRINT
 *
 * Per-conjugacy-class breakdown of the invariant.
 * The 24-element Clifford group has 5 conjugacy classes:
 *   [0]: {I}             — 1 element
 *   [1]: {±Z, ±X, ±Y}   — 6 elements (half-turns)
 *   [2]: {±S, ±S†, ...}  — 6 elements (quarter-turns)
 *   [3]: {H, SH, ...}    — 8 elements (third-turns, Hadamard-like)
 *   [4]: {HS, ...}       — 3 elements (sixth-turns)
 *
 * Returns 5 values (one per class).
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define CLIFF_NUM_CLASSES 5

void cliff_exotic_fingerprint(const double *re, const double *im,
                                double *class_deltas);

/* ═══════════════════════════════════════════════════════════════════════════════
 * TRIALITY BRIDGE — Connect Clifford rotation to triality views
 *
 * Given a Clifford element, returns which triality view it maps the
 * standard Z-basis to. This enables O(1) "exotic view" access.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int cliff_target_view(int clifford_idx);  /* 0=Edge, 1=Vertex, 2=Diagonal */

#endif /* CLIFFORD_EXOTIC_H */
