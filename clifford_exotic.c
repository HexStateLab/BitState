/*
 * clifford_exotic.c — Clifford Group Implementation
 *
 * Generates all 24 single-qubit Clifford elements and their 2×2 matrices.
 *
 * The 24 elements are constructed from:
 *   6 axis permutations of {Z, X, Y} (= S₃)
 *   × 4 valid sign combinations (constrained by chirality)
 *   = 24 total
 *
 * Each element C maps:
 *   Z → sign[0] · P_{axes[0]}
 *   X → sign[1] · P_{axes[1]}
 *   Y → sign[2] · P_{axes[2]}
 *
 * The 2×2 gate matrix is derived from the axis mapping.
 */

#include "clifford_exotic.h"
#include <string.h>
#include <stdio.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * STATIC DATA — Generated on cliff_init()
 * ═══════════════════════════════════════════════════════════════════════════════ */

static CliffordElement cliff_elements[CLIFF_ORDER];
static CliffordGate    cliff_gates[CLIFF_ORDER];
static int             cliff_ready = 0;
static int             cliff_n_elements = 0;

/* Conjugacy class assignment for each element */
static int cliff_conj_class[CLIFF_ORDER];

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Compute parity of a 3-element permutation
 * ═══════════════════════════════════════════════════════════════════════════════ */

static int perm_parity(const int p[3])
{
    int inv = 0;
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            if (p[i] > p[j]) inv++;
    return (inv % 2 == 0) ? 1 : -1;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Build 2×2 gate matrix from Clifford element
 *
 * Strategy: the element tells us how Paulis transform.
 * We find the unitary U such that U·σ_P·U† = (target Pauli).
 *
 * Key gates that generate all 24 Clifford elements:
 *   I, H, S, H·S, S·H, H·S·H, S·H·S, etc.
 *
 * Instead of deriving analytically, we construct from the 6 generators:
 *   I, X, Y, Z, H, S
 * and their products, then verify the axis mapping.
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Multiply two 2×2 complex matrices: C = A · B */
static void mat_mul(const double Ar[2][2], const double Ai[2][2],
                    const double Br[2][2], const double Bi[2][2],
                    double Cr[2][2], double Ci[2][2])
{
    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
        Cr[i][j] = 0; Ci[i][j] = 0;
        for (int k = 0; k < 2; k++) {
            Cr[i][j] += Ar[i][k] * Br[k][j] - Ai[i][k] * Bi[k][j];
            Ci[i][j] += Ar[i][k] * Bi[k][j] + Ai[i][k] * Br[k][j];
        }
    }
}

/* Compute U·σ·U† and check which Pauli it maps to */
static int identify_pauli_image(const double Ur[2][2], const double Ui[2][2],
                                 const double Pr[2][2], const double Pi[2][2],
                                 int *out_sign)
{
    /* Compute U† */
    double Udr[2][2], Udi[2][2];
    Udr[0][0] =  Ur[0][0]; Udi[0][0] = -Ui[0][0];
    Udr[0][1] =  Ur[1][0]; Udi[0][1] = -Ui[1][0];
    Udr[1][0] =  Ur[0][1]; Udi[1][0] = -Ui[0][1];
    Udr[1][1] =  Ur[1][1]; Udi[1][1] = -Ui[1][1];

    /* Compute U · P · U† */
    double tmp_r[2][2], tmp_i[2][2];
    double res_r[2][2], res_i[2][2];
    mat_mul(Ur, Ui, Pr, Pi, tmp_r, tmp_i);
    mat_mul(tmp_r, tmp_i, Udr, Udi, res_r, res_i);

    /* Pauli matrices for comparison */
    static const double PZ_r[2][2] = {{1,0},{0,-1}};
    static const double PX_r[2][2] = {{0,1},{1,0}};
    static const double PY_r[2][2] = {{0,0},{0,0}};
    static const double PY_i[2][2] = {{0,-1},{1,0}};
    static const double zero[2][2] = {{0,0},{0,0}};

    const double *paulis_r[3] = { &PZ_r[0][0], &PX_r[0][0], &PY_r[0][0] };
    const double *paulis_i[3] = { &zero[0][0],  &zero[0][0],  &PY_i[0][0] };

    for (int p = 0; p < 3; p++) {
        for (int sgn = 1; sgn >= -1; sgn -= 2) {
            double diff = 0;
            for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) {
                double dr = res_r[i][j] - sgn * ((const double(*)[2])paulis_r[p])[i][j];
                double di = res_i[i][j] - sgn * ((const double(*)[2])paulis_i[p])[i][j];
                diff += dr * dr + di * di;
            }
            if (diff < 1e-8) {
                *out_sign = sgn;
                return p;
            }
        }
    }

    *out_sign = 1;
    return -1;  /* Not a Pauli — shouldn't happen for Clifford */
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * INITIALIZATION — Generate all 24 Clifford elements
 *
 * Strategy: compose products of {I, H, S} up to length 4,
 * identify the axis mapping, and deduplicate.
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Generator matrices */
static const double GEN_I_R[2][2] = {{1,0},{0,1}};
static const double GEN_I_I[2][2] = {{0,0},{0,0}};

static const double GEN_H_R[2][2] = {{0.7071067811865475, 0.7071067811865475},
                                      {0.7071067811865475,-0.7071067811865475}};
static const double GEN_H_I[2][2] = {{0,0},{0,0}};

static const double GEN_S_R[2][2] = {{1,0},{0,0}};
static const double GEN_S_I[2][2] = {{0,0},{0,1}};

void cliff_init(void)
{
    if (cliff_ready) return;

    const double (*gens_r[3])[2] = {
        (const double(*)[2])GEN_I_R,
        (const double(*)[2])GEN_H_R,
        (const double(*)[2])GEN_S_R
    };
    const double (*gens_i[3])[2] = {
        (const double(*)[2])GEN_I_I,
        (const double(*)[2])GEN_H_I,
        (const double(*)[2])GEN_S_I
    };

    /* Pauli Z/X/Y matrices */
    static const double PZ_r[2][2] = {{1,0},{0,-1}};
    static const double PX_r[2][2] = {{0,1},{1,0}};
    static const double PY_r[2][2] = {{0,0},{0,0}};
    static const double PY_i[2][2] = {{0,-1},{1,0}};
    static const double zero[2][2] = {{0,0},{0,0}};

    const double (*pauli_r[3])[2] = {
        (const double(*)[2])PZ_r,
        (const double(*)[2])PX_r,
        (const double(*)[2])PY_r
    };
    const double (*pauli_i[3])[2] = {
        (const double(*)[2])zero,
        (const double(*)[2])zero,
        (const double(*)[2])PY_i
    };

    cliff_n_elements = 0;

    /* Generate products of up to 4 generators */
    for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
    for (int c = 0; c < 3; c++)
    for (int d = 0; d < 3; d++) {
        /* U = gen[a] · gen[b] · gen[c] · gen[d] */
        double t1_r[2][2], t1_i[2][2];
        double t2_r[2][2], t2_i[2][2];
        double U_r[2][2], U_i[2][2];

        mat_mul(gens_r[a], gens_i[a], gens_r[b], gens_i[b], t1_r, t1_i);
        mat_mul(t1_r, t1_i, gens_r[c], gens_i[c], t2_r, t2_i);
        mat_mul(t2_r, t2_i, gens_r[d], gens_i[d], U_r, U_i);

        /* Identify axis mapping */
        CliffordElement elem;
        int valid = 1;
        for (int p = 0; p < 3 && valid; p++) {
            int target = identify_pauli_image(U_r, U_i,
                                              pauli_r[p], pauli_i[p],
                                              &elem.sign[p]);
            if (target < 0) { valid = 0; break; }
            elem.axes[p] = target;
        }
        if (!valid) continue;

        /* Check uniqueness */
        int duplicate = 0;
        for (int e = 0; e < cliff_n_elements; e++) {
            if (cliff_elements[e].axes[0] == elem.axes[0] &&
                cliff_elements[e].axes[1] == elem.axes[1] &&
                cliff_elements[e].axes[2] == elem.axes[2] &&
                cliff_elements[e].sign[0] == elem.sign[0] &&
                cliff_elements[e].sign[1] == elem.sign[1] &&
                cliff_elements[e].sign[2] == elem.sign[2]) {
                duplicate = 1;
                break;
            }
        }

        if (!duplicate && cliff_n_elements < CLIFF_ORDER) {
            cliff_elements[cliff_n_elements] = elem;
            memcpy(cliff_gates[cliff_n_elements].re, U_r, sizeof(U_r));
            memcpy(cliff_gates[cliff_n_elements].im, U_i, sizeof(U_i));

            /* Classify conjugacy class by trace */
            double tr_re = U_r[0][0] + U_r[1][1];
            double tr_im = U_i[0][0] + U_i[1][1];
            double tr_mag = sqrt(tr_re*tr_re + tr_im*tr_im);

            if (tr_mag > 1.99)       cliff_conj_class[cliff_n_elements] = 0; /* I */
            else if (tr_mag < 0.01)  cliff_conj_class[cliff_n_elements] = 1; /* half-turn */
            else if (tr_mag > 1.39 && tr_mag < 1.43)
                                     cliff_conj_class[cliff_n_elements] = 2; /* quarter */
            else if (tr_mag > 0.99 && tr_mag < 1.01)
                                     cliff_conj_class[cliff_n_elements] = 3; /* third */
            else                     cliff_conj_class[cliff_n_elements] = 4; /* sixth */

            cliff_n_elements++;
        }

        if (cliff_n_elements >= CLIFF_ORDER) goto done;
    }

done:
    cliff_ready = 1;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ACCESS
 * ═══════════════════════════════════════════════════════════════════════════════ */

const CliffordElement *cliff_get_element(int idx)
{
    if (!cliff_ready) cliff_init();
    if (idx < 0 || idx >= cliff_n_elements) return NULL;
    return &cliff_elements[idx];
}

const CliffordGate *cliff_get_gate(int idx)
{
    if (!cliff_ready) cliff_init();
    if (idx < 0 || idx >= cliff_n_elements) return NULL;
    return &cliff_gates[idx];
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GROUP OPERATIONS
 * ═══════════════════════════════════════════════════════════════════════════════ */

int cliff_identity(void)
{
    if (!cliff_ready) cliff_init();
    for (int i = 0; i < cliff_n_elements; i++) {
        if (cliff_elements[i].axes[0] == 0 &&
            cliff_elements[i].axes[1] == 1 &&
            cliff_elements[i].axes[2] == 2 &&
            cliff_elements[i].sign[0] == 1 &&
            cliff_elements[i].sign[1] == 1 &&
            cliff_elements[i].sign[2] == 1)
            return i;
    }
    return 0;
}

int cliff_compose(int a, int b)
{
    if (!cliff_ready) cliff_init();
    if (a < 0 || a >= cliff_n_elements) return 0;
    if (b < 0 || b >= cliff_n_elements) return 0;

    /* Compose: (a∘b)(P) = a(b(P)) */
    CliffordElement result;
    for (int p = 0; p < 3; p++) {
        int b_target = cliff_elements[b].axes[p];
        int b_sign = cliff_elements[b].sign[p];
        result.axes[p] = cliff_elements[a].axes[b_target];
        result.sign[p] = b_sign * cliff_elements[a].sign[b_target];
    }

    /* Find matching element */
    for (int i = 0; i < cliff_n_elements; i++) {
        if (cliff_elements[i].axes[0] == result.axes[0] &&
            cliff_elements[i].axes[1] == result.axes[1] &&
            cliff_elements[i].axes[2] == result.axes[2] &&
            cliff_elements[i].sign[0] == result.sign[0] &&
            cliff_elements[i].sign[1] == result.sign[1] &&
            cliff_elements[i].sign[2] == result.sign[2])
            return i;
    }
    return 0;
}

int cliff_inverse(int idx)
{
    if (!cliff_ready) cliff_init();
    int id = cliff_identity();
    for (int i = 0; i < cliff_n_elements; i++) {
        if (cliff_compose(idx, i) == id)
            return i;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC GATE APPLICATION — Apply C·|ψ⟩ using Clifford gate matrix
 * ═══════════════════════════════════════════════════════════════════════════════ */

void cliff_apply_exotic_gate(const double *in_re, const double *in_im,
                              double *out_re, double *out_im,
                              int clifford_idx)
{
    if (!cliff_ready) cliff_init();
    if (clifford_idx < 0 || clifford_idx >= cliff_n_elements) {
        out_re[0] = in_re[0]; out_re[1] = in_re[1];
        out_im[0] = in_im[0]; out_im[1] = in_im[1];
        return;
    }

    const CliffordGate *g = &cliff_gates[clifford_idx];

    for (int k = 0; k < 2; k++) {
        out_re[k] = 0; out_im[k] = 0;
        for (int j = 0; j < 2; j++) {
            out_re[k] += g->re[k][j] * in_re[j] - g->im[k][j] * in_im[j];
            out_im[k] += g->re[k][j] * in_im[j] + g->im[k][j] * in_re[j];
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DUAL MEASUREMENT — Standard vs Clifford-conjugate probabilities
 * ═══════════════════════════════════════════════════════════════════════════════ */

void cliff_dual_probabilities(const double *re, const double *im,
                               double *probs_std, double *probs_exo,
                               int clifford_idx)
{
    /* Standard */
    probs_std[0] = re[0]*re[0] + im[0]*im[0];
    probs_std[1] = re[1]*re[1] + im[1]*im[1];

    /* Exotic: P_C(k) = |⟨k|C|ψ⟩|² */
    double out_re[2], out_im[2];
    cliff_apply_exotic_gate(re, im, out_re, out_im, clifford_idx);
    probs_exo[0] = out_re[0]*out_re[0] + out_im[0]*out_im[0];
    probs_exo[1] = out_re[1]*out_re[1] + out_im[1]*out_im[1];
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CLIFFORD INVARIANT Δ_C
 *
 * Δ_C(ψ) = (1/24) Σ_{C ∈ C₁} Σ_k (P_std(k) - P_C(k))²
 *
 * Measures how much the Born probabilities change under Clifford rotations.
 * ═══════════════════════════════════════════════════════════════════════════════ */

double cliff_exotic_invariant(const double *re, const double *im)
{
    if (!cliff_ready) cliff_init();

    double p_std[2];
    p_std[0] = re[0]*re[0] + im[0]*im[0];
    p_std[1] = re[1]*re[1] + im[1]*im[1];

    double total_delta = 0.0;

    for (int c = 0; c < cliff_n_elements; c++) {
        double out_re[2], out_im[2];
        cliff_apply_exotic_gate(re, im, out_re, out_im, c);

        double p_c[2];
        p_c[0] = out_re[0]*out_re[0] + out_im[0]*out_im[0];
        p_c[1] = out_re[1]*out_re[1] + out_im[1]*out_im[1];

        double d0 = p_std[0] - p_c[0];
        double d1 = p_std[1] - p_c[1];
        total_delta += d0*d0 + d1*d1;
    }

    return total_delta / cliff_n_elements;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC ENTROPY
 * ═══════════════════════════════════════════════════════════════════════════════ */

static double binary_entropy(double p)
{
    if (p < 1e-14 || p > 1.0 - 1e-14) return 0.0;
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p);
}

double cliff_exotic_entropy(const double *re, const double *im,
                              int clifford_idx)
{
    double p_std[2], p_exo[2];
    cliff_dual_probabilities(re, im, p_std, p_exo, clifford_idx);

    double s_std = binary_entropy(p_std[0]);
    double s_exo = binary_entropy(p_exo[0]);

    return s_std - s_exo;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * EXOTIC FINGERPRINT — Per-conjugacy-class breakdown
 * ═══════════════════════════════════════════════════════════════════════════════ */

void cliff_exotic_fingerprint(const double *re, const double *im,
                                double *class_deltas)
{
    if (!cliff_ready) cliff_init();

    double p_std[2];
    p_std[0] = re[0]*re[0] + im[0]*im[0];
    p_std[1] = re[1]*re[1] + im[1]*im[1];

    double class_sum[CLIFF_NUM_CLASSES] = {0};
    int    class_count[CLIFF_NUM_CLASSES] = {0};

    for (int c = 0; c < cliff_n_elements; c++) {
        double out_re[2], out_im[2];
        cliff_apply_exotic_gate(re, im, out_re, out_im, c);

        double p_c[2];
        p_c[0] = out_re[0]*out_re[0] + out_im[0]*out_im[0];
        p_c[1] = out_re[1]*out_re[1] + out_im[1]*out_im[1];

        double d0 = p_std[0] - p_c[0];
        double d1 = p_std[1] - p_c[1];
        double delta = d0*d0 + d1*d1;

        int cls = cliff_conj_class[c];
        if (cls >= 0 && cls < CLIFF_NUM_CLASSES) {
            class_sum[cls] += delta;
            class_count[cls]++;
        }
    }

    for (int cls = 0; cls < CLIFF_NUM_CLASSES; cls++) {
        class_deltas[cls] = (class_count[cls] > 0) ?
                            class_sum[cls] / class_count[cls] : 0.0;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TRIALITY BRIDGE — Which view does this Clifford element map Z to?
 * ═══════════════════════════════════════════════════════════════════════════════ */

int cliff_target_view(int clifford_idx)
{
    if (!cliff_ready) cliff_init();
    if (clifford_idx < 0 || clifford_idx >= cliff_n_elements) return 0;

    /* axes[0] tells us where Z maps to: 0=Z (Edge), 1=X (Vertex), 2=Y (Diag) */
    return cliff_elements[clifford_idx].axes[0];
}
