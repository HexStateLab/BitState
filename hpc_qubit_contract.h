/*
 * hpc_qubit_contract.h — Pauli-Aware Bond Encoding for D=2
 *
 * The D=2 adaptation of HexState's syntheme-based contraction.
 *
 * D=6 has 15 synthemes (pairings of 6 elements).
 * D=2 has 3 "Pauli projections" — project a 2×2 interaction onto the
 * Pauli basis {I, Z, X, Y}. Each Pauli defines a "lens" through
 * which the entangling interaction can be viewed.
 *
 * HexState's vesica fold (sum/diff of antipodal pairs):
 *   vesica[c] = (state[c] + state[c+3]) / √2
 *   wave[c]   = (state[c] - state[c+3]) / √2
 *
 * D=2 "Hadamard fold" (sum/diff):
 *   plus  = (state[0] + state[1]) / √2 = ⟨+|ψ⟩
 *   minus = (state[0] - state[1]) / √2 = ⟨-|ψ⟩
 *
 * This IS the Hadamard transform viewed as a fold.
 */

#ifndef HPC_QUBIT_CONTRACT_H
#define HPC_QUBIT_CONTRACT_H

#include "hpc_qubit_graph.h"
#include <math.h>
#include <string.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI PROJECTIONS — The D=2 "Synthemes"
 *
 * For D=6: 15 synthemes (pairings of {0..5}) define natural decomposition bases.
 * For D=2: 3 Pauli projections define natural "lenses":
 *
 *   Pauli-Z: diagonal → {|0⟩, |1⟩}
 *   Pauli-X: Hadamard → {|+⟩, |−⟩}
 *   Pauli-Y: S†·H     → {|+i⟩, |−i⟩}
 *
 * These correspond exactly to the 3 triality views.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef enum {
    PAULI_I = 0,    /* Identity  */
    PAULI_Z = 1,    /* Z-basis   */
    PAULI_X = 2,    /* X-basis   */
    PAULI_Y = 3     /* Y-basis   */
} PauliChannel;

/* Pauli matrices (real part, imaginary part) — stored row-major as 2×2 */
/* σ_I = [[1,0],[0,1]] */
static const double PAULI_I_RE[4] = {1,0, 0,1};
static const double PAULI_I_IM[4] = {0,0, 0,0};

/* σ_Z = [[1,0],[0,-1]] */
static const double PAULI_Z_RE[4] = {1,0, 0,-1};
static const double PAULI_Z_IM[4] = {0,0, 0,0};

/* σ_X = [[0,1],[1,0]] */
static const double PAULI_X_RE[4] = {0,1, 1,0};
static const double PAULI_X_IM[4] = {0,0, 0,0};

/* σ_Y = [[0,-i],[i,0]] */
static const double PAULI_Y_RE[4] = {0,0, 0,0};
static const double PAULI_Y_IM[4] = {0,-1, 1,0};

/* ═══════════════════════════════════════════════════════════════════════════════
 * HADAMARD FOLD — The D=2 Vesica Fold
 *
 * Decomposes a 2-vector into sum/diff channels:
 *   plus  = (ψ[0] + ψ[1]) / √2    — symmetric (vesica)
 *   minus = (ψ[0] - ψ[1]) / √2    — antisymmetric (wave)
 *
 * This IS the Hadamard transform viewed as a channel decomposition.
 * Cost: O(2), zero multiplies.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double plus_re,  plus_im;     /* Symmetric channel  (|+⟩)  */
    double minus_re, minus_im;    /* Antisymmetric channel (|−⟩) */
} HadamardFold;

static inline HadamardFold hpcq_hadamard_fold(const double re[2],
                                                const double im[2])
{
    static const double S = 0.7071067811865475244;  /* 1/√2 */
    HadamardFold hf;
    hf.plus_re  = S * (re[0] + re[1]);
    hf.plus_im  = S * (im[0] + im[1]);
    hf.minus_re = S * (re[0] - re[1]);
    hf.minus_im = S * (im[0] - im[1]);
    return hf;
}

static inline void hpcq_hadamard_unfold(const HadamardFold *hf,
                                          double re[2], double im[2])
{
    static const double S = 0.7071067811865475244;
    re[0] = S * (hf->plus_re + hf->minus_re);
    im[0] = S * (hf->plus_im + hf->minus_im);
    re[1] = S * (hf->plus_re - hf->minus_re);
    im[1] = S * (hf->plus_im - hf->minus_im);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI DECOMPOSITION OF A 2×2 INTERACTION
 *
 * Any 2×2 matrix M can be decomposed as:
 *   M = c_I·I + c_Z·Z + c_X·X + c_Y·Y
 * where c_P = Tr(P·M) / 2.
 *
 * This is the D=2 analog of syntheme energy computation:
 * find which Pauli channel captures the most structure.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double coeff_re[4];   /* Pauli coefficients: I, Z, X, Y (real part) */
    double coeff_im[4];   /* Pauli coefficients: I, Z, X, Y (imag part) */
    double energy[4];     /* |c_P|² for each Pauli                      */
    int    dominant;      /* Which Pauli has the most energy             */
} PauliDecomposition;

static inline PauliDecomposition hpcq_pauli_decompose(
    const double M_re[2][2], const double M_im[2][2])
{
    PauliDecomposition pd;
    memset(&pd, 0, sizeof(pd));

    /* c_I = (M[0][0] + M[1][1]) / 2 */
    pd.coeff_re[PAULI_I] = 0.5 * (M_re[0][0] + M_re[1][1]);
    pd.coeff_im[PAULI_I] = 0.5 * (M_im[0][0] + M_im[1][1]);

    /* c_Z = (M[0][0] - M[1][1]) / 2 */
    pd.coeff_re[PAULI_Z] = 0.5 * (M_re[0][0] - M_re[1][1]);
    pd.coeff_im[PAULI_Z] = 0.5 * (M_im[0][0] - M_im[1][1]);

    /* c_X = (M[0][1] + M[1][0]) / 2 */
    pd.coeff_re[PAULI_X] = 0.5 * (M_re[0][1] + M_re[1][0]);
    pd.coeff_im[PAULI_X] = 0.5 * (M_im[0][1] + M_im[1][0]);

    /* c_Y = i(M[0][1] - M[1][0]) / 2 = (im part trick)
     * Actually: Tr(σ_Y · M) / 2
     * σ_Y · M = [[0,-i],[i,0]] · M = [[-i·M[1][·]], [i·M[0][·]]]
     * Tr = -i·M[1][0] + i·M[0][1]
     * c_Y = Tr / 2 = i(M[0][1] - M[1][0]) / 2 */
    pd.coeff_re[PAULI_Y] = 0.5 * (M_im[1][0] - M_im[0][1]);
    pd.coeff_im[PAULI_Y] = 0.5 * (M_re[0][1] - M_re[1][0]);

    /* Compute energies and find dominant */
    double max_energy = -1;
    pd.dominant = PAULI_I;
    for (int p = 0; p < 4; p++) {
        pd.energy[p] = pd.coeff_re[p] * pd.coeff_re[p]
                     + pd.coeff_im[p] * pd.coeff_im[p];
        if (pd.energy[p] > max_energy) {
            max_energy = pd.energy[p];
            pd.dominant = p;
        }
    }

    return pd;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI PROJECTION — The D=2 "Syntheme Projection"
 *
 * Retains only the dominant Pauli component of an interaction.
 * This is the D=2 analog of HexState's syntheme projection.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_pauli_project(const double M_re[2][2],
                                        const double M_im[2][2],
                                        PauliChannel channel,
                                        double out_re[2][2],
                                        double out_im[2][2])
{
    PauliDecomposition pd = hpcq_pauli_decompose(M_re, M_im);
    double cr = pd.coeff_re[channel];
    double ci = pd.coeff_im[channel];

    const double *P_re = NULL, *P_im = NULL;
    switch (channel) {
        case PAULI_I: P_re = PAULI_I_RE; P_im = PAULI_I_IM; break;
        case PAULI_Z: P_re = PAULI_Z_RE; P_im = PAULI_Z_IM; break;
        case PAULI_X: P_re = PAULI_X_RE; P_im = PAULI_X_IM; break;
        case PAULI_Y: P_re = PAULI_Y_RE; P_im = PAULI_Y_IM; break;
    }

    /* out = c · P */
    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
        int k = i * 2 + j;
        out_re[i][j] = cr * P_re[k] - ci * P_im[k];
        out_im[i][j] = cr * P_im[k] + ci * P_re[k];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * FIDELITY — How much of the gate was captured?
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double hpcq_compute_fidelity(const double orig_re[2][2],
                                             const double orig_im[2][2],
                                             const double proj_re[2][2],
                                             const double proj_im[2][2])
{
    double norm_orig = 0.0, norm_proj = 0.0;
    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
        norm_orig += orig_re[i][j] * orig_re[i][j] +
                     orig_im[i][j] * orig_im[i][j];
        norm_proj += proj_re[i][j] * proj_re[i][j] +
                     proj_im[i][j] * proj_im[i][j];
    }
    return (norm_orig > 1e-30) ? norm_proj / norm_orig : 0.0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * OPTIMAL PAULI SELECTION — O(4) lookup (vs O(45) for D=6 synthemes)
 *
 * Selects the Pauli channel that captures the most structure.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline PauliChannel hpcq_select_pauli(const double M_re[2][2],
                                               const double M_im[2][2])
{
    PauliDecomposition pd = hpcq_pauli_decompose(M_re, M_im);
    return (PauliChannel)pd.dominant;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENCODE GATE AS CLIFFORD EDGE — The D=2 syntheme contraction
 *
 * 1. Decompose into Pauli basis — O(4)
 * 2. Select dominant Pauli — O(1)
 * 3. Project onto dominant — O(4)
 * 4. Compute fidelity — O(4)
 * 5. Store as Clifford edge — O(1)
 *
 * Total: O(16). (vs O(36) for D=6 syntheme encoding)
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_encode_clifford(HPCQGraph *g,
                                          uint64_t site_a, uint64_t site_b,
                                          const double phase_re[2][2],
                                          const double phase_im[2][2])
{
    PauliChannel best = hpcq_select_pauli(phase_re, phase_im);

    double proj_re[2][2], proj_im[2][2];
    hpcq_pauli_project(phase_re, phase_im, best, proj_re, proj_im);

    double fidelity = hpcq_compute_fidelity(phase_re, phase_im,
                                             proj_re, proj_im);

    hpcq_grow_edges(g);
    HPCQEdge *e = &g->edges[g->n_edges];
    memset(e, 0, sizeof(HPCQEdge));
    e->type = HPCQ_EDGE_CLIFFORD;
    e->site_a = site_a;
    e->site_b = site_b;
    e->pauli_channel = (uint8_t)best;
    e->fidelity = fidelity;

    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
        double mag = sqrt(proj_re[i][j] * proj_re[i][j] +
                          proj_im[i][j] * proj_im[i][j]);
        if (mag > 1e-15) {
            e->w_re[i][j] = proj_re[i][j] / mag;
            e->w_im[i][j] = proj_im[i][j] / mag;
        } else {
            e->w_re[i][j] = 1.0;
            e->w_im[i][j] = 0.0;
        }
    }

    g->n_edges++;
    g->clifford_edges++;
    hpcq_update_fidelity_stats(g);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HIGH-LEVEL ENCODE — Auto-selects encoding strategy
 *
 * 1. If CZ: exact (fidelity=1.0)
 * 2. If Pauli fidelity ≥ threshold: Clifford edge
 * 3. Otherwise: general phase edge
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define HPCQ_CLIFFORD_THRESHOLD 0.80

static inline void hpcq_encode_2site(HPCQGraph *g,
                                       uint64_t site_a, uint64_t site_b,
                                       const double *G_re, const double *G_im)
{
    /* Extract phase matrix from 4×4 gate */
    double phase_re[2][2], phase_im[2][2];
    for (int j = 0; j < HPCQ_D; j++)
    for (int k = 0; k < HPCQ_D; k++) {
        int idx = (j * HPCQ_D + k) * HPCQ_D * HPCQ_D + (j * HPCQ_D + k);
        double g_re = G_re[idx], g_im = G_im[idx];
        double mag = sqrt(g_re * g_re + g_im * g_im);
        if (mag > 1e-15) {
            phase_re[j][k] = g_re / mag;
            phase_im[j][k] = g_im / mag;
        } else {
            phase_re[j][k] = 1.0;
            phase_im[j][k] = 0.0;
        }
    }

    /* Test for CZ: w(j,k) = (-1)^(j·k) */
    int is_cz = 1;
    for (int j = 0; j < HPCQ_D && is_cz; j++)
    for (int k = 0; k < HPCQ_D && is_cz; k++) {
        uint32_t pi = (j * k) % HPCQ_D;
        double diff_re = phase_re[j][k] - HPCQ_W2_RE[pi];
        double diff_im = phase_im[j][k] - HPCQ_W2_IM[pi];
        if (diff_re * diff_re + diff_im * diff_im > 1e-10)
            is_cz = 0;
    }

    if (is_cz) { hpcq_cz(g, site_a, site_b); return; }

    /* Try Clifford encoding */
    PauliChannel best = hpcq_select_pauli(phase_re, phase_im);
    double proj_re[2][2], proj_im[2][2];
    hpcq_pauli_project(phase_re, phase_im, best, proj_re, proj_im);
    double fidelity = hpcq_compute_fidelity(phase_re, phase_im,
                                             proj_re, proj_im);

    if (fidelity >= HPCQ_CLIFFORD_THRESHOLD) {
        hpcq_encode_clifford(g, site_a, site_b, phase_re, phase_im);
    } else {
        hpcq_general_2site(g, site_a, site_b, G_re, G_im);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HADAMARD FOLD ANALYSIS — Channel decomposition of entanglement
 *
 * For D=2, the vesica analysis becomes the Hadamard fold analysis:
 * how much entanglement is in the symmetric (|+⟩) vs antisymmetric (|−⟩)
 * channel.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double plus_fidelity;      /* Probability in |+⟩ channel */
    double minus_fidelity;     /* Probability in |−⟩ channel */
    double fold_entropy;       /* Binary entropy of the split */
} HPCQFoldAnalysis;

static inline HPCQFoldAnalysis hpcq_analyze_fold(const HPCQGraph *g,
                                                    uint64_t site)
{
    HPCQFoldAnalysis fa;
    memset(&fa, 0, sizeof(fa));

    double re[2], im[2];
    tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE, re, im);

    HadamardFold hf = hpcq_hadamard_fold(re, im);

    double plus_prob = hf.plus_re * hf.plus_re + hf.plus_im * hf.plus_im;
    double minus_prob = hf.minus_re * hf.minus_re + hf.minus_im * hf.minus_im;
    double total = plus_prob + minus_prob;

    if (total > 1e-15) {
        fa.plus_fidelity = plus_prob / total;
        fa.minus_fidelity = minus_prob / total;

        double p = fa.plus_fidelity;
        if (p > 1e-14 && p < 1.0 - 1e-14)
            fa.fold_entropy = -p * log2(p) - (1.0 - p) * log2(1.0 - p);
    }

    return fa;
}

#endif /* HPC_QUBIT_CONTRACT_H */
