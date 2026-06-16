/*
 * hpc_qubit_graph.h — Holographic Phase Graph for D=2
 *
 * The D=2 adaptation of HexState's HPC module.
 * State vector is never materialized. Entanglement is a graph.
 * Amplitudes computed on demand via O(N+E) traversal.
 *
 * Core formula:
 *   ψ(i₁,...,iₙ) = [Π_k a_k(i_k)] × [Π_edges w_e(i_a, i_b)]
 *
 * For CZ edges: w_e(a,b) = (-1)^(a·b) — EXACT, fidelity = 1.0
 * For general edges: w_e(a,b) = arbitrary 2×2 phase matrix
 * For Clifford edges: w_e determined by Clifford-group projector
 *
 * D=2 simplifications vs D=6:
 *   ω = -1 (vs 6th root of unity)
 *   CZ phase table: only 2 entries {1, -1}
 *   Phase edges: 2×2 = 4 entries (vs 6×6 = 36)
 *   Memory per edge: 64 bytes (vs 576 bytes)
 */

#ifndef HPC_QUBIT_GRAPH_H
#define HPC_QUBIT_GRAPH_H

#include "qubit_triality.h"
#include "born_rule.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * CONSTANTS
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define HPCQ_D          2       /* Physical dimension per site            */
#define HPCQ_INIT_EDGES 4096    /* Initial edge capacity (grows)          */
#define HPCQ_INIT_LOG   8192    /* Initial gate log capacity (grows)      */

/* ω = exp(2πi/2) roots of unity — D=2 */
static const double HPCQ_W2_RE[2] = { 1.0, -1.0 };
static const double HPCQ_W2_IM[2] = { 0.0,  0.0 };

/* Evaluate CZ weight with per-edge X parity: (-1)^((a^xp_a)·(b^xp_b)) */
#define HPCQ_CZ_W(a, b, xpa, xpb) \
    ((((a) ^ (xpa)) * ((b) ^ (xpb)) & 1) ? -1.0 : 1.0)

/* ═══════════════════════════════════════════════════════════════════════════════
 * EDGE TYPES
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef enum {
    HPCQ_EDGE_CZ,        /* Exact CZ: w(a,b) = (-1)^(a·b), fidelity=1.0  */
    HPCQ_EDGE_PHASE,     /* General phase: w(a,b) = arbitrary 2×2 matrix   */
    HPCQ_EDGE_CLIFFORD,  /* Clifford-projected: from Pauli decomposition   */
    HPCQ_EDGE_ABSORBED   /* Edge consumed by multi-edge H absorption       */
} HPCQEdgeType;

/* ═══════════════════════════════════════════════════════════════════════════════
 * WEIGHTED PHASE EDGE — One entangling interaction
 *
 * For CZ edges: only type + site indices are used (4 bytes implicit).
 * For general/Clifford edges: full 2×2 phase matrix stored (32 bytes).
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    HPCQEdgeType type;
    uint64_t     site_a;
    uint64_t     site_b;

    /* Phase matrix: w(a,b) — stored for PHASE and CLIFFORD types.
     * For CZ: implicitly (-1)^(a·b), never stored.
     * 2×2 = 4 complex entries = 64 bytes. */
    double       w_re[HPCQ_D][HPCQ_D];
    double       w_im[HPCQ_D][HPCQ_D];

    /* Clifford metadata (for CLIFFORD type) */
    uint8_t      pauli_channel;   /* Which Pauli basis (0=I, 1=Z, 2=X, 3=Y) */

    /* Per-edge X parity: when X is applied to an endpoint of this edge,
     * the corresponding xp toggles.  The CZ weight becomes
     * (-1)^((a^xp_a)·(b^xp_b)) instead of (-1)^(a·b). */
    uint8_t      xp_a;
    uint8_t      xp_b;

    /* Quality metric */
    double       fidelity;        /* 1.0 = lossless                         */

} HPCQEdge;

/* ═══════════════════════════════════════════════════════════════════════════════
 * GATE LOG ENTRY
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef enum {
    HPCQ_GATE_LOCAL_H,
    HPCQ_GATE_LOCAL_PHASE,
    HPCQ_GATE_LOCAL_X,
    HPCQ_GATE_LOCAL_UNITARY,
    HPCQ_GATE_CZ,
    HPCQ_GATE_GENERAL_2SITE,
    HPCQ_GATE_INIT
} HPCQGateType;

typedef struct {
    HPCQGateType type;
    uint64_t     site_a;
    uint64_t     site_b;
    double       params[4];
    double       fidelity;
} HPCQGateEntry;

/* Absorb entry: site whose H gate was absorbed with >1 incident edges.
 * Each neighbor k stores a 2×2 edge matrix in y-z form:
 *   w_re[k*4 + y*2 + z_k], w_im[k*4 + y*2 + z_k]
 * where y is the pre-H (summed-over) index of the center site.
 *
 * Single-layer (n_layers=1):
 *   ψ_v(x_v, zs) = Σ_y H_{x_v,y} · a_re(y) · Π_{k} w_k(y, z_k)
 * with all neighbors in the "inner" group (k < n_inner).
 *
 * Two-layer (n_layers=2): re-absorption — H applied again on an
 * already-absorbed site.  Edges from the FIRST absorption remain
 * in the inner group; edges added AFTER that go to the outer group
 * (k >= n_inner).  The evaluation uses two nested sums:
 *   ψ_v(x_v, zs) = Σ_q H[x_v][q] · a_cur(q) · Π_{outer} w_k(q, z_k)
 *                   · Σ_y  H[q][y]  · a_re(y)   · Π_{inner} w_k(y, z_k)
 * Multi-layer extension: for L layers, variable y_ℓ (ℓ=0..L-1, ℓ=0 innermost).
 * Neighbor k at layer[k] couples y_{layer[k]} of this center to
 * y_{layer[k]} of the partner center.  State before layer ℓ's H:
 *   ℓ = 0 → a_re/a_im
 *   ℓ >= 1 → a_layer_re[(ℓ-1)*2 .. (ℓ-1)*2+1] / a_layer_im
 * Evaluation for L layers:
 *   ψ_v(x_v) = Σ_{y_0..y_{L-1}} H[x_v][y_{L-1}] · state(L-1,y_{L-1})
 *               · Π_{k: layer[k]=L-1} w_k(y_{L-1}, partner_y)
 *               · H[y_{L-1}][y_{L-2}] · state(L-2,y_{L-2})
 *               · Π_{k: layer[k]=L-2} w_k(y_{L-2}, partner_y) · ...
 *               · H[y_1][y_0] · state(0,y_0) · Π_{k: layer[k]=0} w_k(y_0, partner_y)
 * Diagonal gates (T, S) between H operations modify per-layer states. */
typedef struct {
    uint64_t  center;
    uint64_t  n_nbrs;
    uint64_t *nbrs;
    uint8_t  *layer;      /* [n_nbrs]: layer index for each neighbor */
    double   *w_re;       /* [n_nbrs * 4] */
    double   *w_im;
    int       n_layers;
    double    a_re[2];    /* state BEFORE layer 0's H */
    double    a_im[2];
    double   *a_layer_re; /* [max(0, n_layers-1) * 2] — state before layer ℓ (ℓ>=1) */
    double   *a_layer_im;
    uint8_t   x_parity;   /* number of X gates mod 2 since last H absorption */
    uint8_t  *layer_x_parity; /* [max(0, n_layers-1)] — X parity stored when each layer was created */
} HPCQAbsorbEntry;

/* ═══════════════════════════════════════════════════════════════════════════════
 * HPC QUBIT GRAPH — The state representation
 *
 * This struct IS the state. The 2^N state vector does not exist.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    /* Sites */
    uint64_t        n_sites;
    TrialityQubit  *locals;          /* Per-site local states              */

    /* Phase graph */
    uint64_t        n_edges;
    uint64_t        edge_cap;
    HPCQEdge       *edges;

    /* Absorbed multi-edge H entries */
    uint64_t        n_absorb;
    uint64_t        absorb_cap;
    HPCQAbsorbEntry *absorb;

    /* Gate log */
    uint64_t        n_log;
    uint64_t        log_cap;
    HPCQGateEntry  *gate_log;

    /* Per-site incident edge lists (for O(1) lookup in hadamard_absorb) */
    uint64_t       *inc_counts;      /* [n_sites] count of incident non-absorbed edges */
    uint64_t       *inc_cap;         /* [n_sites] capacity */
    uint64_t      **inc_edges;       /* [n_sites] arrays of edge indices */

    /* Per-site absorb entry index (-1 = none) */
    int64_t        *absorb_idx;      /* [n_sites] index into absorb[], or -1 */

    /* Statistics */
    uint64_t        amp_evals;
    uint64_t        prob_evals;
    uint64_t        measurements;
    uint64_t        cz_edges;
    uint64_t        phase_edges;
    uint64_t        clifford_edges;
    double          min_fidelity;
    double          avg_fidelity;

    /* Global phase parity: toggled when a CZ edge is cancelled with
     * both xp_a and xp_b non-zero.  This accounts for the anti-commutation
     * Z·X = -X·Z that arises when both endpoints have X gates between the
     * two CZ applications.  The net effect is a global -1 phase for each
     * such cancellation, applied as (-1)^global_phase_parity during amplitude
     * evaluation. */
    uint64_t        global_phase_parity;
} HPCQGraph;

/* ═══════════════════════════════════════════════════════════════════════════════
 * LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_destroy(HPCQGraph *g); /* forward decl for create */

static inline HPCQGraph *hpcq_create(uint64_t n_sites)
{
    HPCQGraph *g = (HPCQGraph *)calloc(1, sizeof(HPCQGraph));
    if (!g) return NULL;

    g->n_sites = n_sites;
    g->locals = (TrialityQubit *)calloc(n_sites, sizeof(TrialityQubit));
    if (!g->locals) { free(g); return NULL; }

    for (uint64_t i = 0; i < n_sites; i++)
        tri_init(&g->locals[i]);

    g->edge_cap = (n_sites < HPCQ_INIT_EDGES) ? n_sites * 2 + 16 : HPCQ_INIT_EDGES;
    g->edges = (HPCQEdge *)calloc(g->edge_cap, sizeof(HPCQEdge));
    g->n_edges = 0;

    g->n_absorb = 0;
    g->absorb_cap = 4;
    g->absorb = (HPCQAbsorbEntry *)calloc(g->absorb_cap, sizeof(HPCQAbsorbEntry));
    if (!g->absorb) { free(g->locals); free(g->edges); free(g); return NULL; }

    g->log_cap = HPCQ_INIT_LOG;
    g->gate_log = (HPCQGateEntry *)calloc(g->log_cap, sizeof(HPCQGateEntry));
    g->n_log = 0;

    g->inc_counts = (uint64_t *)calloc(n_sites, sizeof(uint64_t));
    g->inc_cap = (uint64_t *)calloc(n_sites, sizeof(uint64_t));
    g->inc_edges = (uint64_t **)calloc(n_sites, sizeof(uint64_t *));
    g->absorb_idx = (int64_t *)malloc(n_sites * sizeof(int64_t));
    if (!g->inc_counts || !g->inc_cap || !g->inc_edges || !g->absorb_idx)
        { hpcq_destroy(g); return NULL; }
    for (uint64_t i = 0; i < n_sites; i++) g->absorb_idx[i] = -1;

    g->min_fidelity = 1.0;
    g->avg_fidelity = 1.0;

    return g;
}

static inline void hpcq_destroy(HPCQGraph *g)
{
    if (!g) return;
    for (uint64_t a = 0; a < g->n_absorb; a++) {
        free(g->absorb[a].nbrs);
        free(g->absorb[a].w_re);
        free(g->absorb[a].w_im);
        free(g->absorb[a].layer);
        free(g->absorb[a].a_layer_re);
        free(g->absorb[a].a_layer_im);
        free(g->absorb[a].layer_x_parity);
    }
    free(g->absorb);
    free(g->locals);
    free(g->edges);
    free(g->gate_log);
    if (g->inc_edges) {
        for (uint64_t i = 0; i < g->n_sites; i++) free(g->inc_edges[i]);
        free(g->inc_edges);
    }
    free(g->absorb_idx);
    free(g->inc_cap);
    free(g->inc_counts);
    free(g);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * INTERNAL: grow arrays
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_grow_edges(HPCQGraph *g)
{
    if (g->n_edges < g->edge_cap) return;
    g->edge_cap *= 2;
    g->edges = (HPCQEdge *)realloc(g->edges, g->edge_cap * sizeof(HPCQEdge));
}

static inline void hpcq_grow_log(HPCQGraph *g)
{
    if (g->n_log < g->log_cap) return;
    g->log_cap *= 2;
    g->gate_log = (HPCQGateEntry *)realloc(g->gate_log,
                                            g->log_cap * sizeof(HPCQGateEntry));
}

static inline void hpcq_log_gate(HPCQGraph *g, HPCQGateEntry entry)
{
    hpcq_grow_log(g);
    g->gate_log[g->n_log++] = entry;
}

/* ── Incident edge list management ── */
static inline void hpcq_inc_add(HPCQGraph *g, uint64_t site, uint64_t edge_idx)
{
    if (g->inc_counts[site] >= g->inc_cap[site]) {
        uint64_t new_cap = g->inc_cap[site] ? g->inc_cap[site] * 2 : 4;
        g->inc_edges[site] = (uint64_t *)realloc(g->inc_edges[site],
                                                    new_cap * sizeof(uint64_t));
        g->inc_cap[site] = new_cap;
    }
    g->inc_edges[site][g->inc_counts[site]++] = edge_idx;
}

static inline void hpcq_inc_remove(HPCQGraph *g, uint64_t site, uint64_t edge_idx)
{
    for (uint64_t i = 0; i < g->inc_counts[site]; i++) {
        if (g->inc_edges[site][i] == edge_idx) {
            g->inc_edges[site][i] = g->inc_edges[site][--g->inc_counts[site]];
            return;
        }
    }
}

static inline void hpcq_update_fidelity_stats(HPCQGraph *g)
{
    if (g->n_edges == 0) {
        g->min_fidelity = 1.0;
        g->avg_fidelity = 1.0;
        return;
    }
    double sum = 0.0, min_f = 1.0;
    for (uint64_t e = 0; e < g->n_edges; e++) {
        double f = g->edges[e].fidelity;
        sum += f;
        if (f < min_f) min_f = f;
    }
    g->min_fidelity = min_f;
    g->avg_fidelity = sum / g->n_edges;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * LOCAL GATES — Absorbed into the local triality qubit state
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_set_local(HPCQGraph *g, uint64_t site,
                                   const double re[2], const double im[2])
{
    tri_init_state(&g->locals[site], VIEW_EDGE, re, im);
    HPCQGateEntry entry = { .type = HPCQ_GATE_INIT, .site_a = site,
                            .fidelity = 1.0 };
    entry.params[0] = re[0]; entry.params[1] = re[1];
    hpcq_log_gate(g, entry);
}

static inline void hpcq_hadamard(HPCQGraph *g, uint64_t site)
{
    tri_apply_hadamard(&g->locals[site]);
    HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                            .fidelity = 1.0 };
    hpcq_log_gate(g, entry);
}

/* ──────────────────────────────────────────────────────────────────────────────
 * hpcq_hadamard_absorb — H gate on a site with incident phase edges.
 *
 * Standard hpcq_hadamard only modifies the local Z-basis amplitude, which
 * gives WRONG results for qubits with incident CZ/PHASE edges (the edge
 * phase must be evaluated under the sum over y created by H).
 *
 * Cases:
 *   0 incident edges → standard hpcq_hadamard (no sum needed).
 *   1 incident edge → absorb into edge matrix (single-edge absorption):
 *        w'(x_k, x_j) = Σ_y H_{x_k,y} · a_k(y) · w_old(y, x_j)
 *   >1 incident edges → multi-edge absorption:
 *        stores edge matrices in absorb entry, evaluated as
 *        ψ_v(x_v, zs) = Σ_y H_{x_v,y} · a_orig(y) · Π_k w_k(y, z_k)
 * Existing absorb entry → re-absorption (H applied again):
 *        old edges stay in inner group (variable y), new incident
 *        edges go to outer group (variable q).  Two-layer evaluation:
 *        ψ_v(x_v) = Σ_q H[x_v][q] · a_cur(q) · Π_outer w_k(q, z_k)
 *                  · Σ_y  H[q][y]  · a_re(y)   · Π_inner w_k(y, z_k)
 * Diagonal gates (Z, S, T) between H operations modify a_cur.
 * Local state becomes uniform representer (1,1) after absorption.
 * ────────────────────────────────────────────────────────────────────────────── */
static inline void hpcq_hadamard_absorb(HPCQGraph *g, uint64_t site)
{
    static const double S = 0.7071067811865475244;  /* 1/√2 */

    int n_inc = (int)g->inc_counts[site];

    /* Check if this site already has an absorb entry (re-apply H) — O(1) via index */
    int existing_absorb = (int)g->absorb_idx[site];

    /* ─── Re-absorption: H gate on an already-absorbed site ─── */
    if (existing_absorb >= 0) {
        HPCQAbsorbEntry *old = &g->absorb[existing_absorb];
        double cur_re[2], cur_im[2];
        tri_get_amplitudes(&g->locals[site], VIEW_EDGE, cur_re, cur_im);

        if (n_inc == 0) {
            /* No new edges: just record a new H layer */
            int L = old->n_layers;
            old->a_layer_re = (double *)realloc(old->a_layer_re, (L) * 2 * sizeof(double));
            old->a_layer_im = (double *)realloc(old->a_layer_im, (L) * 2 * sizeof(double));
            old->a_layer_re[(L-1)*2] = cur_re[0];
            old->a_layer_re[(L-1)*2+1] = cur_re[1];
            old->a_layer_im[(L-1)*2] = cur_im[0];
            old->a_layer_im[(L-1)*2+1] = cur_im[1];
            old->layer_x_parity = (uint8_t *)realloc(old->layer_x_parity,
                                    (L) * sizeof(uint8_t));
            old->layer_x_parity[L - 1] = old->x_parity;
            old->x_parity = 0;
            old->n_layers = L + 1;
            double uni[2] = {1.0, 1.0}, zero[2] = {0.0, 0.0};
            tri_init_state(&g->locals[site], VIEW_EDGE, uni, zero);
            HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                                    .fidelity = 1.0 };
            hpcq_log_gate(g, entry);
            return;
        }

        /* Merge: old neighbors keep their layers; new edges get layer = old->n_layers
         * (the next layer index).  NEVER combine weights — each layer has its own
         * independent variable, so same partner at different layers is correct. */
        int new_layer = old->n_layers;  /* the layer index being added now */
        uint64_t total_max = old->n_nbrs + (uint64_t)n_inc;
        uint64_t *new_nbrs = (uint64_t *)calloc(total_max, sizeof(uint64_t));
        uint8_t *new_layer_arr = (uint8_t *)calloc(total_max, sizeof(uint8_t));
        double *new_w_re = (double *)calloc(total_max * 4, sizeof(double));
        double *new_w_im = (double *)calloc(total_max * 4, sizeof(double));
        uint64_t n_nbrs = old->n_nbrs;

        memcpy(new_nbrs, old->nbrs, old->n_nbrs * sizeof(uint64_t));
        memcpy(new_layer_arr, old->layer, old->n_nbrs * sizeof(uint8_t));
        memcpy(new_w_re, old->w_re, old->n_nbrs * 4 * sizeof(double));
        memcpy(new_w_im, old->w_im, old->n_nbrs * 4 * sizeof(double));

        /* Iterate IN REVERSE to avoid swap-remove corruption */
        for (int64_t ii = (int64_t)n_inc - 1; ii >= 0; ii--) {
            uint64_t e = g->inc_edges[site][ii];
            HPCQEdge *edge = &g->edges[e];
            uint64_t p = (edge->site_a == site) ? edge->site_b : edge->site_a;
            double e_re[2][2], e_im[2][2];
            if (edge->type == HPCQ_EDGE_CZ) {
                uint8_t xp_q = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                uint8_t xp_z = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                for (int q = 0; q < 2; q++)
                    for (int z = 0; z < 2; z++) {
                        e_re[q][z] = HPCQ_CZ_W(q, z, xp_q, xp_z);
                        e_im[q][z] = 0.0;
                    }
            } else {
                for (int q = 0; q < 2; q++)
                    for (int z = 0; z < 2; z++) {
                        if (edge->site_a == site) {
                            e_re[q][z] = edge->w_re[q][z];
                            e_im[q][z] = edge->w_im[q][z];
                        } else {
                            e_re[q][z] = edge->w_re[z][q];
                            e_im[q][z] = edge->w_im[z][q];
                        }
                    }
            }

            /* Always add as NEW neighbor with current layer — never combine */
            new_nbrs[n_nbrs] = p;
            new_layer_arr[n_nbrs] = (uint8_t)new_layer;
            for (int q = 0; q < 2; q++)
                for (int z = 0; z < 2; z++) {
                    int idx = n_nbrs * 4 + q * 2 + z;
                    new_w_re[idx] = e_re[q][z];
                    new_w_im[idx] = e_im[q][z];
                }
            n_nbrs++;

            hpcq_inc_remove(g, site, e);
            /* Keep in partner's incident list for bilayer center approach */
            if (edge->type == HPCQ_EDGE_CZ) {
                uint8_t xp_y = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                uint8_t xp_z = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        edge->w_re[y][z] = HPCQ_CZ_W(y, z, xp_y, xp_z);
                        edge->w_im[y][z] = 0.0;
                    }
            }
            edge->type = HPCQ_EDGE_ABSORBED;
        }

        free(old->nbrs);
        free(old->layer);
        free(old->w_re);
        free(old->w_im);
        old->n_nbrs = n_nbrs;
        old->nbrs = new_nbrs;
        old->layer = new_layer_arr;
        old->w_re = new_w_re;
        old->w_im = new_w_im;

        /* Store state before this H as a_layer entry for the new layer */
        int L_new = old->n_layers + 1;  /* old->n_layers is the count BEFORE this re-absorption */
        old->a_layer_re = (double *)realloc(old->a_layer_re, (L_new - 1) * 2 * sizeof(double));
        old->a_layer_im = (double *)realloc(old->a_layer_im, (L_new - 1) * 2 * sizeof(double));
        old->a_layer_re[(L_new - 2) * 2] = cur_re[0];
        old->a_layer_re[(L_new - 2) * 2 + 1] = cur_re[1];
        old->a_layer_im[(L_new - 2) * 2] = cur_im[0];
        old->a_layer_im[(L_new - 2) * 2 + 1] = cur_im[1];
        old->layer_x_parity = (uint8_t *)realloc(old->layer_x_parity, (L_new - 1) * sizeof(uint8_t));
        old->layer_x_parity[L_new - 2] = old->x_parity;
        old->x_parity = 0;
        old->n_layers = L_new;

        double uni[2] = {1.0, 1.0}, zero[2] = {0.0, 0.0};
        tri_init_state(&g->locals[site], VIEW_EDGE, uni, zero);
        HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                                .fidelity = 1.0 };
        hpcq_log_gate(g, entry);
        return;
    }

    /* ─── First-time absorption (no prior absorb entry) ─── */
    if (n_inc == 0) {
        hpcq_hadamard(g, site);
        return;
    }

    if (n_inc >= 1) {
        /* Multi-edge absorption (n_inc >= 1) */
        uint64_t *nbrs = (uint64_t *)calloc(n_inc, sizeof(uint64_t));
        double *w_re = (double *)calloc(n_inc * 4, sizeof(double));
        double *w_im = (double *)calloc(n_inc * 4, sizeof(double));
        uint64_t n_nbrs = 0;

        for (int64_t ii = (int64_t)n_inc - 1; ii >= 0; ii--) {
            uint64_t e = g->inc_edges[site][ii];
            HPCQEdge *edge = &g->edges[e];
            uint64_t p = (edge->site_a == site) ? edge->site_b : edge->site_a;

            double e_re[2][2], e_im[2][2];
            if (edge->type == HPCQ_EDGE_CZ) {
                uint8_t xp_y = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                uint8_t xp_z = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        e_re[y][z] = HPCQ_CZ_W(y, z, xp_y, xp_z);
                        e_im[y][z] = 0.0;
                    }
            } else {
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        if (edge->site_a == site) {
                            e_re[y][z] = edge->w_re[y][z];
                            e_im[y][z] = edge->w_im[y][z];
                        } else {
                            e_re[y][z] = edge->w_re[z][y];
                            e_im[y][z] = edge->w_im[z][y];
                        }
                    }
            }

            int found = -1;
            for (uint64_t k = 0; k < n_nbrs; k++)
                if (nbrs[k] == p) { found = (int)k; break; }

            if (found >= 0) {
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        int idx = found * 4 + y * 2 + z;
                        double wr = w_re[idx], wi = w_im[idx];
                        w_re[idx] = wr * e_re[y][z] - wi * e_im[y][z];
                        w_im[idx] = wr * e_im[y][z] + wi * e_re[y][z];
                    }
            } else {
                nbrs[n_nbrs] = p;
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        int idx = n_nbrs * 4 + y * 2 + z;
                        w_re[idx] = e_re[y][z];
                        w_im[idx] = e_im[y][z];
                    }
                n_nbrs++;
            }

            hpcq_inc_remove(g, site, e);
            /* Keep edge in PARTNER's incident list so partner can also create
             * an absorb entry for this edge — both qubits become centers and
             * the CZ weight is evaluated in the component joint sum using
             * both qubits' y variables.
             * Store the weight matrix for the partner to read later. */
            if (edge->type == HPCQ_EDGE_CZ) {
                uint8_t xp_y = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                uint8_t xp_z = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                for (int y = 0; y < 2; y++)
                    for (int z = 0; z < 2; z++) {
                        edge->w_re[y][z] = HPCQ_CZ_W(y, z, xp_y, xp_z);
                        edge->w_im[y][z] = 0.0;
                    }
            }
            edge->type = HPCQ_EDGE_ABSORBED;
        }

        /* Grow absorb array if needed */
        if (g->n_absorb >= g->absorb_cap) {
            g->absorb_cap = g->absorb_cap ? g->absorb_cap * 2 : 4;
            g->absorb = (HPCQAbsorbEntry *)realloc(g->absorb,
                        g->absorb_cap * sizeof(HPCQAbsorbEntry));
        }

        double orig_re[2], orig_im[2];
        tri_get_amplitudes(&g->locals[site], VIEW_EDGE, orig_re, orig_im);
        double uni[2] = {1.0, 1.0}, zero[2] = {0.0, 0.0};
        tri_init_state(&g->locals[site], VIEW_EDGE, uni, zero);

        uint64_t ai = g->n_absorb;
        g->absorb[ai].center = site;
        g->absorb[ai].n_nbrs = n_nbrs;
        g->absorb[ai].nbrs = nbrs;
        g->absorb[ai].layer = (uint8_t *)calloc(n_nbrs, sizeof(uint8_t));
        memset(g->absorb[ai].layer, 0, n_nbrs * sizeof(uint8_t));
        g->absorb[ai].w_re = w_re;
        g->absorb[ai].w_im = w_im;
        g->absorb[ai].a_re[0] = orig_re[0];
        g->absorb[ai].a_re[1] = orig_re[1];
        g->absorb[ai].a_im[0] = orig_im[0];
        g->absorb[ai].a_im[1] = orig_im[1];
        g->absorb[ai].n_layers = 1;
        g->absorb[ai].a_layer_re = NULL;
        g->absorb[ai].a_layer_im = NULL;
        g->absorb[ai].x_parity = 0;
        g->absorb[ai].layer_x_parity = NULL;
        g->absorb_idx[site] = (int64_t)ai;
        g->n_absorb++;

        HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                                .fidelity = 1.0 };
        hpcq_log_gate(g, entry);
        return;
    }

    /* Single-edge absorption (first time, exactly 1 incident edge) */
    uint64_t e_idx = g->inc_edges[site][0];
    double a_re[2], a_im[2];
    tri_get_amplitudes(&g->locals[site], VIEW_EDGE, a_re, a_im);

    HPCQEdge *edge = &g->edges[e_idx];

    double wr[2][2] = {{0}}, wi[2][2] = {{0}};
    for (int xk = 0; xk < 2; xk++) {
        for (int xj = 0; xj < 2; xj++) {
            double sum_re = 0.0, sum_im = 0.0;
            for (int y = 0; y < 2; y++) {
                double H_re = (xk == 0) ? S : (y == 0 ? S : -S);
                double ha_re = H_re * a_re[y];
                double ha_im = H_re * a_im[y];
                double w_re_yj, w_im_yj;
                if (edge->type == HPCQ_EDGE_CZ) {
                    uint8_t xp_y = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                    uint8_t xp_xj = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                    w_re_yj = HPCQ_CZ_W(y, xj, xp_y, xp_xj);
                    w_im_yj = 0.0;
                } else {
                    /* Edge stores w[site_a_val][site_b_val].
                     * We need w(absorbing_preH_val, partner_val).
                     * When site_a is absorbing: y = absorbing_preH, xj = partner → w_re[y][xj] ✓
                     * When site_b is absorbing: y = absorbing_preH, xj = partner but
                     *   w_re is w[site_a=partner][site_b=absorbing] → transpose: w_re[xj][y] */
                    if (edge->site_a == site) {
                        w_re_yj = edge->w_re[y][xj];
                        w_im_yj = edge->w_im[y][xj];
                    } else {
                        w_re_yj = edge->w_re[xj][y];
                        w_im_yj = edge->w_im[xj][y];
                    }
                }
                sum_re += ha_re * w_re_yj - ha_im * w_im_yj;
                sum_im += ha_re * w_im_yj + ha_im * w_re_yj;
            }
            wr[xk][xj] = sum_re;
            wi[xk][xj] = sum_im;
        }
    }

    double uni[2] = {1.0, 1.0}, zero[2] = {0.0, 0.0};
    tri_init_state(&g->locals[site], VIEW_EDGE, uni, zero);

    /* Store wr at edge->w_re[site_a_val][site_b_val].
     * wr[xk][xj] where xk = absorbing_postH_val, xj = partner_val.
     * When site_a is absorbing: w_re[xk][xj] = wr[xk][xj] ✓ (xk=site_a, xj=site_b)
     * When site_b is absorbing: w_re[xj][xk] = wr[xk][xj] (xj=site_a, xk=site_b) */
    edge->type = HPCQ_EDGE_PHASE;
    edge->fidelity = 1.0;
    if (edge->site_a == site) {
        memcpy(edge->w_re, wr, sizeof(wr));
        memcpy(edge->w_im, wi, sizeof(wi));
    } else {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) {
                edge->w_re[j][i] = wr[i][j];
                edge->w_im[j][i] = wi[i][j];
            }
    }

    HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                            .fidelity = 1.0 };
    hpcq_log_gate(g, entry);
}

static inline void hpcq_phase(HPCQGraph *g, uint64_t site, double theta)
{
    tri_apply_z(&g->locals[site], theta);
    HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_PHASE, .site_a = site,
                            .fidelity = 1.0 };
    entry.params[0] = theta;
    hpcq_log_gate(g, entry);
}

/* T gate: |1⟩ → e^{iπ/4}|1⟩ — Z(π/4), exact, works everywhere */
static inline void hpcq_t(HPCQGraph *g, uint64_t site)
{
    hpcq_phase(g, site, atan(1.0));  /* π/4 */
}

/* T† gate: |1⟩ → e^{-iπ/4}|1⟩ */
static inline void hpcq_td(HPCQGraph *g, uint64_t site)
{
    hpcq_phase(g, site, -atan(1.0));
}

static inline void hpcq_x(HPCQGraph *g, uint64_t site)
{
    int64_t ai = g->absorb_idx[site];
    if (ai >= 0) {
        /* Center qubit: toggle X parity */
        g->absorb[ai].x_parity ^= 1;
    } else {
        /* Regular site: X and CZ don't commute. When X is applied to a
         * site with incident CZ edges, the CZ weight transforms:
         *   (-1)^(a·b) → (-1)^((a⊕1)·b) = (-1)^(a·b)·(-1)^b
         * Instead of applying Z(π) to site b's local state (which would
         * leak to other edges involving site b), toggle the per-edge X
         * parity flag on this edge's endpoint. */
        tri_apply_x(&g->locals[site]);
        for (uint64_t ii = 0; ii < g->inc_counts[site]; ii++) {
            uint64_t e = g->inc_edges[site][ii];
            HPCQEdge *edge = &g->edges[e];
            if (edge->type == HPCQ_EDGE_CZ || edge->type == HPCQ_EDGE_ABSORBED) {
                if (edge->site_a == site) edge->xp_a ^= 1;
                else edge->xp_b ^= 1;
                if (edge->type == HPCQ_EDGE_ABSORBED) {
                    uint64_t center = (edge->site_a == site) ? edge->site_b : edge->site_a;
                    int64_t ai = g->absorb_idx[center];
                    if (ai >= 0) {
                        HPCQAbsorbEntry *entry = &g->absorb[ai];
                        for (uint64_t k = 0; k < entry->n_nbrs; k++) {
                            if (entry->nbrs[k] == site) {
                                uint8_t xp_y = (edge->site_a == center) ? edge->xp_a : edge->xp_b;
                                uint8_t xp_z = (edge->site_a == center) ? edge->xp_b : edge->xp_a;
                                for (int y = 0; y < 2; y++)
                                    for (int z = 0; z < 2; z++) {
                                        int idx = (int)k * 4 + y * 2 + z;
                                        entry->w_re[idx] = HPCQ_CZ_W(y, z, xp_y, xp_z);
                                        entry->w_im[idx] = 0.0;
                                    }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_X, .site_a = site,
                            .fidelity = 1.0 };
    hpcq_log_gate(g, entry);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CZ GATE — EXACT in HPC (INVOLUTION: CZ² = I)
 *
 * CZ is EXACT: w(a,b) = (-1)^(a·b).
 * Only |1,1⟩ picks up a phase of -1. All others unchanged.
 * Fidelity = 1.0. Always.
 *
 * KEY PROPERTY: CZ · CZ = I. If a CZ edge already exists between
 * (site_a, site_b), applying CZ again CANCELS it (swap-remove).
 * This keeps the edge count bounded by the actual entanglement
 * structure rather than growing linearly with circuit depth.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_cz(HPCQGraph *g, uint64_t site_a, uint64_t site_b)
{
    /* Check for existing CZ edge between this pair — CZ² = I cancellation
     * Uses incident edge list of site_a for O(deg) lookup instead of O(E). */
    for (uint64_t ii = 0; ii < g->inc_counts[site_a]; ii++) {
        uint64_t e = g->inc_edges[site_a][ii];
        HPCQEdge *edge = &g->edges[e];
        if (edge->type == HPCQ_EDGE_CZ &&
            ((edge->site_a == site_a && edge->site_b == site_b) ||
             (edge->site_a == site_b && edge->site_b == site_a))) {
            /* Save xp values BEFORE swap-remove invalidates edge pointer */
            uint8_t xp_a = edge->xp_a;
            uint8_t xp_b = edge->xp_b;
            /* Cancel: swap-remove this edge */
            hpcq_inc_remove(g, edge->site_a, e);
            hpcq_inc_remove(g, edge->site_b, e);
            uint64_t moved_idx = --g->n_edges;
            if (moved_idx != e) {
                uint64_t moved_sa = g->edges[moved_idx].site_a;
                uint64_t moved_sb = g->edges[moved_idx].site_b;
                int moved_type = g->edges[moved_idx].type;
                hpcq_inc_remove(g, moved_sa, moved_idx);
                hpcq_inc_remove(g, moved_sb, moved_idx);
                g->edges[e] = g->edges[moved_idx];
                if (moved_type != HPCQ_EDGE_ABSORBED) {
                    hpcq_inc_add(g, moved_sa, e);
                    hpcq_inc_add(g, moved_sb, e);
                }
            }
            g->cz_edges--;
            /* Non-zero xp means X gate(s) were applied between CZ applications.
             * From CZ·X·CZ = X·Z·CZ·CZ = X·Z (when X is between two CZ gates),
             * a Z phase on the neighbor survives the CZ cancellation.
             * Apply Z(π) to the opposite site for each non-zero xp. */
            if (xp_a) tri_apply_z(&g->locals[site_b], M_PI);
            if (xp_b) tri_apply_z(&g->locals[site_a], M_PI);
            /* Both xp non-zero → Z applied to both endpoints after X.
             * Z·X = -X·Z on each qubit, giving an overall -1 phase on
             * the state.  Track as global phase parity. */
            if (xp_a && xp_b) g->global_phase_parity ^= 1;

            HPCQGateEntry entry = {
                .type = HPCQ_GATE_CZ,
                .site_a = site_a, .site_b = site_b,
                .fidelity = 1.0
            };
            hpcq_log_gate(g, entry);
            return;
        }
    }

    /* NOTE: ABSORBED edges (from H absorption) should NOT be cancelled here.
     * CZ·H·CZ ≠ H·CZ² = H — the second CZ must act on the OUTER variable xv,
     * not the inner variable y.  The absorbed edge handles the first CZ via the
     * inner sum; the second CZ needs a new edge that the late-edge handler will
     * evaluate with xv.  So ABSORBED edges are intentionally NOT checked here —
     * the code falls through to create a new CZ edge below. */

    /* No existing edge — add new one */
    hpcq_grow_edges(g);

    uint64_t idx = g->n_edges;
    HPCQEdge *e = &g->edges[idx];
    memset(e, 0, sizeof(HPCQEdge));
    e->type = HPCQ_EDGE_CZ;
    e->site_a = site_a;
    e->site_b = site_b;
    e->fidelity = 1.0;

    g->n_edges++;
    g->cz_edges++;
    hpcq_inc_add(g, site_a, idx);
    hpcq_inc_add(g, site_b, idx);

    HPCQGateEntry entry = {
        .type = HPCQ_GATE_CZ,
        .site_a = site_a, .site_b = site_b,
        .fidelity = 1.0
    };
    hpcq_log_gate(g, entry);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GENERAL 2-SITE GATE — Encoded as weighted phase edge
 *
 * For a general 4×4 gate G on 2 qubits:
 * Extract diagonal phases w(j,k) = G_{(j,k),(j,k)} / |G_{(j,k),(j,k)}|
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_general_2site(HPCQGraph *g, uint64_t site_a,
                                       uint64_t site_b,
                                       const double *G_re, const double *G_im)
{
    hpcq_grow_edges(g);

    HPCQEdge *e = &g->edges[g->n_edges];
    memset(e, 0, sizeof(HPCQEdge));
    e->type = HPCQ_EDGE_PHASE;
    e->site_a = site_a;
    e->site_b = site_b;

    double fidelity_sum = 0.0;
    int fidelity_count = 0;

    for (int j = 0; j < HPCQ_D; j++) {
        for (int k = 0; k < HPCQ_D; k++) {
            int idx = (j * HPCQ_D + k) * HPCQ_D * HPCQ_D + (j * HPCQ_D + k);
            double g_re = G_re[idx];
            double g_im = G_im[idx];
            double mag = sqrt(g_re * g_re + g_im * g_im);

            if (mag > 1e-15) {
                e->w_re[j][k] = g_re / mag;
                e->w_im[j][k] = g_im / mag;
            } else {
                e->w_re[j][k] = 1.0;
                e->w_im[j][k] = 0.0;
            }

            /* Row fidelity */
            double row_norm2 = 0.0;
            for (int m = 0; m < HPCQ_D; m++)
                for (int n = 0; n < HPCQ_D; n++) {
                    int ridx = (j * HPCQ_D + k) * HPCQ_D * HPCQ_D + (m * HPCQ_D + n);
                    row_norm2 += G_re[ridx] * G_re[ridx] + G_im[ridx] * G_im[ridx];
                }
            if (row_norm2 > 1e-30) {
                fidelity_sum += (g_re * g_re + g_im * g_im) / row_norm2;
                fidelity_count++;
            }
        }
    }

    e->fidelity = (fidelity_count > 0) ? fidelity_sum / fidelity_count : 0.0;

    uint64_t e_idx2 = g->n_edges;
    hpcq_inc_add(g, site_a, e_idx2);
    hpcq_inc_add(g, site_b, e_idx2);
    g->n_edges++;
    g->phase_edges++;
    hpcq_update_fidelity_stats(g);

    HPCQGateEntry entry = {
        .type = HPCQ_GATE_GENERAL_2SITE,
        .site_a = site_a, .site_b = site_b,
        .fidelity = e->fidelity
    };
    hpcq_log_gate(g, entry);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * AMPLITUDE EVALUATION — O(N + E)
 *
 * ψ(i₁,...,iₙ) = [Π_k a_k(i_k)] × [Π_edges w_e(i_a, i_b)]
 *
 * For CZ edges: w_e(a,b) = (-1)^(a·b) — branch-free bit-AND
 * For PHASE/CLIFFORD edges: lookup from stored 2×2 matrix
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_amplitude(const HPCQGraph *g,
                                   const uint32_t *indices,
                                   double *out_re, double *out_im)
{
    double re = 1.0, im = 0.0;

    /* Step 1: Product of local amplitudes — O(N)
     * Skip absorbed centers (their contribution comes from Step 3). */
    for (uint64_t k = 0; k < g->n_sites; k++) {
        if (g->absorb_idx[k] >= 0) continue;
        uint32_t idx = indices[k];
        double a_re, a_im;
        double re_buf[HPCQ_D], im_buf[HPCQ_D];
        tri_get_amplitudes((TrialityQubit *)&g->locals[k], VIEW_EDGE,
                          re_buf, im_buf);
        a_re = re_buf[idx];
        a_im = im_buf[idx];
        double new_re = re * a_re - im * a_im;
        double new_im = re * a_im + im * a_re;
        re = new_re;
        im = new_im;
    }

    /* Step 2: Phase edge accumulation — O(E)
     * Skip ABSORBED edges and any CZ edge involving a center
     * (center-regular handled in Step 3a, center-center in component sum). */
    for (uint64_t e = 0; e < g->n_edges; e++) {
        const HPCQEdge *edge = &g->edges[e];

        /* Skip edges consumed by multi-edge absorption */
        if (edge->type == HPCQ_EDGE_ABSORBED) continue;

        /* Skip CZ edges where either endpoint is a center — handled in
         * Step 3a (center-regular) or component sum (center-center). */
        if (edge->type == HPCQ_EDGE_CZ) {
            int ca = (g->absorb_idx[edge->site_a] >= 0) ? 1 : 0;
            int cb = (g->absorb_idx[edge->site_b] >= 0) ? 1 : 0;
            if (ca || cb) continue;
        }

        uint32_t ia = indices[edge->site_a];
        uint32_t ib = indices[edge->site_b];

        double w_re, w_im;

        if (edge->type == HPCQ_EDGE_CZ) {
            /* CZ: (-1)^((ia^xp_a)·(ib^xp_b)) */
            double w_cz = HPCQ_CZ_W(ia, ib, edge->xp_a, edge->xp_b);
            w_re = w_cz; w_im = 0.0;
        } else {
            w_re = edge->w_re[ia][ib];
            w_im = edge->w_im[ia][ib];
        }

        double new_re = re * w_re - im * w_im;
        double new_im = re * w_im + im * w_re;
        re = new_re;
        im = new_im;
    }

    /* Step 3: Handle absorbed centers via connected-component joint sums.
     *
     * Each absorbed center has an inner variable y_i ∈ {0,1}.  The
     * correct amplitude for absorbed centers is:
     *
     *   ψ = Σ_{y_0..y_{n-1}} Π_i H[x_i][y_i] · a_re(i,y_i) · Π_{non-center k} w_k(y_i, z_k)
     *                          · Π_{center-center edges (i,j)} w_{ij}(y_i, y_j)
     *
     * where the sum is over ALL inner variables jointly.  Centers that share
     * a center-center edge must be summed together as a connected component.
     *
     * For n_layers==2, the inner variable y_i couples through H[q_i][y_i]
     * with an outer variable q_i.  The joint sum then includes both.
     */
    static const double SQ = 0.7071067811865475244;
    uint64_t n_absorb = g->n_absorb;

    if (n_absorb == 0) {
        if (g->global_phase_parity & 1) { re = -re; im = -im; }
        *out_re = re;
        *out_im = im;
        ((HPCQGraph *)g)->amp_evals++;
        return;
    }

    /* ── 3a: Per-layer self-factors (excluding center-center edges) ── */
    double (*sf_re)[2] = (double (*)[2])calloc(n_absorb, 2 * sizeof(double));
    double (*sf_im)[2] = (double (*)[2])calloc(n_absorb, 2 * sizeof(double));
    /* per-layer outer factors: so_layer[a][(li-1)*2 + y] for li=1..L-1 */
    double **so_layer_re = (double **)calloc(n_absorb, sizeof(double *));
    double **so_layer_im = (double **)calloc(n_absorb, sizeof(double *));
    double *na_cz_re = (double *)calloc(n_absorb, sizeof(double));
    double *na_cz_im = (double *)calloc(n_absorb, sizeof(double));
    int *nl_arr = (int *)calloc(n_absorb, sizeof(int));
    for (uint64_t a = 0; a < n_absorb; a++) { na_cz_re[a] = 1.0; na_cz_im[a] = 0.0; }

    for (uint64_t a = 0; a < n_absorb; a++) {
        uint64_t nn = g->absorb[a].n_nbrs;
        int L = (int)g->absorb[a].n_layers;
        nl_arr[a] = L;

        double pi_re[2] = {1.0, 1.0}, pi_im[2] = {0.0, 0.0};

        /* Allocate per-layer outer factors for layers 1..L-1 */
        if (L > 1) {
            so_layer_re[a] = (double *)calloc((L-1) * 2, sizeof(double));
            so_layer_im[a] = (double *)calloc((L-1) * 2, sizeof(double));
            for (int li = 0; li < L-1; li++) {
                so_layer_re[a][li*2] = 1.0; so_layer_re[a][li*2+1] = 1.0;
            }
        }

        /* Group neighbors by layer: layer 0 → pi, layers 1..L-1 → so_layer[li] */
        for (uint64_t k = 0; k < nn; k++) {
            uint64_t nb = g->absorb[a].nbrs[k];
            if (g->absorb_idx[nb] >= 0) continue; /* center-center → component sum */
            int li = (int)g->absorb[a].layer[k];
            int z = (int)indices[nb];
            for (int y = 0; y < 2; y++) {
                int idx = (int)k * 4 + y * 2 + z;
                double wr = g->absorb[a].w_re[idx], wi = g->absorb[a].w_im[idx];
                if (li == 0) {
                    double pr = pi_re[y], pim = pi_im[y];
                    pi_re[y] = pr * wr - pim * wi; pi_im[y] = pr * wi + pim * wr;
                } else {
                    int sli = (li - 1) * 2 + y;
                    double pr = so_layer_re[a][sli], pim = so_layer_im[a][sli];
                    so_layer_re[a][sli] = pr * wr - pim * wi;
                    so_layer_im[a][sli] = pr * wi + pim * wr;
                }
            }
        }
        /* Add non-absorbed CZ edges to outermost layer */
        if (L >= 1) {
            uint64_t q = g->absorb[a].center;
            for (uint64_t e = 0; e < g->n_edges; e++) {
                const HPCQEdge *edge = &g->edges[e];
                if (edge->type != HPCQ_EDGE_CZ) continue;
                uint64_t other;
                if (edge->site_a == q && g->absorb_idx[edge->site_b] < 0)
                    other = edge->site_b;
                else if (edge->site_b == q && g->absorb_idx[edge->site_a] < 0)
                    other = edge->site_a;
                else
                    continue;
                if (L == 1) {
                    int z = (int)indices[other];
                    int xv = (int)(indices[q] ^ g->absorb[a].x_parity);
                    uint8_t xp_q = (edge->site_a == q) ? edge->xp_a : edge->xp_b;
                    uint8_t xp_o = (edge->site_a == q) ? edge->xp_b : edge->xp_a;
                    double wr = HPCQ_CZ_W(xv, z, xp_q, xp_o);
                    na_cz_re[a] *= wr; na_cz_im[a] *= wr;
                } else {
                    int z = (int)indices[other];
                    uint8_t xp_q = (edge->site_a == q) ? edge->xp_a : edge->xp_b;
                    uint8_t xp_o = (edge->site_a == q) ? edge->xp_b : edge->xp_a;
                    int osli = (L - 2) * 2; /* outermost layer index in so_layer */
                    for (int yo = 0; yo < 2; yo++) {
                        double wr = HPCQ_CZ_W(yo, z, xp_q, xp_o);
                        so_layer_re[a][osli + yo] *= wr;
                        so_layer_im[a][osli + yo] *= wr;
                    }
                }
            }
        }
        for (int y = 0; y < 2; y++) {
            sf_re[a][y] = g->absorb[a].a_re[y] * pi_re[y] - g->absorb[a].a_im[y] * pi_im[y];
            sf_im[a][y] = g->absorb[a].a_re[y] * pi_im[y] + g->absorb[a].a_im[y] * pi_re[y];
        }
    }

    /* ── 3b: Find connected components in center-center graph ── */
    uint64_t *comp_of = (uint64_t *)calloc(n_absorb, sizeof(uint64_t));
    uint64_t n_comp = 0;

    /* Stack-based DFS */
    uint64_t *stack = (uint64_t *)calloc(n_absorb, sizeof(uint64_t));
    for (uint64_t a = 0; a < n_absorb; a++) {
        if (comp_of[a]) continue;
        n_comp++;
        comp_of[a] = n_comp;
        uint64_t sp = 0;
        stack[sp++] = a;
        while (sp) {
            uint64_t cur = stack[--sp];
            for (uint64_t k = 0; k < g->absorb[cur].n_nbrs; k++) {
                uint64_t nb = g->absorb[cur].nbrs[k];
                int64_t aj = g->absorb_idx[nb];
                if (aj < 0) continue;
                uint64_t aju = (uint64_t)aj;
                if (comp_of[aju]) continue;
                comp_of[aju] = n_comp;
                stack[sp++] = aju;
            }
        }
    }
    free(stack);

    /* ── 3c: Per-component joint evaluation ── */
    /* Build component member lists */
    uint64_t *comp_size = (uint64_t *)calloc(n_comp + 1, sizeof(uint64_t));
    for (uint64_t a = 0; a < n_absorb; a++) comp_size[comp_of[a]]++;
    uint64_t **comp_members = (uint64_t **)calloc(n_comp + 1, sizeof(uint64_t *));
    uint64_t *comp_cursor = (uint64_t *)calloc(n_comp + 1, sizeof(uint64_t));
    for (uint64_t c = 1; c <= n_comp; c++)
        comp_members[c] = (uint64_t *)calloc(comp_size[c], sizeof(uint64_t));
    for (uint64_t a = 0; a < n_absorb; a++) {
        uint64_t c = comp_of[a];
        comp_members[c][comp_cursor[c]++] = a;
    }
    free(comp_cursor);

    for (uint64_t c = 1; c <= n_comp; c++) {
        uint64_t sz = comp_size[c];
        uint64_t *mems = comp_members[c];

        if (sz == 1) {
            /* Single center — no other centers; evaluate H-chain directly */
            uint64_t a = mems[0];
            uint64_t xv = indices[g->absorb[a].center] ^ g->absorb[a].x_parity;
            int L = nl_arr[a];

            double cur_re[2] = {sf_re[a][0], sf_re[a][1]};
            double cur_im[2] = {sf_im[a][0], sf_im[a][1]};

            for (int li = 1; li < L; li++) {
                double nxt_re[2] = {0,0}, nxt_im[2] = {0,0};
                int use_z = (g->absorb[a].layer_x_parity &&
                             g->absorb[a].layer_x_parity[li-1]);
                for (int yo = 0; yo < 2; yo++) {
                    int sli = (li - 1) * 2 + yo;
                    double ar = g->absorb[a].a_layer_re[sli];
                    double ai = g->absorb[a].a_layer_im[sli];
                    double or = so_layer_re[a][sli];
                    double oi = so_layer_im[a][sli];
                    double sfr = ar * or - ai * oi;
                    double sfi = ar * oi + ai * or;
                    for (int yi = 0; yi < 2; yi++) {
                        /* H·X·H = Z = diag(1, -1): replace H with Z */
                        double M;
                        if (use_z) {
                            M = (yo == yi) ? (yo == 0 ? 1.0 : -1.0) : 0.0;
                        } else {
                            M = (yo == 0) ? SQ : (yi == 0 ? SQ : -SQ);
                        }
                        nxt_re[yo] += M * (sfr * cur_re[yi] - sfi * cur_im[yi]);
                        nxt_im[yo] += M * (sfr * cur_im[yi] + sfi * cur_re[yi]);
                    }
                }
                cur_re[0] = nxt_re[0]; cur_re[1] = nxt_re[1];
                cur_im[0] = nxt_im[0]; cur_im[1] = nxt_im[1];
            }
            double sum_re = 0, sum_im = 0;
            for (int y = 0; y < 2; y++) {
                double Hxy = (xv == 0) ? SQ : (y == 0 ? SQ : -SQ);
                sum_re += Hxy * cur_re[y];
                sum_im += Hxy * cur_im[y];
            }
            /* Multiply by local state at xv (captures T/S/Z after last H) */
            double lst_re[2], lst_im[2];
            tri_get_amplitudes((TrialityQubit *)&g->locals[g->absorb[a].center],
                              VIEW_EDGE, lst_re, lst_im);
            double lr = lst_re[xv], li = lst_im[xv];
            double sr = sum_re * lr - sum_im * li;
            double si = sum_re * li + sum_im * lr;
            sum_re = sr; sum_im = si;
            /* Apply non-absorbed CZ factors (L=1 center-regular edges) */
            double nre = sum_re * na_cz_re[a] - sum_im * na_cz_im[a];
            double nim = sum_re * na_cz_im[a] + sum_im * na_cz_re[a];
            double new_re = re * nre - im * nim;
            double new_im = re * nim + im * nre;
            re = new_re; im = new_im;
            continue;
        }

        /* sz > 1: Multi-center component exhaustive joint sum over all layers.
         * Each center i has L_i = nl_arr[mems[i]] variables y_0..y_{L_i-1}.
         * Total variables = sum(L_i).  Enumerate all 2^total_vars assignments. */
        uint64_t total_vars = 0;
        uint64_t *var_start = (uint64_t *)calloc(sz, sizeof(uint64_t));
        uint64_t *var_count = (uint64_t *)calloc(sz, sizeof(uint64_t));
        for (uint64_t mi = 0; mi < sz; mi++) {
            var_start[mi] = total_vars;
            var_count[mi] = (uint64_t)nl_arr[mems[mi]];
            total_vars += var_count[mi];
        }

        /* Collect non-absorbed CZ edges within this component */
        #define MAX_NA_CE 128
        uint64_t na_ce_a[MAX_NA_CE], na_ce_b[MAX_NA_CE], na_ce_e[MAX_NA_CE], n_na_ce = 0;
        for (uint64_t e = 0; e < g->n_edges && n_na_ce < MAX_NA_CE; e++) {
            const HPCQEdge *edge = &g->edges[e];
            if (edge->type != HPCQ_EDGE_CZ) continue;
            int64_t aia = g->absorb_idx[edge->site_a];
            int64_t aib = g->absorb_idx[edge->site_b];
            if (aia < 0 || aib < 0) continue;
            uint64_t aa = (uint64_t)aia, ab = (uint64_t)aib;
            if (comp_of[aa] != c || comp_of[ab] != c) continue;
            if (aa > ab) { uint64_t tmp = aa; aa = ab; ab = tmp; }
            int dup = 0;
            for (uint64_t ei = 0; ei < n_na_ce; ei++)
                if (na_ce_a[ei] == aa && na_ce_b[ei] == ab) { dup = 1; break; }
            if (!dup) { na_ce_a[n_na_ce] = aa; na_ce_b[n_na_ce] = ab; na_ce_e[n_na_ce] = e; n_na_ce++; }
        }

        /* ── Variable elimination for center-center component sum ──
         * Instead of enumerating 2^total_vars assignments, contract the factor
         * graph one variable at a time.  For a grid with treewidth W = min(R,C),
         * max intermediate factor size is 2^{W} ≪ 2^{total_vars}. */
        double comp_re = 0.0, comp_im = 0.0;
        if (total_vars <= 20) {
            /* Small components: fast-path with existing exhaustive enumeration */
            #define MAX_LAYERS 16
            uint64_t (*y_val)[MAX_LAYERS] = (uint64_t (*)[MAX_LAYERS])calloc(sz, MAX_LAYERS * sizeof(uint64_t));
            uint64_t n_assign = (uint64_t)1 << total_vars;
            for (uint64_t assign = 0; assign < n_assign; assign++) {
                double term_re = 1.0, term_im = 0.0;
                for (uint64_t mi = 0; mi < sz; mi++) {
                    uint64_t vs = var_start[mi];
                    uint64_t vc = var_count[mi];
                    for (uint64_t li = 0; li < vc; li++)
                        y_val[mi][li] = (assign >> (vs + li)) & 1;
                }
                for (uint64_t mi = 0; mi < sz; mi++) {
                    uint64_t a = mems[mi];
                    uint64_t xv = indices[g->absorb[a].center] ^ g->absorb[a].x_parity;
                    int L = nl_arr[a];
                    double factor_re, factor_im;
                    if (L >= 1) { factor_re = sf_re[a][y_val[mi][0]]; factor_im = sf_im[a][y_val[mi][0]]; }
                    else { factor_re = 1.0; factor_im = 0.0; }
                    for (int li = 1; li < L; li++) {
                        int use_z = (g->absorb[a].layer_x_parity &&
                                     g->absorb[a].layer_x_parity[li-1]);
                        double M;
                        if (use_z) {
                            M = (y_val[mi][li] == y_val[mi][li-1])
                                ? (y_val[mi][li] == 0 ? 1.0 : -1.0) : 0.0;
                        } else {
                            M = (y_val[mi][li] == 0)
                                ? SQ : (y_val[mi][li-1] == 0 ? SQ : -SQ);
                        }
                        double sfr, sfi;
                        {
                            int vv = (int)y_val[mi][li];
                            int sli = (li - 1) * 2 + vv;
                            double ar = g->absorb[a].a_layer_re[sli];
                            double ai = g->absorb[a].a_layer_im[sli];
                            double or = so_layer_re[a][sli];
                            double oi = so_layer_im[a][sli];
                            sfr = ar * or - ai * oi;
                            sfi = ar * oi + ai * or;
                        }
                        double new_re = factor_re * (M * sfr) - factor_im * (M * sfi);
                        double new_im = factor_re * (M * sfi) + factor_im * (M * sfr);
                        factor_re = new_re; factor_im = new_im;
                    }
                    if (L >= 1) {
                        double H_outer = (xv == 0) ? SQ : (y_val[mi][L-1] == 0 ? SQ : -SQ);
                        factor_re *= H_outer; factor_im *= H_outer;
                    }
                    double lst_re[2], lst_im[2];
                    tri_get_amplitudes((TrialityQubit *)&g->locals[g->absorb[a].center], VIEW_EDGE, lst_re, lst_im);
                    double lr = lst_re[xv], li = lst_im[xv];
                    double fr = factor_re * lr - factor_im * li;
                    double fi = factor_re * li + factor_im * lr;
                    factor_re = fr; factor_im = fi;
                    if (na_cz_re[a] != 1.0 || na_cz_im[a] != 0.0) {
                        double nr = factor_re * na_cz_re[a] - factor_im * na_cz_im[a];
                        double ni = factor_re * na_cz_im[a] + factor_im * na_cz_re[a];
                        factor_re = nr; factor_im = ni;
                    }
                    double new_re = term_re * factor_re - term_im * factor_im;
                    double new_im = term_re * factor_im + term_im * factor_re;
                    term_re = new_re; term_im = new_im;
                }
                for (uint64_t mi = 0; mi < sz; mi++) {
                    uint64_t a = mems[mi];
                    for (uint64_t k = 0; k < g->absorb[a].n_nbrs; k++) {
                        uint64_t nb = g->absorb[a].nbrs[k];
                        int64_t aj_idx = g->absorb_idx[nb];
                        if (aj_idx < 0) continue;
                        uint64_t aj = (uint64_t)aj_idx;
                        uint64_t mi2 = 0;
                        for (; mi2 < sz && mems[mi2] != aj; mi2++);
                        if (mi2 == sz) continue;
                        if (a >= aj) continue;
                        int li = (int)g->absorb[a].layer[k];
                        if (li >= (int)var_count[mi] || li >= (int)var_count[mi2]) continue;
                        uint64_t ca = g->absorb[a].center, cb = nb;
                        uint8_t xp_va = 0, xp_vb = 0;
                        for (uint64_t ee = 0; ee < g->n_edges; ee++) {
                            const HPCQEdge *ce = &g->edges[ee];
                            if ((ce->site_a == ca && ce->site_b == cb) || (ce->site_a == cb && ce->site_b == ca)) {
                                xp_va = (ce->site_a == ca) ? ce->xp_a : ce->xp_b;
                                xp_vb = (ce->site_a == ca) ? ce->xp_b : ce->xp_a;
                                break;
                            }
                        }
                        double wr = HPCQ_CZ_W(y_val[mi][li], y_val[mi2][li], xp_va, xp_vb);
                        term_re *= wr; term_im *= wr;
                    }
                }
                for (uint64_t ei = 0; ei < n_na_ce; ei++) {
                    uint64_t aa = na_ce_a[ei], ab = na_ce_b[ei], ee = na_ce_e[ei];
                    const HPCQEdge *edge = &g->edges[ee];
                    uint64_t mi_a, mi_b;
                    for (mi_a = 0; mi_a < sz && mems[mi_a] != aa; mi_a++);
                    for (mi_b = 0; mi_b < sz && mems[mi_b] != ab; mi_b++);
                    if (mi_a >= sz || mi_b >= sz) continue;
                    int L_a = nl_arr[aa], L_b = nl_arr[ab];
                    if (L_a < 1 || L_b < 1) continue;
                    uint64_t va = (L_a == 1) ? (indices[g->absorb[aa].center] ^ g->absorb[aa].x_parity) : y_val[mi_a][L_a - 1];
                    uint64_t vb = (L_b == 1) ? (indices[g->absorb[ab].center] ^ g->absorb[ab].x_parity) : y_val[mi_b][L_b - 1];
                    uint8_t xp_a = (edge->site_a == g->absorb[aa].center) ? edge->xp_a : edge->xp_b;
                    uint8_t xp_b = (edge->site_a == g->absorb[aa].center) ? edge->xp_b : edge->xp_a;
                    double wr = HPCQ_CZ_W(va, vb, xp_a, xp_b);
                    term_re *= wr; term_im *= wr;
                }
                comp_re += term_re;
                comp_im += term_im;
            }
            free(y_val);
        } else {
            /* Large components: variable elimination using min-degree ordering.
             * Build factors for each center's chain and each CZ edge, then
             * eliminate variables one at a time. */

            #define VE_MAX_VARS 128
            /* Maximum factor scope (variables in a factor). Factor storage
             * is 2^MAX_SCOPE entries → 4096 × 2 × 8 B = 64 KB per factor.
             * For grid-like center graphs, treewidth ≤ smaller grid dimension.
             * Runtime checks prevent scope overflow from memory corruption. */
            #define VE_MAX_SCOPE 12
            #define VE_MAX_NVALS (1 << VE_MAX_SCOPE)

            typedef struct { uint64_t vars[VE_MAX_SCOPE]; int n_vars; int n_vals; double re[VE_MAX_NVALS]; double im[VE_MAX_NVALS]; } VE_F;
            int ve_max_factors = 4096;
            VE_F *vf = (VE_F *)calloc(ve_max_factors, sizeof(VE_F)); int nvf = 0;
            #define VE_CHECK() do { if (nvf >= ve_max_factors) { \
                ve_max_factors *= 2; \
                vf = (VE_F *)realloc(vf, ve_max_factors * sizeof(VE_F)); \
                memset(&vf[nvf], 0, (ve_max_factors - nvf) * sizeof(VE_F)); \
            } } while(0)

            /* Helper: add a 1-variable factor over variable v with values v0, v1 */
            #define ve_add1(v, v0r, v0i, v1r, v1i) do { \
                VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F)); \
                f->vars[0]=v; f->n_vars=1; f->n_vals=2; \
                f->re[0]=(v0r); f->im[0]=(v0i); f->re[1]=(v1r); f->im[1]=(v1i); \
            } while(0)

            /* Helper: multiply two factors, store in vf[nvf] */
            #define ve_mul(fi, fj) do { \
                VE_CHECK(); VE_F *fa = &vf[(fi)], *fb = &vf[(fj)]; \
                VE_F *fc = &vf[nvf]; memset(fc,0,sizeof(VE_F)); nvf++; \
                int ni=0, ia=0, ib=0; \
                while (ia < fa->n_vars || ib < fb->n_vars) { \
                    if (ib>=fb->n_vars || (ia<fa->n_vars && fa->vars[ia]<fb->vars[ib])) \
                        fc->vars[ni++] = fa->vars[ia++]; \
                    else if (ia>=fa->n_vars || (ib<fb->n_vars && fb->vars[ib]<fa->vars[ia])) \
                        fc->vars[ni++] = fb->vars[ib++]; \
                    else { fc->vars[ni++] = fa->vars[ia]; ia++; ib++; } \
                } \
                fc->n_vars = ni; fc->n_vals = 1 << ni; \
                if (ni > VE_MAX_SCOPE) fc->n_vals = 0; \
                memset(fc->re, 0, fc->n_vals * sizeof(double)); \
                memset(fc->im, 0, fc->n_vals * sizeof(double)); \
                for (int a2 = 0; a2 < fb->n_vals; a2++) { \
                    uint64_t fb_assign = 0; \
                    for (int k = 0; k < fb->n_vars; k++) \
                        fb_assign |= ((a2 >> k) & 1) << fb->vars[k]; \
                    for (int a1 = 0; a1 < fa->n_vals; a1++) { \
                        int ok = 1; uint64_t fa_assign = 0; \
                        for (int k = 0; k < fa->n_vars; k++) { \
                            int bit = (a1 >> k) & 1; \
                            fa_assign |= bit << fa->vars[k]; \
                            for (int j = 0; j < fb->n_vars; j++) \
                                if (fa->vars[k] == fb->vars[j]) { \
                                    if (bit != ((fb_assign >> fb->vars[j]) & 1)) ok = 0; \
                                    break; \
                                } \
                        } \
                        if (!ok) continue; \
                        uint64_t full = fa_assign | fb_assign; \
                        uint64_t idx = 0; \
                        for (int k = 0; k < ni; k++) \
                            idx |= ((full >> fc->vars[k]) & 1) << k; \
                        fc->re[idx] += fa->re[a1] * fb->re[a2] - fa->im[a1] * fb->im[a2]; \
                        fc->im[idx] += fa->re[a1] * fb->im[a2] + fa->im[a1] * fb->re[a2]; \
                    } \
                } \
            } while(0)

            /* Build center chain factors */
            uint64_t *active = (uint64_t *)calloc(total_vars, sizeof(uint64_t));
            for (uint64_t mi = 0; mi < sz; mi++) {
                uint64_t a = mems[mi];
                uint64_t xv = indices[g->absorb[a].center] ^ g->absorb[a].x_parity;
                int L = nl_arr[a];
                uint64_t vs = var_start[mi];
                for (int li = 0; li < L; li++) active[vs + li] = 1;

                /* Factor over the last L-1 chain transitions plus outer Hadamard */
                /* Start with the last variable's factor */
                /* Instead of building one big factor, build chain of pair factors */
                /* Layer 0: a_re[y0] * pi_re[y0] (sf_re) */
                ve_add1(vs, sf_re[a][0], sf_im[a][0], sf_re[a][1], sf_im[a][1]);
                for (int li = 1; li < L; li++) {
                    uint64_t vp = vs + li - 1, vc = vs + li;
                    double sr0, si0, sr1, si1;
                    {
                        int sl0 = (li - 1) * 2;
                        int sl1 = sl0 + 1;
                        double ar0 = g->absorb[a].a_layer_re[sl0];
                        double ai0 = g->absorb[a].a_layer_im[sl0];
                        double ar1 = g->absorb[a].a_layer_re[sl1];
                        double ai1 = g->absorb[a].a_layer_im[sl1];
                        double or0 = so_layer_re[a][sl0];
                        double oi0 = so_layer_im[a][sl0];
                        double or1 = so_layer_re[a][sl1];
                        double oi1 = so_layer_im[a][sl1];
                        sr0 = ar0 * or0 - ai0 * oi0;
                        si0 = ar0 * oi0 + ai0 * or0;
                        sr1 = ar1 * or1 - ai1 * oi1;
                        si1 = ar1 * oi1 + ai1 * or1;
                    }
                    /* Build 2-variable factor for H[yc][yp] * sf[yc] * Z_parity[yp] */
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->vars[0] = vp; f->vars[1] = vc; f->n_vars = 2; f->n_vals = 4;
                    for (int yp = 0; yp < 2; yp++) {
                        for (int yc = 0; yc < 2; yc++) {
                            int use_z = (g->absorb[a].layer_x_parity &&
                                         g->absorb[a].layer_x_parity[li-1]);
                            double M;
                            if (use_z) {
                                M = (yc == yp) ? (yc == 0 ? 1.0 : -1.0) : 0.0;
                            } else {
                                M = (yc == 0) ? SQ : (yp == 0 ? SQ : -SQ);
                            }
                            double sf = (yc == 0) ? sr0 : sr1;
                            double si = (yc == 0) ? si0 : si1;
                            int idx = yp * 2 + yc;
                            f->re[idx] = M * sf;
                            f->im[idx] = M * si;
                        }
                    }
                }
                /* Outer Hadamard + local state + na_cz: H[xv][yL-1] * lst[xv] * na_cz */
                if (L >= 1) {
                    int li_last = L - 1;
                    uint64_t vl = vs + li_last;
                    double lst_re[2], lst_im[2];
                    tri_get_amplitudes((TrialityQubit *)&g->locals[g->absorb[a].center], VIEW_EDGE, lst_re, lst_im);
                    double lr0 = lst_re[0], li0 = lst_im[0], lr1 = lst_re[1], li1 = lst_im[1];
                    double ncr = na_cz_re[a], nci = na_cz_im[a];
                    double p0_re = lr0 * ncr - li0 * nci, p0_im = lr0 * nci + li0 * ncr;
                    double p1_re = lr1 * ncr - li1 * nci, p1_im = lr1 * nci + li1 * ncr;
                    if (xv == 0) {
                        ve_add1(vl,
                            SQ * p0_re, SQ * p0_im,
                            SQ * p0_re, SQ * p0_im);
                    } else {
                        ve_add1(vl,
                            SQ * p1_re, SQ * p1_im,
                            -SQ * p1_re, -SQ * p1_im);
                    }
                }
            }

            /* Build absorbed CZ pair factors */
            for (uint64_t mi = 0; mi < sz; mi++) {
                uint64_t a = mems[mi];
                for (uint64_t k = 0; k < g->absorb[a].n_nbrs; k++) {
                    uint64_t nb = g->absorb[a].nbrs[k];
                    int64_t aj_idx = g->absorb_idx[nb];
                    if (aj_idx < 0) continue;
                    uint64_t aj = (uint64_t)aj_idx;
                    uint64_t mi2 = 0;
                    for (; mi2 < sz && mems[mi2] != aj; mi2++);
                    if (mi2 == sz) continue;
                    if (a >= aj) continue;
                    int li = (int)g->absorb[a].layer[k];
                    if (li >= (int)var_count[mi] || li >= (int)var_count[mi2]) continue;
                    uint64_t va = var_start[mi] + li;
                    uint64_t vb = var_start[mi2] + li;
                    uint64_t ca = g->absorb[a].center, cb = nb;
                    uint8_t xp_va = 0, xp_vb = 0;
                    for (uint64_t ee = 0; ee < g->n_edges; ee++) {
                        const HPCQEdge *ce = &g->edges[ee];
                        if ((ce->site_a == ca && ce->site_b == cb) || (ce->site_a == cb && ce->site_b == ca)) {
                            xp_va = (ce->site_a == ca) ? ce->xp_a : ce->xp_b;
                            xp_vb = (ce->site_a == ca) ? ce->xp_b : ce->xp_a;
                            break;
                        }
                    }
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->vars[0] = va; f->vars[1] = vb; f->n_vars = 2; f->n_vals = 4;
                    for (int aa2 = 0; aa2 < 2; aa2++)
                        for (int bb2 = 0; bb2 < 2; bb2++) {
                            int idx = aa2 * 2 + bb2;
                            f->re[idx] = HPCQ_CZ_W(aa2, bb2, xp_va, xp_vb);
                            f->im[idx] = 0.0;
                        }
                }
            }

            /* Build non-absorbed CZ pair factors (na_ce) */
            for (uint64_t ei = 0; ei < n_na_ce; ei++) {
                uint64_t aa = na_ce_a[ei], ab = na_ce_b[ei], ee = na_ce_e[ei];
                const HPCQEdge *edge = &g->edges[ee];
                uint64_t mi_a, mi_b;
                for (mi_a = 0; mi_a < sz && mems[mi_a] != aa; mi_a++);
                for (mi_b = 0; mi_b < sz && mems[mi_b] != ab; mi_b++);
                if (mi_a >= sz || mi_b >= sz) continue;
                int L_a = nl_arr[aa], L_b = nl_arr[ab];
                if (L_a < 1 || L_b < 1) continue;
                uint64_t va = (L_a == 1) ? (uint64_t)-1 : var_start[mi_a] + L_a - 1;
                uint64_t vb = (L_b == 1) ? (uint64_t)-1 : var_start[mi_b] + L_b - 1;
                uint8_t xp_a = (edge->site_a == g->absorb[aa].center) ? edge->xp_a : edge->xp_b;
                uint8_t xp_b = (edge->site_a == g->absorb[aa].center) ? edge->xp_b : edge->xp_a;
                if (L_a == 1 && L_b == 1) {
                    /* Both L=1: depends on indices, not inner vars — constant factor */
                    uint64_t va_fixed = indices[g->absorb[aa].center] ^ g->absorb[aa].x_parity;
                    uint64_t vb_fixed = indices[g->absorb[ab].center] ^ g->absorb[ab].x_parity;
                    double wr = HPCQ_CZ_W(va_fixed, vb_fixed, xp_a, xp_b);
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->n_vars = 0; f->n_vals = 1;
                    f->re[0] = wr; f->im[0] = 0.0;
                } else if (L_a == 1) {
                    uint64_t va_fixed = indices[g->absorb[aa].center] ^ g->absorb[aa].x_parity;
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->vars[0] = vb; f->n_vars = 1; f->n_vals = 2;
                    for (int b = 0; b < 2; b++) {
                        f->re[b] = HPCQ_CZ_W(va_fixed, b, xp_a, xp_b);
                        f->im[b] = 0.0;
                    }
                } else if (L_b == 1) {
                    uint64_t vb_fixed = indices[g->absorb[ab].center] ^ g->absorb[ab].x_parity;
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->vars[0] = va; f->n_vars = 1; f->n_vals = 2;
                    for (int a = 0; a < 2; a++) {
                        f->re[a] = HPCQ_CZ_W(a, vb_fixed, xp_a, xp_b);
                        f->im[a] = 0.0;
                    }
                } else {
                    VE_CHECK(); VE_F *f = &vf[nvf++]; memset(f,0,sizeof(VE_F));
                    f->vars[0] = va; f->vars[1] = vb; f->n_vars = 2; f->n_vals = 4;
                    for (int aa3 = 0; aa3 < 2; aa3++)
                        for (int bb3 = 0; bb3 < 2; bb3++) {
                            int idx = aa3 * 2 + bb3;
                            f->re[idx] = HPCQ_CZ_W(aa3, bb3, xp_a, xp_b);
                            f->im[idx] = 0.0;
                        }
                }
            }

            /* Eliminate variables using min-degree ordering */
            int *elim_order = (int *)calloc(total_vars, sizeof(int));
            int *var_active = (int *)calloc(total_vars, sizeof(int));
            for (uint64_t i = 0; i < total_vars; i++) var_active[i] = (int)active[i];
            free(active);

            int n_elim = 0;
            for (uint64_t i = 0; i < total_vars; i++) {
                if (!var_active[i]) continue;
                /* Find variable with minimum factor scope product (min-degree heuristic) */
                int best_v = -1; uint64_t best_cost = (uint64_t)-1;
                for (uint64_t v = 0; v < total_vars; v++) {
                    if (!var_active[v]) continue;
                    uint64_t cost = 1;
                    for (int fi = 0; fi < nvf; fi++) {
                        for (int k = 0; k < vf[fi].n_vars; k++)
                            if (vf[fi].vars[k] == v) {
                                cost *= (uint64_t)vf[fi].n_vals;
                                break;
                            }
                    }
                    if (cost < best_cost) { best_cost = cost; best_v = (int)v; }
                }
                if (best_v < 0) break;
                elim_order[n_elim++] = best_v;
                int v = best_v;
                var_active[v] = 0;

                /* Find all factors containing v, multiply them, sum out v */
                int first = -1;
                for (int fi = 0; fi < nvf; fi++) {
                    int found = 0;
                    for (int k = 0; k < vf[fi].n_vars; k++)
                        if (vf[fi].vars[k] == (uint64_t)v) { found = 1; break; }
                    if (!found) continue;
                    if (first < 0) { first = fi; continue; }
                    /* Multiply vf[fi] into vf[first] */
                    int sav_nvf = nvf;
                    /* Use temporary factor */
                    VE_F tmp = vf[first];
                    /* Recompute first = tmp * vf[fi] */
                    VE_F *fa = &tmp, *fb = &vf[fi];
                    VE_F *fc = &vf[first]; memset(fc,0,sizeof(VE_F));
                    int ni=0, ia=0, ib=0;
                    while (ia < fa->n_vars || ib < fb->n_vars) {
                        if (ib>=fb->n_vars || (ia<fa->n_vars && fa->vars[ia]<fb->vars[ib]))
                            fc->vars[ni++] = fa->vars[ia++];
                        else if (ia>=fa->n_vars || (ib<fb->n_vars && fb->vars[ib]<fa->vars[ia]))
                            fc->vars[ni++] = fb->vars[ib++];
                        else { fc->vars[ni++] = fa->vars[ia]; ia++; ib++; }
                    }
                    fc->n_vars = ni; fc->n_vals = 1 << ni;
                    memset(fc->re, 0, fc->n_vals * sizeof(double));
                    memset(fc->im, 0, fc->n_vals * sizeof(double));
                    for (int a2 = 0; a2 < fb->n_vals; a2++) {
                        uint64_t fb_assign = 0;
                        for (int k = 0; k < fb->n_vars; k++)
                            fb_assign |= ((a2 >> k) & 1) << fb->vars[k];
                        for (int a1 = 0; a1 < fa->n_vals; a1++) {
                            int ok = 1; uint64_t fa_assign = 0;
                            for (int k = 0; k < fa->n_vars; k++) {
                                int bit = (a1 >> k) & 1;
                                fa_assign |= bit << fa->vars[k];
                                for (int j = 0; j < fb->n_vars; j++)
                                    if (fa->vars[k] == fb->vars[j]) {
                                        if (bit != ((fb_assign >> fb->vars[j]) & 1)) ok = 0;
                                        break;
                                    }
                            }
                            if (!ok) continue;
                            uint64_t full = fa_assign | fb_assign;
                            uint64_t idx = 0;
                            for (int k = 0; k < ni; k++)
                                idx |= ((full >> fc->vars[k]) & 1) << k;
                            fc->re[idx] += fa->re[a1] * fb->re[a2] - fa->im[a1] * fb->im[a2];
                            fc->im[idx] += fa->re[a1] * fb->im[a2] + fa->im[a1] * fb->re[a2];
                        }
                    }
                    /* Mark vf[fi] as inactive by swapping with last */
                    /* For simplicity, just mark it; we'll compact later if needed */
                    vf[fi].n_vars = 0; /* invalid */
                }
                /* Now sum out v from vf[first] */
                if (first >= 0) {
                    VE_F *f = &vf[first];
                    VE_F tmp = *f;
                    int pos = -1;
                    for (int k = 0; k < tmp.n_vars; k++)
                        if (tmp.vars[k] == (uint64_t)v) { pos = k; break; }
                    if (pos >= 0) {
                        int new_n = tmp.n_vars - 1;
                        int ni = 0;
                        for (int k = 0; k < tmp.n_vars; k++)
                            if (k != pos) f->vars[ni++] = tmp.vars[k];
                        f->n_vars = new_n;
                        f->n_vals = new_n > 0 ? (1 << new_n) : 1;
                        memset(f->re, 0, f->n_vals * sizeof(double));
                        memset(f->im, 0, f->n_vals * sizeof(double));
                        int full_vals = 1 << tmp.n_vars;
                        for (int fi = 0; fi < full_vals; fi++) {
                            uint64_t out_idx = 0;
                            for (int k = 0; k < tmp.n_vars; k++) {
                                if (k == pos) continue;
                                int bit = (fi >> k) & 1;
                                int dest = k < pos ? k : k - 1;
                                out_idx |= bit << dest;
                            }
                            int out_ni = new_n > 0 ? new_n : 1;
                            if (out_ni == 1) {
                                f->re[0] += tmp.re[fi];
                                f->im[0] += tmp.im[fi];
                            } else {
                                f->re[out_idx] += tmp.re[fi];
                                f->im[out_idx] += tmp.im[fi];
                            }
                        }
                        if (new_n == 0) { f->n_vals = 1; }
                    }
                }
                /* Compact: remove invalid factors (keep constants) */
                int w = 0;
                for (int r = 0; r < nvf; r++)
                    if (vf[r].n_vars > 0 || r == first || (vf[r].n_vars == 0 && vf[r].n_vals == 1)) {
                        if (w != r) vf[w] = vf[r];
                        w++;
                    }
                nvf = w;
            }
            free(elim_order);
            free(var_active);

            /* The result scalar is in the remaining factor(s) */
            comp_re = 1.0; comp_im = 0.0;
            for (int fi = 0; fi < nvf; fi++) {
                double nr = comp_re * vf[fi].re[0] - comp_im * vf[fi].im[0];
                double ni = comp_re * vf[fi].im[0] + comp_im * vf[fi].re[0];
                comp_re = nr; comp_im = ni;
            }
            free(vf);
            #undef ve_add1
            #undef ve_mul
        }
        free(var_start);
        free(var_count);

        double new_re = re * comp_re - im * comp_im;
        double new_im = re * comp_im + im * comp_re;
        re = new_re; im = new_im;
    }

    /* Free temp arrays */
    for (uint64_t c = 1; c <= n_comp; c++) free(comp_members[c]);
    free(comp_members);
    free(comp_size);
    free(comp_of);
    free(sf_re); free(sf_im);
    for (uint64_t a = 0; a < n_absorb; a++) { free(so_layer_re[a]); free(so_layer_im[a]); }
    free(so_layer_re); free(so_layer_im);
    free(na_cz_re); free(na_cz_im);
    free(nl_arr);

    if (g->global_phase_parity & 1) { re = -re; im = -im; }
    *out_re = re;
    *out_im = im;
    ((HPCQGraph *)g)->amp_evals++;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PROBABILITY — |ψ(i₁,...,iₙ)|²
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double hpcq_probability(const HPCQGraph *g,
                                       const uint32_t *indices)
{
    double re, im;
    hpcq_amplitude(g, indices, &re, &im);
    ((HPCQGraph *)g)->prob_evals++;
    return re * re + im * im;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MARGINAL PROBABILITY — P(site_k = v)
 *
 * Sums |ψ(..., i_k=v, ...)|² over connected partner configurations.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double hpcq_marginal(const HPCQGraph *g,
                                    uint64_t site, uint32_t value)
{
    /* Find connected sites */
    uint64_t connected[128];
    uint64_t n_connected = 0;

    for (uint64_t e = 0; e < g->n_edges; e++) {
        uint64_t sa = g->edges[e].site_a;
        uint64_t sb = g->edges[e].site_b;
        if (sa == site || sb == site) {
            uint64_t partner = (sa == site) ? sb : sa;
            int found = 0;
            for (uint64_t c = 0; c < n_connected; c++)
                if (connected[c] == partner) { found = 1; break; }
            if (!found && n_connected < 128)
                connected[n_connected++] = partner;
        }
    }

    /* Product state: no edges touching this site */
    if (n_connected == 0) {
        double re_buf[HPCQ_D], im_buf[HPCQ_D];
        tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE,
                          re_buf, im_buf);
        return re_buf[value] * re_buf[value] + im_buf[value] * im_buf[value];
    }

    /* Entangled: enumerate 2^n_connected configurations */
    double total_prob = 0.0;
    uint64_t n_configs = 1ULL << n_connected;  /* 2^n_connected */

    for (uint64_t cfg = 0; cfg < n_configs; cfg++) {
        uint32_t partner_vals[128];
        for (uint64_t c = 0; c < n_connected; c++)
            partner_vals[c] = (cfg >> c) & 1;

        /* Local amplitude for measured site */
        double re_buf[HPCQ_D], im_buf[HPCQ_D];
        tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE,
                          re_buf, im_buf);
        double amp_re = re_buf[value];
        double amp_im = im_buf[value];

        /* Partner amplitudes */
        for (uint64_t c = 0; c < n_connected; c++) {
            double p_re_buf[HPCQ_D], p_im_buf[HPCQ_D];
            tri_get_amplitudes((TrialityQubit *)&g->locals[connected[c]],
                              VIEW_EDGE, p_re_buf, p_im_buf);
            double p_re = p_re_buf[partner_vals[c]];
            double p_im = p_im_buf[partner_vals[c]];
            double new_re = amp_re * p_re - amp_im * p_im;
            double new_im = amp_re * p_im + amp_im * p_re;
            amp_re = new_re;
            amp_im = new_im;
        }

        /* Phase contributions from edges */
        for (uint64_t e = 0; e < g->n_edges; e++) {
            uint64_t sa = g->edges[e].site_a;
            uint64_t sb = g->edges[e].site_b;

            uint32_t va = 0, vb = 0;
            int involves = 0;

            if (sa == site) {
                va = value;
                for (uint64_t c = 0; c < n_connected; c++)
                    if (connected[c] == sb) { vb = partner_vals[c]; involves = 1; break; }
            } else if (sb == site) {
                vb = value;
                for (uint64_t c = 0; c < n_connected; c++)
                    if (connected[c] == sa) { va = partner_vals[c]; involves = 1; break; }
            } else {
                int found_a = 0, found_b = 0;
                for (uint64_t c = 0; c < n_connected; c++) {
                    if (connected[c] == sa) { va = partner_vals[c]; found_a = 1; }
                    if (connected[c] == sb) { vb = partner_vals[c]; found_b = 1; }
                }
                if (found_a && found_b) involves = 1;
            }

            if (!involves) continue;

            double w_re, w_im;
            if (g->edges[e].type == HPCQ_EDGE_CZ) {
                w_re = HPCQ_CZ_W(va, vb, g->edges[e].xp_a, g->edges[e].xp_b);
                w_im = 0.0;
            } else {
                w_re = g->edges[e].w_re[va][vb];
                w_im = g->edges[e].w_im[va][vb];
            }

            double new_re = amp_re * w_re - amp_im * w_im;
            double new_im = amp_re * w_im + amp_im * w_re;
            amp_re = new_re;
            amp_im = new_im;
        }

        total_prob += amp_re * amp_re + amp_im * amp_im;
    }

    return total_prob;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BORN SAMPLING — Measure site k
 *
 * For D=2: compute P(0), compare against random.
 * Collapses local state to |outcome⟩, absorbs phases into partners,
 * removes resolved edges. This IS measurement-induced disentanglement.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline uint32_t hpcq_measure(HPCQGraph *g, uint64_t site,
                                     double random_01)
{
    /* Marginals */
    double probs[HPCQ_D];
    double total = 0.0;
    for (int v = 0; v < HPCQ_D; v++) {
        probs[v] = hpcq_marginal(g, site, v);
        total += probs[v];
    }
    if (total > 0) {
        probs[0] /= total;
        probs[1] /= total;
    }

    /* Sample */
    uint32_t outcome = (random_01 < probs[0]) ? 0 : 1;

    /* Collapse local state */
    double init_re[2] = {0, 0}, init_im[2] = {0, 0};
    init_re[outcome] = 1.0;
    tri_init_state(&g->locals[site], VIEW_EDGE, init_re, init_im);

    /* Absorb edge phases into partners and remove resolved edges */
    for (uint64_t e = 0; e < g->n_edges; ) {
        HPCQEdge *edge = &g->edges[e];
        if (edge->site_a == site || edge->site_b == site) {
            uint64_t partner = (edge->site_a == site) ?
                                edge->site_b : edge->site_a;
            TrialityQubit *p = &g->locals[partner];

            /* Absorb the phase */
            double p_re[HPCQ_D], p_im[HPCQ_D];
            tri_get_amplitudes(p, VIEW_EDGE, p_re, p_im);

            for (int k = 0; k < HPCQ_D; k++) {
                double w_re, w_im;
                if (edge->type == HPCQ_EDGE_CZ) {
                    uint8_t xp_s = (edge->site_a == site) ? edge->xp_a : edge->xp_b;
                    uint8_t xp_p = (edge->site_a == site) ? edge->xp_b : edge->xp_a;
                    w_re = HPCQ_CZ_W(outcome, k, xp_s, xp_p);
                    w_im = 0.0;
                } else if (edge->site_a == site) {
                    w_re = edge->w_re[outcome][k];
                    w_im = edge->w_im[outcome][k];
                } else {
                    w_re = edge->w_re[k][outcome];
                    w_im = edge->w_im[k][outcome];
                }

                double old_re = p_re[k], old_im = p_im[k];
                p_re[k] = old_re * w_re - old_im * w_im;
                p_im[k] = old_re * w_im + old_im * w_re;
            }

            tri_init_state(p, VIEW_EDGE, p_re, p_im);

            /* Track edge type removal */
            if (edge->type == HPCQ_EDGE_CZ) g->cz_edges--;
            else if (edge->type == HPCQ_EDGE_PHASE) g->phase_edges--;
            else g->clifford_edges--;

            /* Swap-remove */
            g->edges[e] = g->edges[--g->n_edges];
        } else {
            e++;
        }
    }

    g->measurements++;
    hpcq_update_fidelity_stats(g);
    return outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * NORM CHECK — Σ|ψ|² over ALL indices (small N only)
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double hpcq_norm_sq(const HPCQGraph *g)
{
    if (g->n_sites > 20) {
        return -1.0;
    }

    uint64_t total_configs = 1ULL << g->n_sites;
    double norm = 0.0;
    uint32_t indices[20];

    for (uint64_t cfg = 0; cfg < total_configs; cfg++) {
        for (uint64_t i = 0; i < g->n_sites; i++)
            indices[i] = (cfg >> i) & 1;
        norm += hpcq_probability(g, indices);
    }
    return norm;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENTROPY ESTIMATE — across a bipartition cut
 *
 * CZ edges contribute exactly 1 bit per crossing edge.
 * General edges contribute fidelity-weighted 1 bit.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline double hpcq_entropy_cut(const HPCQGraph *g, uint64_t cut_after)
{
    double entropy = 0.0;
    for (uint64_t e = 0; e < g->n_edges; e++) {
        uint64_t sa = g->edges[e].site_a;
        uint64_t sb = g->edges[e].site_b;
        if ((sa <= cut_after && sb > cut_after) ||
            (sb <= cut_after && sa > cut_after)) {
            entropy += g->edges[e].fidelity * 1.0; /* log₂(2) = 1 bit */
        }
    }
    return entropy;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DIAGNOSTICS
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_print_stats(const HPCQGraph *g)
{
    printf("╔═════════════════════════════════════════════════════╗\n");
    printf("║  Qubit Phase Graph Statistics                      ║\n");
    printf("╠═════════════════════════════════════════════════════╣\n");
    printf("║  Sites:           %10lu                       ║\n", g->n_sites);
    printf("║  Total edges:     %10lu                       ║\n", g->n_edges);
    printf("║    CZ (exact):    %10lu                       ║\n", g->cz_edges);
    printf("║    Phase (lossy): %10lu                       ║\n", g->phase_edges);
    printf("║    Clifford:      %10lu                       ║\n", g->clifford_edges);
    printf("║  Amp evals:       %10lu                       ║\n", g->amp_evals);
    printf("║  Measurements:    %10lu                       ║\n", g->measurements);
    printf("║  Min fidelity:    %10.6f                       ║\n", g->min_fidelity);
    printf("║  Avg fidelity:    %10.6f                       ║\n", g->avg_fidelity);

    uint64_t mem_bytes = g->n_sites * sizeof(TrialityQubit) +
                         g->n_edges * sizeof(HPCQEdge) +
                         g->n_log * sizeof(HPCQGateEntry) +
                         sizeof(HPCQGraph);
    printf("║  Memory:          %10lu bytes                ║\n", mem_bytes);

    double full_sv_log = g->n_sites * log10(2.0) + log10(16.0);
    printf("║  Full SV:         10^%.1f bytes                   ║\n", full_sv_log);
    printf("╚═════════════════════════════════════════════════════╝\n");
}

#endif /* HPC_QUBIT_GRAPH_H */
