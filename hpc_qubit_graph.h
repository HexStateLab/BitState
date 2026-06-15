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
 * This correctly handles new CZ edges (which use the outer variable q)
 * alongside the original absorbed edges (which use the inner variable y).
 * Diagonal gates (T, S) between the two H operations modify a_cur. */
typedef struct {
    uint64_t  center;
    uint64_t  n_nbrs;
    uint64_t  n_inner;   /* first n_inner neighbors are inner group */
    uint64_t *nbrs;
    double   *w_re;       /* [n_nbrs * 4] — inner then outer */
    double   *w_im;
    double    a_re[2];    /* Z-basis amplitude at FIRST absorption time */
    double    a_im[2];
    int       n_layers;   /* 1 or 2 */
    double    a_cur_re[2];/* pre-H local state for second layer */
    double    a_cur_im[2];
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

    /* Statistics */
    uint64_t        amp_evals;
    uint64_t        prob_evals;
    uint64_t        measurements;
    uint64_t        cz_edges;
    uint64_t        phase_edges;
    uint64_t        clifford_edges;
    double          min_fidelity;
    double          avg_fidelity;
} HPCQGraph;

/* ═══════════════════════════════════════════════════════════════════════════════
 * LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

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
    }
    free(g->absorb);
    free(g->locals);
    free(g->edges);
    free(g->gate_log);
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

    /* Count incident non-absorbed edges */
    int n_inc = 0;
    uint64_t e_idx = 0;

    for (uint64_t e = 0; e < g->n_edges; e++) {
        HPCQEdge *edge = &g->edges[e];
        if (edge->type == HPCQ_EDGE_ABSORBED) continue;
        if (edge->site_a == site || edge->site_b == site) {
            n_inc++;
            e_idx = e;
        }
    }

    /* Check if this site already has an absorb entry (re-apply H) */
    int existing_absorb = -1;
    for (uint64_t a = 0; a < g->n_absorb; a++) {
        if (g->absorb[a].center == site) { existing_absorb = (int)a; break; }
    }

    /* ─── Re-absorption: H gate applied again on an already-absorbed site ─── */
    if (existing_absorb >= 0) {
        HPCQAbsorbEntry *old = &g->absorb[existing_absorb];
        double cur_re[2], cur_im[2];
        tri_get_amplitudes(&g->locals[site], VIEW_EDGE, cur_re, cur_im);

        if (n_inc == 0) {
            /* No new edges: just record the second H layer */
            old->n_layers = 2;
            old->a_cur_re[0] = cur_re[0];
            old->a_cur_re[1] = cur_re[1];
            old->a_cur_im[0] = cur_im[0];
            old->a_cur_im[1] = cur_im[1];
            double uni[2] = {1.0, 1.0}, zero[2] = {0.0, 0.0};
            tri_init_state(&g->locals[site], VIEW_EDGE, uni, zero);
            HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                                    .fidelity = 1.0 };
            hpcq_log_gate(g, entry);
            return;
        }

        /* Merge: old edges stay inner, new edges go to outer group.
         * Dedup only within the outer group (same outer variable q).
         * Duplicates across inner/outer are correct — different variables. */
        uint64_t inner_start = old->n_inner;
        for (uint64_t e = 0; e < g->n_edges; e++) {
            HPCQEdge *edge = &g->edges[e];
            if (edge->type == HPCQ_EDGE_ABSORBED) continue;
            if (edge->site_a == site || edge->site_b == site) n_inc++;
        }
        uint64_t total_max = old->n_nbrs + (uint64_t)n_inc;
        uint64_t *new_nbrs = (uint64_t *)calloc(total_max, sizeof(uint64_t));
        double *new_w_re = (double *)calloc(total_max * 4, sizeof(double));
        double *new_w_im = (double *)calloc(total_max * 4, sizeof(double));
        uint64_t n_nbrs = old->n_nbrs;

        memcpy(new_nbrs, old->nbrs, old->n_nbrs * sizeof(uint64_t));
        memcpy(new_w_re, old->w_re, old->n_nbrs * 4 * sizeof(double));
        memcpy(new_w_im, old->w_im, old->n_nbrs * 4 * sizeof(double));

        for (uint64_t e = 0; e < g->n_edges; e++) {
            HPCQEdge *edge = &g->edges[e];
            if (edge->type == HPCQ_EDGE_ABSORBED) continue;
            if (edge->site_a == site || edge->site_b == site) {
                uint64_t p = (edge->site_a == site) ? edge->site_b : edge->site_a;
                double e_re[2][2], e_im[2][2];
                if (edge->type == HPCQ_EDGE_CZ) {
                    for (int q = 0; q < 2; q++)
                        for (int z = 0; z < 2; z++) {
                            e_re[q][z] = ((q * z) % 2 == 0) ? 1.0 : -1.0;
                            e_im[q][z] = 0.0;
                        }
                } else {
                    for (int q = 0; q < 2; q++)
                        for (int z = 0; z < 2; z++) {
                            e_re[q][z] = edge->w_re[q][z];
                            e_im[q][z] = edge->w_im[q][z];
                        }
                }
                /* Dedup only against outer group (uses same variable q) */
                int found = -1;
                for (uint64_t k = inner_start; k < n_nbrs; k++)
                    if (new_nbrs[k] == p) { found = (int)k; break; }
                if (found >= 0) {
                    for (int q = 0; q < 2; q++)
                        for (int z = 0; z < 2; z++) {
                            int idx = found * 4 + q * 2 + z;
                            double wr = new_w_re[idx], wi = new_w_im[idx];
                            new_w_re[idx] = wr * e_re[q][z] - wi * e_im[q][z];
                            new_w_im[idx] = wr * e_im[q][z] + wi * e_re[q][z];
                        }
                } else {
                    new_nbrs[n_nbrs] = p;
                    for (int q = 0; q < 2; q++)
                        for (int z = 0; z < 2; z++) {
                            int idx = n_nbrs * 4 + q * 2 + z;
                            new_w_re[idx] = e_re[q][z];
                            new_w_im[idx] = e_im[q][z];
                        }
                    n_nbrs++;
                }
                edge->type = HPCQ_EDGE_ABSORBED;
            }
        }

        free(old->nbrs);
        free(old->w_re);
        free(old->w_im);
        old->n_nbrs = n_nbrs;
        old->n_inner = old->n_inner;  /* preserved */
        old->nbrs = new_nbrs;
        old->w_re = new_w_re;
        old->w_im = new_w_im;
        old->n_layers = 2;
        old->a_cur_re[0] = cur_re[0];
        old->a_cur_re[1] = cur_re[1];
        old->a_cur_im[0] = cur_im[0];
        old->a_cur_im[1] = cur_im[1];

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

    if (n_inc > 1) {
        /* Multi-edge absorption: for each incident edge, store the edge matrix
         * in y-z form where y is the summed-over pre-H index of the center site
         * and z is the neighbor's bit value.
         *   - CZ edges: w(y, z) = (-1)^{y·z}
         *   - PHASE edges: use the stored w_re[y][z] directly (the first index
         *     IS the site's own value — post-absorb = pre-next-absorb) */
        uint64_t *nbrs = (uint64_t *)calloc(n_inc, sizeof(uint64_t));
        double *w_re = (double *)calloc(n_inc * 4, sizeof(double));
        double *w_im = (double *)calloc(n_inc * 4, sizeof(double));
        uint64_t n_nbrs = 0;

        for (uint64_t e = 0; e < g->n_edges; e++) {
            HPCQEdge *edge = &g->edges[e];
            if (edge->type == HPCQ_EDGE_ABSORBED) continue;
            if (edge->site_a == site || edge->site_b == site) {
                uint64_t p = (edge->site_a == site) ? edge->site_b : edge->site_a;

                double e_re[2][2], e_im[2][2];
                if (edge->type == HPCQ_EDGE_CZ) {
                    for (int y = 0; y < 2; y++)
                        for (int z = 0; z < 2; z++) {
                            e_re[y][z] = ((y * z) % 2 == 0) ? 1.0 : -1.0;
                            e_im[y][z] = 0.0;
                        }
                } else {
                    for (int y = 0; y < 2; y++)
                        for (int z = 0; z < 2; z++) {
                            e_re[y][z] = edge->w_re[y][z];
                            e_im[y][z] = edge->w_im[y][z];
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

                edge->type = HPCQ_EDGE_ABSORBED;
            }
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

        g->absorb[g->n_absorb].center = site;
        g->absorb[g->n_absorb].n_nbrs = n_nbrs;
        g->absorb[g->n_absorb].n_inner = n_nbrs;
        g->absorb[g->n_absorb].nbrs = nbrs;
        g->absorb[g->n_absorb].w_re = w_re;
        g->absorb[g->n_absorb].w_im = w_im;
        g->absorb[g->n_absorb].a_re[0] = orig_re[0];
        g->absorb[g->n_absorb].a_re[1] = orig_re[1];
        g->absorb[g->n_absorb].a_im[0] = orig_im[0];
        g->absorb[g->n_absorb].a_im[1] = orig_im[1];
        g->absorb[g->n_absorb].n_layers = 1;
        g->n_absorb++;

        HPCQGateEntry entry = { .type = HPCQ_GATE_LOCAL_H, .site_a = site,
                                .fidelity = 1.0 };
        hpcq_log_gate(g, entry);
        return;
    }

    /* Single-edge absorption (first time, exactly 1 incident edge) */
    double a_re[2], a_im[2];
    tri_get_amplitudes(&g->locals[site], VIEW_EDGE, a_re, a_im);

    HPCQEdge *edge = &g->edges[e_idx];
    uint64_t partner = (edge->site_a == site) ? edge->site_b : edge->site_a;

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
                    uint32_t pi = (y * xj) % 2;
                    w_re_yj = (pi == 0) ? 1.0 : -1.0;
                    w_im_yj = 0.0;
                } else {
                    w_re_yj = edge->w_re[y][xj];
                    w_im_yj = edge->w_im[y][xj];
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

    edge->type = HPCQ_EDGE_PHASE;
    edge->fidelity = 1.0;
    memcpy(edge->w_re, wr, sizeof(wr));
    memcpy(edge->w_im, wi, sizeof(wi));

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
    tri_apply_x(&g->locals[site]);
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
    /* Check for existing CZ edge between this pair — CZ² = I cancellation */
    for (uint64_t e = 0; e < g->n_edges; e++) {
        HPCQEdge *edge = &g->edges[e];
        if (edge->type == HPCQ_EDGE_CZ &&
            ((edge->site_a == site_a && edge->site_b == site_b) ||
             (edge->site_a == site_b && edge->site_b == site_a))) {
            /* Cancel: swap-remove this edge */
            g->edges[e] = g->edges[--g->n_edges];
            g->cz_edges--;

            HPCQGateEntry entry = {
                .type = HPCQ_GATE_CZ,
                .site_a = site_a, .site_b = site_b,
                .fidelity = 1.0
            };
            hpcq_log_gate(g, entry);
            return;
        }
    }

    /* No existing edge — add new one */
    hpcq_grow_edges(g);

    HPCQEdge *e = &g->edges[g->n_edges];
    memset(e, 0, sizeof(HPCQEdge));
    e->type = HPCQ_EDGE_CZ;
    e->site_a = site_a;
    e->site_b = site_b;
    e->fidelity = 1.0;

    g->n_edges++;
    g->cz_edges++;

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

    /* Step 1: Product of local amplitudes — O(N) */
    for (uint64_t k = 0; k < g->n_sites; k++) {
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

    /* Step 2: Phase edge accumulation — O(E) */
    for (uint64_t e = 0; e < g->n_edges; e++) {
        const HPCQEdge *edge = &g->edges[e];

        /* Skip edges consumed by multi-edge absorption */
        if (edge->type == HPCQ_EDGE_ABSORBED) continue;

        uint32_t ia = indices[edge->site_a];
        uint32_t ib = indices[edge->site_b];

        double w_re, w_im;

        if (edge->type == HPCQ_EDGE_CZ) {
            /* CZ: (-1)^(ia·ib) — only -1 when both are 1 */
            uint32_t phase_idx = (ia * ib) % HPCQ_D;
            w_re = HPCQ_W2_RE[phase_idx];
            w_im = HPCQ_W2_IM[phase_idx];
        } else {
            w_re = edge->w_re[ia][ib];
            w_im = edge->w_im[ia][ib];
        }

        double new_re = re * w_re - im * w_im;
        double new_im = re * w_im + im * w_re;
        re = new_re;
        im = new_im;
    }

    /* Step 3: Multi-edge absorb correction
     * For each absorbed center:
     *   ψ_v(x_v, z_1..z_m) = Σ_y H_{x_v,y} · a_original(y) · Π_k w_k(y, z_k)
     * where w_k(y, z_k) is the y-form edge matrix for neighbor k.
     * The local state already contributes a'_v(x_v) (post-absorb * diagonal gates).
     * Total center factor = a'_v(x_v) × ψ_v(x_v, ...). */
    static const double SQ = 0.7071067811865475244;
    for (uint64_t a = 0; a < g->n_absorb; a++) {
        uint64_t center = g->absorb[a].center;
        int xv = (int)indices[center];
        uint64_t n_inner = g->absorb[a].n_inner;

        /* Inner product: Π w_k(y, z_k) for y = 0, 1 (k < n_inner) */
        double pi_re[2] = {1.0, 1.0};
        double pi_im[2] = {0.0, 0.0};
        for (uint64_t k = 0; k < n_inner; k++) {
            int z = (int)indices[g->absorb[a].nbrs[k]];
            for (int y = 0; y < 2; y++) {
                int idx = (int)k * 4 + y * 2 + z;
                double wr = g->absorb[a].w_re[idx];
                double wi = g->absorb[a].w_im[idx];
                double pr = pi_re[y], pim = pi_im[y];
                pi_re[y] = pr * wr - pim * wi;
                pi_im[y] = pr * wi + pim * wr;
            }
        }

        if (g->absorb[a].n_layers >= 2) {
            /* Outer product: Π w_k(q, z_k) for q = 0, 1 (k >= n_inner) */
            double po_re[2] = {1.0, 1.0};
            double po_im[2] = {0.0, 0.0};
            for (uint64_t k = n_inner; k < g->absorb[a].n_nbrs; k++) {
                int z = (int)indices[g->absorb[a].nbrs[k]];
                for (int q = 0; q < 2; q++) {
                    int idx = (int)k * 4 + q * 2 + z;
                    double wr = g->absorb[a].w_re[idx];
                    double wi = g->absorb[a].w_im[idx];
                    double pr = po_re[q], pim = po_im[q];
                    po_re[q] = pr * wr - pim * wi;
                    po_im[q] = pr * wi + pim * wr;
                }
            }

            /* mid[q] = Σ_y H[q][y] · a_re[y] · pi[y] */
            double mid_re[2] = {0.0, 0.0}, mid_im[2] = {0.0, 0.0};
            for (int q = 0; q < 2; q++) {
                for (int y = 0; y < 2; y++) {
                    double Hqy = (q == 0) ? SQ : (y == 0 ? SQ : -SQ);
                    double ar = g->absorb[a].a_re[y];
                    double ai = g->absorb[a].a_im[y];
                    double pr = pi_re[y], pim = pi_im[y];
                    double hr = Hqy * ar, hi = Hqy * ai;
                    mid_re[q] += hr * pr - hi * pim;
                    mid_im[q] += hr * pim + hi * pr;
                }
            }

            /* Σ_q H[xv][q] · a_cur(q) · po(q) · mid[q] */
            double sum_re = 0.0, sum_im = 0.0;
            for (int q = 0; q < 2; q++) {
                double Hxq = (xv == 0) ? SQ : (q == 0 ? SQ : -SQ);
                double cr = g->absorb[a].a_cur_re[q];
                double ci = g->absorb[a].a_cur_im[q];
                double por = po_re[q], poi = po_im[q];
                double mr = mid_re[q], mi = mid_im[q];
                double cpor = cr * por - ci * poi;
                double cpoi = cr * poi + ci * por;
                double comb_re = cpor * mr - cpoi * mi;
                double comb_im = cpor * mi + cpoi * mr;
                sum_re += Hxq * comb_re;
                sum_im += Hxq * comb_im;
            }
            double new_re = re * sum_re - im * sum_im;
            double new_im = re * sum_im + im * sum_re;
            re = new_re;
            im = new_im;
        } else {
            /* Single-layer: Σ_y H_{xv,y} · a_re(y) · pi(y) */
            double sum_re = 0.0, sum_im = 0.0;
            for (int y = 0; y < 2; y++) {
                double Hy = (xv == 0) ? SQ : (y == 0 ? SQ : -SQ);
                double ar = g->absorb[a].a_re[y];
                double ai = g->absorb[a].a_im[y];
                double pr = pi_re[y], pim = pi_im[y];
                double hr = Hy * ar, hi = Hy * ai;
                sum_re += hr * pr - hi * pim;
                sum_im += hr * pim + hi * pr;
            }
            double new_re = re * sum_re - im * sum_im;
            double new_im = re * sum_im + im * sum_re;
            re = new_re;
            im = new_im;
        }
    }

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
                uint32_t phase_idx = (va * vb) % HPCQ_D;
                w_re = HPCQ_W2_RE[phase_idx];
                w_im = HPCQ_W2_IM[phase_idx];
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
                    uint32_t phase_idx = (outcome * (uint32_t)k) % HPCQ_D;
                    w_re = HPCQ_W2_RE[phase_idx];
                    w_im = HPCQ_W2_IM[phase_idx];
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
        fprintf(stderr, "hpcq_norm_sq: N=%lu too large\n", g->n_sites);
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
