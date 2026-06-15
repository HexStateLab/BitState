/*
 * hpc_qubit_amplitude.h — On-Demand State Vector for D=2
 *
 * Sparse state vector reconstruction and Monte Carlo expectation values.
 * Never materializes the full 2^N state vector.
 *
 * Three modes:
 *   1. POINT QUERY:    ψ(i₁,...,iₙ) → O(N+E)
 *   2. SPARSE RECON:   All |ψ| > threshold → tree-pruned
 *   3. EXPECTATION:    ⟨ψ|O|ψ⟩ → importance sampling
 *
 * D=2 simplifications:
 *   - Configurations are bitstrings → use uint64_t for up to 64 qubits
 *   - Enumeration: 2^n_connected (vs 6^n for D=6)
 *   - Tree pruning: binary tree (2 children vs 6)
 */

#ifndef HPC_QUBIT_AMPLITUDE_H
#define HPC_QUBIT_AMPLITUDE_H

#include "hpc_qubit_graph.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * SPARSE STATE VECTOR ENTRY
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    uint32_t *indices;
    double    re, im;
    double    prob;
} HPCQSparseEntry;

typedef struct {
    HPCQSparseEntry *entries;
    uint64_t         count;
    uint64_t         capacity;
    uint64_t         n_sites;
    double           total_prob;
    double           threshold;
} HPCQSparseVector;

/* ═══════════════════════════════════════════════════════════════════════════════
 * SPARSE VECTOR LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline HPCQSparseVector *hpcq_sv_create(uint64_t n_sites,
                                                 uint64_t initial_cap)
{
    HPCQSparseVector *sv = (HPCQSparseVector *)calloc(1, sizeof(HPCQSparseVector));
    if (!sv) return NULL;
    sv->n_sites = n_sites;
    sv->capacity = initial_cap;
    sv->entries = (HPCQSparseEntry *)calloc(initial_cap, sizeof(HPCQSparseEntry));
    for (uint64_t i = 0; i < initial_cap; i++)
        sv->entries[i].indices = (uint32_t *)calloc(n_sites, sizeof(uint32_t));
    return sv;
}

static inline void hpcq_sv_destroy(HPCQSparseVector *sv)
{
    if (!sv) return;
    for (uint64_t i = 0; i < sv->capacity; i++)
        free(sv->entries[i].indices);
    free(sv->entries);
    free(sv);
}

static inline void hpcq_sv_grow(HPCQSparseVector *sv)
{
    if (sv->count < sv->capacity) return;
    uint64_t new_cap = sv->capacity * 2;
    sv->entries = (HPCQSparseEntry *)realloc(sv->entries,
                                               new_cap * sizeof(HPCQSparseEntry));
    for (uint64_t i = sv->capacity; i < new_cap; i++) {
        sv->entries[i].indices = (uint32_t *)calloc(sv->n_sites, sizeof(uint32_t));
        sv->entries[i].re = 0; sv->entries[i].im = 0; sv->entries[i].prob = 0;
    }
    sv->capacity = new_cap;
}

static inline void hpcq_sv_add(HPCQSparseVector *sv,
                                 const uint32_t *indices,
                                 double re, double im)
{
    hpcq_sv_grow(sv);
    HPCQSparseEntry *e = &sv->entries[sv->count];
    memcpy(e->indices, indices, sv->n_sites * sizeof(uint32_t));
    e->re = re; e->im = im;
    e->prob = re * re + im * im;
    sv->total_prob += e->prob;
    sv->count++;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BRUTE-FORCE SPARSE RECONSTRUCTION — O(2^N × (N+E))
 *
 * For small N: enumerate all 2^N configurations.
 * D=2 advantage: bitstring enumeration, up to N=20 practical.
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline HPCQSparseVector *hpcq_sparse_brute(const HPCQGraph *g,
                                                     double threshold,
                                                     uint64_t max_entries)
{
    if (g->n_sites > 20) {
        fprintf(stderr, "hpcq_sparse_brute: N=%lu too large\n", g->n_sites);
        return NULL;
    }

    HPCQSparseVector *sv = hpcq_sv_create(g->n_sites, 256);
    if (!sv) return NULL;
    sv->threshold = threshold;

    uint64_t total_configs = 1ULL << g->n_sites;
    uint32_t indices[20];

    for (uint64_t cfg = 0; cfg < total_configs && sv->count < max_entries; cfg++) {
        for (uint64_t i = 0; i < g->n_sites; i++)
            indices[i] = (cfg >> i) & 1;

        double re, im;
        hpcq_amplitude(g, indices, &re, &im);
        double prob = re * re + im * im;

        if (prob >= threshold)
            hpcq_sv_add(sv, indices, re, im);
    }

    return sv;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TREE-PRUNED SPARSE RECONSTRUCTION
 *
 * For larger N: binary tree pruning.
 * At each site, extend each live branch to value 0 and 1.
 * Prune branches whose cumulative prob falls below threshold.
 *
 * D=2 advantage: only 2 children per node (vs 6 for HexState).
 * Cost: O(active_branches × 2 × E_local) per site.
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    uint32_t *indices;
    double    re, im;
} HPCQTreeNode;

static inline HPCQSparseVector *hpcq_sparse_tree(const HPCQGraph *g,
                                                    double threshold,
                                                    uint64_t max_branches)
{
    HPCQSparseVector *sv = hpcq_sv_create(g->n_sites, 256);
    if (!sv) return NULL;
    sv->threshold = threshold;

    uint64_t pool_cap = max_branches * 2 + 16;
    HPCQTreeNode *current = (HPCQTreeNode *)calloc(pool_cap, sizeof(HPCQTreeNode));
    HPCQTreeNode *next    = (HPCQTreeNode *)calloc(pool_cap, sizeof(HPCQTreeNode));
    for (uint64_t i = 0; i < pool_cap; i++) {
        current[i].indices = (uint32_t *)calloc(g->n_sites, sizeof(uint32_t));
        next[i].indices    = (uint32_t *)calloc(g->n_sites, sizeof(uint32_t));
    }

    /* Seed: one root */
    uint64_t n_current = 1;
    current[0].re = 1.0;
    current[0].im = 0.0;

    /* Grow site by site */
    for (uint64_t site = 0; site < g->n_sites; site++) {
        uint64_t n_next = 0;

        double local_re[HPCQ_D], local_im[HPCQ_D];
        tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE,
                          local_re, local_im);

        for (uint64_t b = 0; b < n_current; b++) {
            for (int v = 0; v < HPCQ_D; v++) {
                double a_re = local_re[v];
                double a_im = local_im[v];

                double new_re = current[b].re * a_re - current[b].im * a_im;
                double new_im = current[b].re * a_im + current[b].im * a_re;

                /* Apply phase from edges connecting this site to prior sites */
                for (uint64_t e = 0; e < g->n_edges; e++) {
                    uint64_t sa = g->edges[e].site_a;
                    uint64_t sb = g->edges[e].site_b;
                    int partner_site = -1;

                    if (sa == site && sb < site) partner_site = (int)sb;
                    else if (sb == site && sa < site) partner_site = (int)sa;

                    if (partner_site >= 0) {
                        /* Skip edges consumed by multi-edge absorption */
                        if (g->edges[e].type == HPCQ_EDGE_ABSORBED) continue;

                        uint32_t pv = current[b].indices[partner_site];
                        double w_re, w_im;

                        if (g->edges[e].type == HPCQ_EDGE_CZ) {
                            uint32_t pi = ((uint32_t)v * pv) % HPCQ_D;
                            w_re = HPCQ_W2_RE[pi];
                            w_im = HPCQ_W2_IM[pi];
                        } else {
                            if (sa == site) {
                                w_re = g->edges[e].w_re[v][pv];
                                w_im = g->edges[e].w_im[v][pv];
                            } else {
                                w_re = g->edges[e].w_re[pv][v];
                                w_im = g->edges[e].w_im[pv][v];
                            }
                        }

                        double tmp_re = new_re * w_re - new_im * w_im;
                        double tmp_im = new_re * w_im + new_im * w_re;
                        new_re = tmp_re;
                        new_im = tmp_im;
                    }
                }

                double prob = new_re * new_re + new_im * new_im;
                if (prob < threshold && site < g->n_sites - 1) continue;

                if (n_next < pool_cap) {
                    memcpy(next[n_next].indices, current[b].indices,
                           g->n_sites * sizeof(uint32_t));
                    next[n_next].indices[site] = (uint32_t)v;
                    next[n_next].re = new_re;
                    next[n_next].im = new_im;
                    n_next++;
                }
            }
        }

        /* Swap pools */
        HPCQTreeNode *tmp = current;
        current = next;
        next = tmp;
        n_current = n_next;

        /* Truncate to max_branches */
        if (n_current > max_branches && site < g->n_sites - 1) {
            /* Keep top by probability */
            for (uint64_t i = max_branches; i < n_current; i++) {
                uint64_t min_idx = 0;
                double min_p = current[0].re * current[0].re +
                               current[0].im * current[0].im;
                for (uint64_t j = 1; j < max_branches; j++) {
                    double p = current[j].re * current[j].re +
                               current[j].im * current[j].im;
                    if (p < min_p) { min_p = p; min_idx = j; }
                }
                double p_i = current[i].re * current[i].re +
                             current[i].im * current[i].im;
                if (p_i > min_p) {
                    HPCQTreeNode swap = current[min_idx];
                    current[min_idx] = current[i];
                    current[i] = swap;
                }
            }
            n_current = max_branches;
        }
    }

    /* Collect results — apply absorb correction, then filter by threshold */
    static const double SQRT2 = 0.7071067811865475244;
    for (uint64_t b = 0; b < n_current; b++) {
        double re = current[b].re, im = current[b].im;
        for (uint64_t a = 0; a < g->n_absorb; a++) {
            uint64_t center = g->absorb[a].center;
            int xv = (int)current[b].indices[center];
            uint64_t n_inner = g->absorb[a].n_inner;

            /* Inner product: Π w_k(y, z_k) for y = 0, 1 (k < n_inner) */
            double pi_re[2] = {1.0, 1.0};
            double pi_im[2] = {0.0, 0.0};
            for (uint64_t k = 0; k < n_inner; k++) {
                int z = (int)current[b].indices[g->absorb[a].nbrs[k]];
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
                    int z = (int)current[b].indices[g->absorb[a].nbrs[k]];
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
                        double Hqy = (q == 0) ? SQRT2 : (y == 0 ? SQRT2 : -SQRT2);
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
                    double Hxq = (xv == 0) ? SQRT2 : (q == 0 ? SQRT2 : -SQRT2);
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
                double tmp_re = re * sum_re - im * sum_im;
                double tmp_im = re * sum_im + im * sum_re;
                re = tmp_re; im = tmp_im;
            } else {
                /* Single-layer: Σ_y H[xv][y] · a_re(y) · pi(y) */
                double sum_re = 0.0, sum_im = 0.0;
                for (int y = 0; y < 2; y++) {
                    double Hy = (xv == 0) ? SQRT2 : (y == 0 ? SQRT2 : -SQRT2);
                    double ar = g->absorb[a].a_re[y];
                    double ai = g->absorb[a].a_im[y];
                    double pr = pi_re[y], pim = pi_im[y];
                    double hr = Hy * ar, hi = Hy * ai;
                    sum_re += hr * pr - hi * pim;
                    sum_im += hr * pim + hi * pr;
                }
                double tmp_re = re * sum_re - im * sum_im;
                double tmp_im = re * sum_im + im * sum_re;
                re = tmp_re; im = tmp_im;
            }
        }
        double prob = re * re + im * im;
        if (prob >= threshold)
            hpcq_sv_add(sv, current[b].indices, re, im);
    }

    for (uint64_t i = 0; i < pool_cap; i++) {
        free(current[i].indices);
        free(next[i].indices);
    }
    free(current);
    free(next);

    return sv;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MONTE CARLO EXPECTATION VALUE
 *
 * ⟨ψ|O|ψ⟩ via importance sampling.
 * Cost: O(n_samples × (N + E))
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef double (*HPCQObservable)(const uint32_t *indices, uint64_t n_sites,
                                  void *ctx);

static inline double hpcq_expectation(const HPCQGraph *g,
                                        HPCQObservable obs, void *obs_ctx,
                                        int n_samples, uint64_t rng_seed)
{
    uint64_t rng = rng_seed;
    #define HPCQ_LCG(r) ((r) = (r) * 6364136223846793005ULL + 1442695040888963407ULL)
    #define HPCQ_RAND(r) (((double)((r) >> 11)) * 0x1.0p-53)

    double sum_obs = 0.0;
    int valid = 0;

    for (int s = 0; s < n_samples; s++) {
        uint32_t config[64];
        if (g->n_sites > 64) break;

        /* Sample each site from local distribution */
        for (uint64_t site = 0; site < g->n_sites; site++) {
            double re[HPCQ_D], im[HPCQ_D];
            tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE,
                              re, im);
            double p0 = re[0]*re[0] + im[0]*im[0];
            double p1 = re[1]*re[1] + im[1]*im[1];
            double total = p0 + p1;

            HPCQ_LCG(rng);
            double r = HPCQ_RAND(rng) * total;
            config[site] = (r <= p0) ? 0 : 1;
        }

        /* Importance weight */
        double prob_psi = hpcq_probability(g, config);
        double prob_q = 1.0;
        for (uint64_t site = 0; site < g->n_sites; site++) {
            double re[HPCQ_D], im[HPCQ_D];
            tri_get_amplitudes((TrialityQubit *)&g->locals[site], VIEW_EDGE,
                              re, im);
            uint32_t v = config[site];
            prob_q *= re[v]*re[v] + im[v]*im[v];
        }

        if (prob_q > 1e-30) {
            double weight = prob_psi / prob_q;
            double obs_val = obs(config, g->n_sites, obs_ctx);
            sum_obs += weight * obs_val;
            valid++;
        }
    }

    #undef HPCQ_LCG
    #undef HPCQ_RAND

    return (valid > 0) ? sum_obs / valid : 0.0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRINT SPARSE VECTOR
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline void hpcq_sv_print(const HPCQSparseVector *sv, int max_show)
{
    printf("── Sparse State Vector (D=2) ──\n");
    printf("  Entries: %lu, Captured prob: %.6f, Threshold: %.2e\n",
           sv->count, sv->total_prob, sv->threshold);

    uint64_t show = sv->count;
    if (max_show > 0 && show > (uint64_t)max_show) show = max_show;

    for (uint64_t i = 0; i < show; i++) {
        printf("  |");
        for (uint64_t s = 0; s < sv->n_sites; s++)
            printf("%u", sv->entries[i].indices[s]);
        printf("⟩ → %.6f%+.6fi  (P=%.6e)\n",
               sv->entries[i].re, sv->entries[i].im, sv->entries[i].prob);
    }
    if (show < sv->count)
        printf("  ... (%lu more entries)\n", sv->count - show);
}

#endif /* HPC_QUBIT_AMPLITUDE_H */
