/*
 * mipt_2d.c — Measurement-Induced Phase Transition Implementation
 *
 * Random circuit model on a 2D square lattice using HPC phase graph.
 *
 * Circuit per timestep:
 *   1. Random single-qubit rotations (H + random phase) on all qubits
 *   2. Checkerboard CZ gates (even/odd parity alternating)
 *   3. Projective measurements with probability p
 *
 * The HPC graph tracks entanglement without materializing 2^N state vector.
 *
 * Build:
 *   gcc -O2 -march=native -o mipt_2d mipt_2d.c qubit_triality.c \
 *       clifford_exotic.c -lm -DMIPT_STANDALONE
 */

#include "mipt_2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRNG — LCG with SplitMix64 output function
 * ═══════════════════════════════════════════════════════════════════════════════ */

static inline uint64_t mipt_rng(uint64_t *state)
{
    *state = (*state) * 6364136223846793005ULL + 1442695040888963407ULL;
    uint64_t z = *state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline double mipt_rand01(uint64_t *state)
{
    return ((double)(mipt_rng(state) >> 11)) * 0x1.0p-53;
}

static inline double mipt_rand_angle(uint64_t *state)
{
    return mipt_rand01(state) * 2.0 * M_PI;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

MIPT2D *mipt_create(uint32_t Lx, uint32_t Ly, uint64_t seed)
{
    MIPT2D *m = (MIPT2D *)calloc(1, sizeof(MIPT2D));
    if (!m) return NULL;

    m->Lx = Lx;
    m->Ly = Ly;
    m->N = Lx * Ly;
    m->rng = seed;
    m->graph = hpcq_create(m->N);

    if (!m->graph) { free(m); return NULL; }

    /* Initialize all qubits to |0⟩ (default from hpcq_create/tri_init) */
    return m;
}

void mipt_destroy(MIPT2D *m)
{
    if (!m) return;
    hpcq_destroy(m->graph);
    free(m);
}

void mipt_reset(MIPT2D *m)
{
    if (!m) return;
    hpcq_destroy(m->graph);
    m->graph = hpcq_create(m->N);
    m->depth = 0;
    m->total_gates = 0;
    m->total_measurements = 0;
    m->total_collapses = 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * UNITARY LAYER — Random CZ gates with random local rotations
 *
 * Each step:
 *   1. Apply random single-qubit rotation (H + random phase) to all qubits
 *   2. For each nearest-neighbor bond in the 2D lattice, apply CZ with
 *      probability 0.5. This breaks the periodic CZ²=I cancellation that
 *      would occur with a fixed checkerboard pattern at even depths.
 *
 * With CZ involution (CZ²=I), an edge toggles on/off each time the same
 * bond is activated. Random per-bond selection ensures edge count saturates
 * at a steady-state value rather than periodically collapsing to zero.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void mipt_unitary_layer(MIPT2D *m, int parity)
{
    (void)parity;  /* No longer used — bonds are randomly selected */

    /* Apply random single-qubit rotations to all sites */
    for (uint32_t i = 0; i < m->N; i++) {
        /* Random Hadamard + phase: creates random basis rotation */
        hpcq_hadamard(m->graph, i);
        double theta = mipt_rand_angle(&m->rng);
        hpcq_phase(m->graph, i, theta);
        m->total_gates++;
    }

    /* Horizontal bonds: (x, y) — (x+1, y), each with prob 0.5 */
    for (uint32_t y = 0; y < m->Ly; y++) {
        for (uint32_t x = 0; x + 1 < m->Lx; x++) {
            if (mipt_rand01(&m->rng) < 0.5) {
                uint32_t a = mipt_site(m, x, y);
                uint32_t b = mipt_site(m, x + 1, y);
                hpcq_cz(m->graph, a, b);
                m->total_gates++;
            }
        }
    }

    /* Vertical bonds: (x, y) — (x, y+1), each with prob 0.5 */
    for (uint32_t x = 0; x < m->Lx; x++) {
        for (uint32_t y = 0; y + 1 < m->Ly; y++) {
            if (mipt_rand01(&m->rng) < 0.5) {
                uint32_t a = mipt_site(m, x, y);
                uint32_t b = mipt_site(m, x, y + 1);
                hpcq_cz(m->graph, a, b);
                m->total_gates++;
            }
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT LAYER — Each qubit measured with probability p
 *
 * The HPC graph makes this cheap:
 *   - Measurement collapses the qubit to |0⟩ or |1⟩
 *   - All edges touching that qubit are absorbed into partners
 *   - Edge removal IS disentanglement
 * ═══════════════════════════════════════════════════════════════════════════════ */

void mipt_measure_layer(MIPT2D *m, double p)
{
    for (uint32_t i = 0; i < m->N; i++) {
        if (mipt_rand01(&m->rng) < p) {
            /* Count edges touching this site before measurement */
            uint64_t edges_before = m->graph->n_edges;

            double r = mipt_rand01(&m->rng);
            hpcq_measure(m->graph, i, r);

            m->total_measurements++;
            if (m->graph->n_edges < edges_before)
                m->total_collapses++;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * STEP — One complete circuit timestep
 * ═══════════════════════════════════════════════════════════════════════════════ */

void mipt_step(MIPT2D *m, double p)
{
    int parity = m->depth % 2;

    /* Unitary layer (entangling) */
    mipt_unitary_layer(m, parity);

    /* Measurement layer (disentangling) */
    mipt_measure_layer(m, p);

    m->depth++;
    m->meas_rate = p;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * RUN — Execute circuit to target depth
 * ═══════════════════════════════════════════════════════════════════════════════ */

void mipt_run(MIPT2D *m, uint32_t depth, double p)
{
    for (uint32_t t = 0; t < depth; t++)
        mipt_step(m, p);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DIAGNOSTICS — Entanglement structure analysis
 *
 * Uses HPC edge crossings to estimate entanglement entropy.
 * For CZ edges: each crossing contributes exactly 1 bit.
 * For general edges: contribution weighted by fidelity.
 * ═══════════════════════════════════════════════════════════════════════════════ */

MIPTDiagnostics mipt_diagnostics(const MIPT2D *m)
{
    MIPTDiagnostics d;
    memset(&d, 0, sizeof(d));

    /* Half-system entropy: cut between rows Ly/2-1 and Ly/2 */
    uint32_t half_cut = (m->Ly / 2 - 1) * m->Lx + m->Lx - 1;
    d.half_entropy = hpcq_entropy_cut(m->graph, half_cut);

    /* Quarter-system entropy: cut at Ly/4 */
    uint32_t quarter_cut = (m->Ly / 4) * m->Lx + m->Lx - 1;
    d.quarter_entropy = hpcq_entropy_cut(m->graph, quarter_cut);

    /* Edge statistics */
    d.n_edges = m->graph->n_edges;
    d.cz_edges = m->graph->cz_edges;
    d.avg_fidelity = m->graph->avg_fidelity;

    /* Max possible edges: each site can have up to 4 neighbors in 2D */
    uint64_t max_edges = 2ULL * m->Lx * m->Ly;
    d.edge_density = (max_edges > 0) ? (double)d.n_edges / max_edges : 0.0;

    /* Per-row entropy contribution */
    uint32_t max_rows = (m->Ly < 64) ? m->Ly : 64;
    for (uint32_t row = 0; row + 1 < max_rows; row++) {
        uint32_t cut = (row + 1) * m->Lx - 1;

        /* Count edges crossing between row and row+1 */
        double row_entropy = 0.0;
        for (uint64_t e = 0; e < m->graph->n_edges; e++) {
            uint64_t sa = m->graph->edges[e].site_a;
            uint64_t sb = m->graph->edges[e].site_b;
            uint32_t ya = (uint32_t)(sa / m->Lx);
            uint32_t yb = (uint32_t)(sb / m->Lx);

            if ((ya <= row && yb > row) || (yb <= row && ya > row)) {
                row_entropy += m->graph->edges[e].fidelity;
            }
        }
        d.entropy_per_row[row] = row_entropy;
    }

    return d;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PHASE SWEEP — Scan measurement rate to find the transition
 *
 * For each p value:
 *   1. Create fresh lattice
 *   2. Run circuit to depth (≈ 2L for thermalization)
 *   3. Record half-cut entropy
 *   4. Average over multiple realizations
 *
 * Output: array of MIPTSweepPoint with S(p) curve data.
 * ═══════════════════════════════════════════════════════════════════════════════ */

MIPTSweepPoint *mipt_sweep(uint32_t Lx, uint32_t Ly,
                            double p_min, double p_max, int n_points,
                            uint32_t depth, int samples_per_p,
                            uint64_t seed)
{
    MIPTSweepPoint *points = (MIPTSweepPoint *)calloc(n_points,
                                                        sizeof(MIPTSweepPoint));
    if (!points) return NULL;

    double dp = (n_points > 1) ? (p_max - p_min) / (n_points - 1) : 0;

    for (int pi = 0; pi < n_points; pi++) {
        double p = p_min + pi * dp;
        points[pi].p = p;
        points[pi].n_samples = samples_per_p;

        double sum_S = 0, sum_S2 = 0, sum_density = 0;

        for (int s = 0; s < samples_per_p; s++) {
            MIPT2D *m = mipt_create(Lx, Ly, seed + (uint64_t)pi * 10000 + s);
            if (!m) continue;

            mipt_run(m, depth, p);

            MIPTDiagnostics diag = mipt_diagnostics(m);
            sum_S += diag.half_entropy;
            sum_S2 += diag.half_entropy * diag.half_entropy;
            sum_density += diag.edge_density;

            mipt_destroy(m);
        }

        points[pi].S_mean = sum_S / samples_per_p;
        points[pi].S_std = sqrt(sum_S2 / samples_per_p -
                                points[pi].S_mean * points[pi].S_mean);
        points[pi].S_over_L = points[pi].S_mean / Lx;
        points[pi].edge_density = sum_density / samples_per_p;
    }

    return points;
}

void mipt_sweep_free(MIPTSweepPoint *points)
{
    free(points);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * OUTPUT — Diagnostics and sweep results
 * ═══════════════════════════════════════════════════════════════════════════════ */

void mipt_print_diagnostics(const MIPTDiagnostics *d, uint32_t Lx, uint32_t Ly)
{
    printf("╔═════════════════════════════════════════════════════╗\n");
    printf("║  MIPT 2D Diagnostics (%ux%u lattice)               ║\n", Lx, Ly);
    printf("╠═════════════════════════════════════════════════════╣\n");
    printf("║  Half-cut entropy:    %10.4f bits              ║\n", d->half_entropy);
    printf("║  Quarter-cut entropy: %10.4f bits              ║\n", d->quarter_entropy);
    printf("║  Total edges:         %10lu                    ║\n", d->n_edges);
    printf("║    CZ edges:          %10lu                    ║\n", d->cz_edges);
    printf("║  Edge density:        %10.4f                    ║\n", d->edge_density);
    printf("║  Avg fidelity:        %10.6f                    ║\n", d->avg_fidelity);
    printf("╠═════════════════════════════════════════════════════╣\n");
    printf("║  S/L = %.4f (%s regime)                       ║\n",
           d->half_entropy / Lx,
           (d->half_entropy / Lx > 0.3) ? "volume-law" :
           (d->half_entropy / Lx > 0.1) ? "critical  " : "area-law  ");
    printf("╚═════════════════════════════════════════════════════╝\n");
}

void mipt_print_sweep(const MIPTSweepPoint *points, int n_points,
                       uint32_t Lx, uint32_t Ly)
{
    printf("\n╔═══════════════════════════════════════════════════════════╗\n");
    printf("║  MIPT Phase Sweep — %ux%u lattice                        ║\n",
           Lx, Ly);
    printf("╠═════════╤══════════╤══════════╤══════════╤═════════════════╣\n");
    printf("║    p    │  S(L/2)  │   S/L    │  σ(S)   │   edge dens    ║\n");
    printf("╠═════════╪══════════╪══════════╪══════════╪═════════════════╣\n");

    for (int i = 0; i < n_points; i++) {
        const MIPTSweepPoint *pt = &points[i];

        /* Visual bar for S/L */
        int bar_len = (int)(pt->S_over_L * 40);
        if (bar_len > 20) bar_len = 20;
        char bar[21];
        memset(bar, ' ', 20);
        bar[20] = '\0';
        for (int b = 0; b < bar_len; b++) bar[b] = '#';

        printf("║  %.3f  │  %6.2f  │  %6.4f  │  %5.3f  │ %s║\n",
               pt->p, pt->S_mean, pt->S_over_L, pt->S_std, bar);
    }

    printf("╚═════════╧══════════╧══════════╧══════════╧═════════════════╝\n");

    /* Estimate p_c from max dS/dp (steepest drop in entropy) */
    double max_drop = 0;
    int max_idx = 0;
    for (int i = 1; i < n_points; i++) {
        double drop = (points[i-1].S_mean - points[i].S_mean) /
                      (points[i].p - points[i-1].p);
        if (drop > max_drop) {
            max_drop = drop;
            max_idx = i;
        }
    }
    if (max_idx > 0) {
        printf("\n  Estimated p_c ≈ %.3f (max |dS/dp| between p=%.3f and p=%.3f)\n",
               (points[max_idx-1].p + points[max_idx].p) / 2,
               points[max_idx-1].p, points[max_idx].p);
    }
}

void mipt_save_sweep_tsv(const MIPTSweepPoint *points, int n_points,
                           uint32_t Lx, uint32_t Ly, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", filename); return; }

    fprintf(f, "# MIPT 2D Phase Sweep: %ux%u lattice\n", Lx, Ly);
    fprintf(f, "# p\tS_mean\tS_std\tS_over_L\tedge_density\tn_samples\n");

    for (int i = 0; i < n_points; i++) {
        fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
                points[i].p, points[i].S_mean, points[i].S_std,
                points[i].S_over_L, points[i].edge_density,
                points[i].n_samples);
    }

    fclose(f);
    printf("  Saved sweep data to %s\n", filename);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * STANDALONE DRIVER — When compiled with -DMIPT_STANDALONE
 *
 * Usage: ./mipt_2d [Lx] [Ly] [depth] [samples]
 *
 * Default: 8×8 lattice, depth=16, 10 samples per p-point, 20 p-values
 * ═══════════════════════════════════════════════════════════════════════════════ */

#ifdef MIPT_STANDALONE
int main(int argc, char **argv)
{
    uint32_t Lx = 8, Ly = 8;
    uint32_t depth = 16;
    int samples = 10;

    if (argc > 1) Lx = (uint32_t)atoi(argv[1]);
    if (argc > 2) Ly = (uint32_t)atoi(argv[2]);
    if (argc > 3) depth = (uint32_t)atoi(argv[3]);
    if (argc > 4) samples = atoi(argv[4]);

    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║  MIPT 2D — Measurement-Induced Phase Transition       ║\n");
    printf("║  Lattice: %ux%u, Depth: %u, Samples: %d              ║\n",
           Lx, Ly, depth, samples);
    printf("╚════════════════════════════════════════════════════════╝\n\n");

    /* Phase sweep: p ∈ [0.0, 0.5] with 21 points */
    int n_points = 21;
    printf("Running phase sweep with %d p-values...\n\n", n_points);

    MIPTSweepPoint *sweep = mipt_sweep(Lx, Ly, 0.0, 0.5, n_points,
                                        depth, samples, 42);

    if (sweep) {
        mipt_print_sweep(sweep, n_points, Lx, Ly);

        /* Save for plotting */
        char fname[256];
        snprintf(fname, sizeof(fname), "mipt_%ux%u_d%u.tsv", Lx, Ly, depth);
        mipt_save_sweep_tsv(sweep, n_points, Lx, Ly, fname);

        mipt_sweep_free(sweep);
    }

    /* Single run at p=0 and p=0.5 for comparison */
    printf("\n━━━ Detailed diagnostics ━━━\n\n");

    printf("── p = 0.0 (no measurement, volume-law expected) ──\n");
    MIPT2D *m0 = mipt_create(Lx, Ly, 123);
    mipt_run(m0, depth, 0.0);
    MIPTDiagnostics d0 = mipt_diagnostics(m0);
    mipt_print_diagnostics(&d0, Lx, Ly);
    printf("  Gates: %lu, Measurements: %lu\n\n",
           m0->total_gates, m0->total_measurements);
    mipt_destroy(m0);

    printf("── p = 0.5 (heavy measurement, area-law expected) ──\n");
    MIPT2D *m5 = mipt_create(Lx, Ly, 456);
    mipt_run(m5, depth, 0.5);
    MIPTDiagnostics d5 = mipt_diagnostics(m5);
    mipt_print_diagnostics(&d5, Lx, Ly);
    printf("  Gates: %lu, Measurements: %lu\n",
           m5->total_gates, m5->total_measurements);
    mipt_destroy(m5);

    return 0;
}
#endif /* MIPT_STANDALONE */
