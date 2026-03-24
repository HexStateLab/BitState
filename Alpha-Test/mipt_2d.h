/*
 * mipt_2d.h — Measurement-Induced Phase Transition on a 2D Square Lattice
 *
 * Uses the HPC phase graph to efficiently track entanglement structure
 * through random circuits with competing unitaries and measurements.
 *
 * The HPC graph IS the entanglement map:
 *   - CZ edges = entangling bonds (fidelity 1.0)
 *   - hpcq_measure() = projective measurement (removes edges)
 *   - Edge count across a cut = entanglement entropy estimate
 *
 * Circuit model:
 *   Each timestep applies two layers:
 *   1. UNITARY LAYER: checkerboard CZ gates preceded by random H+Phase
 *   2. MEASUREMENT LAYER: each qubit measured with probability p
 *
 * Phase transition at p_c ≈ 0.17 for 2D:
 *   p < p_c: volume-law entanglement (S ∝ L²)
 *   p > p_c: area-law entanglement (S ∝ L)
 */

#ifndef MIPT_2D_H
#define MIPT_2D_H

#include "hpc_qubit_graph.h"
#include "hpc_qubit_contract.h"
#include <stdint.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * LATTICE GEOMETRY
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    uint32_t Lx, Ly;         /* Lattice dimensions                        */
    uint32_t N;              /* Total qubits = Lx × Ly                    */
    HPCQGraph *graph;        /* HPC phase graph                           */
    uint64_t  rng;           /* PRNG state                                */

    /* ── Circuit parameters ── */
    uint32_t  depth;         /* Current circuit depth                     */
    double    meas_rate;     /* Measurement probability p                 */

    /* ── Diagnostics ── */
    uint64_t  total_gates;
    uint64_t  total_measurements;
    uint64_t  total_collapses;  /* Measurements that actually collapsed   */
} MIPT2D;

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENTANGLEMENT DIAGNOSTICS
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double   half_entropy;       /* S(L/2): entropy across middle cut     */
    double   quarter_entropy;    /* S(L/4): entropy across quarter cut    */
    double   edge_density;       /* Edges / max_possible_edges            */
    uint64_t n_edges;            /* Total edges in graph                  */
    uint64_t cz_edges;           /* CZ edges                              */
    double   avg_fidelity;       /* Average edge fidelity                 */
    double   entropy_per_row[64];/* Per-row contribution (max 64 rows)    */
} MIPTDiagnostics;

/* ═══════════════════════════════════════════════════════════════════════════════
 * SWEEP RESULT — One data point in the p-scan
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double p;                    /* Measurement rate                      */
    double S_mean;               /* Mean half-cut entropy                 */
    double S_std;                /* Standard deviation                    */
    double S_over_L;             /* S / L — scaling indicator             */
    double edge_density;         /* Mean edge density                     */
    int    n_samples;            /* Number of circuit realizations        */
} MIPTSweepPoint;

/* ═══════════════════════════════════════════════════════════════════════════════
 * API
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Lifecycle */
MIPT2D *mipt_create(uint32_t Lx, uint32_t Ly, uint64_t seed);
void    mipt_destroy(MIPT2D *m);
void    mipt_reset(MIPT2D *m);  /* Reset to |0...0⟩ product state */

/* Circuit operations */
void mipt_unitary_layer(MIPT2D *m, int parity);  /* 0=even, 1=odd */
void mipt_measure_layer(MIPT2D *m, double p);
void mipt_step(MIPT2D *m, double p);             /* One full timestep */
void mipt_run(MIPT2D *m, uint32_t depth, double p);  /* Run to depth */

/* Diagnostics */
MIPTDiagnostics mipt_diagnostics(const MIPT2D *m);

/* Phase sweep */
MIPTSweepPoint *mipt_sweep(uint32_t Lx, uint32_t Ly,
                            double p_min, double p_max, int n_points,
                            uint32_t depth, int samples_per_p,
                            uint64_t seed);
void mipt_sweep_free(MIPTSweepPoint *points);

/* Output */
void mipt_print_diagnostics(const MIPTDiagnostics *d, uint32_t Lx, uint32_t Ly);
void mipt_print_sweep(const MIPTSweepPoint *points, int n_points,
                       uint32_t Lx, uint32_t Ly);
void mipt_save_sweep_tsv(const MIPTSweepPoint *points, int n_points,
                           uint32_t Lx, uint32_t Ly, const char *filename);

/* Coordinate helpers */
static inline uint32_t mipt_site(const MIPT2D *m, uint32_t x, uint32_t y)
{
    return y * m->Lx + x;
}

static inline void mipt_coords(const MIPT2D *m, uint32_t site,
                                 uint32_t *x, uint32_t *y)
{
    *x = site % m->Lx;
    *y = site / m->Lx;
}

#endif /* MIPT_2D_H */
