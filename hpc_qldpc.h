#ifndef HPC_QLDPC_H
#define HPC_QLDPC_H
#include "hpc_qubit_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

static uint64_t qlcg(uint64_t *s){return *s=*s*6364136223846793005ULL+1442695040888963407ULL;}

/* ── Generate random (w_c, w_r)-regular LDPC matrix H[r×n]
 * w_c = column weight, w_r = row weight.  Sparse: only w_r entries per row. ── */
static uint8_t *qldpc_gen_H(int r, int n, int w_c, int w_r, uint64_t *seed) {
    uint8_t *H = (uint8_t *)calloc((size_t)r * n, sizeof(uint8_t));
    if (!H) return NULL;
    /* For each row, randomly select w_r columns */
    for (int i = 0; i < r; i++) {
        int placed = 0;
        while (placed < w_r) {
            int j = (int)((qlcg(seed) >> 11) % (uint64_t)n);
            if (!H[i * n + j]) { H[i * n + j] = 1; placed++; }
        }
    }
    return H;
}

/* ── Hypergraph product code: H₁[r₁×n₁] ⊗ H₂[r₂×n₂]
 * Physical data qubits: n₁·n₂ (indexed by (i,j) Cartesian product)
 * Z-stabilizers: for each check in H₁ (row a) and each column in H₂ (col j):
 *    connects to data qubits (a',j) where H₁[a][a'] = 1
 * X-stabilizers: for each row in H₁^T and each check in H₂
 *
 * Total qubits in simulation:
 *   data qubits = n₁·n₂
 *   Z-ancilla = r₁·n₂  (one per Z-stabilizer)
 *   X-ancilla = n₁·r₂  (one per X-stabilizer)
 *   Total = n₁·n₂ + r₁·n₂ + n₁·r₂
 *
 * Code rate: (n₁-r₁)(n₂-r₂) / (n₁·n₂)  ── */

typedef struct {
    uint8_t *H1, *H2;
    int r1, n1, r2, n2;
    int w_c, w_r;
    uint64_t *qubits;      /* flat array of qubit indices */
    int n_data, n_z_anc, n_x_anc, n_total;
} QLDPCCode;

static QLDPCCode *qldpc_create(int n1, int n2, int w_c, int w_r, uint64_t *seed) {
    QLDPCCode *c = (QLDPCCode *)calloc(1, sizeof(QLDPCCode));
    c->r1 = n1 / 2; /* H1 is (n1/2) × n1, rate 1/2 */
    c->r2 = n2 / 2; /* H2 is (n2/2) × n2, rate 1/2 */
    c->n1 = n1; c->n2 = n2; c->w_c = w_c; c->w_r = w_r;

    c->H1 = qldpc_gen_H(c->r1, c->n1, w_c, w_r, seed);
    c->H2 = qldpc_gen_H(c->r2, c->n2, w_c, w_r, seed);

    c->n_data  = c->n1 * c->n2;
    c->n_z_anc = c->r1 * c->n2;
    c->n_x_anc = c->n1 * c->r2;
    c->n_total = c->n_data + c->n_z_anc + c->n_x_anc;

    c->qubits = (uint64_t *)calloc(c->n_total, sizeof(uint64_t));
    for (int i = 0; i < c->n_total; i++) c->qubits[i] = (uint64_t)i;
    return c;
}

static void qldpc_destroy(QLDPCCode *c) {
    if (!c) return;
    free(c->H1); free(c->H2); free(c->qubits); free(c);
}

/* ── Build qLDPC syndrome extraction circuit on HPC graph ── */
static void qldpc_build_circuit(HPCQGraph *g, QLDPCCode *c) {
    uint64_t *q = c->qubits;

    /* Data qubits: placed at indices 0 .. n_data-1 */
    int data_off = 0;
    int z_off    = c->n_data;
    int x_off    = c->n_data + c->n_z_anc;

    /* Z-stabilizers: for each row a of H1 and each column j of H2 */
    for (int a = 0; a < c->r1; a++) {
        for (int j = 0; j < c->n2; j++) {
            int anc_idx = z_off + a * c->n2 + j;
            /* Connect ancilla to data qubits (a', j) where H1[a][a']=1 */
            for (int ap = 0; ap < c->n1; ap++) {
                if (c->H1[a * c->n1 + ap]) {
                    int data_idx = data_off + ap * c->n2 + j;
                    hpcq_cz(g, q[anc_idx], q[data_idx]);
                }
            }
        }
    }

    /* X-stabilizers: for each row i of H1 and each row b of H2 */
    for (int i = 0; i < c->n1; i++) {
        for (int b = 0; b < c->r2; b++) {
            int anc_idx = x_off + i * c->r2 + b;
            /* Connect ancilla to data qubits (i, b') where H2[b][b']=1 */
            for (int bp = 0; bp < c->n2; bp++) {
                if (c->H2[b * c->n2 + bp]) {
                    int data_idx = data_off + i * c->n2 + bp;
                    hpcq_cz(g, q[anc_idx], q[data_idx]);
                }
            }
        }
    }
}

/* ── Print code statistics ── */
static void qldpc_print(QLDPCCode *c) {
    int logical = (c->n1 - c->r1) * (c->n2 - c->r2);
    double rate = (double)logical / (double)c->n_data;
    printf("  Physical data qubits: %d (%d×%d)\n", c->n_data, c->n1, c->n2);
    printf("  Logical qubits:       %d (rate=%.4f)\n", logical, rate);
    printf("  Z-stabilizers:        %d\n", c->n_z_anc);
    printf("  X-stabilizers:        %d\n", c->n_x_anc);
    printf("  Total simulation q:   %d\n", c->n_total);
    printf("  (w_c=%d, w_r=%d)\n", c->w_c, c->w_r);
}

#endif
