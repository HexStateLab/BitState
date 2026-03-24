/*
 * qubit_engine.h — QubitEngine Master Header
 *
 * D=2 qubit quantum engine, adapted from HexState D=6 architecture.
 * Uses square-geometry triality (Z/X/Y ↔ Edge/Vertex/Diagonal).
 *
 * Memory model (from headers):
 *   N qubits    = N × 32 bytes
 *   P pairs     = P × 64 bytes
 *   Total       = O(N + P), never O(2^N)
 */

#ifndef QUBIT_ENGINE_H
#define QUBIT_ENGINE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════════
 * CONSTANTS
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define QUBIT_D          2          /* Dimension                                */
#define QUBIT_D2         4          /* D² for joint states                      */
#define MAX_QUBITS       262144     /* 256K qubits                              */
#define MAX_PAIRS        262144     /* 256K entangled pairs                     */
#define MAX_REGISTERS    16384      /* 16K registers                            */

/* Magic pointer: encodes (chunk_id, offset) for sparse addressing */
#define MAGIC_PTR(chunk, offset)  (((uint64_t)(chunk) << 32) | (uint64_t)(offset))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ═══════════════════════════════════════════════════════════════════════════════
 * BASIS TYPE — uint64_t for qubit registers (up to 64 qubits natively)
 * For larger registers, use BigInt from bigint.h
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef uint64_t basis_t;

/* ═══════════════════════════════════════════════════════════════════════════════
 * INCLUDE HEADER-ONLY PRIMITIVES
 * ═══════════════════════════════════════════════════════════════════════════════ */

#include "arithmetic.h"
#include "statevector.h"
#include "superposition.h"
#include "born_rule.h"
#include "entanglement.h"
#include "qubit_management.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT — Single qubit with metadata
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    QubitState state;           /* 2 complex amplitudes (32 bytes)            */
    int32_t    pair_id;         /* Index into pairs[] (-1 if not entangled)  */
    int32_t    pair_side;       /* 0 = side A, 1 = side B                    */
    uint32_t   id;              /* This qubit's ID within the engine         */
} Qubit;

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT PAIR — Two entangled qubits
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    JointState joint;           /* 4 complex amplitudes (64 bytes)           */
    uint32_t   id_a;            /* First qubit ID                            */
    uint32_t   id_b;            /* Second qubit ID                           */
    uint8_t    active;          /* 1 = pair is live                          */
} QubitPair;

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT REGISTER — Sparse state vector for multi-qubit states
 *
 * Stores up to 4096 nonzero amplitudes in the 2^N Hilbert space.
 * Basis states are uint64_t (supports up to 64 qubits natively).
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    basis_t  basis_state;       /* Basis state index                         */
    double   amp_re;            /* Real amplitude                            */
    double   amp_im;            /* Imaginary amplitude                       */
} RegisterEntry;

typedef struct {
    uint64_t       chunk_id;                /* Chunk identifier              */
    uint64_t       n_qubits;                /* Number of logical qubits      */
    uint32_t       dim;                     /* Dimension per qubit (always 2)*/
    uint8_t        collapsed;               /* 1 if register has collapsed   */
    uint32_t       collapse_outcome;        /* Outcome if collapsed          */
    uint64_t       magic_base;              /* Magic pointer base address    */
    uint8_t        bulk_rule;               /* 0=general, 1=GHZ, 2=circuit  */

    /* Sparse amplitude storage */
    RegisterEntry  entries[4096];           /* Nonzero amplitudes            */
    uint32_t       num_nonzero;             /* Count of active entries       */

    /* GHZ / circuit mode metadata */
    uint16_t       gauss_n_dft;
    uint16_t       gauss_n_cz;
    uint8_t        gauss_ready;
    uint16_t       gauss_cz_a[256];
    uint16_t       gauss_cz_b[256];
} QubitRegister;

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT ENGINE — Top-level container
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    Qubit          qubits[MAX_QUBITS];
    QubitPair      pairs[MAX_PAIRS];
    QubitRegister  registers[MAX_REGISTERS];

    uint32_t       num_qubits;
    uint32_t       num_pairs;
    uint32_t       num_registers;

    /* PRNG state (LCG) */
    uint64_t       rng_state;
} QubitEngine;

/* ═══════════════════════════════════════════════════════════════════════════════
 * CORE API — qubit_core.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

void     qubit_engine_init(QubitEngine *eng);
void     qubit_engine_destroy(QubitEngine *eng);
uint32_t qubit_init(QubitEngine *eng);                /* Allocate |0⟩       */
uint32_t qubit_init_plus(QubitEngine *eng);            /* Allocate |+⟩       */
uint32_t qubit_init_basis(QubitEngine *eng, int k);    /* Allocate |k⟩       */
double   qubit_prng_double(QubitEngine *eng);           /* Random [0,1)       */

/* ═══════════════════════════════════════════════════════════════════════════════
 * GATE API — qubit_gates.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_hadamard(QubitEngine *eng, uint32_t id);
void qubit_apply_x(QubitEngine *eng, uint32_t id);
void qubit_apply_y(QubitEngine *eng, uint32_t id);
void qubit_apply_z(QubitEngine *eng, uint32_t id);
void qubit_apply_s(QubitEngine *eng, uint32_t id);
void qubit_apply_t(QubitEngine *eng, uint32_t id);
void qubit_apply_phase(QubitEngine *eng, uint32_t id, double theta);
void qubit_apply_cz(QubitEngine *eng, uint32_t id_a, uint32_t id_b);
void qubit_apply_cx(QubitEngine *eng, uint32_t ctrl, uint32_t target);
void qubit_apply_unitary(QubitEngine *eng, uint32_t id,
                         const double *U_re, const double *U_im);
void qubit_apply_unitary_pair(QubitEngine *eng, uint32_t id_a, uint32_t id_b,
                              const double *U_re, const double *U_im);

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT API — qubit_measure.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

int  qubit_measure(QubitEngine *eng, uint32_t id);
int  qubit_measure_in_pair(QubitEngine *eng, uint32_t id);

typedef struct {
    double prob[QUBIT_D];
    double entropy;
    double purity;
    int    schmidt_rank;
} QubitInspect;

QubitInspect qubit_inspect(QubitEngine *eng, uint32_t id);

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENTANGLEMENT API — qubit_entangle.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

int  qubit_entangle_bell(QubitEngine *eng, uint32_t id_a, uint32_t id_b);
int  qubit_entangle_product(QubitEngine *eng, uint32_t id_a, uint32_t id_b);
void qubit_disentangle(QubitEngine *eng, uint32_t id_a, uint32_t id_b);

/* ═══════════════════════════════════════════════════════════════════════════════
 * REGISTER API — qubit_register.c
 * ═══════════════════════════════════════════════════════════════════════════════ */

int         qubit_reg_init(QubitEngine *eng, uint64_t chunk_id,
                           uint64_t n_qubits, uint32_t dim);
void        qubit_reg_entangle_all(QubitEngine *eng, int reg_idx);
void        qubit_reg_apply_hadamard(QubitEngine *eng, int reg_idx,
                                     uint64_t qubit_idx);
void        qubit_reg_apply_cz(QubitEngine *eng, int reg_idx,
                                uint64_t idx_a, uint64_t idx_b);
void        qubit_reg_apply_unitary_pos(QubitEngine *eng, int reg_idx,
                                         uint64_t pos,
                                         const double *U_re, const double *U_im);
uint64_t    qubit_reg_measure(QubitEngine *eng, int reg_idx, uint64_t qubit_idx);
SV_Amplitude qubit_reg_sv_get(QubitEngine *eng, int reg_idx, basis_t basis_k);
void        qubit_reg_sv_set(QubitEngine *eng, int reg_idx,
                              basis_t basis_k, double re, double im);
void        qubit_reg_sv_stream(QubitEngine *eng, int reg_idx,
                                 sv_stream_fn callback, void *user_data);
double      qubit_reg_sv_total_prob(QubitEngine *eng, int reg_idx);
SV_Amplitude qubit_reg_sv_inner(QubitEngine *eng, int reg_a, int reg_b);
QubitState  qubit_reg_local_sv(QubitEngine *eng, int reg_idx, uint64_t qubit_pos);

#endif /* QUBIT_ENGINE_H */
