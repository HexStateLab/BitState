/*
 * qubit_core.c — QubitEngine Lifecycle, PRNG, and Qubit Initialization
 *
 * Engine init/destroy, PRNG (LCG with Java's constants), and
 * qubit allocation to |0⟩, |+⟩, or |k⟩.
 */

#include "qubit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRNG — Linear Congruential Generator
 *
 * Same constants as java.util.Random:
 *   state = state * 6364136223846793005 + 1442695040888963407
 * ═══════════════════════════════════════════════════════════════════════════════ */

#define PRNG_A  6364136223846793005ULL
#define PRNG_C  1442695040888963407ULL

static inline uint64_t prng_next(uint64_t *state)
{
    *state = *state * PRNG_A + PRNG_C;
    return *state;
}

double qubit_prng_double(QubitEngine *eng)
{
    uint64_t bits = prng_next(&eng->rng_state);
    /* Use top 53 bits for double precision */
    return (double)(bits >> 11) * (1.0 / 9007199254740992.0);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENGINE LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_engine_init(QubitEngine *eng)
{
    eng->num_qubits    = 0;
    eng->num_pairs     = 0;
    eng->num_registers = 0;
    eng->rng_state     = 0x5DEECE66DULL;  /* Java seed */
}

void qubit_engine_destroy(QubitEngine *eng)
{
    eng->num_qubits    = 0;
    eng->num_pairs     = 0;
    eng->num_registers = 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT ALLOCATION
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Allocate a qubit in state |0⟩ */
uint32_t qubit_init(QubitEngine *eng)
{
    if (eng->num_qubits >= MAX_QUBITS) {
        fprintf(stderr, "[QUBIT] ERROR: max qubits (%d) reached\n", MAX_QUBITS);
        return UINT32_MAX;
    }

    uint32_t id = eng->num_qubits++;
    Qubit *q = &eng->qubits[id];
    qm_init_zero(&q->state);
    q->pair_id   = -1;
    q->pair_side = -1;
    q->id        = id;
    return id;
}

/* Allocate a qubit in state |+⟩ = (|0⟩+|1⟩)/√2 */
uint32_t qubit_init_plus(QubitEngine *eng)
{
    if (eng->num_qubits >= MAX_QUBITS) return UINT32_MAX;

    uint32_t id = eng->num_qubits++;
    Qubit *q = &eng->qubits[id];
    qm_init_plus(&q->state);
    q->pair_id   = -1;
    q->pair_side = -1;
    q->id        = id;
    return id;
}

/* Allocate a qubit in state |k⟩ */
uint32_t qubit_init_basis(QubitEngine *eng, int k)
{
    if (eng->num_qubits >= MAX_QUBITS) return UINT32_MAX;
    if (k < 0 || k >= QUBIT_D) return UINT32_MAX;

    uint32_t id = eng->num_qubits++;
    Qubit *q = &eng->qubits[id];
    qm_init_basis(&q->state, k);
    q->pair_id   = -1;
    q->pair_side = -1;
    q->id        = id;
    return id;
}
