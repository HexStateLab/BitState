/*
 * qubit_entangle.c — Entanglement: Bell, Product, Disentangle
 *
 * Two modes of entanglement creation:
 *   Bell:    (|00⟩+|11⟩)/√2 — maximally entangled (64 bytes)
 *   Product: |ψ_a⟩ ⊗ |ψ_b⟩ — separable (64 bytes, becomes entangled via CZ)
 *
 * Disentangle extracts marginals back into local qubit states.
 */

#include "qubit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Allocate a pair slot
 * ═══════════════════════════════════════════════════════════════════════════════ */

static int alloc_pair(QubitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (eng->num_pairs >= MAX_PAIRS) {
        fprintf(stderr, "[QUBIT] ERROR: max pairs (%d) reached\n", MAX_PAIRS);
        return -1;
    }

    int slot = (int)eng->num_pairs++;
    QubitPair *p = &eng->pairs[slot];
    memset(p, 0, sizeof(*p));
    p->id_a   = id_a;
    p->id_b   = id_b;
    p->active = 1;

    Qubit *qa = &eng->qubits[id_a];
    Qubit *qb = &eng->qubits[id_b];
    qa->pair_id   = slot;
    qa->pair_side = 0;
    qb->pair_id   = slot;
    qb->pair_side = 1;

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BELL PAIR — (|00⟩+|11⟩)/√2
 *
 * Maximally entangled state. 64 bytes.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int qubit_entangle_bell(QubitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_qubits || id_b >= eng->num_qubits) return -1;

    int slot = alloc_pair(eng, id_a, id_b);
    if (slot < 0) return -1;

    qm_entangle_bell(&eng->pairs[slot].joint);

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRODUCT PAIR — |ψ_a⟩ ⊗ |ψ_b⟩
 *
 * Creates joint state from tensor product of current local states.
 * Initially separable. CZ gate then creates genuine entanglement.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int qubit_entangle_product(QubitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_qubits || id_b >= eng->num_qubits) return -1;

    int slot = alloc_pair(eng, id_a, id_b);
    if (slot < 0) return -1;

    qm_entangle_product(&eng->pairs[slot].joint,
                        &eng->qubits[id_a].state,
                        &eng->qubits[id_b].state);

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DISENTANGLE — Break pair, extract marginals to local states
 *
 * For side A: ψ_A(k) = √(Σ_b |ψ(k,b)|²) × phase_of_dominant_component
 * This preserves probabilities but not full phase coherence
 * (which is physically correct — tracing out the partner destroys coherence).
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_disentangle(QubitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_qubits || id_b >= eng->num_qubits) return;

    Qubit *qa = &eng->qubits[id_a];
    Qubit *qb = &eng->qubits[id_b];

    /* Verify they share a pair */
    if (qa->pair_id < 0 || qa->pair_id != qb->pair_id) return;

    QubitPair *p = &eng->pairs[qa->pair_id];

    /* Extract marginal for A (rows) */
    for (int k = 0; k < QUBIT_D; k++) {
        double prob = 0;
        double dom_re = 0, dom_im = 0;
        double dom_p = 0;
        for (int b = 0; b < QUBIT_D; b++) {
            int idx = k * QUBIT_D + b;
            double r = p->joint.re[idx], i = p->joint.im[idx];
            double p2 = r * r + i * i;
            prob += p2;
            if (p2 > dom_p) { dom_p = p2; dom_re = r; dom_im = i; }
        }
        if (prob > 1e-30) {
            double amp = sqrt(prob);
            double phase = atan2(dom_im, dom_re);
            qa->state.re[k] = amp * cos(phase);
            qa->state.im[k] = amp * sin(phase);
        } else {
            qa->state.re[k] = 0;
            qa->state.im[k] = 0;
        }
    }

    /* Extract marginal for B (columns) */
    for (int k = 0; k < QUBIT_D; k++) {
        double prob = 0;
        double dom_re = 0, dom_im = 0;
        double dom_p = 0;
        for (int a = 0; a < QUBIT_D; a++) {
            int idx = a * QUBIT_D + k;
            double r = p->joint.re[idx], i = p->joint.im[idx];
            double p2 = r * r + i * i;
            prob += p2;
            if (p2 > dom_p) { dom_p = p2; dom_re = r; dom_im = i; }
        }
        if (prob > 1e-30) {
            double amp = sqrt(prob);
            double phase = atan2(dom_im, dom_re);
            qb->state.re[k] = amp * cos(phase);
            qb->state.im[k] = amp * sin(phase);
        } else {
            qb->state.re[k] = 0;
            qb->state.im[k] = 0;
        }
    }

    /* Renormalize both */
    sup_renormalize(qa->state.re, qa->state.im, QUBIT_D);
    sup_renormalize(qb->state.re, qb->state.im, QUBIT_D);

    /* Deactivate pair */
    p->active = 0;
    qa->pair_id = -1;
    qb->pair_id = -1;
}
