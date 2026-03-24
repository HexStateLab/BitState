/*
 * qubit_measure.c — Measurement, Collapse, and Non-Destructive Inspection
 *
 * Born-rule sampling for D=2.
 * For D=2, sampling is trivial: compute P(0), compare against random.
 */

#include "qubit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASURE LOCAL QUBIT — Born-rule sampling + collapse
 *
 * Returns outcome 0 or 1.
 * Post-measurement: qubit collapses to |outcome⟩ with original phase.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int qubit_measure(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return -1;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        return qubit_measure_in_pair(eng, id);
    }

    /* Born sampling */
    double rand_01 = qubit_prng_double(eng);
    int outcome = born_sample(q->state.re, q->state.im, QUBIT_D, rand_01);

    /* Collapse */
    born_collapse(q->state.re, q->state.im, QUBIT_D, outcome);

    return outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASURE IN PAIR — Measure one qubit of an entangled pair
 *
 * Computes marginal probabilities from joint state, samples,
 * then collapses both qubits according to the measurement outcome.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int qubit_measure_in_pair(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return -1;
    Qubit *q = &eng->qubits[id];
    if (q->pair_id < 0) return qubit_measure(eng, id);

    QubitPair *p = &eng->pairs[q->pair_id];
    int side = q->pair_side;

    /* Compute marginal probabilities */
    double probs[QUBIT_D] = {0};
    if (side == 0) {
        ent_marginal_a(&p->joint, probs);
    } else {
        ent_marginal_b(&p->joint, probs);
    }

    /* Born sampling */
    double rand_01 = qubit_prng_double(eng);
    int outcome = (rand_01 < probs[0]) ? 0 : 1;

    /* Collapse joint state: zero out entries where measured qubit ≠ outcome */
    for (int a = 0; a < QUBIT_D; a++)
    for (int b = 0; b < QUBIT_D; b++) {
        int s = (side == 0) ? a : b;
        if (s != outcome) {
            int idx = a * QUBIT_D + b;
            p->joint.re[idx] = 0;
            p->joint.im[idx] = 0;
        }
    }

    /* Renormalize joint state */
    ent_renormalize(&p->joint);

    return outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * INSPECT — Non-destructive qubit analysis
 *
 * Returns probabilities, entropy, purity, and Schmidt rank (if entangled).
 * Does NOT collapse the state.
 * ═══════════════════════════════════════════════════════════════════════════════ */

QubitInspect qubit_inspect(QubitEngine *eng, uint32_t id)
{
    QubitInspect info;
    memset(&info, 0, sizeof(info));

    if (id >= eng->num_qubits) return info;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        /* Entangled: compute marginals from joint state */
        QubitPair *p = &eng->pairs[q->pair_id];

        if (q->pair_side == 0) {
            ent_marginal_a(&p->joint, info.prob);
        } else {
            ent_marginal_b(&p->joint, info.prob);
        }

        info.schmidt_rank = ent_schmidt_rank(&p->joint);
    } else {
        /* Local: direct probability computation */
        for (int k = 0; k < QUBIT_D; k++) {
            info.prob[k] = q->state.re[k] * q->state.re[k]
                         + q->state.im[k] * q->state.im[k];
        }
        info.schmidt_rank = 1;
    }

    /* Entropy: S = -Σ P(k) log₂ P(k) */
    info.entropy = 0;
    for (int k = 0; k < QUBIT_D; k++) {
        if (info.prob[k] > 1e-14)
            info.entropy -= info.prob[k] * log2(info.prob[k]);
    }

    /* Purity: Tr(ρ²) = Σ P(k)² for pure states */
    info.purity = 0;
    for (int k = 0; k < QUBIT_D; k++) {
        info.purity += info.prob[k] * info.prob[k];
    }

    return info;
}
