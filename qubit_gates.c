/*
 * qubit_gates.c — Quantum Gate Operations for D=2
 *
 * Standard qubit gates: H, X, Y, Z, S, T, Phase, CZ, CX, arbitrary unitary.
 * Each gate handles both local qubits and entangled pairs.
 *
 * CZ phase table for D=2:
 *   ω = e^(2πi/2) = -1
 *   CZ|a,b⟩ = (-1)^(a·b)|a,b⟩
 *   Only |1,1⟩ picks up a phase of -1.
 */

#include "qubit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * HADAMARD (DFT₂) — The square's fundamental symmetry
 *
 *   H = (1/√2) [[1,  1],
 *                [1, -1]]
 *
 * Local: apply to 2 amplitudes.
 * Entangled: apply to 4 joint amplitudes (on one side).
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_hadamard(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        /* Entangled: apply H to one side of the joint state */
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        double new_re[QUBIT_D2], new_im[QUBIT_D2];
        memset(new_re, 0, sizeof(new_re));
        memset(new_im, 0, sizeof(new_im));

        static const double INV_SQRT2 = 0.7071067811865475244;

        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int idx = a * QUBIT_D + b;
            /* Apply H to the 'side' qubit index */
            for (int k = 0; k < QUBIT_D; k++) {
                int new_idx;
                double h_val;

                if (side == 0) {
                    new_idx = k * QUBIT_D + b;
                    h_val = DFT2_RE[k][a];  /* H[k,a] */
                } else {
                    new_idx = a * QUBIT_D + k;
                    h_val = DFT2_RE[k][b];  /* H[k,b] */
                }

                new_re[new_idx] += h_val * p->joint.re[idx];
                new_im[new_idx] += h_val * p->joint.im[idx];
            }
        }

        memcpy(p->joint.re, new_re, sizeof(new_re));
        memcpy(p->joint.im, new_im, sizeof(new_im));
    } else {
        /* Local: direct Hadamard */
        sup_apply_hadamard(q->state.re, q->state.im);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI X — Bit flip: |0⟩↔|1⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_x(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        double new_re[QUBIT_D2], new_im[QUBIT_D2];
        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int idx_old = a * QUBIT_D + b;
            int idx_new;
            if (side == 0)
                idx_new = (1 - a) * QUBIT_D + b;
            else
                idx_new = a * QUBIT_D + (1 - b);
            new_re[idx_new] = p->joint.re[idx_old];
            new_im[idx_new] = p->joint.im[idx_old];
        }
        memcpy(p->joint.re, new_re, sizeof(new_re));
        memcpy(p->joint.im, new_im, sizeof(new_im));
    } else {
        double tr = q->state.re[0], ti = q->state.im[0];
        q->state.re[0] = q->state.re[1]; q->state.im[0] = q->state.im[1];
        q->state.re[1] = tr;              q->state.im[1] = ti;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI Y — |0⟩→i|1⟩, |1⟩→-i|0⟩
 * Y = iXZ
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_y(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        double new_re[QUBIT_D2], new_im[QUBIT_D2];
        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int idx = a * QUBIT_D + b;
            int s = (side == 0) ? a : b;  /* the index being acted on */
            int idx_new;
            if (side == 0)
                idx_new = (1 - a) * QUBIT_D + b;
            else
                idx_new = a * QUBIT_D + (1 - b);

            /* Y[0,1] = -i,  Y[1,0] = i */
            double sign_re, sign_im;
            if (s == 0) {
                /* Y|0⟩ = i|1⟩:  multiply by i = (0 + i) */
                sign_re = 0.0; sign_im = 1.0;
            } else {
                /* Y|1⟩ = -i|0⟩: multiply by -i = (0 - i) */
                sign_re = 0.0; sign_im = -1.0;
            }

            double r = p->joint.re[idx], m = p->joint.im[idx];
            new_re[idx_new] = sign_re * r - sign_im * m;
            new_im[idx_new] = sign_re * m + sign_im * r;
        }
        memcpy(p->joint.re, new_re, sizeof(new_re));
        memcpy(p->joint.im, new_im, sizeof(new_im));
    } else {
        /* Y|ψ⟩: swap + phase */
        double r0 = q->state.re[0], i0 = q->state.im[0];
        double r1 = q->state.re[1], i1 = q->state.im[1];
        /* |0⟩ → i|1⟩ : new[1] += i * old[0] */
        /* |1⟩ → -i|0⟩: new[0] += -i * old[1] */
        q->state.re[0] =  i1;  /* -i * (r1 + i*i1) → i1 - i*r1 */
        q->state.im[0] = -r1;
        q->state.re[1] = -i0;  /* i * (r0 + i*i0) → -i0 + i*r0 */
        q->state.im[1] =  r0;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PAULI Z — Phase flip: |0⟩→|0⟩, |1⟩→-|1⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_z(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int s = (side == 0) ? a : b;
            if (s == 1) {
                int idx = a * QUBIT_D + b;
                p->joint.re[idx] = -p->joint.re[idx];
                p->joint.im[idx] = -p->joint.im[idx];
            }
        }
    } else {
        q->state.re[1] = -q->state.re[1];
        q->state.im[1] = -q->state.im[1];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * S GATE — |0⟩→|0⟩, |1⟩→i|1⟩    (√Z)
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_s(QubitEngine *eng, uint32_t id)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int s = (side == 0) ? a : b;
            if (s == 1) {
                int idx = a * QUBIT_D + b;
                /* multiply by i: (r + im*i) * i = -im + r*i */
                double r = p->joint.re[idx], m = p->joint.im[idx];
                p->joint.re[idx] = -m;
                p->joint.im[idx] =  r;
            }
        }
    } else {
        double r = q->state.re[1], m = q->state.im[1];
        q->state.re[1] = -m;
        q->state.im[1] =  r;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * T GATE — |0⟩→|0⟩, |1⟩→e^{iπ/4}|1⟩    (√S)
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_t(QubitEngine *eng, uint32_t id)
{
    qubit_apply_phase(eng, id, M_PI / 4.0);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PHASE GATE — |0⟩→|0⟩, |1⟩→e^{iθ}|1⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_phase(QubitEngine *eng, uint32_t id, double theta)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    double cs = cos(theta), sn = sin(theta);

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int s = (side == 0) ? a : b;
            if (s == 1) {
                int idx = a * QUBIT_D + b;
                double r = p->joint.re[idx], m = p->joint.im[idx];
                p->joint.re[idx] = r * cs - m * sn;
                p->joint.im[idx] = r * sn + m * cs;
            }
        }
    } else {
        double r = q->state.re[1], m = q->state.im[1];
        q->state.re[1] = r * cs - m * sn;
        q->state.im[1] = r * sn + m * cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CZ GATE — Controlled-Z: |1,1⟩ → -|1,1⟩
 *
 * For D=2, ω = e^(2πi/2) = -1.
 * CZ|a,b⟩ = (-1)^(a·b)|a,b⟩
 * Only |1,1⟩ gets a phase of -1; all others unchanged.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_cz(QubitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_qubits || id_b >= eng->num_qubits) return;

    Qubit *qa = &eng->qubits[id_a];
    Qubit *qb = &eng->qubits[id_b];

    if (qa->pair_id >= 0 && qa->pair_id == qb->pair_id) {
        /* Already entangled as a pair — apply CZ to joint state */
        QubitPair *p = &eng->pairs[qa->pair_id];
        /* |1,1⟩ → -|1,1⟩ */
        p->joint.re[3] = -p->joint.re[3];
        p->joint.im[3] = -p->joint.im[3];
    } else if (qa->pair_id < 0 && qb->pair_id < 0) {
        /* Both local: create product pair, then apply CZ */
        int slot = qubit_entangle_product(eng, id_a, id_b);
        if (slot < 0) return;
        QubitPair *p = &eng->pairs[slot];
        p->joint.re[3] = -p->joint.re[3];
        p->joint.im[3] = -p->joint.im[3];
    } else {
        fprintf(stderr, "[QUBIT] CZ: qubits in different pairs not supported\n");
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CX (CNOT) GATE — Controlled-X: |c,t⟩ → |c, c⊕t⟩
 *
 * CX = (I⊗H)·CZ·(I⊗H)
 * Or directly: swap |1,0⟩ ↔ |1,1⟩ in joint state.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_cx(QubitEngine *eng, uint32_t ctrl, uint32_t target)
{
    if (ctrl >= eng->num_qubits || target >= eng->num_qubits) return;

    Qubit *qc = &eng->qubits[ctrl];
    Qubit *qt = &eng->qubits[target];

    /* Ensure they share a pair */
    if (qc->pair_id < 0 || qt->pair_id < 0) {
        /* Create product pair first */
        int slot = qubit_entangle_product(eng, ctrl, target);
        if (slot < 0) return;
    }

    if (qc->pair_id != qt->pair_id) {
        fprintf(stderr, "[QUBIT] CX: qubits in different pairs not supported\n");
        return;
    }

    QubitPair *p = &eng->pairs[qc->pair_id];
    int sc = qc->pair_side;

    if (sc == 0) {
        /* ctrl=A, target=B: swap |1,0⟩ ↔ |1,1⟩ (indices 2 and 3) */
        double tr, ti;
        tr = p->joint.re[2]; ti = p->joint.im[2];
        p->joint.re[2] = p->joint.re[3]; p->joint.im[2] = p->joint.im[3];
        p->joint.re[3] = tr;              p->joint.im[3] = ti;
    } else {
        /* ctrl=B, target=A: swap |0,1⟩ ↔ |1,1⟩ (indices 1 and 3) */
        double tr, ti;
        tr = p->joint.re[1]; ti = p->joint.im[1];
        p->joint.re[1] = p->joint.re[3]; p->joint.im[1] = p->joint.im[3];
        p->joint.re[3] = tr;              p->joint.im[3] = ti;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ARBITRARY 2×2 UNITARY — General single-qubit gate
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_unitary(QubitEngine *eng, uint32_t id,
                         const double *U_re, const double *U_im)
{
    if (id >= eng->num_qubits) return;
    Qubit *q = &eng->qubits[id];

    if (q->pair_id >= 0) {
        QubitPair *p = &eng->pairs[q->pair_id];
        int side = q->pair_side;

        double new_re[QUBIT_D2], new_im[QUBIT_D2];
        memset(new_re, 0, sizeof(new_re));
        memset(new_im, 0, sizeof(new_im));

        for (int a = 0; a < QUBIT_D; a++)
        for (int b = 0; b < QUBIT_D; b++) {
            int idx = a * QUBIT_D + b;
            for (int k = 0; k < QUBIT_D; k++) {
                int new_idx, s_old;
                if (side == 0) {
                    new_idx = k * QUBIT_D + b;
                    s_old = a;
                } else {
                    new_idx = a * QUBIT_D + k;
                    s_old = b;
                }

                int u = k * QUBIT_D + s_old;
                double ur = U_re[u], ui = U_im[u];
                double r = p->joint.re[idx], m = p->joint.im[idx];
                new_re[new_idx] += ur * r - ui * m;
                new_im[new_idx] += ur * m + ui * r;
            }
        }
        memcpy(p->joint.re, new_re, sizeof(new_re));
        memcpy(p->joint.im, new_im, sizeof(new_im));
    } else {
        /* Local 2×2 unitary */
        double r0 = q->state.re[0], i0 = q->state.im[0];
        double r1 = q->state.re[1], i1 = q->state.im[1];

        for (int k = 0; k < QUBIT_D; k++) {
            double sum_re = 0, sum_im = 0;
            double src_re[2] = {r0, r1};
            double src_im[2] = {i0, i1};
            for (int j = 0; j < QUBIT_D; j++) {
                int u = k * QUBIT_D + j;
                sum_re += U_re[u] * src_re[j] - U_im[u] * src_im[j];
                sum_im += U_re[u] * src_im[j] + U_im[u] * src_re[j];
            }
            q->state.re[k] = sum_re;
            q->state.im[k] = sum_im;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ARBITRARY 4×4 UNITARY ON PAIR — General two-qubit gate
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_apply_unitary_pair(QubitEngine *eng, uint32_t id_a, uint32_t id_b,
                              const double *U_re, const double *U_im)
{
    if (id_a >= eng->num_qubits || id_b >= eng->num_qubits) return;

    Qubit *qa = &eng->qubits[id_a];
    Qubit *qb = &eng->qubits[id_b];

    /* Ensure they share a pair */
    if (qa->pair_id < 0 || qa->pair_id != qb->pair_id) {
        if (qa->pair_id < 0 && qb->pair_id < 0) {
            int slot = qubit_entangle_product(eng, id_a, id_b);
            if (slot < 0) return;
        } else {
            fprintf(stderr, "[QUBIT] unitary_pair: qubits not in same pair\n");
            return;
        }
    }

    QubitPair *p = &eng->pairs[qa->pair_id];

    /* Apply 4×4 unitary to joint state */
    double new_re[QUBIT_D2] = {0}, new_im[QUBIT_D2] = {0};
    for (int i = 0; i < QUBIT_D2; i++) {
        for (int j = 0; j < QUBIT_D2; j++) {
            int u = i * QUBIT_D2 + j;
            new_re[i] += U_re[u] * p->joint.re[j] - U_im[u] * p->joint.im[j];
            new_im[i] += U_re[u] * p->joint.im[j] + U_im[u] * p->joint.re[j];
        }
    }
    memcpy(p->joint.re, new_re, sizeof(new_re));
    memcpy(p->joint.im, new_im, sizeof(new_im));
}
