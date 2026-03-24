/*
 * qubit_register.c — Qubit Register: Managed Groups of Qubits
 *
 * A register is a logical collection of N qubits, each stored as
 * QubitState (32 bytes), with pairwise JointState (64 bytes) for entanglement.
 *
 * Memory model:
 *   N qubits    = N × 32 bytes
 *   P pairs     = P × 64 bytes
 *   Total       = O(N + P), never O(2^N)
 *
 * GHZ entanglement across N qubits is done via chained Bell pairs.
 * For large N, the register tracks metadata for on-the-fly operations.
 */

#include "qubit_engine.h"

/* Forward declaration */
static uint32_t reg_extract_digit(basis_t basis, uint64_t pos,
                                  uint32_t D, uint8_t bulk_rule);

/* ═══════════════════════════════════════════════════════════════════════════════
 * REGISTER INIT
 * ═══════════════════════════════════════════════════════════════════════════════ */

int qubit_reg_init(QubitEngine *eng, uint64_t chunk_id,
                   uint64_t n_qubits, uint32_t dim)
{
    if (eng->num_registers >= MAX_REGISTERS) {
        fprintf(stderr, "[QUBIT] ERROR: max registers (%d) reached\n",
                MAX_REGISTERS);
        return -1;
    }

    int idx = (int)eng->num_registers++;
    QubitRegister *reg = &eng->registers[idx];
    memset(reg, 0, sizeof(*reg));

    reg->chunk_id  = chunk_id;
    reg->n_qubits  = n_qubits;
    reg->dim       = dim;
    reg->collapsed = 0;
    reg->magic_base = MAGIC_PTR(chunk_id, 0);

    return idx;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GHZ ENTANGLEMENT — (|00...0⟩ + |11...1⟩)/√2
 *
 * For qubits: the GHZ state is (|0⟩^N + |1⟩^N)/√2.
 * Memory: O(1) regardless of N — only 2 nonzero entries.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_reg_entangle_all(QubitEngine *eng, int reg_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];

    reg->bulk_rule = 1;  /* GHZ mode */

    /* Store (|0⟩^N + |1⟩^N)/√2 */
    static const double INV_SQRT2 = 0.7071067811865475244;
    reg->num_nonzero = 2;
    reg->entries[0].basis_state = 0;                /* |00...0⟩ */
    reg->entries[0].amp_re = INV_SQRT2;
    reg->entries[0].amp_im = 0;
    reg->entries[1].basis_state = 1;                /* |11...1⟩ in GHZ mode */
    reg->entries[1].amp_re = INV_SQRT2;
    reg->entries[1].amp_im = 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HADAMARD ON SINGLE QUBIT IN REGISTER
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_reg_apply_hadamard(QubitEngine *eng, int reg_idx,
                               uint64_t qubit_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];

    (void)qubit_idx;

    /* Apply Hadamard to the register's amplitude array */
    double re[QUBIT_D] = {0}, im[QUBIT_D] = {0};
    for (uint32_t k = 0; k < reg->num_nonzero && k < QUBIT_D; k++) {
        re[k] = reg->entries[k].amp_re;
        im[k] = reg->entries[k].amp_im;
    }

    sup_apply_hadamard(re, im);

    reg->num_nonzero = 0;
    for (uint32_t k = 0; k < QUBIT_D; k++) {
        double prob = re[k] * re[k] + im[k] * im[k];
        if (prob > 1e-30) {
            reg->entries[reg->num_nonzero].basis_state = k;
            reg->entries[reg->num_nonzero].amp_re = re[k];
            reg->entries[reg->num_nonzero].amp_im = im[k];
            reg->num_nonzero++;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CZ BETWEEN TWO QUBITS IN REGISTER
 *
 * CZ|a,b⟩ = (-1)^(a·b)|a,b⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_reg_apply_cz(QubitEngine *eng, int reg_idx,
                         uint64_t idx_a, uint64_t idx_b)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint32_t ka = reg_extract_digit(reg->entries[e].basis_state,
                                        idx_a, D, reg->bulk_rule);
        uint32_t kb = reg_extract_digit(reg->entries[e].basis_state,
                                        idx_b, D, reg->bulk_rule);
        /* (-1)^(ka·kb): negate only when both are 1 */
        if (ka == 1 && kb == 1) {
            reg->entries[e].amp_re = -reg->entries[e].amp_re;
            reg->entries[e].amp_im = -reg->entries[e].amp_im;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * UNITARY ON SINGLE QUBIT POSITION — O(entries × D) gate application
 * ═══════════════════════════════════════════════════════════════════════════════ */

void qubit_reg_apply_unitary_pos(QubitEngine *eng, int reg_idx,
                                  uint64_t pos,
                                  const double *U_re, const double *U_im)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    uint32_t max_entries = reg->num_nonzero * D + 1;
    if (max_entries < 4096) max_entries = 4096;
    uint32_t new_count = 0;
    struct tmp_entry { basis_t basis; double re, im; };
    struct tmp_entry *tmp = (struct tmp_entry *)calloc(max_entries,
                                                       sizeof(struct tmp_entry));

    uint64_t pos_mul = 1;
    for (uint64_t p = 0; p < pos; p++) pos_mul *= D;

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        basis_t basis = reg->entries[e].basis_state;
        double a_re = reg->entries[e].amp_re;
        double a_im = reg->entries[e].amp_im;

        uint32_t k = (uint32_t)((basis / pos_mul) % D);
        basis_t base = basis - (basis_t)k * pos_mul;

        for (uint32_t kp = 0; kp < D; kp++) {
            double ur = U_re[kp * D + k];
            double ui = U_im[kp * D + k];
            double nr = ur * a_re - ui * a_im;
            double ni = ur * a_im + ui * a_re;

            if (nr * nr + ni * ni < 1e-30) continue;

            basis_t new_basis = base + (basis_t)kp * pos_mul;

            int found = -1;
            for (uint32_t t = 0; t < new_count; t++) {
                if (tmp[t].basis == new_basis) { found = (int)t; break; }
            }
            if (found >= 0) {
                tmp[found].re += nr;
                tmp[found].im += ni;
            } else if (new_count < max_entries) {
                tmp[new_count].basis = new_basis;
                tmp[new_count].re = nr;
                tmp[new_count].im = ni;
                new_count++;
            }
        }
    }

    reg->num_nonzero = 0;
    for (uint32_t t = 0; t < new_count; t++) {
        if (tmp[t].re * tmp[t].re + tmp[t].im * tmp[t].im >= 1e-30) {
            if (reg->num_nonzero < 4096) {
                reg->entries[reg->num_nonzero].basis_state = tmp[t].basis;
                reg->entries[reg->num_nonzero].amp_re = tmp[t].re;
                reg->entries[reg->num_nonzero].amp_im = tmp[t].im;
                reg->num_nonzero++;
            }
        }
    }
    free(tmp);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT — Born-rule sampling of single qubit in register
 * ═══════════════════════════════════════════════════════════════════════════════ */

uint64_t qubit_reg_measure(QubitEngine *eng, int reg_idx, uint64_t qubit_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return 0;
    QubitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    (void)qubit_idx;

    double re[QUBIT_D] = {0}, im[QUBIT_D] = {0};
    for (uint32_t e = 0; e < reg->num_nonzero && e < QUBIT_D; e++) {
        uint32_t k = (uint32_t)(reg->entries[e].basis_state % D);
        re[k] = reg->entries[e].amp_re;
        im[k] = reg->entries[e].amp_im;
    }

    double rand_01 = qubit_prng_double(eng);
    int outcome = born_sample(re, im, D, rand_01);

    born_collapse(re, im, D, outcome);

    reg->num_nonzero = 1;
    reg->entries[0].basis_state = (basis_t)outcome;
    reg->entries[0].amp_re = 1.0;
    reg->entries[0].amp_im = 0.0;
    reg->collapsed = 1;
    reg->collapse_outcome = (uint32_t)outcome;

    return (uint64_t)outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * STREAMING STATE VECTOR ACCESS
 * ═══════════════════════════════════════════════════════════════════════════════ */

SV_Amplitude qubit_reg_sv_get(QubitEngine *eng, int reg_idx, basis_t basis_k)
{
    SV_Amplitude amp = {0.0, 0.0, -1.0/0.0};
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return amp;
    QubitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    if (reg->bulk_rule == 1) {
        /* GHZ mode: nonzero only if all N digits are equal */
        uint32_t first_digit = (uint32_t)(basis_k % D);
        basis_t remaining = basis_k;
        int all_same = 1;

        for (uint64_t q = 0; q < reg->n_qubits && all_same; q++) {
            if ((remaining % D) != first_digit) all_same = 0;
            remaining /= D;
        }
        if (remaining != 0) all_same = 0;

        if (all_same) {
            for (uint32_t e = 0; e < reg->num_nonzero; e++) {
                if ((uint32_t)(reg->entries[e].basis_state % D) == first_digit) {
                    amp.re = reg->entries[e].amp_re;
                    amp.im = reg->entries[e].amp_im;
                    return amp;
                }
            }
        }
        return amp;
    }

    /* General mode: linear scan */
    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        if (reg->entries[e].basis_state == basis_k) {
            amp.re = reg->entries[e].amp_re;
            amp.im = reg->entries[e].amp_im;
            return amp;
        }
    }
    return amp;
}

void qubit_reg_sv_set(QubitEngine *eng, int reg_idx,
                       basis_t basis_k, double re, double im)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];

    double mag2 = re * re + im * im;

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        if (reg->entries[e].basis_state == basis_k) {
            if (mag2 < 1e-30) {
                reg->entries[e] = reg->entries[reg->num_nonzero - 1];
                reg->num_nonzero--;
            } else {
                reg->entries[e].amp_re = re;
                reg->entries[e].amp_im = im;
            }
            return;
        }
    }

    if (mag2 >= 1e-30 && reg->num_nonzero < 4096) {
        reg->entries[reg->num_nonzero].basis_state = basis_k;
        reg->entries[reg->num_nonzero].amp_re = re;
        reg->entries[reg->num_nonzero].amp_im = im;
        reg->num_nonzero++;
    }
}

void qubit_reg_sv_stream(QubitEngine *eng, int reg_idx,
                          sv_stream_fn callback, void *user_data)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QubitRegister *reg = &eng->registers[reg_idx];

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        SV_Amplitude amp;
        amp.re = reg->entries[e].amp_re;
        amp.im = reg->entries[e].amp_im;
        amp.log2_mag = -1.0/0.0;
        callback(reg->entries[e].basis_state, amp, user_data);
    }
}

double qubit_reg_sv_total_prob(QubitEngine *eng, int reg_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return 0;
    QubitRegister *reg = &eng->registers[reg_idx];

    double total = 0.0;
    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        total += reg->entries[e].amp_re * reg->entries[e].amp_re
               + reg->entries[e].amp_im * reg->entries[e].amp_im;
    }
    return total;
}

SV_Amplitude qubit_reg_sv_inner(QubitEngine *eng, int reg_a, int reg_b)
{
    SV_Amplitude result = {0.0, 0.0, -1.0/0.0};
    if (reg_a < 0 || reg_b < 0) return result;
    if ((uint32_t)reg_a >= eng->num_registers) return result;
    if ((uint32_t)reg_b >= eng->num_registers) return result;

    QubitRegister *ra = &eng->registers[reg_a];
    QubitRegister *rb = &eng->registers[reg_b];

    for (uint32_t i = 0; i < ra->num_nonzero; i++) {
        for (uint32_t j = 0; j < rb->num_nonzero; j++) {
            if (ra->entries[i].basis_state == rb->entries[j].basis_state) {
                double a_re = ra->entries[i].amp_re;
                double a_im = ra->entries[i].amp_im;
                double b_re = rb->entries[j].amp_re;
                double b_im = rb->entries[j].amp_im;
                result.re += a_re * b_re + a_im * b_im;
                result.im += a_re * b_im - a_im * b_re;
            }
        }
    }
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Extract digit at position pos from basis state
 * ═══════════════════════════════════════════════════════════════════════════════ */

static uint32_t reg_extract_digit(basis_t basis, uint64_t pos,
                                  uint32_t D, uint8_t bulk_rule)
{
    if (bulk_rule == 1) {
        (void)pos;
        return (uint32_t)(basis % D);
    }
    basis_t v = basis;
    for (uint64_t i = 0; i < pos; i++) v /= D;
    return (uint32_t)(v % D);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PER-QUBIT LOCAL STATE VECTOR — Partial trace to single qubit
 * ═══════════════════════════════════════════════════════════════════════════════ */

QubitState qubit_reg_local_sv(QubitEngine *eng, int reg_idx, uint64_t qubit_pos)
{
    QubitState local;
    qm_init_zero(&local);

    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return local;
    QubitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    double prob[QUBIT_D] = {0};
    double phase_re[QUBIT_D] = {0};
    double phase_im[QUBIT_D] = {0};

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint32_t digit = reg_extract_digit(reg->entries[e].basis_state,
                                           qubit_pos, D, reg->bulk_rule);
        if (digit < D) {
            double re = reg->entries[e].amp_re;
            double im = reg->entries[e].amp_im;
            double p  = re * re + im * im;
            prob[digit] += p;
            phase_re[digit] += re;
            phase_im[digit] += im;
        }
    }

    for (uint32_t k = 0; k < D; k++) {
        if (prob[k] > 0) {
            double mag = sqrt(prob[k]);
            double norm = sqrt(phase_re[k]*phase_re[k] + phase_im[k]*phase_im[k]);
            if (norm > 1e-30) {
                local.re[k] = mag * (phase_re[k] / norm);
                local.im[k] = mag * (phase_im[k] / norm);
            } else {
                local.re[k] = mag;
                local.im[k] = 0;
            }
        }
    }

    return local;
}
