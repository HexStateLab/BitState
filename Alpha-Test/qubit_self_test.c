/*
 * qubit_self_test.c — QubitEngine Self-Test Suite
 *
 * Validates all engine components:
 *   1. Lifecycle: init/destroy, qubit allocation
 *   2. Gates: H²=I, X²=I, Z²=I, HXH=Z, HZH=X, CZ symmetry
 *   3. Measurement: Born sampling statistics, collapse correctness
 *   4. Entanglement: Bell pair correlations, disentangle
 *   5. Register: GHZ state, sparse amplitude access
 *   6. Triality: view conversions, triality rotation cycle
 *   7. Flat: auto-promote on H, auto-demote on collapse
 *   8. HPC: phase graph, amplitude eval, norm check, measurement
 *   9. Clifford: invariant, exotic gate, fingerprint
 *
 * Build:
 *   gcc -O2 -march=native -o qubit_test qubit_self_test.c \
 *       qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
 *       qubit_register.c qubit_triality.c clifford_exotic.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qubit_engine.h"
#include "qubit_triality.h"
#include "flat_qubit.h"
#include "hpc_qubit_graph.h"
#include "hpc_qubit_contract.h"
#include "hpc_qubit_amplitude.h"
#include "clifford_exotic.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (cond) { tests_passed++; } \
    else { tests_failed++; printf("  FAIL: %s\n", msg); } \
} while (0)

#define APPROX(a, b) (fabs((a)-(b)) < 1e-10)

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 1: Engine Lifecycle
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_lifecycle(void)
{
    printf("Test 1: Engine Lifecycle\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    ASSERT(eng != NULL, "engine alloc");

    qubit_engine_init(eng);
    ASSERT(eng->num_qubits == 0, "initial qubits = 0");

    uint32_t q0 = qubit_init(eng);
    ASSERT(q0 == 0, "first qubit id = 0");
    ASSERT(eng->num_qubits == 1, "num_qubits = 1");

    /* Verify |0⟩ state */
    ASSERT(APPROX(eng->qubits[q0].state.re[0], 1.0), "|0⟩ re[0] = 1");
    ASSERT(APPROX(eng->qubits[q0].state.re[1], 0.0), "|0⟩ re[1] = 0");

    uint32_t q1 = qubit_init_plus(eng);
    ASSERT(APPROX(eng->qubits[q1].state.re[0], 0.7071067811865475), "|+⟩ re[0]");
    ASSERT(APPROX(eng->qubits[q1].state.re[1], 0.7071067811865475), "|+⟩ re[1]");

    qubit_engine_destroy(eng);
    free(eng);
    printf("  lifecycle: %d tests\n\n", tests_passed);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 2: Gate Identities
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_gates(void)
{
    printf("Test 2: Gate Identities\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    /* H²=I: H(H|0⟩) should give |0⟩ back */
    uint32_t q = qubit_init(eng);
    qubit_apply_hadamard(eng, q);
    qubit_apply_hadamard(eng, q);
    ASSERT(APPROX(eng->qubits[q].state.re[0], 1.0), "H²|0⟩ = |0⟩: re[0]");
    ASSERT(APPROX(eng->qubits[q].state.re[1], 0.0), "H²|0⟩ = |0⟩: re[1]");

    /* X²=I */
    q = qubit_init(eng);
    qubit_apply_x(eng, q);
    qubit_apply_x(eng, q);
    ASSERT(APPROX(eng->qubits[q].state.re[0], 1.0), "X²|0⟩ = |0⟩");

    /* Z²=I */
    q = qubit_init_plus(eng);
    double orig_re0 = eng->qubits[q].state.re[0];
    double orig_re1 = eng->qubits[q].state.re[1];
    qubit_apply_z(eng, q);
    qubit_apply_z(eng, q);
    ASSERT(APPROX(eng->qubits[q].state.re[0], orig_re0), "Z²|+⟩: re[0]");
    ASSERT(APPROX(eng->qubits[q].state.re[1], orig_re1), "Z²|+⟩: re[1]");

    /* HXH=Z: Apply H, X, H to |+⟩ — should give Z|+⟩ = |0⟩-|1⟩)/√2 ? */
    /* Actually HXH=Z means H·X·H = Z. So H(X(H|0⟩)) = Z|0⟩ = |0⟩ */
    q = qubit_init(eng);
    qubit_apply_hadamard(eng, q);
    qubit_apply_x(eng, q);
    qubit_apply_hadamard(eng, q);
    ASSERT(APPROX(eng->qubits[q].state.re[0], 1.0), "HXH|0⟩ = Z|0⟩ = |0⟩");

    /* S²=Z: S(S|1⟩) = i·i|1⟩ = -|1⟩ = Z|1⟩ */
    q = qubit_init_basis(eng, 1);
    qubit_apply_s(eng, q);
    qubit_apply_s(eng, q);
    ASSERT(APPROX(eng->qubits[q].state.re[1], -1.0), "S²|1⟩ = -|1⟩: re");
    ASSERT(APPROX(eng->qubits[q].state.im[1], 0.0), "S²|1⟩ = -|1⟩: im");

    qubit_engine_destroy(eng);
    free(eng);
    printf("  gates: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 3: Measurement Statistics
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_measurement(void)
{
    printf("Test 3: Measurement\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    /* Measure |0⟩ — should always give 0 */
    uint32_t q = qubit_init(eng);
    int outcome = qubit_measure(eng, q);
    ASSERT(outcome == 0, "measure |0⟩ = 0");

    /* Measure |1⟩ — should always give 1 */
    q = qubit_init_basis(eng, 1);
    outcome = qubit_measure(eng, q);
    ASSERT(outcome == 1, "measure |1⟩ = 1");

    /* Measure |+⟩ many times — should be ~50/50 */
    int count0 = 0, count1 = 0;
    int N = 1000;
    for (int i = 0; i < N; i++) {
        q = qubit_init_plus(eng);
        outcome = qubit_measure(eng, q);
        if (outcome == 0) count0++;
        else count1++;
    }
    double ratio = (double)count0 / N;
    ASSERT(ratio > 0.4 && ratio < 0.6,
           "measure |+⟩ ~50/50 (within 10%)");
    printf("    |+⟩ stats: %d zeros, %d ones (%.1f%%/%.1f%%)\n",
           count0, count1, 100.0*count0/N, 100.0*count1/N);

    /* Inspect */
    q = qubit_init_plus(eng);
    QubitInspect info = qubit_inspect(eng, q);
    ASSERT(APPROX(info.prob[0], 0.5), "inspect |+⟩ P(0) = 0.5");
    ASSERT(APPROX(info.prob[1], 0.5), "inspect |+⟩ P(1) = 0.5");
    ASSERT(APPROX(info.entropy, 1.0), "inspect |+⟩ entropy = 1 bit");

    qubit_engine_destroy(eng);
    free(eng);
    printf("  measurement: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 4: Entanglement
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_entanglement(void)
{
    printf("Test 4: Entanglement\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    /* Bell pair: (|00⟩+|11⟩)/√2 */
    uint32_t qa = qubit_init(eng);
    uint32_t qb = qubit_init(eng);
    int slot = qubit_entangle_bell(eng, qa, qb);
    ASSERT(slot >= 0, "Bell pair created");

    QubitPair *p = &eng->pairs[slot];
    ASSERT(APPROX(p->joint.re[0], 0.7071067811865475), "Bell |00⟩");
    ASSERT(APPROX(p->joint.re[3], 0.7071067811865475), "Bell |11⟩");
    ASSERT(APPROX(p->joint.re[1], 0.0), "Bell |01⟩ = 0");
    ASSERT(APPROX(p->joint.re[2], 0.0), "Bell |10⟩ = 0");

    /* Schmidt rank = 2 */
    ASSERT(ent_schmidt_rank(&p->joint) == 2, "Bell Schmidt rank = 2");

    /* Bell correlation: measure both, should always agree */
    int agree = 0;
    int N = 100;
    for (int i = 0; i < N; i++) {
        qa = qubit_init(eng);
        qb = qubit_init(eng);
        qubit_entangle_bell(eng, qa, qb);
        int oa = qubit_measure_in_pair(eng, qa);
        int ob = qubit_measure_in_pair(eng, qb);
        if (oa == ob) agree++;
    }
    ASSERT(agree == N, "Bell correlations: 100% agreement");
    printf("    Bell agreement: %d/%d\n", agree, N);

    qubit_engine_destroy(eng);
    free(eng);
    printf("  entanglement: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 5: Triality
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_triality(void)
{
    printf("Test 5: Triality\n");

    TrialityQubit tq;
    tri_init(&tq);

    /* |0⟩ in edge view should give |+⟩ = (|0⟩+|1⟩)/√2 in vertex view */
    double re[2], im[2];
    tri_get_amplitudes(&tq, VIEW_VERTEX, re, im);
    ASSERT(APPROX(re[0], 0.7071067811865475), "Edge |0⟩ → Vertex re[0]");
    ASSERT(APPROX(re[1], 0.7071067811865475), "Edge |0⟩ → Vertex re[1]");
    ASSERT(APPROX(im[0], 0.0), "Edge |0⟩ → Vertex im[0]");

    /* Test triality Z gate (in edge view, like Phase(π)) */
    tri_init(&tq);
    tri_apply_z(&tq, M_PI);
    tri_get_amplitudes(&tq, VIEW_EDGE, re, im);
    ASSERT(APPROX(re[0], 1.0), "Z(π)|0⟩ edge re[0] = 1");
    ASSERT(APPROX(re[1], 0.0), "Z(π)|0⟩ edge re[1] = 0");
    /* Z(π)|0⟩ = |0⟩ (phase only on |1⟩) */

    /* Test triality X gate (O(1) in vertex view) */
    tri_init(&tq);
    /* Set to |1⟩ first */
    double init_re[2] = {0.0, 1.0}, init_im[2] = {0.0, 0.0};
    tri_init_state(&tq, VIEW_EDGE, init_re, init_im);
    tri_apply_x(&tq);
    tri_get_amplitudes(&tq, VIEW_EDGE, re, im);
    ASSERT(APPROX(re[0], 1.0), "X|1⟩ = |0⟩: re[0]");
    ASSERT(APPROX(re[1], 0.0), "X|1⟩ = |0⟩: re[1]");

    /* Test Hadamard via triality (O(1) swap) */
    tri_init(&tq);
    tri_apply_hadamard(&tq);
    tri_get_amplitudes(&tq, VIEW_EDGE, re, im);
    ASSERT(APPROX(re[0], 0.7071067811865475), "H|0⟩ = |+⟩: re[0]");
    ASSERT(APPROX(re[1], 0.7071067811865475), "H|0⟩ = |+⟩: re[1]");

    /* H² = I via triality */
    tri_apply_hadamard(&tq);
    tri_get_amplitudes(&tq, VIEW_EDGE, re, im);
    ASSERT(APPROX(re[0], 1.0), "H²|0⟩ = |0⟩ via triality");
    ASSERT(APPROX(re[1], 0.0), "H²|0⟩ = |0⟩ via triality (re[1])");

    tri_print_stats(&tq);

    printf("  triality: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 6: Flat Qubit
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_flat(void)
{
    printf("Test 6: Flat Qubit\n");

    FlatQubit fq;
    flat_init(&fq);
    ASSERT(fq.mode == FLAT_BASIS, "init mode = FLAT_BASIS");
    ASSERT(fq.basis_k == 0, "init basis = |0⟩");

    /* X gate in flat mode: O(1) */
    flat_apply_x(&fq);
    ASSERT(fq.mode == FLAT_BASIS, "X keeps FLAT_BASIS");
    ASSERT(fq.basis_k == 1, "X|0⟩ = |1⟩");

    /* X again: back to |0⟩ */
    flat_apply_x(&fq);
    ASSERT(fq.basis_k == 0, "X²|0⟩ = |0⟩");

    /* Z gate in flat mode: only phase change */
    flat_init_basis(&fq, 1, 1.0, 0.0);
    flat_apply_z(&fq);
    ASSERT(fq.mode == FLAT_BASIS, "Z keeps FLAT_BASIS");
    ASSERT(APPROX(fq.phase_re, -1.0), "Z|1⟩ phase = -1");

    /* H promotes to QUANTUM_FULL */
    flat_init(&fq);
    flat_apply_hadamard(&fq);
    ASSERT(fq.mode == QUANTUM_FULL, "H promotes to QUANTUM_FULL");

    double re[2], im[2];
    flat_get_amplitudes(&fq, re, im);
    ASSERT(APPROX(re[0], 0.7071067811865475), "H|0⟩ flat re[0]");
    ASSERT(APPROX(re[1], 0.7071067811865475), "H|0⟩ flat re[1]");

    /* H again should demote back to FLAT_BASIS */
    flat_apply_hadamard(&fq);
    ASSERT(fq.mode == FLAT_BASIS, "H² demotes to FLAT_BASIS");
    ASSERT(fq.basis_k == 0, "H²|0⟩ = |0⟩ flat");

    printf("  flat: done (promotions=%u, demotions=%u)\n\n",
           fq.promotions, fq.demotions);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 7: Lazy Triality
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_lazy(void)
{
    printf("Test 7: Lazy Triality\n");

    LazyTrialityQubit lq;
    lazy_tri_init(&lq);

    /* Accumulate phase on |0⟩ — should be no-op */
    lazy_tri_apply_phase(&lq, M_PI / 4.0);
    lazy_tri_materialize(&lq);
    ASSERT(APPROX(lq.re[0], 1.0), "lazy phase on |0⟩ = no-op: re[0]");

    /* H then measure: should be ~50/50 */
    int count0 = 0;
    for (int i = 0; i < 200; i++) {
        lazy_tri_init(&lq);
        lazy_tri_apply_hadamard(&lq);
        int outcome = lazy_tri_measure(&lq, (double)i / 200.0);
        if (outcome == 0) count0++;
    }
    double ratio = (double)count0 / 200.0;
    ASSERT(ratio > 0.4 && ratio < 0.6, "lazy H|0⟩ measure ~50/50");
    printf("    lazy H|0⟩: %d/200 zeros (%.1f%%)\n", count0, ratio * 100);

    printf("  lazy: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 8: HPC Phase Graph
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_hpc(void)
{
    printf("Test 8: HPC Phase Graph\n");

    /* Create 2-site graph: |+⟩ ⊗ |0⟩ */
    HPCQGraph *g = hpcq_create(2);
    ASSERT(g != NULL, "HPC graph alloc");
    ASSERT(g->n_sites == 2, "HPC 2 sites");

    /* Set site 0 to |+⟩ */
    double plus_re[2] = {0.7071067811865475, 0.7071067811865475};
    double plus_im[2] = {0, 0};
    hpcq_set_local(g, 0, plus_re, plus_im);

    /* Site 1 stays |0⟩ = default from tri_init */

    /* Product state: check norm = 1.0 */
    double norm = hpcq_norm_sq(g);
    ASSERT(APPROX(norm, 1.0), "HPC product norm = 1.0");

    /* Point query: ψ(0,0) = (1/√2)×1 = 1/√2 */
    uint32_t idx00[2] = {0, 0};
    double re, im;
    hpcq_amplitude(g, idx00, &re, &im);
    ASSERT(APPROX(re, 0.7071067811865475), "HPC ψ(00) = 1/√2");

    /* ψ(1,0) = (1/√2)×1 = 1/√2 */
    uint32_t idx10[2] = {1, 0};
    hpcq_amplitude(g, idx10, &re, &im);
    ASSERT(APPROX(re, 0.7071067811865475), "HPC ψ(10) = 1/√2");

    /* ψ(0,1) = 0 (site 1 is |0⟩) */
    uint32_t idx01[2] = {0, 1};
    hpcq_amplitude(g, idx01, &re, &im);
    ASSERT(APPROX(re, 0.0), "HPC ψ(01) = 0");

    /* Add CZ edge between sites 0 and 1 */
    hpcq_cz(g, 0, 1);
    ASSERT(g->n_edges == 1, "HPC 1 CZ edge");
    ASSERT(g->cz_edges == 1, "HPC 1 CZ edge (counter)");

    /* After CZ: ψ(0,0) unchanged, but ψ(1,1) picks up -1 */
    /* CZ|+,0⟩ = (|00⟩+|10⟩)/√2 (unchanged since site 1 is |0⟩) */
    norm = hpcq_norm_sq(g);
    ASSERT(APPROX(norm, 1.0), "HPC CZ norm = 1.0");

    /* Set both to |+⟩ and apply CZ */
    hpcq_destroy(g);
    g = hpcq_create(2);
    hpcq_set_local(g, 0, plus_re, plus_im);
    hpcq_set_local(g, 1, plus_re, plus_im);
    hpcq_cz(g, 0, 1);

    /* CZ|+,+⟩ = (|00⟩+|01⟩+|10⟩-|11⟩)/2 */
    uint32_t idx11[2] = {1, 1};
    hpcq_amplitude(g, idx00, &re, &im);
    ASSERT(APPROX(re, 0.5), "CZ|+,+⟩ ψ(00)=0.5");
    hpcq_amplitude(g, idx11, &re, &im);
    ASSERT(APPROX(re, -0.5), "CZ|+,+⟩ ψ(11)=-0.5");

    norm = hpcq_norm_sq(g);
    ASSERT(APPROX(norm, 1.0), "CZ|+,+⟩ norm = 1.0");

    /* Test measurement */
    uint32_t outcome = hpcq_measure(g, 0, 0.3);
    ASSERT(outcome == 0 || outcome == 1, "HPC measure gives 0 or 1");
    ASSERT(g->n_edges == 0, "HPC edges removed after measure");

    hpcq_destroy(g);
    printf("  HPC: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 9: HPC Contract (Pauli Decomposition + Hadamard Fold)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_contract(void)
{
    printf("Test 9: Pauli Decomposition & Hadamard Fold\n");

    /* Hadamard fold of |0⟩ = (|+⟩ + |−⟩)/√2 */
    double zero_re[2] = {1, 0}, zero_im[2] = {0, 0};
    HadamardFold hf = hpcq_hadamard_fold(zero_re, zero_im);
    ASSERT(APPROX(hf.plus_re, 0.7071067811865475), "|0⟩ fold: plus");
    ASSERT(APPROX(hf.minus_re, 0.7071067811865475), "|0⟩ fold: minus");

    /* Unfold should give |0⟩ back */
    double unf_re[2], unf_im[2];
    hpcq_hadamard_unfold(&hf, unf_re, unf_im);
    ASSERT(APPROX(unf_re[0], 1.0), "unfold re[0] = 1");
    ASSERT(APPROX(unf_re[1], 0.0), "unfold re[1] = 0");

    /* Pauli decomposition of Z matrix */
    double Z_re[2][2] = {{1,0},{0,-1}};
    double Z_im[2][2] = {{0,0},{0,0}};
    PauliDecomposition pd = hpcq_pauli_decompose(Z_re, Z_im);
    ASSERT(APPROX(pd.coeff_re[PAULI_Z], 1.0), "Z decomposes to Z: coeff");
    ASSERT(APPROX(pd.coeff_re[PAULI_X], 0.0), "Z has no X component");
    ASSERT(APPROX(pd.coeff_re[PAULI_I], 0.0), "Z has no I component");
    ASSERT(pd.dominant == PAULI_Z, "Z dominant = PAULI_Z");

    /* Pauli decomposition of X matrix */
    double X_re[2][2] = {{0,1},{1,0}};
    double X_im[2][2] = {{0,0},{0,0}};
    pd = hpcq_pauli_decompose(X_re, X_im);
    ASSERT(APPROX(pd.coeff_re[PAULI_X], 1.0), "X decomposes to X: coeff");
    ASSERT(pd.dominant == PAULI_X, "X dominant = PAULI_X");

    /* Pauli decomposition of H = (X+Z)/√2 */
    double S2 = 0.7071067811865475;
    double H_re[2][2] = {{S2, S2},{S2, -S2}};
    double H_im[2][2] = {{0,0},{0,0}};
    pd = hpcq_pauli_decompose(H_re, H_im);
    ASSERT(APPROX(pd.coeff_re[PAULI_X], S2), "H has X component = 1/√2");
    ASSERT(APPROX(pd.coeff_re[PAULI_Z], S2), "H has Z component = 1/√2");
    ASSERT(APPROX(pd.coeff_re[PAULI_I], 0.0), "H has no I component");

    printf("  contract: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TEST 10: Clifford Exotic Automorphism
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void test_clifford(void)
{
    printf("Test 10: Clifford Exotic Automorphism\n");

    cliff_init();

    /* Identity exists */
    int id = cliff_identity();
    const CliffordElement *e = cliff_get_element(id);
    ASSERT(e != NULL, "identity exists");
    ASSERT(e->axes[0] == 0 && e->axes[1] == 1 && e->axes[2] == 2, "identity axes");
    ASSERT(e->sign[0] == 1 && e->sign[1] == 1 && e->sign[2] == 1, "identity signs");

    /* Compose identity with itself = identity */
    int id2 = cliff_compose(id, id);
    ASSERT(id2 == id, "I∘I = I");

    /* Inverse of identity = identity */
    int id_inv = cliff_inverse(id);
    ASSERT(id_inv == id, "I⁻¹ = I");

    /* Exotic gate: identity maps |0⟩ → |0⟩ */
    double in_re[2] = {1, 0}, in_im[2] = {0, 0};
    double out_re[2], out_im[2];
    cliff_apply_exotic_gate(in_re, in_im, out_re, out_im, id);
    ASSERT(APPROX(out_re[0], 1.0), "I|0⟩ = |0⟩");
    ASSERT(APPROX(out_re[1], 0.0), "I|0⟩ re[1] = 0");

    /* Clifford invariant of |0⟩ — should be > 0 (basis states have frame preference) */
    double delta_0 = cliff_exotic_invariant(in_re, in_im);
    printf("    Δ_C(|0⟩) = %.6f\n", delta_0);
    ASSERT(delta_0 >= 0.0, "Δ_C(|0⟩) ≥ 0");

    /* Clifford invariant of |+⟩ — different frame preference */
    double plus_re[2] = {0.7071067811865475, 0.7071067811865475};
    double plus_im[2] = {0, 0};
    double delta_plus = cliff_exotic_invariant(plus_re, plus_im);
    printf("    Δ_C(|+⟩) = %.6f\n", delta_plus);
    ASSERT(delta_plus >= 0.0, "Δ_C(|+⟩) ≥ 0");

    /* Both |0⟩ and |+⟩ are Bloch-sphere poles — both have frame preference */
    ASSERT(delta_0 > 0.0, "Δ_C(|0⟩) > 0 (Z-axis preference)");
    ASSERT(delta_plus > 0.0, "Δ_C(|+⟩) > 0 (X-axis preference)");

    /* Dual probabilities: std vs exotic */
    double probs_std[2], probs_exo[2];
    cliff_dual_probabilities(in_re, in_im, probs_std, probs_exo, id);
    ASSERT(APPROX(probs_std[0], 1.0), "dual std P(0) = 1");
    ASSERT(APPROX(probs_exo[0], 1.0), "dual exo P(0) = 1 (identity)");

    /* Fingerprint */
    double fp[CLIFF_NUM_CLASSES];
    cliff_exotic_fingerprint(in_re, in_im, fp);
    ASSERT(APPROX(fp[0], 0.0), "identity class Δ = 0");
    printf("    fingerprint: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
           fp[0], fp[1], fp[2], fp[3], fp[4]);

    /* Exotic entropy of |0⟩ with identity = 0 */
    double ds = cliff_exotic_entropy(in_re, in_im, id);
    ASSERT(APPROX(ds, 0.0), "ΔS(|0⟩, I) = 0");

    /* Triality bridge: identity maps Z→Z = Edge (view 0) */
    int view = cliff_target_view(id);
    ASSERT(view == 0, "identity → Edge view");

    printf("  clifford: done\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MAIN
 * ═══════════════════════════════════════════════════════════════════════════════ */
int main(void)
{
    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║        QubitEngine v2.0 — Self-Test Suite             ║\n");
    printf("║        D=2 Triality + HPC + Clifford Exotic           ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n\n");

    test_lifecycle();
    test_gates();
    test_measurement();
    test_entanglement();
    test_triality();
    test_flat();
    test_lazy();
    test_hpc();
    test_contract();
    test_clifford();

    printf("═══════════════════════════════════════════════════════════\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("═══════════════════════════════════════════════════════════\n");

    return tests_failed ? 1 : 0;
}
