/*
 * template.c — BitState Engine API Usage Template
 *
 * Build:
 *   gcc -O2 -march=native -o template template.c \
 *       qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
 *       qubit_register.c qubit_triality.c clifford_exotic.c -lm
 *
 * Run:
 *   ./template
 *
 * This file demonstrates every public API in the BitState engine.
 * BitState is a graph-native quantum simulator that never materializes
 * the state vector.  Entanglement is tracked as weighted edges in an
 * HPC (Holographic Phase Contraction) graph, giving O(N+E) amplitude
 * evaluation for Clifford-dominated circuits.
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
#include "bigint.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 1: Engine Lifecycle & Qubit Allocation
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_lifecycle(void)
{
    printf("═══ 1. Engine Lifecycle & Qubit Allocation ═══\n");

    /* Allocate engine on the heap or stack */
    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);
    /* eng->rng_state is seeded; eng->num_qubits/pairs/registers = 0 */

    /* Allocate |0⟩ — returns uint32_t id (0, 1, 2, …) */
    uint32_t q0 = qubit_init(eng);           /* |0⟩ */
    uint32_t q1 = qubit_init_plus(eng);      /* |+⟩ = (|0⟩+|1⟩)/√2 */
    uint32_t q2 = qubit_init_basis(eng, 1);  /* |1⟩ */

    /* Access raw amplitudes via eng->qubits[id].state  */
    printf("  q0 (|0⟩): re[0]=%f re[1]=%f\n",
           eng->qubits[q0].state.re[0], eng->qubits[q0].state.re[1]);
    printf("  q1 (|+⟩): re[0]=%f re[1]=%f\n",
           eng->qubits[q1].state.re[0], eng->qubits[q1].state.re[1]);
    printf("  q2 (|1⟩): re[0]=%f re[1]=%f\n",
           eng->qubits[q2].state.re[0], eng->qubits[q2].state.re[1]);

    /* PRNG: double in [0,1) */
    double r = qubit_prng_double(eng);

    qubit_engine_destroy(eng);
    free(eng);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 2: Single-Qubit Gates
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_gates(void)
{
    printf("═══ 2. Single-Qubit Gates ═══\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    uint32_t q = qubit_init(eng);  /* |0⟩ */

    qubit_apply_hadamard(eng, q);  /* → |+⟩ */
    qubit_apply_x(eng, q);         /* X|+⟩ = |+⟩ (swap components) */
    qubit_apply_y(eng, q);         /* Y gate */
    qubit_apply_z(eng, q);         /* Z gate: flips phase of |1⟩ */
    qubit_apply_s(eng, q);         /* S = √Z = diag(1, i) */
    qubit_apply_t(eng, q);         /* T = √S = diag(1, e^{iπ/4}) */
    qubit_apply_phase(eng, q, 0.3);/* Arbitrary phase: diag(1, e^{iθ}) on |1⟩ */

    /* Arbitrary 2×2 unitary: flat array [re00, re01, re10, re11] */
    double U_re[4] = {1, 0, 0, 1};
    double U_im[4] = {0, 0, 0, 0};
    qubit_apply_unitary(eng, q, U_re, U_im);

    printf("  Gates applied. Final state amplitudes: re[0]=%f re[1]=%f\n",
           eng->qubits[q].state.re[0], eng->qubits[q].state.re[1]);

    qubit_engine_destroy(eng);
    free(eng);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 3: Two-Qubit Gates (CZ, CX)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_2q_gates(void)
{
    printf("═══ 3. Two-Qubit Gates ═══\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    uint32_t qa = qubit_init(eng);  /* |0⟩ */
    uint32_t qb = qubit_init(eng);  /* |0⟩ */

    /* CZ: controlled-Z — adds a phase of -1 to |11⟩ */
    qubit_apply_cz(eng, qa, qb);

    /* CX: controlled-X (CNOT) — auto-creates product pair internally */
    uint32_t qc = qubit_init(eng);
    uint32_t qd = qubit_init(eng);
    qubit_apply_cx(eng, qc, qd);    /* control=qc, target=qd */

    /* Two-qubit unitary on an entangled pair */
    uint32_t qe = qubit_init(eng);
    uint32_t qf = qubit_init(eng);
    qubit_entangle_product(eng, qe, qf);
    double U2_re[16] = {0};
    double U2_im[16] = {0};
    U2_re[0] = 1.0;
    U2_re[5] = 1.0;
    U2_re[10] = 1.0;
    U2_re[15] = 1.0;
    qubit_apply_unitary_pair(eng, qe, qf, U2_re, U2_im);

    qubit_engine_destroy(eng);
    free(eng);
    printf("  2-qubit gates applied.\n\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 4: Measurement & Inspection
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_measurement(void)
{
    printf("═══ 4. Measurement & Inspection ═══\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    uint32_t q = qubit_init_plus(eng);  /* |+⟩ */

    /* Inspect without collapsing */
    QubitInspect info = qubit_inspect(eng, q);
    printf("  Inspect: P(0)=%.4f P(1)=%.4f entropy=%.4f purity=%.4f rank=%d\n",
           info.prob[0], info.prob[1], info.entropy, info.purity,
           info.schmidt_rank);

    /* Measure: collapses to |0⟩ or |1⟩, returns outcome */
    int outcome = qubit_measure(eng, q);
    printf("  Measure outcome: %d\n", outcome);

    qubit_engine_destroy(eng);
    free(eng);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 5: Entanglement (Bell Pairs, Product States)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_entanglement(void)
{
    printf("═══ 5. Entanglement ═══\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    uint32_t qa = qubit_init(eng);
    uint32_t qb = qubit_init(eng);

    /* Create a Bell pair (|00⟩+|11⟩)/√2 */
    int pair_id = qubit_entangle_bell(eng, qa, qb);
    printf("  Bell pair: id=%d\n", pair_id);

    /* Measure one qubit — the other collapses instantly (correlated) */
    int outcome_a = qubit_measure_in_pair(eng, qa);
    printf("  Measured qa=%d (partner qb collapses too)\n", outcome_a);

    /* Product state: |ψ_a⟩ ⊗ |ψ_b⟩ (initially separable) */
    uint32_t qc = qubit_init_plus(eng);
    uint32_t qd = qubit_init_basis(eng, 1);
    int prod_id = qubit_entangle_product(eng, qc, qd);
    printf("  Product pair: id=%d\n", prod_id);

    /* Disentangle: extract marginals back into local states */
    qubit_disentangle(eng, qc, qd);

    qubit_engine_destroy(eng);
    free(eng);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 6: Qubit Register (Sparse Multi-Qubit)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_register(void)
{
    printf("═══ 6. Qubit Register ═══\n");

    QubitEngine *eng = (QubitEngine *)calloc(1, sizeof(QubitEngine));
    qubit_engine_init(eng);

    /* Create a register of 10 qubits */
    int reg = qubit_reg_init(eng, 0, 10, 2);

    /* GHZ entanglement: (|00…0⟩ + |11…1⟩)/√2 */
    qubit_reg_entangle_all(eng, reg);

    /* Apply Hadamard to first qubit in register */
    qubit_reg_apply_hadamard(eng, reg, 0);

    /* Apply CZ between qubits 2 and 3 */
    qubit_reg_apply_cz(eng, reg, 2, 3);

    /* Arbitrary unitary at position 5 */
    double U_re[4] = {1, 0, 0, 1};
    double U_im[4] = {0, 0, 0, 0};
    qubit_reg_apply_unitary_pos(eng, reg, 5, U_re, U_im);

    /* Get amplitude of |1010101010⟩ */
    SV_Amplitude amp = qubit_reg_sv_get(eng, reg, 0x2AA);
    printf("  SV amplitude: re=%f im=%f\n", amp.re, amp.im);

    /* Set amplitude manually */
    qubit_reg_sv_set(eng, reg, 0x000, 0.7, 0.0);

    /* Inner product between two registers */
    double prob = qubit_reg_sv_total_prob(eng, reg);
    printf("  Total prob: %f\n", prob);

    qubit_engine_destroy(eng);
    free(eng);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 7: Triality (3-View Qubit Representation)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_triality(void)
{
    printf("═══ 7. Triality (Z/X/Y Views) ═══\n");

    /* Triality Qubit — standalone (not part of engine) */
    TrialityQubit tq;
    tri_init(&tq);               /* |0⟩ in Edge (Z) view */

    /* Gate application routes to optimal view automatically */
    tri_apply_hadamard(&tq);     /* → |+⟩, swaps Edge↔Vertex */
    tri_apply_x(&tq);            /* O(1) in Vertex view */
    tri_apply_z(&tq, 0.5);      /* O(1) in Edge view — auto-converts */

    /* Force a specific view */
    tri_ensure_view(&tq, VIEW_DIAGONAL);

    /* Get amplitudes in any view */
    double re[2], im[2];
    tri_get_amplitudes(&tq, VIEW_EDGE, re, im);
    printf("  Edge view: re[0]=%f re[1]=%f\n", re[0], re[1]);

    /* Triality rotation: Edge→Vertex→Diagonal→Edge (O(1) relabeling) */
    tri_rotate(&tq);

    /* Measure in a specific view */
    int outcome = tri_measure(&tq, VIEW_EDGE, 0.5);
    printf("  Measure outcome: %d\n", outcome);

    /* Lazy Triality — accumulates phase, defers materialization */
    LazyTrialityQubit lq;
    lazy_tri_init(&lq);
    lazy_tri_apply_phase(&lq, 1.23);
    lazy_tri_apply_hadamard(&lq);
    lazy_tri_materialize(&lq);   /* compute amplitudes on demand */
    int lo = lazy_tri_measure(&lq, 0.3);
    printf("  Lazy measure: %d\n", lo);

    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 8: Flat Qubit (Auto Basis↔Full Mode)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_flat_qubit(void)
{
    printf("═══ 8. Flat Qubit (Optimized Representation) ═══\n");

    /* Starts in FLAT_BASIS mode — all gates O(1) */
    FlatQubit fq;
    flat_init(&fq);              /* |0⟩ in basis mode */

    flat_apply_x(&fq);           /* O(1): just flip basis_k */
    flat_apply_z(&fq);           /* O(1): just flip phase on |1⟩ */
    flat_apply_phase(&fq, 1.0);  /* O(1): multiply phase on |1⟩ */

    /* Hadamard auto-promotes to QUANTUM_FULL (creates superposition) */
    flat_apply_hadamard(&fq);
    printf("  Mode after H: %s\n",
           fq.mode == QUANTUM_FULL ? "QUANTUM_FULL" : "FLAT_BASIS");

    /* Try to demote back to basis if state collapsed */
    int demoted = flat_try_demote(&fq);
    printf("  Demoted: %s\n", demoted ? "yes" : "no");

    /* Read out amplitudes regardless of mode */
    double re[2], im[2];
    flat_get_amplitudes(&fq, re, im);
    printf("  Amplitudes: re[0]=%f re[1]=%f\n", re[0], re[1]);

    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 9: HPC Phase Graph (Graph-Native Simulation)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_hpc_graph(void)
{
    printf("═══ 9. HPC Phase Graph ═══\n");

    /* Create graph with 4 sites */
    HPCQGraph *g = hpcq_create(4);

    /* Set custom local states */
    double init_re[2] = {1.0, 0.0};
    double init_im[2] = {0.0, 0.0};
    hpcq_set_local(g, 0, init_re, init_im);
    hpcq_hadamard(g, 0);              /* site 0 → |+⟩ */
    hpcq_phase(g, 1, 3.14159);        /* site 1: Z(π) = -|1⟩⟨1| */
    hpcq_x(g, 2);                     /* site 2: |0⟩ → |1⟩ */

    /* CZ edge: exact, fidelity=1.0. CZ²=I — same edge cancels. */
    hpcq_cz(g, 0, 1);
    hpcq_cz(g, 1, 2);

    /* General 2-site gate (4×4 unitary → phase edge) */
    double G_re[16] = {0}, G_im[16] = {0};
    G_re[0] = 1; G_re[5] = 1; G_re[10] = 1; G_re[15] = 1;
    hpcq_general_2site(g, 2, 3, G_re, G_im);

    /* Amplitude evaluation — O(N+E), never materializes the vector */
    uint32_t config[4] = {0, 0, 0, 0};
    double amp_re, amp_im;
    hpcq_amplitude(g, config, &amp_re, &amp_im);
    printf("  ψ(0000) = %f + %fi\n", amp_re, amp_im);

    /* Probability */
    double prob = hpcq_probability(g, config);
    printf("  P(0000) = %f\n", prob);

    /* Marginal probability for a single site */
    double p0 = hpcq_marginal(g, 0, 0);
    printf("  P(site_0=0) = %f\n", p0);

    /* Born measurement — collapses site, absorbs phases into partners */
    uint32_t outcome = hpcq_measure(g, 0, 0.7);
    printf("  Measure site 0 → %u\n", outcome);

    /* Norm check (small N only) */
    double norm = hpcq_norm_sq(g);
    printf("  Norm² = %f\n", norm);

    /* Entropy across a cut */
    double entropy = hpcq_entropy_cut(g, 1);
    printf("  Entropy across cut after site 1: %f\n", entropy);

    /* Print statistics */
    hpcq_print_stats(g);

    hpcq_destroy(g);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 10: HPC Contract (Pauli Decomposition, Clifford Encoding)
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_hpc_contract(void)
{
    printf("═══ 10. HPC Pauli Decomposition & Clifford Encoding ═══\n");

    HPCQGraph *g = hpcq_create(3);

    /* Hadamard fold — the D=2 "vesica fold" (sum/diff channels) */
    double re[2] = {1.0, 0.5};
    double im[2] = {0.0, 0.0};
    HadamardFold hf = hpcq_hadamard_fold(re, im);
    printf("  Hadamard fold: plus=%f%+fi minus=%f%+fi\n",
           hf.plus_re, hf.plus_im, hf.minus_re, hf.minus_im);

    /* Unfold back */
    double re2[2], im2[2];
    hpcq_hadamard_unfold(&hf, re2, im2);

    /* Pauli decomposition of a 2×2 interaction */
    double M_re[2][2] = {{1, 0}, {0, -1}};  /* Pauli Z */
    double M_im[2][2] = {{0, 0}, {0, 0}};
    PauliDecomposition pd = hpcq_pauli_decompose(M_re, M_im);
    printf("  Pauli decomposition: I=%f Z=%f X=%f Y=%f → dominant=%d\n",
           pd.energy[0], pd.energy[1], pd.energy[2], pd.energy[3],
           pd.dominant);

    /* Encode a gate as a Clifford edge (auto-selects best Pauli) */
    hpcq_encode_clifford(g, 0, 1, M_re, M_im);

    /* Auto-encode: CZ→exact, high-fidelity→Clifford, else→general */
    double gate_re[16] = {0}, gate_im[16] = {0};
    gate_re[0] = 1; gate_re[5] = 1; gate_re[10] = 1; gate_re[15] = 1;
    hpcq_encode_2site(g, 1, 2, gate_re, gate_im);

    /* Fold analysis */
    HPCQFoldAnalysis fa = hpcq_analyze_fold(g, 0);
    printf("  Fold analysis: |+⟩=%.3f |−⟩=%.3f entropy=%.3f\n",
           fa.plus_fidelity, fa.minus_fidelity, fa.fold_entropy);

    hpcq_destroy(g);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 11: HPC Sparse Amplitude Reconstruction
 * ═══════════════════════════════════════════════════════════════════════════════ */
/* Sample observable: count number of 1s */
static double demo_obs_count_ones(const uint32_t *indices, uint64_t n_sites, void *ctx)
{
    (void)ctx;
    double sum = 0;
    for (uint64_t i = 0; i < n_sites; i++) sum += indices[i];
    return sum;
}

static void demo_hpc_sparse(void)
{
    printf("═══ 11. HPC Sparse State Reconstruction ═══\n");

    HPCQGraph *g = hpcq_create(4);
    hpcq_hadamard(g, 0);
    hpcq_cz(g, 0, 1);
    hpcq_cz(g, 1, 2);

    /* Brute-force sparse enumeration (N ≤ 20) */
    HPCQSparseVector *sv = hpcq_sparse_brute(g, 0.01, 100);
    if (sv) {
        hpcq_sv_print(sv, 5);
        hpcq_sv_destroy(sv);
    }

    /* Tree-pruned sparse reconstruction (for larger N) */
    sv = hpcq_sparse_tree(g, 0.01, 200);
    if (sv) {
        printf("  Tree-pruned entries: %lu, total prob: %f\n",
               sv->count, sv->total_prob);
        hpcq_sv_destroy(sv);
    }

    /* Monte Carlo expectation value */
    double avg = hpcq_expectation(g, demo_obs_count_ones, NULL, 1000, 42);
    printf("  MC expectation ⟨#ones⟩: %f\n", avg);

    hpcq_destroy(g);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 12: Clifford Exotic Automorphism
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_clifford_exotic(void)
{
    printf("═══ 12. Clifford Exotic Automorphism ═══\n");

    /* Generate the 24 Clifford elements (call once) */
    cliff_init();

    int idx_id = cliff_identity();
    printf("  Identity index: %d\n", idx_id);

    /* Composition and inverse */
    int a = 3, b = 7;
    int c = cliff_compose(a, b);
    int inv_a = cliff_inverse(a);
    printf("  %d ∘ %d = %d,  %d⁻¹ = %d\n", a, b, c, a, inv_a);

    /* Get element and gate matrix */
    const CliffordElement *elem = cliff_get_element(0);
    const CliffordGate *gate = cliff_get_gate(0);
    printf("  Element 0: axes=(%d,%d,%d) signs=(%d,%d,%d)\n",
           elem->axes[0], elem->axes[1], elem->axes[2],
           elem->sign[0], elem->sign[1], elem->sign[2]);

    /* Exotic gate: C·G·C† — gate in Clifford-conjugated basis */
    double in_re[2] = {1.0, 0.0}, in_im[2] = {0.0, 0.0};
    double out_re[2], out_im[2];
    cliff_apply_exotic_gate(in_re, in_im, out_re, out_im, 5);

    /* Dual probabilities */
    double std_p[2], exo_p[2];
    cliff_dual_probabilities(in_re, in_im, std_p, exo_p, 5);
    printf("  Dual probs: std=(%.3f,%.3f) exo=(%.3f,%.3f)\n",
           std_p[0], std_p[1], exo_p[0], exo_p[1]);

    /* Clifford invariant Δ_C — measures Pauli-frame dependence */
    double delta = cliff_exotic_invariant(in_re, in_im);
    printf("  Δ_C = %f\n", delta);

    /* Exotic entropy */
    double ds = cliff_exotic_entropy(in_re, in_im, 5);
    printf("  ΔS = %f\n", ds);

    /* Fingerprint: per-conjugacy-class breakdown */
    double fingerprint[CLIFF_NUM_CLASSES];
    cliff_exotic_fingerprint(in_re, in_im, fingerprint);
    printf("  Fingerprint: [");
    for (int i = 0; i < CLIFF_NUM_CLASSES; i++)
        printf("%s%.3f", i ? ", " : "", fingerprint[i]);
    printf("]\n");

    /* Triality bridge: which view does this Clifford map Z-basis to? */
    int view = cliff_target_view(5);
    printf("  Clifford 5 maps Z-basis to view %d (%s)\n",
           view, TRIALITY_VIEW_NAMES[view]);

    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * SECTION 13: Header-Only Primitives
 * ═══════════════════════════════════════════════════════════════════════════════ */
static void demo_primitives(void)
{
    printf("═══ 13. Header-Only Primitives ═══\n");

    /* State vector primitives (statevector.h) */
    QubitState qs = { .re = {1.0, 0.0}, .im = {0.0, 0.0} };
    SV_Amplitude a = sv_get_local(&qs, 0);
    printf("  sv_get_local: re=%f im=%f\n", a.re, a.im);

    JointState js;
    memset(&js, 0, sizeof(js));
    js.re[0] = 0.707; js.re[3] = 0.707;
    a = sv_get_joint(&js, 0, 1);

    /* Born rule (born_rule.h) */
    double p = born_prob(qs.re[0], qs.im[0]);
    printf("  P(0) = %f\n", p);

    int sample = born_sample(qs.re, qs.im, 2, 0.3);
    born_collapse(qs.re, qs.im, 2, sample);

    /* Superposition / Hadamard (superposition.h) */
    double h_re[2] = {1.0, 0.0}, h_im[2] = {0.0, 0.0};
    sup_apply_hadamard(h_re, h_im);
    sup_renormalize(h_re, h_im, 2);

    /* Entanglement analysis (entanglement.h) */
    ent_bell(&js);
    int rank = ent_schmidt_rank(&js);
    double entropy = ent_entropy(&js);
    double marg[2];
    ent_marginal_a(&js, marg);
    printf("  Bell: rank=%d entropy=%f P_a(0)=%f\n", rank, entropy, marg[0]);
    ent_renormalize(&js);

    /* Qubit management (qubit_management.h) */
    QubitState sa, sb;
    qm_init_zero(&sa);
    qm_init_plus(&sb);
    qm_renormalize(&sa);
    qm_copy(&sa, &sb);

    /* BigInt (bigint.h) — 4096-bit arbitrary precision */
    BigInt bi;
    bigint_set_u64(&bi, 42);
    char buf[1240];
    bigint_to_decimal(buf, sizeof(buf), &bi);
    printf("  BigInt: %s\n", buf);

    /* Fast inverse sqrt / reciprocal (arithmetic.h) */
    double isqrt = arith_fast_isqrt(2.0);
    printf("  Fast 1/√2 ≈ %f (exact: %f)\n", isqrt, 1.0/sqrt(2.0));

    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MAIN — Run all demos
 * ═══════════════════════════════════════════════════════════════════════════════ */
int main(void)
{
    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║      BitState Engine — API Usage Template             ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n\n");

    demo_lifecycle();
    demo_gates();
    demo_2q_gates();
    demo_measurement();
    demo_entanglement();
    demo_register();
    demo_triality();
    demo_flat_qubit();
    demo_hpc_graph();
    demo_hpc_contract();
    demo_hpc_sparse();
    demo_clifford_exotic();
    demo_primitives();

    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║  All API sections demonstrated.                       ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n");
    return 0;
}
