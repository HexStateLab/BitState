# BitState API Reference — D=2 Holographic Phase Graph Engine

**Version:** D=2 (qubit) adaptation of HexState D=6 architecture  
**Memory Model:** O(N + E) — never materializes the full 2^N state vector  
**Dimension:** 2 (each site is a qubit: |0⟩, |1⟩)

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Data Structures](#2-data-structures)
    - 2.1 Core State Types (statevector.h)
    - 2.2 TrialityQubit (qubit_triality.h)
    - 2.3 FlatQubit (flat_qubit.h)
    - 2.4 HPC Graph Types (hpc_qubit_graph.h)
    - 2.5 Sparse Vector Types (hpc_qubit_amplitude.h)
    - 2.6 QubitEngine Types (qubit_engine.h)
3. [HPC Graph API](#3-hpc-graph-api)
    - 3.1 Lifecycle
    - 3.2 Local Gates
    - 3.3 Two-Qubit Gates
    - 3.4 Amplitude & Probability
    - 3.5 Marginal & Measurement
    - 3.6 Norm & Entropy
    - 3.7 Diagnostics
4. [HPC Sparse Vector & Sampling](#4-hpc-sparse-vector--sampling)
    - 4.1 Sparse Vector Lifecycle
    - 4.2 Reconstruction
    - 4.3 Monte Carlo Expectation
5. [HPC Bond Contraction](#5-hpc-bond-contraction)
    - 5.1 Hadamard Fold / Unfold
    - 5.2 Pauli Decomposition
    - 5.3 Clifford Edge Encoding
    - 5.4 2-Site Gate Encoder
6. [QubitEngine API](#6-qubitengine-api)
    - 6.1 Lifecycle & Initialization
    - 6.2 Gates
    - 6.3 Measurement & Inspection
    - 6.4 Entanglement
    - 6.5 Registers
7. [Qubit Management Primitives](#7-qubit-management-primitives)
8. [Entanglement Utilities](#8-entanglement-utilities)
9. [Superposition (DFT₂)](#9-superposition-dft₂)
10. [Born Rule](#10-born-rule)
11. [Floating-Point Primitives](#11-floating-point-primitives)
12. [BigInt Library](#12-bigint-library)
13. [Example Programs](#13-example-programs)

---

## 1. Architecture Overview

BitState provides **two parallel simulation engines**:

| Engine | File | Approach | Strengths |
|--------|------|----------|-----------|
| **HPC Graph** | `hpc_qubit_graph.h` | Phase-graph: state is a product of local amplitudes × edge phases. No state vector. | O(N+E) amplitude queries; scales to thousands of qubits; absorb mechanism handles H gates on entangled qubits |
| **QubitEngine** | `qubit_engine.h` + `.c` | Pairwise entanglement: each entangled pair is a 4-amplitude JointState (64 bytes). Tracks local qubits + pair bonds. | Full state-vector access on registers; GHZ bulk mode; standard gate set |

Both engines share lower-level modules:

| Module | File | Purpose |
|--------|------|---------|
| Triality | `qubit_triality.h/.c` | Three-basis representation (Z/X/Y) with lazy conversion |
| Flat qubit | `flat_qubit.h` | Two-tier optim. (basis vs full) for single qubits |
| State vector | `statevector.h` | QubitState (32 B), JointState (64 B) structs |
| Superposition | `superposition.h` | DFT₂ (Hadamard) on 2 amplitudes |
| Born rule | `born_rule.h` | Probability & sampling utilities |
| Entanglement | `entanglement.h` | JointState operations: Bell states, marginals, Schmidt rank |
| Contract | `hpc_qubit_contract.h` | Pauli decomposition, Clifford edge encoding |
| Amplitude | `hpc_qubit_amplitude.h` | Sparse state-vector recon and Monte Carlo |
| Arithmetic | `arithmetic.h` | Fast float ops (inverse sqrt, reciprocal) |
| BigInt | `bigint.h/.c` | 4096-bit arbitrary-precision arithmetic |

---

## 2. Data Structures

### 2.1 Core State Types — `statevector.h`

```c
/* Single qubit: 2 complex amplitudes = 32 bytes */
typedef struct {
    double re[2];   /* Real parts: re[0]=|0⟩, re[1]=|1⟩ */
    double im[2];   /* Imag parts: im[0]=|0⟩, im[1]=|1⟩ */
} QubitState;

/* Entangled pair: 4 complex amplitudes = 64 bytes */
typedef struct {
    double re[4];   /* Index order: 0=|00⟩, 1=|01⟩, 2=|10⟩, 3=|11⟩ */
    double im[4];
} JointState;

/* Generalized amplitude (for registers on large-N states) */
typedef struct {
    double re, im;
    double log2_mag;   /* For states where magnitude underflows double */
} SV_Amplitude;

/* Stream callback for register state-vector iteration */
typedef void (*sv_stream_fn)(uint64_t basis_state, SV_Amplitude amp, void *user_data);

SV_Amplitude sv_get_local(const QubitState *s, int k);       /* Get |k⟩ amplitude */
SV_Amplitude sv_get_joint(const JointState *j, int a, int b); /* Get |ab⟩ amplitude */
```

### 2.2 TrialityQubit — `qubit_triality.h`

Three mutually-defining views of a single qubit's state:

```c
typedef enum {
    VIEW_EDGE      = 0,   /* Z-basis (computational): |0⟩, |1⟩ */
    VIEW_VERTEX    = 1,   /* X-basis (Hadamard):     |+⟩, |−⟩ */
    VIEW_DIAGONAL  = 2,   /* Y-basis (S†·H):          |i⟩, |−i⟩ */
    VIEW_COUNT     = 3
} TrialityView;

typedef struct {
    double re[3][2], im[3][2];   /* Amplitudes per view (lazy-computed) */
    uint8_t dirty[3];            /* 1 = view needs recomputation */
    TrialityView current_view;   /* Last-written view */
    uint8_t is_eigenstate;       /* 1 = Pauli eigenstate (O(1) convert) */
    uint8_t is_real;             /* 1 = all imag parts zero */
    uint8_t active_mask;         /* Bitmask: bit 0=|0⟩, bit 1=|1⟩ */
    uint32_t conversions;        /* Total view conversions performed */
    uint32_t gates_applied;      /* Total gates applied */
    uint32_t conversion_savings; /* Gates that avoided conversion */
} TrialityQubit;

/* Lazy evaluation variant (chains H and phase gates without materializing) */
typedef struct {
    struct { double phase; } segments[64];
    uint32_t num_segments;
    uint32_t dft_count;     /* H count mod 2 (H²=I) */
    int materialized;
    double re[2], im[2];    /* Materialized amplitudes */
} LazyTrialityQubit;
```

### 2.3 FlatQubit — `flat_qubit.h`

Two-tier optimization: when a qubit is in a known basis state, all gates are O(1).

```c
typedef enum {
    FLAT_BASIS   = 0,   /* Single |k⟩ with phase — cheapest ops */
    QUANTUM_FULL = 1    /* Full 2-amplitude state — general gates */
} FlatMode;

typedef struct {
    FlatMode mode;
    uint8_t  basis_k;       /* FLAT_BASIS: which basis state (0/1) */
    double   phase_re;      /* FLAT_BASIS: phase factor (unit magnitude) */
    double   phase_im;
    double   re[2], im[2];  /* QUANTUM_FULL: amplitudes */
    uint32_t promotions;    /* Times promoted from basis to full */
    uint32_t demotions;     /* Times demoted from full to basis */
} FlatQubit;
```

### 2.4 HPC Graph Types — `hpc_qubit_graph.h`

```c
/* ── Edge types ── */
typedef enum {
    HPCQ_EDGE_CZ        = 0,   /* Exact CZ: w(a,b)=(-1)^(a·b), fidelity=1.0 */
    HPCQ_EDGE_PHASE     = 1,   /* General 2×2 phase matrix */
    HPCQ_EDGE_CLIFFORD  = 2,   /* Clifford-projected (from Pauli decomp.) */
    HPCQ_EDGE_ABSORBED  = 3    /* Consumed by multi-edge H absorption */
} HPCQEdgeType;

/* ── A single weighted edge between two sites ── */
typedef struct {
    HPCQEdgeType type;
    uint64_t site_a, site_b;
    double w_re[2][2], w_im[2][2];   /* 2×2 complex phase matrix */
    uint8_t pauli_channel;            /* For Clifford edges: 0=I,1=Z,2=X,3=Y */
    double fidelity;                  /* 1.0 = lossless */
} HPCQEdge;

/* ── Absorb entry: H-gate absorption of incident edges ── */
typedef struct {
    uint64_t center;            /* Site where H was applied */
    uint64_t n_nbrs;            /* Total neighbors in this entry */
    uint64_t n_inner;           /* First n_inner are inner group (first H) */
    uint64_t *nbrs;             /* Neighbor site indices */
    double *w_re, *w_im;        /* Edge matrices [n_nbrs * 4] */
    double a_re[2], a_im[2];    /* Local amplitude at FIRST absorption */
    int n_layers;               /* 1 or 2 (re-absorption count) */
    double a_cur_re[2], a_cur_im[2];  /* pre-H state for second layer */
} HPCQAbsorbEntry;

/* ── The graph: state = graph + local sites ── */
typedef struct {
    uint64_t n_sites;
    TrialityQubit *locals;       /* Per-site local states */
    uint64_t n_edges, edge_cap;
    HPCQEdge *edges;
    uint64_t n_absorb, absorb_cap;
    HPCQAbsorbEntry *absorb;    /* Absorbed H entries */
    uint64_t n_log, log_cap;
    HPCQGateEntry *gate_log;    /* Full gate history */
    uint64_t *inc_counts;        /* Per-site incident edge count */
    uint64_t *inc_cap;
    uint64_t **inc_edges;        /* Per-site incident edge lists */
    int64_t *absorb_idx;         /* Per-site absorb index (-1=none) */
    /* Statistics */
    uint64_t amp_evals, prob_evals, measurements;
    uint64_t cz_edges, phase_edges, clifford_edges;
    double min_fidelity, avg_fidelity;
} HPCQGraph;
```

### 2.5 Sparse Vector Types — `hpc_qubit_amplitude.h`

```c
typedef struct {
    uint32_t *indices;   /* Basis state (0/1 per site) */
    double re, im;       /* Complex amplitude */
    double prob;         /* = re²+im² */
} HPCQSparseEntry;

typedef struct {
    HPCQSparseEntry *entries;
    uint64_t count, capacity;
    uint64_t n_sites;
    double total_prob;   /* Σ prob across all entries */
    double threshold;    /* Minimum prob for inclusion */
} HPCQSparseVector;

/* Observable callback for Monte Carlo */
typedef double (*HPCQObservable)(const uint32_t *indices,
                                  uint64_t n_sites, void *ctx);
```

### 2.6 QubitEngine Types — `qubit_engine.h`

```c
#define MAX_QUBITS     262144   /* 256K qubits */
#define MAX_PAIRS      262144   /* 256K entangled pairs */
#define MAX_REGISTERS  16384    /* 16K registers */
typedef uint64_t basis_t;       /* up to 64 qubits natively */

typedef struct {
    QubitState state;          /* 32 bytes: 2 complex amplitudes */
    int32_t pair_id;           /* Index into pairs[] (-1 if not entangled) */
    int32_t pair_side;         /* 0=A, 1=B */
    uint32_t id;
} Qubit;

typedef struct {
    JointState joint;          /* 64 bytes: 4 complex amplitudes */
    uint32_t id_a, id_b;       /* Two qubit IDs */
    uint8_t active;
} QubitPair;

typedef struct {
    uint64_t chunk_id;
    uint64_t n_qubits;
    uint32_t dim;               /* Always 2 for qubits */
    uint8_t collapsed;
    uint32_t collapse_outcome;
    uint64_t magic_base;
    uint8_t bulk_rule;          /* 0=general, 1=GHZ, 2=circuit */
    RegisterEntry entries[4096];
    uint32_t num_nonzero;
    /* GHZ/circuit mode metadata: */
    uint16_t gauss_n_dft, gauss_n_cz;
    uint8_t gauss_ready;
    uint16_t gauss_cz_a[256], gauss_cz_b[256];
} QubitRegister;

typedef struct {
    Qubit qubits[MAX_QUBITS];
    QubitPair pairs[MAX_PAIRS];
    QubitRegister registers[MAX_REGISTERS];
    uint32_t num_qubits, num_pairs, num_registers;
    uint64_t rng_state;          /* LCG PRNG state */
} QubitEngine;

/* Inspection result: probabilities, entropy, purity, Schmidt rank */
typedef struct {
    double prob[2];
    double entropy;
    double purity;
    int schmidt_rank;
} QubitInspect;
```

---

## 3. HPC Graph API — `hpc_qubit_graph.h`

### 3.1 Lifecycle

```c
/* ── Create a new HPC graph with n_sites (all initialized to |0⟩) ── */
HPCQGraph *hpcq_create(uint64_t n_sites);

/* ── Destroy graph and free all memory ── */
void hpcq_destroy(HPCQGraph *g);
```

**Constants:**
| Macro | Value | Description |
|-------|-------|-------------|
| `HPCQ_D` | `2` | Physical dimension |
| `HPCQ_INIT_EDGES` | `4096` | Initial edge capacity |
| `HPCQ_INIT_LOG` | `8192` | Initial gate log capacity |
| `HPCQ_W2_RE[2]` | `{1.0, -1.0}` | Real parts of ω = e^{2πi/2} |
| `HPCQ_W2_IM[2]` | `{0.0, 0.0}` | Imaginary parts |

### 3.2 Local Gates

```c
/* ── Set local state to arbitrary 2-amplitude vector ── */
void hpcq_set_local(HPCQGraph *g, uint64_t site,
                    const double re[2], const double im[2]);

/* ── Hadamard gate (⚠ use hpcq_hadamard_absorb if site has incident edges) ── */
void hpcq_hadamard(HPCQGraph *g, uint64_t site);

/* ── Hadamard with edge absorption ──
 * Correctly handles H on qubits that have incident CZ/phase edges.
 * Cases:
 *   0 incident edges → standard hpcq_hadamard (fast path)
 *   1 incident edge  → absorb into edge matrix
 *   >1 incident edges → multi-edge absorption (stored in absorb entry)
 *   Re-absorption    → two-layer evaluation (H²=I composition)
 */
void hpcq_hadamard_absorb(HPCQGraph *g, uint64_t site);

/* ── Phase gate Z(θ): |1⟩ → e^{iθ}|1⟩ ── */
void hpcq_phase(HPCQGraph *g, uint64_t site, double theta);

/* ── T gate:  |1⟩ → e^{iπ/4}|1⟩ ── */
void hpcq_t(HPCQGraph *g, uint64_t site);

/* ── T† gate: |1⟩ → e^{-iπ/4}|1⟩ ── */
void hpcq_td(HPCQGraph *g, uint64_t site);

/* ── Pauli-X (NOT): |0⟩↔|1⟩ ── */
void hpcq_x(HPCQGraph *g, uint64_t site);
```

### 3.3 Two-Qubit Gates

```c
/* ── CZ gate: w(a,b) = (-1)^(a·b) ──
 * Key property: CZ² = I, so applying CZ twice cancels.
 * If a CZ edge already exists between (a,b), it is removed (swap-removed).
 */
void hpcq_cz(HPCQGraph *g, uint64_t site_a, uint64_t site_b);

/* ── General 2-site phase gate ──
 * Encodes a 4×4 unitary as diagonal phase edge w(j,k).
 * Only the diagonal phases G[(j,k),(j,k)] are extracted.
 * The normalization divides by |G[(j,k),(j,k)]|.
 *
 * @param G_re, G_im: 16-element arrays (row-major 4×4 complex matrix).
 */
void hpcq_general_2site(HPCQGraph *g, uint64_t site_a, uint64_t site_b,
                        const double *G_re, const double *G_im);

/* Convenience: CNOT = H(target)·CZ(control,target)·H(target)
 *   hpcq_hadamard_absorb(g, t);
 *   hpcq_cz(g, c, t);
 *   hpcq_hadamard_absorb(g, t);
 */
```

### 3.4 Amplitude & Probability

```c
/* ── Compute ψ(i₁,...,iₙ) = [∏ₖ aₖ(iₖ)] × [∏_{edges} wₑ(iₐ,i_b)] ──
 * O(N + E). Step 1: product of local amplitudes. Step 2: edge phases.
 * Step 3: absorb correction (evaluates absorb entries).
 *
 * @param indices  Array of n_sites uint32 values (0 or 1).
 * @param out_re   Output real amplitude.
 * @param out_im   Output imaginary amplitude.
 */
void hpcq_amplitude(const HPCQGraph *g, const uint32_t *indices,
                    double *out_re, double *out_im);

/* ── Probability = |ψ(indices)|² ── */
double hpcq_probability(const HPCQGraph *g, const uint32_t *indices);
```

### 3.5 Marginal & Measurement

```c
/* ── Marginal probability P(site = value) ──
 * Enumerates 2^n_connected configurations of direct neighbors.
 * Only sums over sites with DIRECT edges to `site` (not the full
 * connected component). Does NOT evaluate absorb entries — for
 * correct results on absorbed graphs, use hpcq_sparse_tree instead.
 *
 * @param site   Site index.
 * @param value  0 or 1.
 * @returns      Marginal probability P(site = value).
 */
double hpcq_marginal(const HPCQGraph *g, uint64_t site, uint32_t value);

/* ── Measure site ──
 * 1. Compute P(0), P(1) via hpcq_marginal.
 * 2. Sample outcome from random_01.
 * 3. Collapse local state to |outcome⟩.
 * 4. Absorb phase from incident edges into partners.
 * 5. Remove resolved edges.
 *
 * @param site       Site to measure.
 * @param random_01  Random double in [0,1] for Born-rule sampling.
 * @returns          Outcome (0 or 1).
 */
uint32_t hpcq_measure(HPCQGraph *g, uint64_t site, double random_01);
```

### 3.6 Norm & Entropy

```c
/* ── Total norm Σ|ψ|² (brute-force, N ≤ 20 only) ── */
double hpcq_norm_sq(const HPCQGraph *g);

/* ── Entropy estimate across bipartition cut ──
 * CZ edges contribute exactly 1 bit per crossing edge.
 * General edges contribute fidelity-weighted 1 bit.
 *
 * @param cut_after  Split point: sites 0..cut_after vs cut_after+1..N-1.
 */
double hpcq_entropy_cut(const HPCQGraph *g, uint64_t cut_after);
```

### 3.7 Diagnostics

```c
/* ── Print formatted statistics ──
 * Site count, edge counts by type, amp/prob/measure eval counts,
 * fidelity stats, memory usage, estimated full state-vector size.
 */
void hpcq_print_stats(const HPCQGraph *g);
```

---

## 4. HPC Sparse Vector & Sampling — `hpc_qubit_amplitude.h`

### 4.1 Sparse Vector Lifecycle

```c
/* ── Create sparse vector ── */
HPCQSparseVector *hpcq_sv_create(uint64_t n_sites, uint64_t initial_cap);

/* ── Destroy sparse vector ── */
void hpcq_sv_destroy(HPCQSparseVector *sv);

/* ── Add entry (grows if needed) ── */
void hpcq_sv_add(HPCQSparseVector *sv, const uint32_t *indices,
                 double re, double im);

/* ── Print sparse vector (up to max_show entries) ── */
void hpcq_sv_print(const HPCQSparseVector *sv, int max_show);
```

### 4.2 Reconstruction

```c
/* ── Brute-force sparse reconstruction (N ≤ 20 only) ──
 * Enumerates all 2^N configurations. D=2 advantage: bitstring iteration.
 *
 * @param threshold   Minimum probability for inclusion.
 * @param max_entries Maximum entries to collect.
 */
HPCQSparseVector *hpcq_sparse_brute(const HPCQGraph *g,
                                     double threshold,
                                     uint64_t max_entries);

/* ── Tree-pruned sparse reconstruction (larger N) ──
 * At each site, extends live branches to both {0,1}. Prunes branches
 * below threshold. Applies absorb correction after full traversal.
 *
 * @param threshold     Minimum probability for inclusion.
 * @param max_branches  Maximum active branches (truncates lowest).
 */
HPCQSparseVector *hpcq_sparse_tree(const HPCQGraph *g,
                                    double threshold,
                                    uint64_t max_branches);
```

### 4.3 Monte Carlo Expectation

```c
/* ── ⟨ψ|O|ψ⟩ via importance sampling ──
 * Samples from local product distribution, importance-weight by ψ-prob/q-prob.
 * Cost: O(n_samples × (N + E)).
 *
 * @param obs        Observable callback.
 * @param obs_ctx    User context passed to callback.
 * @param n_samples  Number of Monte Carlo samples.
 * @param rng_seed   Seed for internal LCG.
 */
double hpcq_expectation(const HPCQGraph *g,
                         HPCQObservable obs, void *obs_ctx,
                         int n_samples, uint64_t rng_seed);

/* Internal LCG macros (available for user PRNG) */
#define HPCQ_LCG(r)   ((r) = (r) * 6364136223846793005ULL + 1442695040888963407ULL)
#define HPCQ_RAND(r)  (((double)((r) >> 11)) * 0x1.0p-53)
```

---

## 5. HPC Bond Contraction — `hpc_qubit_contract.h`

### 5.1 Hadamard Fold / Unfold

```c
/* ── Decompose 2-vector into symmetric/antisymmetric channels ──
 * plus = (ψ₀ + ψ₁)/√2   (|+⟩ channel)
 * minus = (ψ₀ - ψ₁)/√2   (|−⟩ channel)
 * Cost: O(2), zero multiplies.
 */
HadamardFold hpcq_hadamard_fold(const double re[2], const double im[2]);

/* ── Reconstruct 2-vector from HadamardFold channels ── */
void hpcq_hadamard_unfold(const HadamardFold *hf, double re[2], double im[2]);
```

### 5.2 Pauli Decomposition

```c
/* ── Decompose 2×2 matrix into Pauli basis ──
 * M = c_I·I + c_Z·Z + c_X·X + c_Y·Y
 * where c_P = Tr(P·M)/2.
 *
 * Returns coefficients (re, im), energies (|c_P|²), and dominant channel.
 */
PauliDecomposition hpcq_pauli_decompose(const double M_re[2][2],
                                         const double M_im[2][2]);

/* ── Project onto dominant Pauli component ──
 * Output = c_dominant · P_dominant
 */
void hpcq_pauli_project(const double M_re[2][2], const double M_im[2][2],
                         PauliChannel channel,
                         double out_re[2][2], double out_im[2][2]);

/* ── Fidelity = norm(projected) / norm(original) ── */
double hpcq_compute_fidelity(const double orig_re[2][2], const double orig_im[2][2],
                              const double proj_re[2][2], const double proj_im[2][2]);

/* ── Select dominant Pauli channel (O(4) lookup) ── */
PauliChannel hpcq_select_pauli(const double M_re[2][2], const double M_im[2][2]);
```

**Pauli Channels:**
```c
typedef enum {
    PAULI_I = 0,   /* Identity */
    PAULI_Z = 1,   /* Z-basis */
    PAULI_X = 2,   /* X-basis */
    PAULI_Y = 3    /* Y-basis */
} PauliChannel;
```

### 5.3 Clifford Edge Encoding

```c
/* ── Encode a 2×2 phase interaction as a Clifford edge ──
 * Steps: decompose → select dominant → project → compute fidelity → store.
 * Total: O(16) operations.
 */
void hpcq_encode_clifford(HPCQGraph *g, uint64_t site_a, uint64_t site_b,
                           const double phase_re[2][2], const double phase_im[2][2]);
```

### 5.4 2-Site Gate Encoder

```c
/* ── Auto-select encoding strategy for a 2-qubit gate ──
 * (1) If exact CZ: exact edge (fidelity=1.0)
 * (2) If Pauli fidelity ≥ HPCQ_CLIFFORD_THRESHOLD (0.80): Clifford edge
 * (3) Otherwise: general phase edge via hpcq_general_2site
 */
void hpcq_encode_2site(HPCQGraph *g, uint64_t site_a, uint64_t site_b,
                        const double *G_re, const double *G_im);

/* ── Fold analysis for a single site ──
 * Measures probability in |+⟩ vs |−⟩ channels.
 */
HPCQFoldAnalysis hpcq_analyze_fold(const HPCQGraph *g, uint64_t site);
```

---

## 6. QubitEngine API — `qubit_engine.h` + `.c` files

### 6.1 Lifecycle & Initialization (`qubit_core.c`)

```c
/* ── Initialize engine (zeros all counters, seeds PRNG) ── */
void qubit_engine_init(QubitEngine *eng);

/* ── Destroy engine (zeros all counters) ── */
void qubit_engine_destroy(QubitEngine *eng);

/* ── Allocate a new qubit initialized to |0⟩ ── */
uint32_t qubit_init(QubitEngine *eng);

/* ── Allocate a new qubit initialized to |+⟩ = (|0⟩+|1⟩)/√2 ── */
uint32_t qubit_init_plus(QubitEngine *eng);

/* ── Allocate a new qubit initialized to |k⟩ (k = 0 or 1) ── */
uint32_t qubit_init_basis(QubitEngine *eng, int k);

/* ── PRNG: uniform random double in [0,1) ── */
double qubit_prng_double(QubitEngine *eng);
```

### 6.2 Gates (`qubit_gates.c`)

All gates handle both local qubits and qubits in entangled pairs.

```c
void qubit_apply_hadamard(QubitEngine *eng, uint32_t id);  /* H */
void qubit_apply_x(QubitEngine *eng, uint32_t id);          /* X */
void qubit_apply_y(QubitEngine *eng, uint32_t id);          /* Y */
void qubit_apply_z(QubitEngine *eng, uint32_t id);          /* Z */
void qubit_apply_s(QubitEngine *eng, uint32_t id);          /* S = √Z */
void qubit_apply_t(QubitEngine *eng, uint32_t id);          /* T = √S */
void qubit_apply_phase(QubitEngine *eng, uint32_t id, double theta); /* Z(θ) */

/* ── CZ: controlled-Z ── */
void qubit_apply_cz(QubitEngine *eng, uint32_t id_a, uint32_t id_b);

/* ── CNOT: controlled-X ── */
void qubit_apply_cx(QubitEngine *eng, uint32_t ctrl, uint32_t target);

/* ── Arbitrary single-qubit unitary (2×2 complex matrix) ── */
void qubit_apply_unitary(QubitEngine *eng, uint32_t id,
                          const double *U_re, const double *U_im);

/* ── Arbitrary two-qubit unitary (4×4 complex matrix) ── */
void qubit_apply_unitary_pair(QubitEngine *eng,
                               uint32_t id_a, uint32_t id_b,
                               const double *U_re, const double *U_im);
```

### 6.3 Measurement & Inspection (`qubit_measure.c`)

```c
/* ── Born-rule measurement (collapses to |0⟩ or |1⟩) ──
 * For entangled qubits, measures within the joint state.
 * Returns 0 or 1, or -1 on error.
 */
int qubit_measure(QubitEngine *eng, uint32_t id);

/* ── Measure one qubit of an entangled pair (marginal + collapse) ── */
int qubit_measure_in_pair(QubitEngine *eng, uint32_t id);

/* ── Non-destructive inspection ──
 * Returns probabilities, von Neumann entropy, purity, Schmidt rank.
 */
QubitInspect qubit_inspect(QubitEngine *eng, uint32_t id);
```

### 6.4 Entanglement (`qubit_entangle.c`)

```c
/* ── Create Bell pair (|00⟩+|11⟩)/√2 between id_a and id_b ── */
int qubit_entangle_bell(QubitEngine *eng, uint32_t id_a, uint32_t id_b);

/* ── Entangle as product state |ψ_a⟩⊗|ψ_b⟩ (becomes entangled via CZ) ── */
int qubit_entangle_product(QubitEngine *eng, uint32_t id_a, uint32_t id_b);

/* ── Disentangle: extract marginals back into local states ── */
void qubit_disentangle(QubitEngine *eng, uint32_t id_a, uint32_t id_b);
```

### 6.5 Registers (`qubit_register.c`)

```c
/* ── Initialize a logical register of n_qubits ── */
int qubit_reg_init(QubitEngine *eng, uint64_t chunk_id,
                    uint64_t n_qubits, uint32_t dim);

/* ── GHZ entanglement: (|0⟩^N + |1⟩^N)/√2 (O(1) memory) ── */
void qubit_reg_entangle_all(QubitEngine *eng, int reg_idx);

/* ── Gates on register qubits ── */
void qubit_reg_apply_hadamard(QubitEngine *eng, int reg_idx, uint64_t qubit_idx);
void qubit_reg_apply_cz(QubitEngine *eng, int reg_idx,
                         uint64_t idx_a, uint64_t idx_b);
void qubit_reg_apply_unitary_pos(QubitEngine *eng, int reg_idx, uint64_t pos,
                                  const double *U_re, const double *U_im);

/* ── Measure a register qubit ── */
uint64_t qubit_reg_measure(QubitEngine *eng, int reg_idx, uint64_t qubit_idx);

/* ── State-vector access (for up to ~20 qubit registers) ── */
SV_Amplitude qubit_reg_sv_get(QubitEngine *eng, int reg_idx, basis_t basis_k);
void qubit_reg_sv_set(QubitEngine *eng, int reg_idx,
                       basis_t basis_k, double re, double im);
void qubit_reg_sv_stream(QubitEngine *eng, int reg_idx,
                          sv_stream_fn callback, void *user_data);

/* ── Total probability and inner product ── */
double qubit_reg_sv_total_prob(QubitEngine *eng, int reg_idx);
SV_Amplitude qubit_reg_sv_inner(QubitEngine *eng, int reg_a, int reg_b);

/* ── Get local state vector of a register qubit ── */
QubitState qubit_reg_local_sv(QubitEngine *eng, int reg_idx, uint64_t qubit_pos);
```

---

## 7. Qubit Management Primitives — `qubit_management.h`

Low-level operations on individual QubitState / JointState.

```c
void qm_init_zero(QubitState *s);        /* |0⟩ */
void qm_init_plus(QubitState *s);        /* |+⟩ = (|0⟩+|1⟩)/√2 */
void qm_init_basis(QubitState *s, int k); /* |k⟩ */

void qm_entangle_bell(JointState *j);                     /* (|00⟩+|11⟩)/√2 */
void qm_entangle_product(JointState *j, const QubitState *sa, const QubitState *sb);

double qm_total_prob(const QubitState *s);  /* Σ|ψₖ|² */
void qm_renormalize(QubitState *s);         /* Normalize to unit norm */
void qm_copy(QubitState *dst, const QubitState *src);  /* memcpy */
```

---

## 8. Entanglement Utilities — `entanglement.h`

JointState creation, analysis, and manipulation.

```c
/* ── Bell states ── */
void ent_bell(JointState *j);              /* (|00⟩+|11⟩)/√2 */
void ent_bell_minus(JointState *j);        /* (|00⟩-|11⟩)/√2 */
void ent_bell_psi_plus(JointState *j);     /* (|01⟩+|10⟩)/√2 */
void ent_bell_psi_minus(JointState *j);    /* (|01⟩-|10⟩)/√2 — singlet */

/* ── Product state: |ψ_a⟩⊗|ψ_b⟩ ── */
void ent_product(JointState *j, const QubitState *sa, const QubitState *sb);

/* ── Marginals from joint state ── */
void ent_marginal_a(const JointState *j, double *probs);  /* P_A(a) */
void ent_marginal_b(const JointState *j, double *probs);  /* P_B(b) */

/* ── Entanglement measures ── */
int ent_schmidt_rank(const JointState *j);       /* 1 or 2 */
double ent_entropy(const JointState *j);          /* -Σλₖlog₂(λₖ) */

/* ── Normalization ── */
double ent_total_prob(const JointState *j);      /* Must equal 1.0 */
void ent_renormalize(JointState *j);              /* Force unit norm */
```

---

## 9. Superposition (DFT₂) — `superposition.h`

In-place Hadamard transform on 2 amplitudes.

```c
static const double DFT2_RE[2][2];   /* Twiddle table: 1/√2 × {{1,1},{1,-1}} */
static const double DFT2_IM[2][2];   /* All zeros (Hadamard is real) */

void sup_apply_hadamard(double *re, double *im);      /* In-place H */
void sup_apply_hadamard_inv(double *re, double *im);  /* H† = H */
void sup_renormalize(double *re, double *im, int D);  /* Force Σ|ψₖ|²=1 */
```

---

## 10. Born Rule — `born_rule.h`

```c
double born_prob(double re, double im);                      /* = re²+im² */
double born_total_prob(const double *re, const double *im, int D);  /* Σ|ψₖ|² */
int born_sample(const double *re, const double *im, int D, double rand_01);
    /* Returns 0 if rand_01 < P(0), else 1 */
void born_collapse(double *re, double *im, int D, int outcome);
    /* Post-measurement: set |outcome⟩ to unit magnitude, other to 0 */
double born_fast_isqrt(double x);  /* Fast 1/√x wrapper */
```

---

## 11. Floating-Point Primitives — `arithmetic.h`

IEEE-754 bit-level fast approximations.

```c
/* Fast inverse sqrt (Quake III style, 2 Newton iterations, ~46-bit precision) */
double arith_fast_isqrt(double x);

/* Fast reciprocal (2 Newton iterations) */
double arith_fast_recip(double x);
```

**Magic constants for direct IEEE-754 bit manipulation:**

| Constant | Value | Purpose |
|----------|-------|---------|
| `IEEE754_SIGN_MASK` | `0x8000000000000000ULL` | Sign bit |
| `IEEE754_EXP_MASK` | `0x7FF0000000000000ULL` | Exponent |
| `IEEE754_MANT_MASK` | `0x000FFFFFFFFFFFFFULL` | Mantissa |
| `MAGIC_ISQRT_F64` | `0x5FE6EB50C7B537A9ULL` | Double inv sqrt |
| `MAGIC_RECIP_F64` | `0x7FDE623822FC16E6ULL` | Double reciprocal |
| `MAGIC_SQRT_F64` | `0x1FF7A7EF9DB22D0EULL` | Double sqrt |
| `ARITH_SQRT2_INV` | `0.7071067811865475244` | 1/√2 |

---

## 12. BigInt Library — `bigint.h` / `bigint.c`

4096-bit (64 × 64-bit limbs) arbitrary-precision integer arithmetic.

### Constants

| Macro | Value | Description |
|-------|-------|-------------|
| `BIGINT_LIMBS` | `64` | Number of 64-bit limbs |
| `BIGINT_BYTES` | `512` | Total bytes |
| `BIGINT_BITS` | `4096` | Total bits |

### Core Types

```c
typedef struct {
    uint64_t limbs[64];   /* Little-endian limb order */
} BigInt;
```

### Arithmetic Operations

```c
void bigint_clear(BigInt *a);                        /* a = 0 */
void bigint_copy(BigInt *dst, const BigInt *src);    /* dst = src */
int  bigint_cmp(const BigInt *a, const BigInt *b);   /* returns 1,0,-1 */
int  bigint_is_zero(const BigInt *a);                 /* a == 0? */

void bigint_add(BigInt *result, const BigInt *a, const BigInt *b);   /* r = a+b */
void bigint_sub(BigInt *result, const BigInt *a, const BigInt *b);   /* r = a-b */
void bigint_mul(BigInt *result, const BigInt *a, const BigInt *b);   /* r = a×b */
void bigint_div_mod(const BigInt *dividend, const BigInt *divisor,
                     BigInt *quotient, BigInt *remainder);            /* q,r = N/D */

void bigint_gcd(BigInt *result, const BigInt *a, const BigInt *b);   /* gcd */
void bigint_pow_mod(BigInt *result, const BigInt *base,
                     const BigInt *exp, const BigInt *mod);           /* b^e mod m */
```

### Bit Operations

```c
void bigint_shl1(BigInt *a);                          /* a <<= 1 */
void bigint_shr1(BigInt *a);                          /* a >>= 1 */
int  bigint_get_bit(const BigInt *a, uint32_t bit_index);   /* get bit n */
void bigint_set_bit(BigInt *a, uint32_t bit_index);         /* set bit n */
void bigint_clr_bit(BigInt *a, uint32_t bit_index);         /* clear bit n */
uint32_t bigint_bitlen(const BigInt *a);                     /* MSB position */
```

### Conversion

```c
void bigint_set_u64(BigInt *a, uint64_t val);          /* From uint64 */
uint64_t bigint_to_u64(const BigInt *a);                /* To uint64 (truncates) */
int bigint_from_decimal(BigInt *a, const char *str);    /* Parse decimal string */
void bigint_to_decimal(char *buf, size_t bufsize, const BigInt *a); /* To string */
```

---

## 13. Example Programs

| File | Description | Qubits | Engine |
|------|-------------|--------|--------|
| `Bitstate_template/shor_bs.c` | Shor factoring (up to ~12-bit N) via HPC graph build + classical measurement | 36-60 | HPC Graph |
| `Bitstate_template/vqf_bs.c` | Variational Quantum Factoring (QAOA) — brute-force enumeration ≤20 qubits | ≤20 | HPC Graph |
| `absorb_test.c` | Absorb mechanism unit tests (single/multi/re-absorb) | 2-4 | HPC Graph |
| `supremacy.c` | Random circuit supremacy benchmark | variable | HPC Graph |
| `universal_supremacy.c` | Universal gate set (H,T,CZ) supremacy | variable | HPC Graph |
| `t_gate_full.c` | T-gate depth analysis with absorb verification | 2 | HPC Graph |
| `stress_deep.c` | Deep-circuit stress test | variable | HPC Graph |
| `stress_density.c` | Density / entropy measurement | variable | HPC Graph |
| `clifford_exotic.c` | Clifford exotic benchmarks | variable | QubitEngine |
| `randomized_depth.c` | Randomized circuit depth benchmarks | variable | HPC Graph |

### Quickstart (QubitEngine)

```c
#include "qubit_engine.h"

int main() {
    QubitEngine eng;
    qubit_engine_init(&eng);

    uint32_t q0 = qubit_init(&eng);          /* |0⟩ */
    uint32_t q1 = qubit_init(&eng);          /* |0⟩ */
    qubit_apply_hadamard(&eng, q0);          /* |+⟩ */
    qubit_apply_cx(&eng, q0, q1);            /* Bell state: (|00⟩+|11⟩)/√2 */

    int outcome = qubit_measure(&eng, q0);   /* 0 or 1 */
    printf("q0 = %d\n", outcome);

    qubit_engine_destroy(&eng);
    return 0;
}
```

### Quickstart (HPC Graph)

```c
#include "hpc_qubit_graph.h"
#include "hpc_qubit_amplitude.h"

int main() {
    HPCQGraph *g = hpcq_create(2);           /* 2 qubits */
    hpcq_hadamard(g, 0);                     /* |+⟩ on q0 */
    hpcq_cz(g, 0, 1);                       /* CZ between q0,q1 */
    hpcq_hadamard_absorb(g, 1);              /* H on q1 with absorb */

    /* Amplitude for |11⟩ */
    uint32_t idx[] = {1, 1};
    double re, im;
    hpcq_amplitude(g, idx, &re, &im);
    printf("ψ(1,1) = %.4f + %.4fi\n", re, im);

    /* Sparse reconstruction */
    HPCQSparseVector *sv = hpcq_sparse_tree(g, 1e-10, 1000);
    hpcq_sv_print(sv, 10);
    hpcq_sv_destroy(sv);

    hpcq_destroy(g);
    return 0;
}
```

---

> **Files:** `hpc_qubit_graph.h`, `hpc_qubit_amplitude.h`, `hpc_qubit_contract.h`,
> `qubit_engine.h`, `qubit_core.c`, `qubit_gates.c`, `qubit_measure.c`,
> `qubit_entangle.c`, `qubit_register.c`, `qubit_triality.h/.c`,
> `qubit_management.h`, `entanglement.h`, `superposition.h`, `born_rule.h`,
> `flat_qubit.h`, `statevector.h`, `arithmetic.h`, `bigint.h/.c`,
> `clifford_exotic.c`, `absorb_test.c`, `Bitstate_template/*.c`
