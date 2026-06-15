# BitState — Graph-Native Quantum Simulation

BitState is a **polynomial-memory quantum simulator** that never materializes the state vector. Entanglement is tracked as weighted edges in an **HPC (Holographic Phase Contraction) graph**, giving O(N+E) amplitude evaluation for universal circuits over {H, CZ, T}.

| Simulation target | State-vector memory | BitState memory |
|---|---|---|
| 100 qubits | ~2 × 10¹⁶ PB | ~0.4 MB |
| 1,024 qubits | ~2 × 10²⁸⁵ yottabytes | ~2 MB |
| 262,144 qubits | ~10⁷⁸⁹⁰⁵ configurations | ~1.7 GB |

## How it works

The state is never written as a 2^N array. Instead it is factorized as:

**ψ(x₁, ..., xₙ) = Π_k a_k(x_k) × Π_{edges} w_e(x_a, x_b)**

- **a_k(x_k)** — per-site local amplitude (2 complex numbers per qubit)
- **w_e(x_a, x_b)** — entangling phase on each edge (implicit for CZ, explicit 2×2 matrix for absorbed gates)
- **CZ edges** — 0 extra memory (phase computed on the fly as (-1)^{x_a·x_b})
- **Amplitude evaluation** — O(N + E) per query, no exponential enumeration

The graph is modified gate-by-gate — **locally**:

| Gate | Effect | Complexity |
|---|---|---|
| **H** (0 incident edges) | Modify local state | O(1) |
| **H** (1 incident edge) | Absorb into edge matrix w' | O(1) |
| **H** (>1 incident edges) | Store multi-edge absorb entry | O(deg) |
| **CZ** | Add edge (or cancel existing CZ² = I) | O(deg) via incident list |
| **T, S, phase** | Diagonal update to local state | O(1) |
| **CNOT** | Decomposed as H·CZ·H | built-in |

After an H gate on a site with incident edges, the local state becomes the uniform representer (1, 1) — the site is "consumed" and its amplitude contribution is merged into the edge matrices. Subsequent H gates on the same site are **re-absorptions** using an inner/outer variable split that correctly layers new CZ edges over old ones.

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│                  BitState Engine                           │
├───────────────────┬──────────────┬────────────────────────┤
│   HPC Graph       │   Triality   │   Clifford Exotic     │
│  (hpc_qubit_graph │  (qubit_     │  (clifford_exotic     │
│   .h)             │   triality   │   .h/.c)              │
│  Phase edges      │   .h/.c)     │  24-element group     │
│  O(N+E) amp eval  │  3 views:    │  Invariant Delta_C    │
│  CZ² = I involution│  Z (Edge)   │  Conjugacy fingerprint│
│  Born sampling    │  X (Vertex)  │                       │
│  H-absorb         │  Y (Diagonal)│                       │
├───────────────────┴──────────────┴────────────────────────┤
│                    Core Engine                             │
│  qubit_core.c    — Engine lifecycle, PRNG, allocation      │
│  qubit_gates.c   — H, X, Y, Z, S, T, Phase, CZ, CX       │
│  qubit_measure.c — Born-rule measurement, collapse        │
│  qubit_entangle.c — Bell pairs, product states            │
│  qubit_register.c — Sparse multi-qubit register           │
├───────────────────────────────────────────────────────────┤
│              Header-only Primitives                        │
│  flat_qubit.h   — Auto basis/full mode optimization       │
│  born_rule.h    — Born probability, sampling, collapse    │
│  statevector.h  — QubitState (32B) / JointState (64B)    │
│  entanglement.h — Schmidt analysis                        │
│  bigint.h       — 4096-bit integer arithmetic             │
│  arithmetic.h   — IEEE-754 constants, fast approx math    │
│  superposition.h — Hadamard (DFT₂) transform              │
└───────────────────────────────────────────────────────────┘
```

## Representations

BitState provides three representations of a qubit, each suited to different phases of simulation:

| Representation | File | Memory | Best for |
|---|---|---|---|
| **TrialityQubit** | qubit_triality.h/.c | ~112 B | Gate-heavy circuits with frequent view changes |
| **FlatQubit** | flat_qubit.h | ~32 B + auto-promote | Lightweight single qubits with auto basis/full mode |
| **HPC Graph** | hpc_qubit_graph.h | O(N+E) total | Large-scale circuits where state vector is impossible |

The **HPC Graph** is the flagship representation — it powers all large-scale tests.

## Build

No external dependencies (except `libgmp` for Shor's algorithm).

```bash
# Minimal build (single .c)
gcc -O2 -march=native -o test my_program.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I.
```

Add `bigint.c` and `-lgmp` for Shor's algorithm.

## Tests

```bash
# H-absorption proof of concept (9 tests)
gcc -O2 -march=native -o absorb_test absorb_test.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I. && ./absorb_test

# Full T-gate verification across all subsystems (30 tests)
gcc -O2 -march=native -o t_gate_fl t_gate_full.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I. && ./t_gate_fl

# Universal {H, CZ, T} at scale
gcc -O2 -march=native -o sup universal_supremacy.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I. && ./sup

# Deep randomized circuits (262K qubits, 100 layers)
gcc -O2 -march=native -o stress_deep stress_deep.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I. && ./stress_deep

# Density scaling (1D/2D/3D topologies)
gcc -O2 -march=native -o stress_density stress_density.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c \
    -lm -I. && ./stress_density
```

## Performance

Measured on commodity x86-64 (single core, -O2 -march=native):

| Circuit | Qubits | Edges | Layers | Gates | Time | Throughput |
|---|---|---|---|---|---|---|
| 1D chain T/H | 262,144 | 262,143 | 100 | 26.2M | 1.09s | 24 M ops/s |
| 2D grid T/H | 262,144 | 523,264 | 50 | 13.1M | 0.60s | 22 M ops/s |
| 3D grid T/H | 262,144 | 774,144 | 50 | 13.1M | 0.63s | 21 M ops/s |

Memory: 1.7 GB for 262K qubits at 100 layers. All gates are exact (fidelity = 1.0) — no truncation or approximation.

## Project structure

```
BitState-main/
├── *.c, *.h                 — Engine source and headers
├── absorb_test.c            — H-after-edge absorption tests
├── t_gate_full.c            — T-gate verification (all subsystems)
├── universal_supremacy.c    — {H, CZ, T} scale tests
├── randomized_depth.c       — Deep randomized CZ+T+H circuits
├── stress_deep.c            — 262K-qubit deep T/H stress test
├── stress_density.c         — 1D/2D/3D density scaling
├── supremacy.c              — Graph state scale demo
├── Bitstate_template/
│   ├── template.c           — Full API usage example (13 sections)
└── README.md                — This file
```

## Public API

### HPC Graph — lifecyle and gates

```c
HPCQGraph *hpcq_create(uint64_t n_sites);
void       hpcq_destroy(HPCQGraph *g);

void hpcq_hadamard(HPCQGraph *g, uint64_t site);
void hpcq_hadamard_absorb(HPCQGraph *g, uint64_t site);
void hpcq_t(HPCQGraph *g, uint64_t site);
void hpcq_td(HPCQGraph *g, uint64_t site);
void hpcq_phase(HPCQGraph *g, uint64_t site, double theta);
void hpcq_x(HPCQGraph *g, uint64_t site);
void hpcq_cz(HPCQGraph *g, uint64_t site_a, uint64_t site_b);

void hpcq_amplitude(const HPCQGraph *g, const uint32_t *indices,
                    double *out_re, double *out_im);
double hpcq_probability(const HPCQGraph *g, const uint32_t *indices);
uint32_t hpcq_measure(HPCQGraph *g, uint64_t site, double random_01);
double hpcq_marginal(const HPCQGraph *g, uint64_t site, uint32_t value);
double hpcq_norm_sq(const HPCQGraph *g);

HPCQSparseVector *hpcq_sparse_brute(const HPCQGraph *g, double threshold, uint64_t max_entries);
HPCQSparseVector *hpcq_sparse_tree(const HPCQGraph *g, double threshold, uint64_t max_branches);
```

### QubitEngine — explicit state-vector simulation (small N)

```c
void     qubit_engine_init(QubitEngine *eng);
uint32_t qubit_init(QubitEngine *eng);           // |0⟩
uint32_t qubit_init_plus(QubitEngine *eng);      // |+⟩
void     qubit_apply_hadamard(QubitEngine *eng, uint32_t id);
void     qubit_apply_cz(QubitEngine *eng, uint32_t id_a, uint32_t id_b);
void     qubit_apply_cx(QubitEngine *eng, uint32_t ctrl, uint32_t target);
void     qubit_apply_t(QubitEngine *eng, uint32_t id);
int      qubit_measure(QubitEngine *eng, uint32_t id);
```

### Triality — three-view qubit

```c
void tri_init(TrialityQubit *tq);    // |0⟩ in edge view
void tri_apply_hadamard(TrialityQubit *tq);
void tri_apply_t(TrialityQubit *tq);
void tri_apply_z(TrialityQubit *tq, double theta);
int  tri_measure(TrialityQubit *tq, TrialityView view, double rand_01);
```

See `Bitstate_template/template.c` for complete API coverage across all subsystems.

## Applications

| Algorithm | File | Qubits | Notes |
|---|---|---|---|
| Random circuit sampling | randomized_depth.c | 1,024 | Google Sycamore-style |
| Blind randomized benchmarking | t_gate_full.c | ≤6 | Full 2^N verification |
| Graph state supremacy | universal_supremacy.c | 1,024 | log2(|ψ|²) matches analytical |
| Deep T/H on 262K-chain | stress_deep.c | 262,144 | 100 layers, norm=1 verified |
| Density scaling | stress_density.c | 262,144 | 1D/2D/3D topology comparison |
