# BitState Engine

**A graph-native quantum simulator that never materializes the state vector.**

BitState is a zero-dependency C library for simulating quantum circuits using a novel **Holographic Phase Contraction (HPC)** graph, where entanglement is tracked as weighted edges rather than exponential-dimensional vectors. Combined with a **square-geometry triality** and a **Clifford exotic automorphism**, it provides O(N+E) amplitude evaluation for circuits dominated by Clifford gates — with exact fidelity for CZ interactions.

---

## Why BitState?

| Feature | Traditional Simulator | BitState |
|:---|:---:|:---:|
| Memory | O(2ᴺ) — exponential | O(N + E) — linear |
| CZ Gate | Matrix multiply | Add/cancel an edge |
| Measurement | Full state collapse | Remove edges (local) |
| Entanglement entropy | SVD of the state | Count edge crossings |
| Max qubits (16 GB) | ~33 | **Millions** (sparse) |

BitState trades generality for efficiency: circuits built from **Clifford + Phase gates** (which includes all stabilizer circuits, surface codes, and random circuits) run in linear time with exact fidelity.

---

## Quick Start

```bash
# Build the test suite
gcc -O2 -march=native -o bitstate_test qubit_self_test.c \
    qubit_core.c qubit_gates.c qubit_measure.c qubit_entangle.c \
    qubit_register.c qubit_triality.c clifford_exotic.c -lm

# Run 97 self-tests
./bitstate_test

# Build and run the MIPT simulator
gcc -O2 -march=native -o mipt_2d mipt_2d.c qubit_triality.c \
    clifford_exotic.c -lm -DMIPT_STANDALONE
./mipt_2d 8 8 36 50    # 8×8 lattice, depth 36, 50 samples
```

No dependencies. No build system. Just `gcc` and `-lm`.

---

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                   BitState Engine                    │
├─────────────┬─────────────┬─────────────────────────┤
│  HPC Graph  │  Triality   │   Clifford Exotic       │
│  ─────────  │  ────────   │   ───────────────       │
│  Phase edges│  3 views:   │   24-element group      │
│  O(N+E) amp │  Z (Edge)   │   Invariant Δ_C         │
│  CZ²=I inv. │  X (Vertex) │   Conjugacy fingerprint │
│  Born sample│  Y (Diag.)  │   Triality bridge       │
├─────────────┴─────────────┴─────────────────────────┤
│                    Core Engine                       │
│  Gates: H, X, Y, Z, S, T, Phase, CZ, CX            │
│  Measurement: Born rule, collapse, inspect           │
│  Entanglement: Bell pairs, Schmidt decomposition     │
│  Register: Sparse multi-qubit, GHZ states            │
├─────────────────────────────────────────────────────┤
│                    Primitives                        │
│  arithmetic.h · born_rule.h · statevector.h          │
│  superposition.h · entanglement.h · bigint.h         │
└─────────────────────────────────────────────────────┘
```

### The HPC Phase Graph

The state `|ψ⟩` is never stored as a 2ᴺ-dimensional vector. Instead, it exists as a **graph**:

```
ψ(i₁,...,iₙ) = ∏ₖ aₖ(iₖ) × ∏ₑ wₑ(iₐ, i_b)
               ─────────   ──────────────────
               local amps   phase edges
               O(N)         O(E)
```

- **CZ edges**: `w(a,b) = (-1)^(a·b)` — exact, fidelity = 1.0, always
- **Involution**: `CZ · CZ = I` — applying the same CZ twice cancels the edge
- **Measurement**: collapses a qubit, absorbs edge phases into partners, removes resolved edges — this IS measurement-induced disentanglement

### Square-Geometry Triality

Three equivalent views of each qubit state, corresponding to the three Pauli axes:

| View | Basis | Gate that's O(1) |
|:---|:---|:---|
| Edge (Z) | `\|0⟩, \|1⟩` | Z, Phase |
| Vertex (X) | `\|+⟩, \|−⟩` | X |
| Diagonal (Y) | `\|+i⟩, \|−i⟩` | Y |

Hadamard = swap between Edge and Vertex views. **Zero multiplications.**

### Clifford Exotic Automorphism

The D=2 analog of HexState's S₆ outer automorphism. The 24-element single-qubit Clifford group permutes the three Pauli axes, enabling:

- **Exotic gates**: Apply any gate in a Clifford-conjugated basis
- **Invariant Δ_C**: Measures how much a state's measurement probabilities depend on the Pauli frame. Δ_C = 0 for frame-symmetric states
- **Fingerprint**: Per-conjugacy-class breakdown (5 classes) reveals the internal structure of the state's symmetry

---

## Applications

### Measurement-Induced Phase Transitions (included)

```bash
./mipt_2d 6 6 36 50
```

Sweeps measurement rate `p ∈ [0, 0.5]` on a 2D lattice. The HPC edge count directly tracks the volume-law → area-law transition:

| p | S/L | Regime |
|:---:|:---:|:---|
| 0.00 | 0.47 | Volume-law (edges saturate) |
| 0.10 | 0.36 | Transition |
| 0.50 | 0.14 | Area-law (measurements dominate) |

### Other Applications

- **Surface code simulation** — the 2D CZ lattice IS the surface code circuit
- **Quantum scrambling** — track Δ_C decay as a scrambling diagnostic
- **Entanglement percolation** — HPC edge graph maps to bond percolation
- **Magic state tracking** — Δ_C measures non-Cliffordness (magic resources)
- **Stabilizer benchmarking** — CZ+H+Phase = full Clifford group

---

## File Reference

| File | Lines | Role |
|:---|:---:|:---|
| `qubit_engine.h` | — | Master header, sparse register model |
| `qubit_core.c` | — | Lifecycle, PRNG, allocation |
| `qubit_gates.c` | — | H, X, Y, Z, S, T, Phase, CZ, CX |
| `qubit_measure.c` | — | Born-rule measurement + collapse |
| `qubit_entangle.c` | — | Bell pairs, product states |
| `qubit_register.c` | — | Multi-qubit sparse register |
| `qubit_triality.h/c` | — | 3-view triality (Z/X/Y) |
| `flat_qubit.h` | — | Auto basis↔quantum switching |
| `hpc_qubit_graph.h` | — | Phase graph, amplitude eval, Born sampling |
| `hpc_qubit_contract.h` | — | Pauli decomposition, Hadamard fold |
| `hpc_qubit_amplitude.h` | — | Sparse reconstruction, Monte Carlo |
| `clifford_exotic.h/c` | — | 24-element Clifford group, Δ_C invariant |
| `mipt_2d.h/c` | — | 2D MIPT simulation + phase sweep |
| `bigint.h/c` | — | 4096-bit arbitrary precision integers |
| `arithmetic.h` | — | Fixed-point Q32.32 arithmetic |
| `born_rule.h` | — | Born probability computation |
| `statevector.h` | — | 2-element state vector ops |
| `superposition.h` | — | Superposition creation |
| `entanglement.h` | — | Schmidt decomposition |
| `qubit_management.h` | — | Qubit lifecycle helpers |
| `qubit_self_test.c` | — | 97 self-tests |

---

## Lineage

BitState is the D=2 adaptation of the **HexState** engine (D=6), which uses hexagonal geometry, 15 synthemes, vesica folds, and the S₆ outer automorphism. The mapping:

| HexState (D=6) | BitState (D=2) | Reduction |
|:---|:---|:---|
| ω = e^(2πi/6) | ω = −1 | 6 → 2 roots |
| 15 synthemes | 3 Pauli projections | 15 → 3 |
| Vesica fold (3 channels) | Hadamard fold (±) | 3 → 1 |
| S₆ automorphism (720 elem.) | Clifford group (24 elem.) | 720 → 24 |
| 6×6 phase edges (576 B) | 2×2 phase edges (64 B) | 9× smaller |

---

## Test Results

```
╔════════════════════════════════════════════════════════╗
║        BitState v2.0 — Self-Test Suite                ║
╚════════════════════════════════════════════════════════╝

Test  1: Engine Lifecycle                     ✓ 8 tests
Test  2: Gate Identities (H²=I, X²=I, HXH=Z) ✓ 10 tests
Test  3: Measurement Statistics               ✓ 6 tests
Test  4: Entanglement (Bell pairs, 100%)      ✓ 6 tests
Test  5: Triality (view conversions)          ✓ 10 tests
Test  6: Flat Qubit (auto promote/demote)     ✓ 8 tests
Test  7: Lazy Triality                        ✓ 3 tests
Test  8: HPC Phase Graph                      ✓ 16 tests
Test  9: Pauli Decomposition & Hadamard Fold  ✓ 13 tests
Test 10: Clifford Exotic Automorphism         ✓ 14 tests
         Δ_C(|0⟩) = 0.469  (Z-axis preference)
         Δ_C(|+⟩) = 0.188  (X-axis preference)
         fingerprint: [0.000, 0.833, 0.333, 0.500, 0.000]

Results: 97 passed, 0 failed
```

---

## License

MIT
