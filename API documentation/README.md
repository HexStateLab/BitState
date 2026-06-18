# BitState API Reference

BitState is a quantum circuit simulation engine using a **Holographic Phase Graph (HPC)** representation. States are stored as a graph of per-qubit local amplitudes and per-edge phase weights — memory is `O(N + E)`, never `O(2^N)`.

---

## Module 1: BigInt (`bigint.h` / `bigint.c`)

4096-bit arbitrary precision integer arithmetic (64 limbs × 64 bits).

### Type

```c
typedef struct { uint64_t limbs[64]; } BigInt;  // little-endian
```

### Core Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_clear` | `void (BigInt *a)` | Zeros all limbs |
| `bigint_copy` | `void (BigInt *dst, const BigInt *src)` | Deep copy |
| `bigint_cmp` | `int (const BigInt *a, const BigInt *b)` | Returns 1 (a>b), 0 (a==b), -1 (a<b) |
| `bigint_is_zero` | `int (const BigInt *a)` | Returns 1 if value is zero |

### Arithmetic

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_add` | `void (BigInt *result, const BigInt *a, const BigInt *b)` | `result = a + b` |
| `bigint_sub` | `void (BigInt *result, const BigInt *a, const BigInt *b)` | `result = a - b` |
| `bigint_mul` | `void (BigInt *result, const BigInt *a, const BigInt *b)` | `result = a * b` (truncated to 4096 bits) |
| `bigint_div_mod` | `void (const BigInt *dividend, const BigInt *divisor, BigInt *quotient, BigInt *remainder)` | Long division, outputs quotient and remainder |

### Bit Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_shl1` | `void (BigInt *a)` | Shift left by 1 bit (in-place) |
| `bigint_shr1` | `void (BigInt *a)` | Shift right by 1 bit (in-place) |
| `bigint_get_bit` | `int (const BigInt *a, uint32_t bit_index)` | Returns bit value (0 or 1) |
| `bigint_set_bit` | `void (BigInt *a, uint32_t bit_index)` | Sets bit to 1 |
| `bigint_clr_bit` | `void (BigInt *a, uint32_t bit_index)` | Clears bit to 0 |
| `bigint_bitlen` | `uint32_t (const BigInt *a)` | Returns number of bits needed to represent the value |

### Conversion

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_set_u64` | `void (BigInt *a, uint64_t val)` | Sets to a 64-bit unsigned integer |
| `bigint_to_u64` | `uint64_t (const BigInt *a)` | Returns low 64 bits (truncation) |

### Higher-Level

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_gcd` | `void (BigInt *result, const BigInt *a, const BigInt *b)` | Euclidean GCD |
| `bigint_pow_mod` | `void (BigInt *result, const BigInt *base, const BigInt *exp, const BigInt *mod)` | Modular exponentiation (left-to-right binary) |

### String Conversion

| Function | Signature | Description |
|----------|-----------|-------------|
| `bigint_from_decimal` | `int (BigInt *a, const char *str)` | Parses decimal string; returns 0 on success, -1 on error |
| `bigint_to_decimal` | `void (char *buf, size_t bufsize, const BigInt *a)` | Writes decimal string (buffer needs ~1240 chars) |

---

## Module 2: Arithmetic (`arithmetic.h`)

IEEE-754 constants and fast floating-point primitives via bit manipulation.

### IEEE-754 Mask Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `IEEE754_SIGN_MASK` | `0x8000000000000000ULL` | Sign bit mask |
| `IEEE754_EXP_MASK` | `0x7FF0000000000000ULL` | Exponent mask |
| `IEEE754_MANT_MASK` | `0x000FFFFFFFFFFFFFULL` | Mantissa mask |
| `IEEE754_EXP_BIAS` | `1023` | Double exponent bias |
| `IEEE754_MANT_BITS` | `52` | Mantissa bit count |

### Magic Constants

| Constant | Description |
|----------|-------------|
| `MAGIC_ISQRT_F32` | Quake III fast inverse sqrt (float) |
| `MAGIC_ISQRT_F64` | Quake III fast inverse sqrt (double) |
| `MAGIC_RECIP_F64` | Fast reciprocal magic |
| `MAGIC_SQRT_F64` | Fast sqrt magic |

### Substrate Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `ARITH_PHI` | `1.6180339887498948482` | Golden ratio φ |
| `ARITH_PHI_INV` | `0.6180339887498948482` | 1/φ |
| `ARITH_SQRT2` | `1.4142135623730950488` | √2 |
| `ARITH_SQRT2_INV` | `0.7071067811865475244` | 1/√2 |
| `ARITH_DOTTIE` | `0.7390851332151606416` | Fixed point of cos (Dottie number) |

### Fast Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| `arith_fast_isqrt` | `double (double x)` | Fast inverse sqrt (2 Newton iterations, ~46 bits precision) |
| `arith_fast_recip` | `double (double x)` | Fast reciprocal (2 Newton iterations) |

---

## Module 3: State Vector (`statevector.h`)

Cache-aligned state vector storage types.

### Types

| Type | Fields | Size | Description |
|------|--------|------|-------------|
| `QubitState` | `double re[2], im[2]` | 32 bytes | Single-qubit amplitudes (|0⟩, |1⟩) |
| `JointState` | `double re[4], im[4]` | 64 bytes | Two-qubit joint amplitudes; index `a*2+b` maps to |a,b⟩ |
| `SV_Amplitude` | `double re, im, log2_mag` | 24 bytes | Streaming amplitude with log-magnitude for overflow avoidance |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `SV_D` | 2 | Qubit dimension |
| `SV_D2` | 4 | Joint dimension (D×D) |
| `SV_ELEMENT_SIZE` | 16 | Bytes per complex amplitude |
| `SV_STATE_SIZE` | 32 | Bytes per QubitState |
| `SV_JOINT_SIZE` | 64 | Bytes per JointState |

### Access Primitives

| Function | Signature | Description |
|----------|-----------|-------------|
| `sv_get_local` | `SV_Amplitude (const QubitState *s, int k)` | Get amplitude at index k from a QubitState |
| `sv_get_joint` | `SV_Amplitude (const JointState *j, int a, int b)` | Get amplitude at |a,b⟩ from a JointState |

---

## Module 4: Superposition (`superposition.h`)

DFT₂ (Hadamard) transform — the D=2 Fourier transform.

### DFT₂ Twiddle Table

| Variable | Description |
|----------|-------------|
| `DFT2_RE[2][2]` | Real part: `{{1/√2, 1/√2}, {1/√2, -1/√2}}` |
| `DFT2_IM[2][2]` | Imaginary part: all zeros |

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `sup_apply_hadamard` | `void (double *re, double *im)` | In-place Hadamard transform: `H|ψ⟩ = (|ψ[0]+ψ[1]⟩/√2, |ψ[0]-ψ[1]⟩/√2)` |
| `sup_apply_hadamard_inv` | `void (double *re, double *im)` | Inverse Hadamard (identical to forward since H²=I) |
| `sup_renormalize` | `void (double *re, double *im, int D)` | Normalize total probability to 1.0 |

---

## Module 5: Born Rule (`born_rule.h`)

Probability computation, sampling, and state collapse for D=2.

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `born_prob` | `double (double re, double im)` | `re² + im²` |
| `born_total_prob` | `double (const double *re, const double *im, int D)` | Total probability across all D basis states |
| `born_sample` | `int (const double *re, const double *im, int D, double rand_01)` | Born-rule sampling: returns 0 if `rand_01 < P(0)`, else 1 |
| `born_collapse` | `void (double *re, double *im, int D, int outcome)` | Post-measurement collapse to |outcome⟩ (preserves phase of measured component) |
| `born_fast_isqrt` | `double (double x)` | Fast inverse sqrt (~46-bit precision) — delegates to `arith_fast_isqrt` |

---

## Module 6: Entanglement (`entanglement.h`)

Two-qubit joint state storage and entanglement analysis.

### Bell States

| Function | Signature | State Produced |
|----------|-----------|----------------|
| `ent_bell` | `void (JointState *j)` | Φ⁺ = (|00⟩+|11⟩)/√2 |
| `ent_bell_minus` | `void (JointState *j)` | Φ⁻ = (|00⟩-|11⟩)/√2 |
| `ent_bell_psi_plus` | `void (JointState *j)` | Ψ⁺ = (|01⟩+|10⟩)/√2 |
| `ent_bell_psi_minus` | `void (JointState *j)` | Ψ⁻ = (|01⟩-|10⟩)/√2 (singlet) |

### Product States

| Function | Signature | Description |
|----------|-----------|-------------|
| `ent_product` | `void (JointState *j, const QubitState *sa, const QubitState *sb)` | Tensor product |ψ_a⟩ ⊗ |ψ_b⟩ |

### Partial Trace / Marginals

| Function | Signature | Description |
|----------|-----------|-------------|
| `ent_marginal_a` | `void (const JointState *j, double *probs)` | P_A(a) = Σ_b |ψ(a,b)|² |
| `ent_marginal_b` | `void (const JointState *j, double *probs)` | P_B(b) = Σ_a |ψ(a,b)|² |

### Entanglement Metrics

| Function | Signature | Description |
|----------|-----------|-------------|
| `ent_schmidt_rank` | `int (const JointState *j)` | 1 = separable, 2 = entangled (via 2×2 complex determinant) |
| `ent_entropy` | `double (const JointState *j)` | Von Neumann entropy of reduced density matrix ρ_A |
| `ent_total_prob` | `double (const JointState *j)` | Sum of |ψ|² over all 4 basis states |
| `ent_renormalize` | `void (JointState *j)` | Normalize total probability to 1.0 |

---

## Module 7: Qubit Engine (`qubit_engine.h` / `qubit_core.c`)

Top-level container: manages qubits, pairs, registers, and PRNG.

### Types

| Type | Description |
|------|-------------|
| `Qubit` | Single qubit with `QubitState`, `pair_id`, `pair_side` (0=A, 1=B), `id` |
| `QubitPair` | Entangled pair: `JointState`, `id_a`, `id_b`, `active` flag |
| `QubitRegister` | Multi-qubit sparse state vector register (up to 4096 nonzero entries) |
| `QubitEngine` | Master container: `qubits[]`, `pairs[]`, `registers[]`, `rng_state` |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `QUBIT_D` | 2 | Qubit dimension |
| `MAX_QUBITS` | 262144 | Maximum allocatable qubits |
| `MAX_PAIRS` | 262144 | Maximum entangled pairs |
| `MAX_REGISTERS` | 16384 | Maximum registers |
| `MAGIC_PTR(chunk, offset)` | `(((chunk)<<32)\|(offset))` | Encodes sparse pointer |

### Core API — Lifecycle

| Function | Signature | Description |
|----------|-----------|-------------|
| `qubit_engine_init` | `void (QubitEngine *eng)` | Initialize engine with Java LCG seed |
| `qubit_engine_destroy` | `void (QubitEngine *eng)` | Reset engine to empty state |

### Core API — Qubit Allocation

| Function | Signature | Description |
|----------|-----------|-------------|
| `qubit_init` | `uint32_t (QubitEngine *eng)` | Allocate new qubit in |0⟩ state |
| `qubit_init_plus` | `uint32_t (QubitEngine *eng)` | Allocate new qubit in |+⟩ = (|0⟩+|1⟩)/√2 |
| `qubit_init_basis` | `uint32_t (QubitEngine *eng, int k)` | Allocate new qubit in |k⟩ |
| `qubit_prng_double` | `double (QubitEngine *eng)` | Random double in [0,1) (LCG, 53-bit precision) |

---

## Module 8: Qubit Management (`qubit_management.h`)

Low-level per-qubit state primitives.

### Initialization

| Function | Signature | Description |
|----------|-----------|-------------|
| `qm_init_zero` | `void (QubitState *s)` | Set to |0⟩ |
| `qm_init_plus` | `void (QubitState *s)` | Set to |+⟩ = (|0⟩+|1⟩)/√2 |
| `qm_init_basis` | `void (QubitState *s, int k)` | Set to |k⟩ |

### Entanglement Primitives

| Function | Signature | Description |
|----------|-----------|-------------|
| `qm_entangle_bell` | `void (JointState *j)` | Create Φ⁺ Bell state in joint state |
| `qm_entangle_product` | `void (JointState *j, const QubitState *sa, const QubitState *sb)` | Create product state |ψ_a⟩⊗|ψ_b⟩ |

### Utility

| Function | Signature | Description |
|----------|-----------|-------------|
| `qm_total_prob` | `double (const QubitState *s)` | Sum of |ψ[k]|² |
| `qm_renormalize` | `void (QubitState *s)` | Normalize to unit total probability |
| `qm_copy` | `void (QubitState *dst, const QubitState *src)` | Deep copy |

---

## Module 9: Qubit Gates (`qubit_gates.c`)

Standard quantum gates — each handles both local qubits and entangled pairs.

### Single-Qubit Gates

| Function | Signature | Gate | Description |
|----------|-----------|------|-------------|
| `qubit_apply_hadamard` | `void (QubitEngine *eng, uint32_t id)` | H | Hadamard: |0⟩↔|+⟩, |1⟩↔|−⟩ |
| `qubit_apply_x` | `void (QubitEngine *eng, uint32_t id)` | X | Pauli X (bit flip): |0⟩↔|1⟩ |
| `qubit_apply_y` | `void (QubitEngine *eng, uint32_t id)` | Y | Pauli Y ≡ iXZ |
| `qubit_apply_z` | `void (QubitEngine *eng, uint32_t id)` | Z | Pauli Z (phase flip): |1⟩→-|1⟩ |
| `qubit_apply_s` | `void (QubitEngine *eng, uint32_t id)` | S | Phase gate: |1⟩→i|1⟩ |
| `qubit_apply_t` | `void (QubitEngine *eng, uint32_t id)` | T | T gate: |1⟩→e^{iπ/4}|1⟩ |
| `qubit_apply_phase` | `void (QubitEngine *eng, uint32_t id, double theta)` | Rz(θ) | Arbitrary Z-rotation: |1⟩→e^{iθ}|1⟩ |
| `qubit_apply_unitary` | `void (QubitEngine *eng, uint32_t id, const double *U_re, const double *U_im)` | U(2) | Arbitrary 2×2 unitary matrix |

### Two-Qubit Gates

| Function | Signature | Gate | Description |
|----------|-----------|------|-------------|
| `qubit_apply_cz` | `void (QubitEngine *eng, uint32_t id_a, uint32_t id_b)` | CZ | Controlled-Z: |1,1⟩→-|1,1⟩. Auto-creates product pair if not entangled. |
| `qubit_apply_cx` | `void (QubitEngine *eng, uint32_t ctrl, uint32_t target)` | CNOT | Controlled-X: |c,t⟩→|c,c⊕t⟩. Auto-creates product pair if not entangled. |
| `qubit_apply_unitary_pair` | `void (QubitEngine *eng, uint32_t id_a, uint32_t id_b, const double *U_re, const double *U_im)` | U(4) | Arbitrary 4×4 unitary on pair |

---

## Module 10: Qubit Measurement (`qubit_measure.c`)

Born-rule measurement, collapse, and non-destructive inspection.

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `qubit_measure` | `int (QubitEngine *eng, uint32_t id)` | Measure qubit, collapse to |outcome⟩. Returns 0 or 1. |
| `qubit_measure_in_pair` | `int (QubitEngine *eng, uint32_t id)` | Measure one qubit of an entangled pair via marginal probabilities |
| `qubit_inspect` | `QubitInspect (QubitEngine *eng, uint32_t id)` | Non-destructive analysis: returns probabilities, entropy, purity, Schmidt rank |

### `QubitInspect` Structure

| Field | Type | Description |
|-------|------|-------------|
| `prob[2]` | `double` | Born probabilities P(0), P(1) |
| `entropy` | `double` | Shannon entropy -Σ P(k) log₂ P(k) |
| `purity` | `double` | Tr(ρ²) = Σ P(k)² |
| `schmidt_rank` | `int` | 1 = separable, 2 = entangled |

---

## Module 11: Qubit Entanglement (`qubit_entangle.c`)

Creation and destruction of entangled qubit pairs.

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `qubit_entangle_bell` | `int (QubitEngine *eng, uint32_t id_a, uint32_t id_b)` | Create Φ⁺ Bell pair. Returns pair slot index, -1 on failure. |
| `qubit_entangle_product` | `int (QubitEngine *eng, uint32_t id_a, uint32_t id_b)` | Create product pair |ψ_a⟩⊗|ψ_b⟩ from current local states |
| `qubit_disentangle` | `void (QubitEngine *eng, uint32_t id_a, uint32_t id_b)` | Break pair: extract marginal probabilities to local states, deactivate pair |

---

## Module 12: Qubit Register (`qubit_register.c`)

Sparse state vector register for multi-qubit states. Stores up to 4096 nonzero amplitudes.

### Types

| Type | Field | Description |
|------|-------|-------------|
| `RegisterEntry` | `basis_t basis_state, double amp_re, amp_im` | One basis → amplitude mapping |
| `QubitRegister` | entries[4096], n_qubits, chunk_id, collapsed, bulk_rule... | Sparse register container |

`bulk_rule`: `0`=general, `1`=GHZ, `2`=circuit.

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `qubit_reg_init` | `int (QubitEngine *eng, uint64_t chunk_id, uint64_t n_qubits, uint32_t dim)` | Create register. Returns index, -1 on failure. |
| `qubit_reg_entangle_all` | `void (QubitEngine *eng, int reg_idx)` | Set register to GHZ state (|0⟩^N + |1⟩^N)/√2. Sets `bulk_rule=1`. |
| `qubit_reg_apply_hadamard` | `void (QubitEngine *eng, int reg_idx, uint64_t qubit_idx)` | Apply Hadamard to one qubit position in the register |
| `qubit_reg_apply_cz` | `void (QubitEngine *eng, int reg_idx, uint64_t idx_a, uint64_t idx_b)` | Apply CZ between two qubit positions |
| `qubit_reg_apply_unitary_pos` | `void (QubitEngine *eng, int reg_idx, uint64_t pos, const double *U_re, const double *U_im)` | Apply D×D unitary to one qubit position |
| `qubit_reg_measure` | `uint64_t (QubitEngine *eng, int reg_idx, uint64_t qubit_idx)` | Born-rule measurement with collapse. Returns outcome. |
| `qubit_reg_sv_get` | `SV_Amplitude (QubitEngine *eng, int reg_idx, basis_t basis_k)` | Get amplitude at specific basis state |
| `qubit_reg_sv_set` | `void (QubitEngine *eng, int reg_idx, basis_t basis_k, double re, double im)` | Set amplitude at specific basis state |
| `qubit_reg_sv_stream` | `void (QubitEngine *eng, int reg_idx, sv_stream_fn callback, void *user_data)` | Stream all nonzero amplitudes through callback |
| `qubit_reg_sv_total_prob` | `double (QubitEngine *eng, int reg_idx)` | Sum of |ψ|² over all entries |
| `qubit_reg_sv_inner` | `SV_Amplitude (QubitEngine *eng, int reg_a, int reg_b)` | Inner product ⟨ψ_a|ψ_b⟩ between two registers |
| `qubit_reg_local_sv` | `QubitState (QubitEngine *eng, int reg_idx, uint64_t qubit_pos)` | Partial trace to get single-qubit reduced state |

---

## Module 13: Triality (`qubit_triality.h` / `qubit_triality.c`)

Three-view Pauli basis representation with lazy conversion.

### Views

| Enum Value | Pauli Basis | Meaning |
|------------|-------------|---------|
| `VIEW_EDGE` (0) | Z-basis | Computational basis: |0⟩, |1⟩ |
| `VIEW_VERTEX` (1) | X-basis | Hadamard basis: |+⟩, |−⟩ |
| `VIEW_DIAGONAL` (2) | Y-basis | S†·H basis: |+i⟩, |−i⟩ |

### Gate Affinity

- **Z, S, T, Phase**: O(1) in EDGE view (only touch |1⟩)
- **X**: O(1) in VERTEX view (diagonal there)
- **Y**: O(1) in DIAGONAL view (diagonal there)
- **Hadamard**: O(1) relabeling (swaps Edge↔Vertex)

### `TrialityQubit` Structure

| Field | Type | Description |
|-------|------|-------------|
| `re[3][2], im[3][2]` | `double` | Amplitudes in each of 3 views |
| `dirty[3]` | `uint8_t` | View needs recomputation flag |
| `current_view` | `TrialityView` | Last written view |
| `active_mask` | `uint8_t` | Which basis states are nonzero (bit 0=|0⟩, bit 1=|1⟩) |
| `is_eigenstate` | `uint8_t` | True if only one nonzero amplitude |
| `conversions` / `conversion_savings` / `gates_applied` | `uint32_t` | Performance counters |

### Triality API

| Function | Signature | Description |
|----------|-----------|-------------|
| `tri_init` | `void (TrialityQubit *tq)` | Initialize to |0⟩ in EDGE view |
| `tri_init_state` | `void (TrialityQubit *tq, TrialityView view, const double *re, const double *im)` | Initialize from given amplitudes in given view |
| `tri_ensure_view` | `void (TrialityQubit *tq, TrialityView target)` | Lazy conversion: recompute target view if dirty |
| `tri_convert` | `void (TrialityQubit *tq, TrialityView from, TrialityView to)` | Explicit view conversion (all 6 routes supported) |
| `tri_get_amplitudes` | `void (TrialityQubit *tq, TrialityView view, double *re, double *im)` | Get amplitudes in given view (lazy recompute if dirty) |
| `tri_apply_z` | `void (TrialityQubit *tq, double theta)` | Z(θ) gate: routes to EDGE view, marks VERTEX/DIAGONAL dirty |
| `tri_apply_x` | `void (TrialityQubit *tq)` | X gate: routes to VERTEX view, applies Z-like diagonal operation |
| `tri_apply_y` | `void (TrialityQubit *tq)` | Y gate: routes to DIAGONAL view |
| `tri_apply_hadamard` | `void (TrialityQubit *tq)` | H gate: O(1) relabeling — swaps EDGE↔VERTEX amplitudes |
| `tri_apply_s` | `void (TrialityQubit *tq)` | S gate: delegates to Z(π/2) |
| `tri_apply_t` | `void (TrialityQubit *tq)` | T gate: delegates to Z(π/4) |
| `tri_rotate` | `void (TrialityQubit *tq)` | Cycle views: Edge→Vertex→Diagonal→Edge |
| `tri_rotate_back` | `void (TrialityQubit *tq)` | Reverse cycle |
| `tri_measure` | `int (TrialityQubit *tq, TrialityView view, double rand_01)` | Born-rule measurement in specified view. Collapses, marks other views dirty. Returns 0 or 1. |
| `tri_print_stats` | `void (const TrialityQubit *tq)` | Print performance counters |

### Lazy Triality API (`LazyTrialityQubit`)

Heisenberg-picture deferred evaluation: accumulates phase gates and DFT counts, materializes only at measurement.

| Function | Signature | Description |
|----------|-----------|-------------|
| `lazy_tri_init` | `void (LazyTrialityQubit *lq)` | Initialize to |0⟩ |
| `lazy_tri_apply_phase` | `void (LazyTrialityQubit *lq, double theta)` | Accumulate phase (fuses with last segment; up to 64 segments) |
| `lazy_tri_apply_hadamard` | `void (LazyTrialityQubit *lq)` | Increment DFT counter (mod 2 since H²=I) |
| `lazy_tri_materialize` | `void (LazyTrialityQubit *lq)` | Apply all accumulated operations to produce actual amplitudes |
| `lazy_tri_measure` | `int (LazyTrialityQubit *lq, double rand_01)` | Materialize then measure. Collapses to |outcome⟩. Returns 0 or 1. |

---

## Module 14: Flat Qubit (`flat_qubit.h`)

Two-tier flat representation optimization that breathes between complexity levels.

### Representation Tiers

| Mode | Value | When Used | Gate Cost |
|------|-------|-----------|-----------|
| `FLAT_BASIS` | 0 | Single |k⟩ with phase (Pauli eigenstate) | O(1) — all gates are trivial label/phase changes |
| `QUANTUM_FULL` | 1 | General 2-amplitude superposition | O(2) — standard matrix multiplication |

### `FlatQubit` Structure

| Field | Type | Description |
|-------|------|-------------|
| `mode` | `FlatMode` | Current representation tier |
| `basis_k` | `uint8_t` | Which basis state in FLAT_BASIS mode |
| `phase_re/im` | `double` | Phase factor (always unit magnitude) |
| `re[2], im[2]` | `double` | Full amplitudes in QUANTUM_FULL mode |
| `promotions` / `demotions` | `uint32_t` | Conversion counters |

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `flat_init` | `void (FlatQubit *fq)` | Init to |0⟩ in FLAT_BASIS mode |
| `flat_init_basis` | `void (FlatQubit *fq, int k, double ph_re, double ph_im)` | Init to |k⟩ with phase |
| `flat_from_amplitudes` | `void (FlatQubit *fq, const double *re, const double *im)` | Auto-detect mode from amplitudes |
| `flat_promote` | `void (FlatQubit *fq)` | FLAT_BASIS → QUANTUM_FULL (expands to full amplitudes) |
| `flat_try_demote` | `int (FlatQubit *fq)` | QUANTUM_FULL → FLAT_BASIS if possible. Returns 1 on success. |
| `flat_apply_x` | `void (FlatQubit *fq)` | X gate (O(1) in FLAT_BASIS) |
| `flat_apply_z` | `void (FlatQubit *fq)` | Z gate (O(1) in FLAT_BASIS) |
| `flat_apply_phase` | `void (FlatQubit *fq, double theta)` | Phase gate (O(1) in FLAT_BASIS) |
| `flat_apply_t` | `void (FlatQubit *fq)` | T gate = Z(π/4) |
| `flat_apply_td` | `void (FlatQubit *fq)` | T† gate = Z(-π/4) |
| `flat_apply_hadamard` | `void (FlatQubit *fq)` | Hadamard (always promotes, tries to demote after) |
| `flat_get_amplitudes` | `void (const FlatQubit *fq, double *re, double *im)` | Read amplitudes regardless of mode |

---

## Module 15: HPC Qubit Graph (`hpc_qubit_graph.h`)

The core holographic phase graph engine. States are represented as a graph:
```
ψ(i₁,...,iₙ) = [∏ aₖ(iₖ)] × [∏ wₑ(i_a, i_b)]
```

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `HPCQ_D` | 2 | Dimension per site |
| `HPCQ_INIT_EDGES` | 4096 | Initial edge capacity |
| `HPCQ_INIT_LOG` | 8192 | Initial gate log capacity |
| `HPCQ_W2_RE[2]` | `{1.0, -1.0}` | Roots of unity (real) |
| `HPCQ_W2_IM[2]` | `{0.0, 0.0}` | Roots of unity (imag) |
| `HPCQ_CZ_W` / `HPCQ_CZ_W_RE` / `HPCQ_CZ_W_IM` | macros | CZ weight with continuous X-parity: `cos/sin(π·ab + xp_a·b + xp_b·a)` |

### Edge Types

| Enum | Description | Fidelity |
|------|-------------|----------|
| `HPCQ_EDGE_CZ` | Exact CZ: w(a,b)=(-1)^(a·b) | 1.0 |
| `HPCQ_EDGE_PHASE` | General phase: arbitrary stored 2×2 matrix | <1.0 |
| `HPCQ_EDGE_CLIFFORD` | Clifford-projected from Pauli decomposition | <1.0 |
| `HPCQ_EDGE_ABSORBED` | Edge consumed by multi-edge H absorption | 1.0 |

### `HPCQEdge` Structure

| Field | Type | Description |
|-------|------|-------------|
| `type` | `HPCQEdgeType` | Edge type |
| `site_a`, `site_b` | `uint64_t` | Endpoint site indices |
| `w_re[2][2]`, `w_im[2][2]` | `double` | 2×2 phase matrix (for PHASE/CLIFFORD types) |
| `pauli_channel` | `uint8_t` | Clifford: which Pauli basis (0=I, 1=Z, 2=X, 3=Y) |
| `xp_a`, `xp_b` | `double` | Continuous X-rotation (radians) on each endpoint |
| `fidelity` | `double` | Quality metric (1.0 = lossless) |

### `HPCQAbsorbEntry` — Multi-layer Absorption

When H gates are applied to qubits with incident edges, edges are "absorbed" into a structured entry supporting multi-layer H-absorption chains for exact H·CZ·H = CNOT simulation.

| Field | Type | Description |
|-------|------|-------------|
| `center` | `uint64_t` | Site index of the absorbed qubit |
| `n_nbrs` | `uint64_t` | Number of incident neighbors |
| `nbrs` | `uint64_t*` | Neighbor site indices |
| `layer` | `uint8_t*` | Layer index per neighbor |
| `w_re/im` | `double*` | Edge weight matrices per neighbor (size `n_nbrs*4`) |
| `n_layers` | `int` | Number of H-absorption layers |
| `a_re[2]/a_im[2]` | `double` | State BEFORE layer 0's H |
| `a_layer_re/im` | `double*` | Per-layer intermediate states |
| `x_parity` | `uint8_t` | X-gate parity since last H absorption |
| `layer_x_parity` | `uint8_t*` | X parity recorded per layer |

### `HPCQGraph` — The State

| Field | Type | Description |
|-------|------|-------------|
| `n_sites` | `uint64_t` | Number of qubits in simulation |
| `locals` | `TrialityQubit*` | Per-site local states (Z/X/Y views) |
| `n_edges`, `edges` | `uint64_t`, `HPCQEdge*` | Entanglement graph edges |
| `n_absorb`, `absorb` | `uint64_t`, `HPCQAbsorbEntry*` | Absorbed multi-edge entries |
| `gate_log` | `HPCQGateEntry*` | Gate operation log |
| `inc_counts/edges` | Per-site incident edge lists | O(1) edge lookup per site |
| `absorb_idx` | `int64_t*` | Per-site absorb entry index (-1 = none) |
| `amp_evals`, `prob_evals`, `measurements` | `uint64_t` | Performance counters |
| `global_phase_parity` | `uint64_t` | Accounts for Z·X = -X·Z anticommutation |

### Lifecycle

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_create` | `HPCQGraph* (uint64_t n_sites)` | Allocate graph with n_sites qubits, initialized to |0⟩ |
| `hpcq_destroy` | `void (HPCQGraph *g)` | Free all memory |

### Local Gates

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_set_local` | `void (HPCQGraph *g, uint64_t site, const double re[2], const double im[2])` | Set local state explicitly |
| `hpcq_hadamard` | `void (HPCQGraph *g, uint64_t site)` | Simple Hadamard (no incident edges handled). **Prefer `hpcq_hadamard_absorb`.** |
| `hpcq_hadamard_absorb` | `void (HPCQGraph *g, uint64_t site)` | **Correct H gate.** Handles: 0 incident edges (simple H), 1 edge (single-edge absorption), >1 edges (multi-edge absorption), or re-absorption on already-absorbed centers (layered chain). |
| `hpcq_phase` | `void (HPCQGraph *g, uint64_t site, double theta)` | Z(θ) phase gate |
| `hpcq_t` | `void (HPCQGraph *g, uint64_t site)` | T gate: phase by π/4 |
| `hpcq_td` | `void (HPCQGraph *g, uint64_t site)` | T† gate: phase by -π/4 |
| `hpcq_x` | `void (HPCQGraph *g, uint64_t site)` | Pauli X: flips bit on local state, adds π to xp on incident CZ edges, or toggles x_parity on absorbed centers |
| `hpcq_rx` | `void (HPCQGraph *g, uint64_t site, double theta)` | Continuous X-rotation Rx(θ). Adds θ to xp on incident CZ edges. |

### Entangling Gates

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_cz` | `void (HPCQGraph *g, uint64_t site_a, uint64_t site_b)` | **Exact CZ.** If edge already exists between pair, CZ²=I cancellation triggers (swap-removes edge, applies residual Z gates for accumulated xp). Fidelity always 1.0. |
| `hpcq_cz_force` | `void (HPCQGraph *g, uint64_t sa, uint64_t sb)` | CZ **without** cancellation (use when CZ·Rx·CZ ≠ I) |
| `hpcq_general_2site` | `void (HPCQGraph *g, uint64_t site_a, uint64_t site_b, const double *G_re, const double *G_im)` | General 4×4 two-qubit gate: extracts diagonal phases and stores as PHASE edge |

### Amplitude Evaluation

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_amplitude` | `void (const HPCQGraph *g, const uint32_t *indices, double *out_re, double *out_im)` | **O(N+E) point query**: compute ψ(i₁,...,iₙ) for a specific configuration. Handles absorbed centers via connected-component joint sums with H-chain evaluation and variable elimination for large components. |
| `hpcq_probability` | `double (const HPCQGraph *g, const uint32_t *indices)` | |ψ(i₁,...,iₙ)|² |
| `hpcq_marginal` | `double (const HPCQGraph *g, uint64_t site, uint32_t value)` | P(site_k = v) — marginal probability for a single site |
| `hpcq_measure` | `uint32_t (HPCQGraph *g, uint64_t site, double random_01)` | Born-rule measurement: computes marginals, samples, collapses local state, absorbs edge phases into partners, removes resolved edges |
| `hpcq_norm_sq` | `double (const HPCQGraph *g)` | Exhaustive Σ|ψ|² over all 2^N configurations (N ≤ 20 only) |
| `hpcq_entropy_cut` | `double (const HPCQGraph *g, uint64_t cut_after)` | Bipartition entanglement entropy estimate: 1 bit per crossing CZ edge × fidelity |

### Diagnostics

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_print_stats` | `void (const HPCQGraph *g)` | Print graph statistics: sites, edges by type, amp evals, measurements, fidelity, memory usage, full SV size comparison |
| `hpcq_update_fidelity_stats` | `void (HPCQGraph *g)` | Recompute min/avg fidelity from current edges |

### Internal

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_grow_edges` | `void (HPCQGraph *g)` | Double edge capacity if full |
| `hpcq_log_gate` | `void (HPCQGraph *g, HPCQGateEntry entry)` | Record gate in operation log |
| `hpcq_inc_add` | `void (HPCQGraph *g, uint64_t site, uint64_t edge_idx)` | Add edge to site's incident list |
| `hpcq_inc_remove` | `void (HPCQGraph *g, uint64_t site, uint64_t edge_idx)` | Remove edge from site's incident list |

---

## Module 16: HPC Qubit Amplitude (`hpc_qubit_amplitude.h`)

On-demand sparse state vector reconstruction and Monte Carlo expectation values. Never materializes the full 2^N state vector.

### Types

| Type | Description |
|------|-------------|
| `HPCQSparseEntry` | `uint32_t *indices`, `double re, im, prob` |
| `HPCQSparseVector` | `HPCQSparseEntry *entries`, `count`, `capacity`, `total_prob`, `threshold` |
| `HPCQTreeNode` | `uint32_t *indices`, `double re, im` — for tree-pruned reconstruction |
| `HPCQObservable` | `double (*)(const uint32_t *indices, uint64_t n_sites, void *ctx)` |

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_sv_create` | `HPCQSparseVector* (uint64_t n_sites, uint64_t initial_cap)` | Allocate sparse vector |
| `hpcq_sv_destroy` | `void (HPCQSparseVector *sv)` | Free sparse vector |
| `hpcq_sv_add` | `void (HPCQSparseVector *sv, const uint32_t *indices, double re, double im)` | Add an entry (auto-grows) |
| `hpcq_sparse_brute` | `HPCQSparseVector* (const HPCQGraph *g, double threshold, uint64_t max_entries)` | **O(2^N × N+E)**. Enumerate all configurations, keep those above threshold. N ≤ 20 only. |
| `hpcq_sparse_tree` | `HPCQSparseVector* (const HPCQGraph *g, double threshold, uint64_t max_branches)` | **Tree-pruned** reconstruction: binary tree with pruning at each site. Grows 2 children per node (D=2). Handles absorb corrections. For N > 20. |
| `hpcq_expectation` | `double (const HPCQGraph *g, HPCQObservable obs, void *obs_ctx, int n_samples, uint64_t rng_seed)` | **Monte Carlo expectation** ⟨ψ|O|ψ⟩ via importance sampling. Cost: O(n_samples × (N+E)). N ≤ 64. |
| `hpcq_sv_print` | `void (const HPCQSparseVector *sv, int max_show)` | Print sparse vector entries |

---

## Module 17: HPC Qubit Contract (`hpc_qubit_contract.h`)

Pauli-aware bond encoding — the D=2 analog of HexState's syntheme-based contraction.

### Pauli Channels

| Enum | Value | Pauli | Description |
|------|-------|-------|-------------|
| `PAULI_I` | 0 | I | Identity |
| `PAULI_Z` | 1 | Z | Computational basis |
| `PAULI_X` | 2 | X | Hadamard basis |
| `PAULI_Y` | 3 | Y | S†·H basis |

### Predefined Pauli Matrices

| Variable | Description |
|----------|-------------|
| `PAULI_I_RE/IM` | Identity matrix |
| `PAULI_Z_RE/IM` | Z = diag(1, -1) |
| `PAULI_X_RE/IM` | X = anti-diag(1, 1) |
| `PAULI_Y_RE/IM` | Y: im part has anti-diag(-1, 1) |

### Types

| Type | Description |
|------|-------------|
| `HadamardFold` | `plus_re/im`, `minus_re/im` — symmetric/antisymmetric channel decomposition |
| `PauliDecomposition` | `coeff_re/im[4]`, `energy[4]`, `dominant` — Pauli basis decomposition of a 2×2 matrix |
| `HPCQFoldAnalysis` | `plus_fidelity`, `minus_fidelity`, `fold_entropy` — entanglement channel analysis |

### Hadamard Fold (D=2 Vesica Fold)

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_hadamard_fold` | `HadamardFold (const double re[2], const double im[2])` | Decompose state into |+⟩ (symmetric) and |−⟩ (antisymmetric) channels |
| `hpcq_hadamard_unfold` | `void (const HadamardFold *hf, double re[2], double im[2])` | Reconstruct state from fold |

### Pauli Decomposition

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_pauli_decompose` | `PauliDecomposition (const double M_re[2][2], const double M_im[2][2])` | Decompose 2×2 matrix into Pauli basis: M = c_I·I + c_Z·Z + c_X·X + c_Y·Y, where c_P = Tr(P·M)/2 |
| `hpcq_pauli_project` | `void (const double M_re[2][2], const double M_im[2][2], PauliChannel channel, double out_re[2][2], double out_im[2][2])` | Project onto single Pauli channel |
| `hpcq_compute_fidelity` | `double (const double orig_re[2][2], const double orig_im[2][2], const double proj_re[2][2], const double proj_im[2][2])` | Compute fidelity: norm(proj)² / norm(orig)² |
| `hpcq_select_pauli` | `PauliChannel (const double M_re[2][2], const double M_im[2][2])` | Select dominant Pauli channel |

### Edge Encoding

| Function | Signature | Description |
|----------|-----------|-------------|
| `hpcq_encode_clifford` | `void (HPCQGraph *g, uint64_t site_a, uint64_t site_b, const double phase_re[2][2], const double phase_im[2][2])` | Encode as Clifford edge: decompose into Pauli, select dominant, project, normalize, add as HPCQ_EDGE_CLIFFORD |
| `hpcq_encode_2site` | `void (HPCQGraph *g, uint64_t site_a, uint64_t site_b, const double *G_re, const double *G_im)` | Auto-encode 4×4 gate: test for CZ → use hpcq_cz, if Pauli fidelity ≥ 0.80 → clifford encoding, else → general phase edge |
| `hpcq_analyze_fold` | `HPCQFoldAnalysis (const HPCQGraph *g, uint64_t site)` | Analyze entanglement in symmetric/antisymmetric channels for a site |

---

## Module 18: HPC qLDPC (`hpc_qldpc.h`)

Quantum Low-Density Parity-Check codes with hypergraph product construction.

### Types

| Type | Description |
|------|-------------|
| `QLDPCCode` | Contains H1[r1×n1], H2[r2×n2], qubit index arrays, ancilla counts, total simulation qubit count |

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `qldpc_create` | `QLDPCCode* (int n1, int n2, int w_c, int w_r, uint64_t *seed)` | Create hypergraph product code: n1×n2 data qubits, r1=r2=n1/2 checks. `w_c`=column weight, `w_r`=row weight. Automatically allocates qubit IDs. |
| `qldpc_destroy` | `void (QLDPCCode *c)` | Free all memory |
| `qldpc_build_circuit` | `void (HPCQGraph *g, QLDPCCode *c)` | Build syndrome extraction circuit on HPC graph: Z-stabilizers (H1 rows × H2 cols) and X-stabilizers (H1 cols × H2 rows) |
| `qldpc_print` | `void (QLDPCCode *c)` | Print code statistics: physical/logical qubits, rate, stabilizer counts |
| `qlcg` | `uint64_t (uint64_t *s)` | LCG: `s = s*6364136223846793005+1442695040888963407` |
| `qldpc_gen_H` | `uint8_t* (int r, int n, int w_c, int w_r, uint64_t *seed)` | Generate sparse random (w_c,w_r)-regular LDPC matrix |

---

## Module 19: HPC Shor's Algorithm (`hpc_shors.h`)

Shor's algorithm factorization on the HPC graph.

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `sgcd` | `uint64_t (uint64_t a, uint64_t b)` | Euclidean GCD |
| `scrz` | `void (HPCQGraph *g, uint64_t c, uint64_t t, double th)` | Controlled-Rz via CNOT = H·CZ·H decomposition |
| `sqft` | `void (HPCQGraph *g, uint64_t *q, int n, int inv)` | Quantum Fourier Transform / Inverse QFT |
| `scmodmul` | `void (HPCQGraph *g, uint64_t ctrl, uint64_t *targ, int tn, uint64_t c, uint64_t N)` | Controlled modulo multiplication: |x⟩ → |c·x mod N⟩ |
| `shors_circuit` | `void (HPCQGraph *g, uint64_t N, uint64_t a, uint64_t *per, int pn, uint64_t *targ, int tn)` | Build full Shor's circuit: initialize period register, apply controlled modular exponentiations a^(2^k) mod N, IQFT |
| `shors_factor` | `int (HPCQGraph *g, uint64_t N, uint64_t a, int pn, int tn, uint64_t *f1, uint64_t *f2)` | Factor N: enumerate period register to find period r, compute gcd(a^(r/2) ± 1, N). Returns 1 on success (f1, f2 = factors). Returns 0 on failure. |

---

## Module 20: Clifford Exotic (`clifford_exotic.h` / `clifford_exotic.c`)

Clifford group exotic automorphism for D=2. The single-qubit Clifford group C₁ has 24 elements (isomorphic to S₄).

### Types

| Type | Description |
|------|-------------|
| `CliffordElement` | `int axes[3], sign[3]` — maps (Z,X,Y)→(sign[0]·P_axes[0], sign[1]·P_axes[1], sign[2]·P_axes[2]) |
| `CliffordGate` | `double re[2][2], im[2][2]` — 2×2 unitary matrix |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `CLIFF_ORDER` | 24 | Number of single-qubit Clifford elements |
| `CLIFF_N_AXES` | 3 | Pauli axes: Z, X, Y |
| `CLIFF_NUM_CLASSES` | 5 | Number of conjugacy classes |

### Conjugacy Classes

| Index | Description | Size |
|-------|-------------|------|
| 0 | Identity {I} | 1 |
| 1 | Half-turns {±Z, ±X, ±Y} | 6 |
| 2 | Quarter-turns {±S, ±S†, ...} | 6 |
| 3 | Third-turns {H, SH, ...} | 8 |
| 4 | Sixth-turns {HS, ...} | 3 |

### Group API

| Function | Signature | Description |
|----------|-----------|-------------|
| `cliff_init` | `void (void)` | Generate all 24 Clifford elements and their 2×2 matrices (idempotent) |
| `cliff_get_element` | `const CliffordElement* (int idx)` | Get Clifford element by index |
| `cliff_get_gate` | `const CliffordGate* (int idx)` | Get 2×2 gate matrix by index |
| `cliff_compose` | `int (int a, int b)` | Compose two elements: returns index of a∘b |
| `cliff_inverse` | `int (int idx)` | Returns index of inverse element |
| `cliff_identity` | `void (void)` | Returns index of identity element |

### Exotic Gate Application

| Function | Signature | Description |
|----------|-----------|-------------|
| `cliff_apply_exotic_gate` | `void (const double *in_re, const double *in_im, double *out_re, double *out_im, int clifford_idx)` | Apply C·|ψ⟩: transform state by Clifford gate matrix |
| `cliff_dual_probabilities` | `void (const double *re, const double *im, double *probs_std, double *probs_exo, int clifford_idx)` | Compute standard AND Clifford-conjugated Born probabilities |

### Invariants

| Function | Signature | Description |
|----------|-----------|-------------|
| `cliff_exotic_invariant` | `double (const double *re, const double *im)` | Δ_C = (1/24) Σ_C ||P_std - P_C||². 0 = symmetric under all Clifford conjugations (no preferred Pauli frame). |
| `cliff_exotic_entropy` | `double (const double *re, const double *im, int clifford_idx)` | ΔS = S_std - S_exotic (Shannon entropy difference) |
| `cliff_exotic_fingerprint` | `void (const double *re, const double *im, double *class_deltas)` | Per-conjugacy-class invariant breakdown (5 values) |

### Triality Bridge

| Function | Signature | Description |
|----------|-----------|-------------|
| `cliff_target_view` | `int (int clifford_idx)` | Which triality view this Clifford element maps Z to: 0=Edge, 1=Vertex, 2=Diagonal |

---

## Architecture Notes

### Memory Model

- **QubitEngine mode**: N qubits × 32 bytes + P pairs × 64 bytes = O(N+P)
- **HPC Graph mode**: N sites × 112 bytes (TrialityQubit) + E edges × ~136 bytes = O(N+E)
- State vector is **never materialized** — amplitudes computed on demand

### Key Design Properties

1. **CZ² = I cancellation**: Re-applying CZ to the same pair removes the edge, preventing unbounded edge growth
2. **Continuous X-parity**: `xp_a, xp_b` track accumulated X-rotations on CZ edges for exact `Rx(θ)·CZ·Rx(θ)` decomposition
3. **Multi-layer absorption**: `hpcq_hadamard_absorb` correctly handles `H·CZ·H = CNOT` via layered absorption chains (verified to depth L=4+)
4. **Topology-agnostic**: Any qubit can connect to any qubit via CZ edges — no grid restriction
5. **Absorption evaluation**: Absorbed centers are evaluated via connected-component joint sums; large components use variable elimination with `O(2^treewidth)` complexity rather than `O(2^total_vars)`

### Performance Benchmarks (from README)

| Metric | Google Willow | BitState |
|--------|---------------|----------|
| Qubits | 105 | **4,000,000** |
| Cycles | ~25 | 100 |
| Memory | N/A | 1.37 GB |
| Build time | N/A | 0.70 s |

### Verification

BitState produces numerically identical amplitudes to brute-force state vector simulation:
- 16 qubits, 16 cycles: 65536 basis states — 0 mismatches, fidelity² = 1.0
- 9 qubits, 12 cycles: 512 basis states — 0 mismatches, fidelity² = 1.0
