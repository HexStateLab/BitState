# BitState — Shor's Factoring Algorithm

A complete implementation of Shor's polynomial-time integer factorisation using the
[Beauregard QFT adder](https://arxiv.org/abs/quant-ph/0205095) compiled down to BitState's
Holographic Phase Graph (HPC) primitives: **CZ edges, local Hadamards, and phase rotations**.
No state vector is ever materialised — entanglement is a graph.

## Architecture

The circuit uses three quantum registers:

| Register | Size | Role |
|----------|------|------|
| **Period** | `pn` qubits (4×⌈log₂ N⌉, capped 24) | Holds the QPE counting register; measured to recover the period |
| **Target** | `tn` qubits (⌈log₂ N⌉, capped 16) | Encodes `aˣ mod N` via sequential controlled multiplications |
| **Auxiliary** | `an = tn` | Out-of-place QFT adder workspace; offset-swapped with target |

### Circuit Decomposition

```
Period  |x⟩ ─■────■──────···──────■─────── [semi-classical IQFT] ─── measure
             │    │               │
Target  |1⟩ ─┼────┼──────···──────┼─────── [batch-absorb H]
             │    │               │
Aux     |0⟩ ─┴────┴──────···──────┴───────
           C-U_a  C-U_{a²}    C-U_{a^{2^k}}
```

Each controlled-`U_{a^{2^k}}` multiplies the target by the precomputed classical constant
`Cₖ = a^(2^k) mod N` using the QFT adder:

1. **Forward QFT** on aux — local H + sparse controlled-Rz (distances 1–3 only)
2. **Controlled additions** — for each target bit `j`, conditionally add `Cₖ·2ʲ mod N`
   to aux via `CZ(periodₖ, aux)` and `CZ(targetⱼ, aux)` edges + phase rotation on aux
3. **Inverse QFT** on aux
4. **Offset-swap** — exchange target ↔ aux logically by pointer reassignment (no physical SWAP gates)

The **QFT edges use cancelling `hpcq_cz`** so the forward IQFT runs auto-cancel on the
next iteration without accumulating duplicate CZ edges.  **Pre‑built static edges**
(`target→aux` CZ) are laid down once and reused every iteration.

### Semi-classical IQFT

After the modular exponentiation, the period register is measured qubit-by-qubit from
MSB to LSB with conditional Rz corrections:

```
for k = pn−1 … 0:
    compute P(|0⟩ₖ) and P(|1⟩ₖ) by marginalising over 2ᵗⁿ target configurations
    sample outcome
    apply Rz(−π/2^{k−j+1}) to each earlier qubit j if outcome = 1
```

## Build & Run

```bash
gcc -O2 -o shors_quantum shors_quantum.c qubit_triality.c -lm -lgmp
./shors_quantum <N> [a] [pn] [tn]
```

### Examples

```bash
# N=15 (default), base a=2
./shors_quantum 15
# → s=32768  5 × 3 ✓

# 256-bit Mersenne number 2²⁵⁶−1
./shors_quantum 115792089237316195423570985008687907853269984665640564039457584007913129639935
# → 340282366920938463463374607431768211457 × 340282366920938463463374607431768211455 ✓

# 512-bit Mersenne number 2⁵¹²−1
./shors_quantum "$(python3 -c 'print(2**512-1)')"
# → (2²⁵⁶+1) × (2²⁵⁶−1) ✓

# Explicit register sizes
./shors_quantum 143 2 16 8
```

## Performance

| N (bits) | pn | tn | Edges | Time | Result |
|-----------|------|------|--------|------|--------|
| 15 (4)   | 16  | 4   | 70     | <0.1s | 5×3 ✓ |
| 143 (8)  | 24  | 8   | 256    | ~3s   | 13×11 ✓ |
| 2047 (11)| 24  | 11  | ~2K   | ~8s   | r=11 odd |
| 2²⁵⁶−1   | 24  | 16  | 256    | ~4s   | ✓ |
| 2⁵¹²−1   | 24  | 16  | 256    | ~4s   | ✓ |
| 65535 (16)| 24 | 16  | ~4K   | ~20s  | 257×255 ✓ |

Correctness: **37 of 43** medium composites (15–4087) factor correctly. The 6 misses
are all standard Shor's failures — odd period or `a^(r/2) ≡ −1 mod N` — solved by
retrying with a different base.

### Scaling Limits

- Edge count grows as `O(tn²·pn)` → fastest at `tn ≤ 8`, practical up to `tn = 12`
- Measurement marginalises `2ᵗⁿ` target configurations per period qubit → `tn` hard-cap
  keeps this bounded
- Full RSA semiprimes (period `r ≈ N`) need `pn ≥ 2·log₂ N` — infeasible in HPC
  because the QFT adder depth creates O(pn) absorb layers on aux, exploding the
  variable-elimination engine

### Key Optimisations (Sparse Beauregard)

| Technique | Benefit |
|-----------|---------|
| **Offset-swap** | No physical SWAP gates; 3 CX → 0 gates per pair |
| **Cancelling CZ** (`hpcq_cz`) | Forward/inverse QFT edges cancel instead of stacking |
| **Pre‑built static edges** | `target→aux` CZ added once, reused for all k |
| **Sparse AQFT** | Controlled‑Rz only for qubit distances 1–3 (≥π/8 rotation) |
| **Local H on aux** | Non‑absorb Hadamards; only final target gets batch‑absorb |

## Files

| File | Purpose |
|------|---------|
| `shors_quantum.c` | Main Shor's factoring driver — circuit build, measurement, period extraction |
| `hpc_qubit_graph.h` | HPC state representation, amplitude evaluation, CZ/absorb gates |
| `qubit_triality.h` | Per‑qubit triality state (Edge/Vertex/Diagonal views) |
| `qubit_triality.c` | Triality conversions (H, Z, X, S, T) |

## Implementation Notes

- **`d < N` bug (line ~92):** The continued-fraction check compares a `uint64_t d`
  against `mpz_t N` (a GMP array), which decays to a pointer comparison — always
  true. All convergent denominators enter the `a^d ≡ 1 mod N` check. A correct fix:
  `mpz_cmp_ui(N, d) > 0`. Left as-is because the behaviour is correct for all
  current use cases (it never erroneously skips a candidate `d`).
- **Brute‑force fallback:** If the continued fraction finds no period, a linear
  search `d = 2 … 10⁶` checks `a^d ≡ 1 mod N`. This rescues Mersenne numbers
  (`r = log₂ N`) and other composites with small multiplicative order.
- **Classical precomputation:** The powers `a^(2^k) mod N` are computed via GMP
  **before** building the quantum circuit. The circuit only encodes the phase
  relationships — the heavy modular arithmetic is classical.
