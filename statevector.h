/*
 * statevector.h — Cache-Aligned D=2 State Vector Storage
 *
 * QubitState: 2 complex amplitudes = 4 doubles = 32 bytes.
 * JointState: 4 complex amplitudes = 8 doubles = 64 bytes.
 *
 * Storage is SoA (Structure of Arrays): separate re[] and im[] arrays
 * for SIMD-friendly access patterns.
 */

#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include <stdint.h>

#define SV_D              2       /* Qubit dimension                    */
#define SV_D2             4       /* Joint dimension (D×D)              */
#define SV_ELEMENT_SIZE   16      /* 1 complex amplitude = 16 bytes     */
#define SV_STATE_SIZE     32      /* 2 amplitudes × 16 bytes            */
#define SV_JOINT_SIZE     64      /* 4 amplitudes × 16 bytes            */

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUBIT STATE — 2 complex amplitudes (32 bytes)
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double re[SV_D];   /* Real parts:      re[0], re[1]         */
    double im[SV_D];   /* Imaginary parts: im[0], im[1]         */
} QubitState;

/* ═══════════════════════════════════════════════════════════════════════════════
 * JOINT STATE — 4 complex amplitudes (64 bytes)
 * For entangled pairs: ψ(a,b) where a,b ∈ {0,1}
 * Layout: idx = a * 2 + b
 *   [0] = |0,0⟩  [1] = |0,1⟩  [2] = |1,0⟩  [3] = |1,1⟩
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double re[SV_D2];   /* Real parts:      re[0..3]             */
    double im[SV_D2];   /* Imaginary parts: im[0..3]             */
} JointState;

/* ═══════════════════════════════════════════════════════════════════════════════
 * STREAMING AMPLITUDE — For register state vector access
 * ═══════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double re;
    double im;
    double log2_mag;   /* For large-N states where direct magnitude underflows */
} SV_Amplitude;

/* Callback for streaming state vector access */
typedef void (*sv_stream_fn)(uint64_t basis_state, SV_Amplitude amp,
                             void *user_data);

/* ═══════════════════════════════════════════════════════════════════════════════
 * ACCESS PRIMITIVES
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Get amplitude at index k from a QubitState */
static inline SV_Amplitude sv_get_local(const QubitState *s, int k)
{
    SV_Amplitude a;
    a.re = s->re[k];
    a.im = s->im[k];
    a.log2_mag = -1.0/0.0;  /* -INFINITY: use re/im directly */
    return a;
}

/* Get amplitude at index (a,b) from a JointState */
static inline SV_Amplitude sv_get_joint(const JointState *j, int a, int b)
{
    int idx = a * SV_D + b;
    SV_Amplitude amp;
    amp.re = j->re[idx];
    amp.im = j->im[idx];
    amp.log2_mag = -1.0/0.0;
    return amp;
}

#endif /* STATEVECTOR_H */
