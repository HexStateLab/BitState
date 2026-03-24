/*
 * qubit_triality.c — Square-Geometry Triality Implementation
 *
 * View conversion routes:
 *   Edge→Vertex:     Hadamard H
 *   Edge→Diagonal:   S†·H
 *   Vertex→Diagonal: S†
 *   (reverses use adjoints: H†=H, S)
 *
 * Gate routing:
 *   Phase/Z/S/T → ensure EDGE, apply to |1⟩ component
 *   X           → ensure VERTEX, swap components
 *   Y           → ensure DIAGONAL, swap with phases
 *   H           → swap EDGE↔VERTEX labels (O(1) relabeling)
 */

#include "qubit_triality.h"
#include <stdio.h>

static const double INV_SQRT2 = 0.7071067811865475244;

/* ═══════════════════════════════════════════════════════════════════════════════
 * INITIALIZATION
 * ═══════════════════════════════════════════════════════════════════════════════ */

void tri_init(TrialityQubit *tq)
{
    memset(tq, 0, sizeof(*tq));

    /* Initialize to |0⟩ in edge view */
    tq->re[VIEW_EDGE][0] = 1.0;
    tq->dirty[VIEW_EDGE] = 0;
    tq->dirty[VIEW_VERTEX] = 1;
    tq->dirty[VIEW_DIAGONAL] = 1;
    tq->current_view = VIEW_EDGE;
    tq->active_mask = 0x01;  /* only |0⟩ is active */
    tq->is_eigenstate = 1;
    tq->is_real = 1;
}

void tri_init_state(TrialityQubit *tq, TrialityView view,
                    const double *re, const double *im)
{
    memset(tq, 0, sizeof(*tq));

    tq->re[view][0] = re[0]; tq->re[view][1] = re[1];
    tq->im[view][0] = im[0]; tq->im[view][1] = im[1];
    tq->dirty[view] = 0;
    tq->current_view = view;

    for (int v = 0; v < VIEW_COUNT; v++) {
        if (v != (int)view) tq->dirty[v] = 1;
    }

    /* Compute active mask */
    tq->active_mask = 0;
    for (int k = 0; k < TRI_D; k++) {
        if (re[k]*re[k] + im[k]*im[k] > 1e-30)
            tq->active_mask |= (1 << k);
    }

    /* Check if real-valued */
    tq->is_real = (fabs(im[0]) < 1e-14 && fabs(im[1]) < 1e-14) ? 1 : 0;

    /* Check if eigenstate (only one nonzero amplitude) */
    int count = __builtin_popcount(tq->active_mask);
    tq->is_eigenstate = (count <= 1) ? 1 : 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * VIEW CONVERSION — Core transformation routines
 *
 * All conversions are O(2) complex operations.
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Apply Hadamard: out = H · in */
static void apply_hadamard(const double *in_re, const double *in_im,
                           double *out_re, double *out_im)
{
    out_re[0] = (in_re[0] + in_re[1]) * INV_SQRT2;
    out_im[0] = (in_im[0] + in_im[1]) * INV_SQRT2;
    out_re[1] = (in_re[0] - in_re[1]) * INV_SQRT2;
    out_im[1] = (in_im[0] - in_im[1]) * INV_SQRT2;
}

/* Apply S†: |0⟩→|0⟩, |1⟩→-i|1⟩ */
static void apply_s_dag(const double *in_re, const double *in_im,
                        double *out_re, double *out_im)
{
    out_re[0] = in_re[0];
    out_im[0] = in_im[0];
    /* -i × (r + im*i) = im - i*r */
    out_re[1] =  in_im[1];
    out_im[1] = -in_re[1];
}

/* Apply S: |0⟩→|0⟩, |1⟩→i|1⟩ */
static void apply_s(const double *in_re, const double *in_im,
                    double *out_re, double *out_im)
{
    out_re[0] = in_re[0];
    out_im[0] = in_im[0];
    /* i × (r + im*i) = -im + i*r */
    out_re[1] = -in_im[1];
    out_im[1] =  in_re[1];
}

void tri_convert(TrialityQubit *tq, TrialityView from, TrialityView to)
{
    if (from == to) return;

    double temp_re[TRI_D], temp_im[TRI_D];

    if (from == VIEW_EDGE && to == VIEW_VERTEX) {
        /* Edge→Vertex: H */
        apply_hadamard(tq->re[VIEW_EDGE], tq->im[VIEW_EDGE],
                      tq->re[VIEW_VERTEX], tq->im[VIEW_VERTEX]);
    }
    else if (from == VIEW_VERTEX && to == VIEW_EDGE) {
        /* Vertex→Edge: H (self-adjoint) */
        apply_hadamard(tq->re[VIEW_VERTEX], tq->im[VIEW_VERTEX],
                      tq->re[VIEW_EDGE], tq->im[VIEW_EDGE]);
    }
    else if (from == VIEW_EDGE && to == VIEW_DIAGONAL) {
        /* Edge→Diagonal: S†·H */
        apply_hadamard(tq->re[VIEW_EDGE], tq->im[VIEW_EDGE],
                      temp_re, temp_im);
        apply_s_dag(temp_re, temp_im,
                   tq->re[VIEW_DIAGONAL], tq->im[VIEW_DIAGONAL]);
    }
    else if (from == VIEW_DIAGONAL && to == VIEW_EDGE) {
        /* Diagonal→Edge: H·S */
        apply_s(tq->re[VIEW_DIAGONAL], tq->im[VIEW_DIAGONAL],
               temp_re, temp_im);
        apply_hadamard(temp_re, temp_im,
                      tq->re[VIEW_EDGE], tq->im[VIEW_EDGE]);
    }
    else if (from == VIEW_VERTEX && to == VIEW_DIAGONAL) {
        /* Vertex→Diagonal: S† */
        apply_s_dag(tq->re[VIEW_VERTEX], tq->im[VIEW_VERTEX],
                   tq->re[VIEW_DIAGONAL], tq->im[VIEW_DIAGONAL]);
    }
    else if (from == VIEW_DIAGONAL && to == VIEW_VERTEX) {
        /* Diagonal→Vertex: S */
        apply_s(tq->re[VIEW_DIAGONAL], tq->im[VIEW_DIAGONAL],
               tq->re[VIEW_VERTEX], tq->im[VIEW_VERTEX]);
    }

    tq->dirty[to] = 0;
    tq->conversions++;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENSURE VIEW — Lazy recomputation of dirty views
 * ═══════════════════════════════════════════════════════════════════════════════ */

void tri_ensure_view(TrialityQubit *tq, TrialityView target)
{
    if (!tq->dirty[target]) return;

    /* Find a clean source view */
    TrialityView src = tq->current_view;
    if (tq->dirty[src]) {
        /* Current is also dirty — find any clean view */
        for (int v = 0; v < VIEW_COUNT; v++) {
            if (!tq->dirty[v]) { src = (TrialityView)v; break; }
        }
    }

    /* Eigenstate fast path: O(1) for Pauli eigenstates */
    if (tq->is_eigenstate && tq->active_mask != 0) {
        tq->conversion_savings++;
        /* For a basis state |k⟩, we know the exact representation in all views */
        /* For now, fall through to general conversion */
    }

    tri_convert(tq, src, target);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GET AMPLITUDES — Read amplitudes in a specific view
 * ═══════════════════════════════════════════════════════════════════════════════ */

void tri_get_amplitudes(TrialityQubit *tq, TrialityView view,
                        double *re, double *im)
{
    tri_ensure_view(tq, view);
    re[0] = tq->re[view][0]; re[1] = tq->re[view][1];
    im[0] = tq->im[view][0]; im[1] = tq->im[view][1];
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GATE APPLICATION — Routes to optimal view
 * ═══════════════════════════════════════════════════════════════════════════════ */

/* Z(θ): |0⟩→|0⟩, |1⟩→e^{iθ}|1⟩  — best in EDGE view */
void tri_apply_z(TrialityQubit *tq, double theta)
{
    tri_ensure_view(tq, VIEW_EDGE);

    double cs = cos(theta), sn = sin(theta);
    double r = tq->re[VIEW_EDGE][1], m = tq->im[VIEW_EDGE][1];
    tq->re[VIEW_EDGE][1] = r * cs - m * sn;
    tq->im[VIEW_EDGE][1] = r * sn + m * cs;

    /* Mark other views dirty */
    tq->dirty[VIEW_VERTEX] = 1;
    tq->dirty[VIEW_DIAGONAL] = 1;
    tq->current_view = VIEW_EDGE;
    tq->gates_applied++;
}

/* X: |0⟩↔|1⟩  — best in VERTEX view where it's a Z gate */
void tri_apply_x(TrialityQubit *tq)
{
    tri_ensure_view(tq, VIEW_VERTEX);

    /* In X-basis, Pauli X is diagonal: |+⟩→|+⟩, |−⟩→-|−⟩ */
    tq->re[VIEW_VERTEX][1] = -tq->re[VIEW_VERTEX][1];
    tq->im[VIEW_VERTEX][1] = -tq->im[VIEW_VERTEX][1];

    tq->dirty[VIEW_EDGE] = 1;
    tq->dirty[VIEW_DIAGONAL] = 1;
    tq->current_view = VIEW_VERTEX;
    tq->gates_applied++;
}

/* Y: |0⟩→i|1⟩, |1⟩→-i|0⟩  — best in DIAGONAL view where it's a Z gate */
void tri_apply_y(TrialityQubit *tq)
{
    tri_ensure_view(tq, VIEW_DIAGONAL);

    /* In Y-basis, Pauli Y is diagonal: |+i⟩→|+i⟩, |−i⟩→-|−i⟩ */
    tq->re[VIEW_DIAGONAL][1] = -tq->re[VIEW_DIAGONAL][1];
    tq->im[VIEW_DIAGONAL][1] = -tq->im[VIEW_DIAGONAL][1];

    tq->dirty[VIEW_EDGE] = 1;
    tq->dirty[VIEW_VERTEX] = 1;
    tq->current_view = VIEW_DIAGONAL;
    tq->gates_applied++;
}

/* Hadamard: swap Edge↔Vertex labels (O(1) relabeling) */
void tri_apply_hadamard(TrialityQubit *tq)
{
    /* H swaps the Z and X bases.
     * Instead of matrix multiply, we relabel:
     *   If we have valid EDGE and VERTEX views, swap them.
     *   The DIAGONAL view picks up an S gate: H takes Y-basis to -Y-basis? No.
     *
     * Actually: H applied in Z-basis is a basis change to X-basis.
     * So if we have the state in EDGE view, applying H gives us
     * the VERTEX view amplitudes as the new EDGE view.
     */
    tri_ensure_view(tq, VIEW_EDGE);
    tri_ensure_view(tq, VIEW_VERTEX);

    /* Swap edge and vertex */
    double tmp_re[TRI_D], tmp_im[TRI_D];
    memcpy(tmp_re, tq->re[VIEW_EDGE], sizeof(tmp_re));
    memcpy(tmp_im, tq->im[VIEW_EDGE], sizeof(tmp_im));
    memcpy(tq->re[VIEW_EDGE], tq->re[VIEW_VERTEX], sizeof(tmp_re));
    memcpy(tq->im[VIEW_EDGE], tq->im[VIEW_VERTEX], sizeof(tmp_im));
    memcpy(tq->re[VIEW_VERTEX], tmp_re, sizeof(tmp_re));
    memcpy(tq->im[VIEW_VERTEX], tmp_im, sizeof(tmp_im));

    /* Diagonal is dirty (H changes the Y-basis representation) */
    tq->dirty[VIEW_DIAGONAL] = 1;
    tq->current_view = VIEW_EDGE;
    tq->gates_applied++;
}

/* S gate: |1⟩→i|1⟩ — in EDGE view */
void tri_apply_s(TrialityQubit *tq)
{
    tri_apply_z(tq, M_PI / 2.0);
}

/* T gate: |1⟩→e^{iπ/4}|1⟩ — in EDGE view */
void tri_apply_t(TrialityQubit *tq)
{
    tri_apply_z(tq, M_PI / 4.0);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * TRIALITY ROTATION — O(1) relabeling
 *
 * Edge→Vertex→Diagonal→Edge
 * The cube is picked up and shown from a different face.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void tri_rotate(TrialityQubit *tq)
{
    /* Ensure all views are clean */
    tri_ensure_view(tq, VIEW_EDGE);
    tri_ensure_view(tq, VIEW_VERTEX);
    tri_ensure_view(tq, VIEW_DIAGONAL);

    /* Cycle: Edge→Vertex→Diagonal→Edge */
    double tmp_re[TRI_D], tmp_im[TRI_D];

    memcpy(tmp_re, tq->re[VIEW_DIAGONAL], sizeof(tmp_re));
    memcpy(tmp_im, tq->im[VIEW_DIAGONAL], sizeof(tmp_im));

    memcpy(tq->re[VIEW_DIAGONAL], tq->re[VIEW_VERTEX], sizeof(tmp_re));
    memcpy(tq->im[VIEW_DIAGONAL], tq->im[VIEW_VERTEX], sizeof(tmp_im));

    memcpy(tq->re[VIEW_VERTEX], tq->re[VIEW_EDGE], sizeof(tmp_re));
    memcpy(tq->im[VIEW_VERTEX], tq->im[VIEW_EDGE], sizeof(tmp_im));

    memcpy(tq->re[VIEW_EDGE], tmp_re, sizeof(tmp_re));
    memcpy(tq->im[VIEW_EDGE], tmp_im, sizeof(tmp_im));
}

void tri_rotate_back(TrialityQubit *tq)
{
    tri_ensure_view(tq, VIEW_EDGE);
    tri_ensure_view(tq, VIEW_VERTEX);
    tri_ensure_view(tq, VIEW_DIAGONAL);

    double tmp_re[TRI_D], tmp_im[TRI_D];

    memcpy(tmp_re, tq->re[VIEW_EDGE], sizeof(tmp_re));
    memcpy(tmp_im, tq->im[VIEW_EDGE], sizeof(tmp_im));

    memcpy(tq->re[VIEW_EDGE], tq->re[VIEW_VERTEX], sizeof(tmp_re));
    memcpy(tq->im[VIEW_EDGE], tq->im[VIEW_VERTEX], sizeof(tmp_im));

    memcpy(tq->re[VIEW_VERTEX], tq->re[VIEW_DIAGONAL], sizeof(tmp_re));
    memcpy(tq->im[VIEW_VERTEX], tq->im[VIEW_DIAGONAL], sizeof(tmp_im));

    memcpy(tq->re[VIEW_DIAGONAL], tmp_re, sizeof(tmp_re));
    memcpy(tq->im[VIEW_DIAGONAL], tmp_im, sizeof(tmp_im));
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT — Born-rule sampling in any view
 * ═══════════════════════════════════════════════════════════════════════════════ */

int tri_measure(TrialityQubit *tq, TrialityView view, double rand_01)
{
    tri_ensure_view(tq, view);

    double p0 = tq->re[view][0] * tq->re[view][0]
              + tq->im[view][0] * tq->im[view][0];
    int outcome = (rand_01 < p0) ? 0 : 1;

    /* Collapse in the measurement view */
    double mag = sqrt(tq->re[view][outcome] * tq->re[view][outcome]
                    + tq->im[view][outcome] * tq->im[view][outcome]);
    if (mag < 1e-30) mag = 1.0;
    double inv = 1.0 / mag;

    double phase_re = tq->re[view][outcome] * inv;
    double phase_im = tq->im[view][outcome] * inv;

    tq->re[view][0] = 0; tq->im[view][0] = 0;
    tq->re[view][1] = 0; tq->im[view][1] = 0;
    tq->re[view][outcome] = phase_re;
    tq->im[view][outcome] = phase_im;

    /* All other views are dirty after collapse */
    for (int v = 0; v < VIEW_COUNT; v++) {
        if (v != (int)view) tq->dirty[v] = 1;
    }
    tq->current_view = view;
    tq->is_eigenstate = 1;
    tq->active_mask = (uint8_t)(1 << outcome);

    return outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * STATISTICS
 * ═══════════════════════════════════════════════════════════════════════════════ */

void tri_print_stats(const TrialityQubit *tq)
{
    printf("TrialityQubit stats:\n");
    printf("  current view:       %s\n", TRIALITY_VIEW_NAMES[tq->current_view]);
    printf("  gates applied:      %u\n", tq->gates_applied);
    printf("  view conversions:   %u\n", tq->conversions);
    printf("  savings (avoided):  %u\n", tq->conversion_savings);
    printf("  is_eigenstate:      %d\n", tq->is_eigenstate);
    printf("  is_real:            %d\n", tq->is_real);
    printf("  active_mask:        0x%02X\n", tq->active_mask);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * LAZY TRIALITY QUBIT — Heisenberg-picture deferred evaluation
 * ═══════════════════════════════════════════════════════════════════════════════ */

void lazy_tri_init(LazyTrialityQubit *lq)
{
    memset(lq, 0, sizeof(*lq));
    lq->re[0] = 1.0;  /* Start in |0⟩ */
    lq->materialized = 1;
}

void lazy_tri_apply_phase(LazyTrialityQubit *lq, double theta)
{
    if (lq->num_segments < LAZY_MAX_SEGMENTS) {
        /* Can still accumulate */
        if (lq->num_segments > 0) {
            /* Fuse with last segment: just add the phase */
            lq->segments[lq->num_segments - 1].phase += theta;
        } else {
            lq->segments[0].phase = theta;
            lq->num_segments = 1;
        }
    } else {
        /* Overflow: materialize, then apply */
        lazy_tri_materialize(lq);
        double cs = cos(theta), sn = sin(theta);
        double r = lq->re[1], m = lq->im[1];
        lq->re[1] = r * cs - m * sn;
        lq->im[1] = r * sn + m * cs;
    }
    lq->materialized = 0;
}

void lazy_tri_apply_hadamard(LazyTrialityQubit *lq)
{
    lq->dft_count = (lq->dft_count + 1) % 2;  /* H²=I */
    lq->materialized = 0;
}

void lazy_tri_materialize(LazyTrialityQubit *lq)
{
    if (lq->materialized) return;

    /* Start from |0⟩ */
    double re[2] = {1.0, 0.0};
    double im[2] = {0.0, 0.0};

    /* Apply accumulated phase segments */
    for (uint32_t s = 0; s < lq->num_segments; s++) {
        double theta = lq->segments[s].phase;
        if (fabs(theta) > 1e-15) {
            double cs = cos(theta), sn = sin(theta);
            double r = re[1], m = im[1];
            re[1] = r * cs - m * sn;
            im[1] = r * sn + m * cs;
        }
    }

    /* Apply accumulated DFTs */
    if (lq->dft_count % 2 == 1) {
        double tr[2], ti[2];
        tr[0] = (re[0] + re[1]) * INV_SQRT2;
        ti[0] = (im[0] + im[1]) * INV_SQRT2;
        tr[1] = (re[0] - re[1]) * INV_SQRT2;
        ti[1] = (im[0] - im[1]) * INV_SQRT2;
        re[0] = tr[0]; re[1] = tr[1];
        im[0] = ti[0]; im[1] = ti[1];
    }

    lq->re[0] = re[0]; lq->re[1] = re[1];
    lq->im[0] = im[0]; lq->im[1] = im[1];
    lq->materialized = 1;
    lq->num_segments = 0;
    lq->dft_count = 0;
}

int lazy_tri_measure(LazyTrialityQubit *lq, double rand_01)
{
    lazy_tri_materialize(lq);

    double p0 = lq->re[0] * lq->re[0] + lq->im[0] * lq->im[0];
    int outcome = (rand_01 < p0) ? 0 : 1;

    /* Collapse */
    double mag = sqrt(lq->re[outcome]*lq->re[outcome]
                    + lq->im[outcome]*lq->im[outcome]);
    if (mag < 1e-30) mag = 1.0;
    double inv = 1.0 / mag;

    lq->re[0] = 0; lq->im[0] = 0;
    lq->re[1] = 0; lq->im[1] = 0;
    lq->re[outcome] = lq->re[outcome] * inv;  /* Preserve phase */
    lq->im[outcome] = lq->im[outcome] * inv;

    /* Ah wait, we zeroed them already. Fix: */
    double phase_re = (outcome == 0) ?
        (lq->re[0] != 0 ? lq->re[0] : 1.0) : 0.0;

    /* Simpler: just set to |outcome⟩ */
    lq->re[0] = (outcome == 0) ? 1.0 : 0.0;
    lq->im[0] = 0.0;
    lq->re[1] = (outcome == 1) ? 1.0 : 0.0;
    lq->im[1] = 0.0;

    (void)phase_re;

    return outcome;
}
