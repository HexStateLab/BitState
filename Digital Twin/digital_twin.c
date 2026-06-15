#include "hpc_qubit_graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <termios.h>
#include <math.h>

#define SYS_Q       8
#define USER_LO     0
#define TWIN_LO     3
#define CTRL_A      6
#define CTRL_B      7
#define N_ROUNDS    6

static struct timeval tv0;
static long timings[N_ROUNDS + 1];
static int n_timings = 0;

static int raw_getchar(void) {
    struct termios old, raw;
    int ch;
    tcgetattr(STDIN_FILENO, &old);
    raw = old;
    raw.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &raw);
    ch = getchar();
    tcsetattr(STDIN_FILENO, TCSANOW, &old);
    return ch;
}

static long record_timing(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    if (tv0.tv_sec == 0) { tv0 = tv; return 0; }
    long dt = (tv.tv_sec - tv0.tv_sec) * 1000000L
            + (tv.tv_usec - tv0.tv_usec);
    tv0 = tv;
    return dt;
}

static void wait_press(const char *msg) {
    printf("  %s", msg); fflush(stdout);
    raw_getchar();
    long dt = record_timing();
    if (n_timings < N_ROUNDS + 1) timings[n_timings++] = dt;
    printf("\r");
}

static int cmp_long(const void *a, const void *b) {
    long la = *(const long *)a, lb = *(const long *)b;
    return (la > lb) - (la < lb);
}

/* 8×8 RDM for 3 qubits via sum over 2^{N-3} complement configurations */
static double trio_entropy(HPCQGraph *g, int qa, int qb, int qc) {
    int comp[8], nc = 0;
    for (int i = 0; i < SYS_Q; i++)
        if (i != qa && i != qb && i != qc)
            comp[nc++] = i;

    uint32_t cfg[SYS_Q];
    double rho[64], re[8], im[8];
    memset(rho, 0, 64 * sizeof(double));
    uint64_t n_cfg = 1ULL << nc;

    for (uint64_t cm = 0; cm < n_cfg; cm++) {
        for (int j = 0; j < nc; j++) cfg[comp[j]] = (cm >> j) & 1;
        for (int va = 0; va <= 1; va++) {
            cfg[qa] = va;
            for (int vb = 0; vb <= 1; vb++) {
                cfg[qb] = vb;
                for (int vc = 0; vc <= 1; vc++) {
                    cfg[qc] = vc;
                    int idx = va*4 + vb*2 + vc;
                    hpcq_amplitude(g, cfg, &re[idx], &im[idx]);
                }
            }
        }
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                rho[i*8+j] += re[i]*re[j] + im[i]*im[j];
    }

    double trace = 0;
    for (int i = 0; i < 8; i++) trace += rho[i*8+i];
    if (trace < 1e-15) return 0;
    double inv = 1.0 / trace;
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            rho[i*8+j] *= inv;
    for (int i = 0; i < 8; i++)
        for (int j = i+1; j < 8; j++) {
            double avg = (rho[i*8+j] + rho[j*8+i]) / 2;
            rho[i*8+j] = rho[j*8+i] = avg;
        }

    for (int iter = 0; iter < 100; iter++) {
        double max_off = 0; int p = 0, q = 0;
        for (int i = 0; i < 8; i++)
            for (int j = i+1; j < 8; j++)
                if (fabs(rho[i*8+j]) > max_off) { max_off = fabs(rho[i*8+j]); p=i; q=j; }
        if (max_off < 1e-14) break;
        double theta = 0.5 * atan2(2*rho[p*8+q], rho[p*8+p]-rho[q*8+q]);
        double c = cos(theta), s = sin(theta);
        for (int i = 0; i < 8; i++) {
            double aip = rho[i*8+p], aiq = rho[i*8+q];
            rho[i*8+p] = c*aip - s*aiq;
            rho[i*8+q] = s*aip + c*aiq;
        }
        for (int j = 0; j < 8; j++) {
            double apj = rho[p*8+j], aqj = rho[q*8+j];
            rho[p*8+j] = c*apj - s*aqj;
            rho[q*8+j] = s*apj + c*aqj;
        }
    }

    double S = 0;
    for (int i = 0; i < 8; i++) {
        double ev = rho[i*8+i];
        if (ev < 0) ev = 0;
        if (ev > 1) ev = 1;
        if (ev > 1e-30) S -= ev * log2(ev);
    }
    return S;
}

/* Compute 3-qubit RDM (rho_out, 8×8) for joint/marginal extraction */
static void trio_rdm(HPCQGraph *g, int qa, int qb, int qc, double *rho) {
    int comp[8], nc = 0;
    for (int i = 0; i < SYS_Q; i++)
        if (i != qa && i != qb && i != qc)
            comp[nc++] = i;

    uint32_t cfg[SYS_Q];
    double re[8], im[8];
    memset(rho, 0, 64 * sizeof(double));
    uint64_t n_cfg = 1ULL << nc;

    for (uint64_t cm = 0; cm < n_cfg; cm++) {
        for (int j = 0; j < nc; j++) cfg[comp[j]] = (cm >> j) & 1;
        for (int va = 0; va <= 1; va++) {
            cfg[qa] = va;
            for (int vb = 0; vb <= 1; vb++) {
                cfg[qb] = vb;
                for (int vc = 0; vc <= 1; vc++) {
                    cfg[qc] = vc;
                    int idx = va*4 + vb*2 + vc;
                    hpcq_amplitude(g, cfg, &re[idx], &im[idx]);
                }
            }
        }
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                rho[i*8+j] += re[i]*re[j] + im[i]*im[j];
    }
    double trace = 0;
    for (int i = 0; i < 8; i++) trace += rho[i*8+i];
    if (trace > 1e-15) {
        double inv = 1.0 / trace;
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                rho[i*8+j] *= inv;
    }
}

/* Marginal P(twin=v) from the 3-qubit RDM (u, t, c) */
static double twin_marginal(double *rho, int user_bit, int twin_val) {
    double p = 0;
    for (int vc = 0; vc <= 1; vc++) {
        int idx = user_bit*4 + twin_val*2 + vc;
        p += rho[idx*8 + idx];
    }
    return p;
}

int main(void) {
    srandom(time(NULL) ^ getpid());

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════╗\n");
    printf("  ║      BitState  Digital Twin                     ║\n");
    printf("  ║  3-qubit entanglement entropy oracle             ║\n");
    printf("  ║  6 Enters → graph,  S(uᵢ,tᵢ,x) = fingerprint   ║\n");
    printf("  ╚══════════════════════════════════════════════════╝\n\n");
    printf("  Each press = CZ(user,twin) + CZ(twin,ctl).  The\n");
    printf("  edge topology encodes timing into the state.\n\n");
    printf("  ──────────────────────────────────────────────\n\n");

    HPCQGraph *g = hpcq_create(SYS_Q);
    for (int i = 0; i < SYS_Q; i++) hpcq_hadamard(g, i);

    printf("  Phase 1:  Press Enter 6×\n\n");

    for (int round = 0; round < N_ROUNDS; round++) {
        char msg[64];
        snprintf(msg, sizeof(msg), "  [%d/6] ⏎  ", round + 1);
        wait_press(msg);

        long dt = timings[n_timings - 1];
        int ui = USER_LO + (round % 3);
        int ti = TWIN_LO + ((round / 3 + round % 3) % 3);

        hpcq_cz(g, ui, ti);

        int e2 = dt & 3;
        const char *desc;
        switch (e2) {
            case 0: {
                int tj = TWIN_LO + ((ti - TWIN_LO + 1) % 3);
                hpcq_cz(g, ti, tj);
                desc = "twin↔twin";
                break;
            }
            case 1: hpcq_cz(g, ti, CTRL_A); desc = "twin↔ctlₐ"; break;
            case 2: hpcq_cz(g, ti, CTRL_B); desc = "twin↔ctl_b"; break;
            case 3:
                hpcq_cz(g, TWIN_LO + ((ti - TWIN_LO + 1) % 3), CTRL_A);
                desc = "twin↔ctlₐ across";
                break;
            default: desc = "";
        }
        printf("∴ CZ(u%d, t%d) + %s   (dt=%ldµs)\n",
               ui, ti - TWIN_LO, desc, dt);
    }

    printf("\n  ────────────────────────────────────────────────\n");
    printf("  Phase 2:  3-qubit entanglement tomography\n\n");

    double S_a[3], S_b[3], S_avg = 0;

    for (int i = 0; i < 3; i++) {
        int u = USER_LO + i;
        int t = TWIN_LO + i;

        S_a[i] = trio_entropy(g, u, t, CTRL_A);
        S_b[i] = trio_entropy(g, u, t, CTRL_B);
        S_avg += (S_a[i] + S_b[i]) / 2;

        printf("  (u%d,t%d,ctlₐ): S = %.4f\n", i, i, S_a[i]);
        printf("  (u%d,t%d,ctl_b): S = %.4f\n", i, i, S_b[i]);
    }
    S_avg /= 3;

    long sorted[N_ROUNDS];
    memcpy(sorted, timings + 1, (N_ROUNDS - 1) * sizeof(long));
    qsort(sorted, N_ROUNDS - 1, sizeof(long), cmp_long);
    long median = sorted[(N_ROUNDS - 1) / 2];

    printf("\n  Timing (µs):");
    for (int i = 1; i < N_ROUNDS; i++)
        printf(" %ld%s", timings[i], i < N_ROUNDS - 1 ? "," : "");
    printf("\n  Median: %ld µs\n\n", median);

    printf("  ────────────────────────────────────────────────\n\n");
    printf("  Phase 3:  Oracle\n\n");

    double deviation = 3.0 - S_avg;

    printf("  Avg 3Q entropy: S = %.4f  (Δ = %.4f below max)\n",
           S_avg, deviation);

    if (deviation > 0.3)
        printf("  🔮  Strong timing signature.\n");
    else if (deviation > 0.1)
        printf("  🔮  Moderate timing signature.\n");
    else
        printf("  🔮  Weak signature — near-random timing.\n");

    printf("\n  Press ⏎ to collapse → ");
    fflush(stdout);
    raw_getchar();
    long dt_test = record_timing();
    printf("\r");

    int user_bit = (int)(dt_test & 1);

    double rho_utc[64];
    trio_rdm(g, USER_LO, TWIN_LO, CTRL_A, rho_utc);

    double p_t0 = twin_marginal(rho_utc, user_bit, 0);
    double p_t1 = twin_marginal(rho_utc, user_bit, 1);
    double r = (double)random() / (double)RAND_MAX;
    int twin_bit = (r < p_t0 / (p_t0 + p_t1 + 1e-30)) ? 0 : 1;
    int agree = (user_bit == twin_bit);

    printf("\n");
    printf("  ╔══════════════════════════════════════════════╗\n");
    printf("  ║  User bᵤ=%d    Twin bₜ=%d    %s  ║\n",
           user_bit, twin_bit, agree ? "AGREE ✓" : "DISAGREE ✗");
    printf("  ╚══════════════════════════════════════════════╝\n");

    printf("\n  Δ = %.4f  |  Outcome: %s\n",
           deviation, agree ? "agree" : "disagree");

    printf("\n  ────────────────────────────────────────────────\n\n");
    hpcq_print_stats(g);

    hpcq_destroy(g);
    return 0;
}
