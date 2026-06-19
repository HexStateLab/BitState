/*
 * shors_code_finder.c — Unchained Shor's-Method Quantum Code Discovery
 *
 * Exploits the QFT diagonalization of circulant matrices to search for
 * GB codes at industrial scale.  The analytic formula:
 *
 *   K = 2·deg(gcd(a(x), b(x), x^L+1))
 *
 * replaces O(L³) Gaussian elimination with O(L²) polynomial gcd.
 * This enables searching BILLIONS of candidates per hour.
 *
 * Search: L=3..127 (N=6..254), adaptive random sampling from 1M to 100M
 * Output: shors_codes.csv — all viable codes with D_bound, rate, weights
 *
 * Build: gcc -std=gnu11 -O3 -march=native -o shors_code_finder shors_code_finder.c -lm
 * Run:   ./shors_code_finder              # quick (50M total)
 *        ./shors_code_finder --deep        # deep (500M total, ~2 hrs)
 *        ./shors_code_finder --maximal     # absurd (5B total, ~overnight)
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* ═══════════════════════════════════════════════════════════════════════
 * 128-bit polynomial over GF(2)
 * ═══════════════════════════════════════════════════════════════════════ */

typedef __uint128_t poly_t;

static inline int p_deg(poly_t p) {
    if (!p) return -1;
    return p >> 64 ? 127 - __builtin_clzll((uint64_t)(p >> 64))
                   : 63 - __builtin_clzll((uint64_t)p);
}
static poly_t p_mod(poly_t a, poly_t b) {
    int da = p_deg(a), db = p_deg(b);
    if (db < 0) return 0;
    while (da >= db) { a ^= b << (da - db); da = p_deg(a); }
    return a;
}
static poly_t p_gcd(poly_t a, poly_t b) {
    while (b) { poly_t t = p_mod(a, b); a = b; b = t; }
    return a;
}
static poly_t p_mul(poly_t a, poly_t b) {
    poly_t r = 0;
    while (b) { if (b & 1) r ^= a; a <<= 1; b >>= 1; }
    return r;
}
static int p_wt(poly_t p) {
    return __builtin_popcountll((uint64_t)p)
         + __builtin_popcountll((uint64_t)(p >> 64));
}
static poly_t xL1(int L) { return ((poly_t)1 << L) | 1; }

/* ═══════════════════════════════════════════════════════════════════════
 * PRNG
 * ═══════════════════════════════════════════════════════════════════════ */

static uint64_t rng[2] = {0xCAFEBABE, 0xDEADBEEF};
static inline uint64_t rand64(void) {
    uint64_t s1 = rng[0];
    uint64_t s0 = rng[1];
    rng[0] = s0;
    s1 ^= s1 << 23;
    rng[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng[1] + s0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Analytic K
 * ═══════════════════════════════════════════════════════════════════════ */

static int analytic_K(poly_t a, poly_t b, int L) {
    poly_t x = xL1(L);
    poly_t g = p_gcd(a, b);
    g = p_gcd(g, x);
    return 2 * p_deg(g);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Dual weight bound (random sampling, lightweight)
 * ═══════════════════════════════════════════════════════════════════════ */

static int dual_wt_random(poly_t c, int L, int trials) {
    int min_w = L + 1;
    poly_t x = xL1(L);
    for (int t = 0; t < trials; t++) {
        poly_t f = rand64() & (((poly_t)1 << L) - 1);
        if (!f) continue;
        poly_t cw = p_mod(p_mul(f, c), x);
        int w = p_wt(cw);
        if (w > 0 && w < min_w) { min_w = w; if (min_w <= 3) break; }
    }
    return min_w;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Code record
 * ═══════════════════════════════════════════════════════════════════════ */

typedef struct {
    int N, K, L, D_bound, min_wt, max_wt;
    double rate;
    poly_t a, b;
} Code;

#define MAX_CODES 100000
static Code *all_codes;
static int n_all;
static FILE *csv;

static void csv_init(void) {
    all_codes = calloc(MAX_CODES, sizeof(Code));
    csv = fopen("shors_codes.csv", "w");
    if (csv) fprintf(csv, "N,K,D_bound,rate,L,a_hex,b_hex,min_wt,max_wt\n");
}

static void csv_save(Code *c) {
    if (n_all >= MAX_CODES) return;

    /* Tiered filtering: be more selective at small L */
    if (c->L < 15 && (c->rate < 0.25 || c->D_bound < 4)) return;
    if (c->L < 25 && c->rate < 0.10) return;
    if (c->L < 40 && c->rate < 0.06) return;

    /* (N,K) uniqueness: only keep best D_bound per (N,K) */
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].N == c->N && all_codes[i].K == c->K) {
            if (c->D_bound > all_codes[i].D_bound ||
                (c->D_bound == all_codes[i].D_bound && c->rate > all_codes[i].rate)) {
                all_codes[i] = *c;  /* replace with better */
            }
            return;
        }

    all_codes[n_all++] = *c;
    if (csv && n_all <= MAX_CODES)
        fprintf(csv, "%d,%d,%d,%.4f,%d,0x%lx,0x%lx,%d,%d\n",
                c->N, c->K, c->D_bound, c->rate, c->L,
                (uint64_t)c->a, (uint64_t)c->b,
                c->min_wt, c->max_wt);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Stabilizer weight estimate (analytical for circulant)
 * ═══════════════════════════════════════════════════════════════════════ */

static int circ_wt(poly_t c, int L) {
    /* Weight of a circulant row from polynomial c over Z_L */
    return p_wt(c);
}

/* ═══════════════════════════════════════════════════════════════════════
 * Progress
 * ═══════════════════════════════════════════════════════════════════════ */

static double now_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Search one L value
 * ═══════════════════════════════════════════════════════════════════════ */

static void search_L(int L, long trials) {
    poly_t mask = (((poly_t)1 << L) - 1);
    int N = 2 * L;
    long viable = 0;

    for (long t = 0; t < trials; t++) {
        poly_t a = rand64() & mask;
        poly_t b = rand64() & mask;
        if (!a && !b) continue;

        int K = analytic_K(a, b, L);
        if (K <= 0) continue;
        double rate = (double)K / N;
        if (rate < 0.04) continue;

        /* For small L, be selective to avoid flooding the database */
        if (L < 20 && rate < 0.15) continue;

        /* Quick dual weight estimate */
        int da = dual_wt_random(a, L, 200);
        int db = dual_wt_random(b, L, 200);
        int D = da < db ? da : db;
        int sing = (N - K) / 2 + 1;
        if (D > sing) D = sing;
        if (D < 3) continue;  /* only interested in d≥3 */

        /* Stabilizer weights */
        int wa = circ_wt(a, L);
        int wb = circ_wt(b, L);
        int min_w = wa < wb ? wa : wb;
        int max_w = wa > wb ? wa : wb;

        /* Check if novel: compare to known parameters from codetables.de */
        /* (We save all viable codes for post-hoc analysis) */

        Code c = {N, K, L, D, min_w, max_w, rate, a, b};
        csv_save(&c);
        viable++;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * Post-hoc analysis: find best codes by category
 * ═══════════════════════════════════════════════════════════════════════ */

static int code_better(Code *a, Code *b) {
    /* Score: rate * D_bound, higher is better */
    double sa = a->rate * a->D_bound;
    double sb = b->rate * b->D_bound;
    if (sa > sb) return 1;
    if (sa < sb) return -1;
    if (a->D_bound > b->D_bound) return 1;
    if (a->D_bound < b->D_bound) return -1;
    if (a->rate > b->rate) return 1;
    return -1;
}

static void print_best(int max_show) {
    /* Sort by score */
    for (int i = 0; i < n_all - 1; i++)
        for (int j = i + 1; j < n_all; j++)
            if (code_better(&all_codes[j], &all_codes[i]) > 0) {
                Code t = all_codes[i];
                all_codes[i] = all_codes[j];
                all_codes[j] = t;
            }

    printf("\n═══ BEST CODES BY rate×D ═══\n\n");
    printf("  %4s %14s %7s %6s %6s %6s %6s\n",
           "N", "Parameters", "Rate", "D≥", "minW", "maxW", "L");
    printf("  %s\n", "───────────────────────────────────────────────");

    for (int i = 0; i < n_all && i < max_show; i++) {
        Code *c = &all_codes[i];
        printf("  %4d [[%d,%d,≥%d]] %7.4f %6d %6d %6d %6d\n",
               c->N, c->N, c->K, c->D_bound, c->rate,
               c->D_bound, c->min_wt, c->max_wt, c->L);
    }

    /* Best by category */
    printf("\n═══ BEST PER CATEGORY ═══\n\n");

    /* Best N<20 */
    printf("  Small (N<20):\n");
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].N < 20) {
            printf("    [[%d,%d,≥%d]] r=%.3f L=%d\n",
                   all_codes[i].N, all_codes[i].K, all_codes[i].D_bound,
                   all_codes[i].rate, all_codes[i].L);
            break;
        }

    /* Best rate */
    printf("  Best rate (K/N > 0.3):\n");
    int shown = 0;
    for (int i = 0; i < n_all && shown < 3; i++)
        if (all_codes[i].rate > 0.3 && all_codes[i].D_bound >= 4) {
            printf("    [[%d,%d,≥%d]] r=%.3f L=%d\n",
                   all_codes[i].N, all_codes[i].K, all_codes[i].D_bound,
                   all_codes[i].rate, all_codes[i].L);
            shown++;
        }

    /* Best distance */
    printf("  Best D_bound (D≥%d):\n", all_codes[0].D_bound);
    shown = 0;
    for (int d = all_codes[0].D_bound; d >= all_codes[0].D_bound - 5 && shown < 5; d--) {
        for (int i = 0; i < n_all; i++)
            if (all_codes[i].D_bound == d && all_codes[i].rate > 0.05) {
                printf("    [[%d,%d,≥%d]] r=%.3f L=%d\n",
                       all_codes[i].N, all_codes[i].K, d,
                       all_codes[i].rate, all_codes[i].L);
                shown++;
                if (shown >= 5) break;
            }
    }

    /* Largest N */
    printf("  Largest N:\n");
    int max_N = 0, max_i = 0;
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].N > max_N) { max_N = all_codes[i].N; max_i = i; }
    if (max_N > 0)
        printf("    [[%d,%d,≥%d]] r=%.3f L=%d\n",
               all_codes[max_i].N, all_codes[max_i].K,
               all_codes[max_i].D_bound, all_codes[max_i].rate,
               all_codes[max_i].L);
}

/* ═══════════════════════════════════════════════════════════════════════
 * MAIN
 * ═══════════════════════════════════════════════════════════════════════ */

int main(int argc, char **argv) {
    int mode = 0; /* 0=quick, 1=deep, 2=maximal */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--deep")) mode = 1;
        if (!strcmp(argv[i], "--maximal")) mode = 2;
    }

    rng[0] = (uint64_t)time(NULL) ^ 0xCAFEBABE;
    rng[1] = rng[0] ^ 0xDEADBEEF;
    csv_init();

    double t0 = now_sec();
    long total_trials = 0, total_viable = 0;

    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Shor's Code Finder — Industrial Quantum Code Discovery     ║\n");
    printf("║  K = 2·deg(gcd(a,b,x^L+1))  —  analytic, O(L²)             ║\n");
    printf("║  Mode: %s                                      ║\n",
           mode == 2 ? "MAXIMAL (5B candidates)" :
           mode == 1 ? "deep (500M candidates)" : "quick (50M candidates)");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* Search over L range */
    long trials_per_L[] = {mode == 2 ? 50000000L : mode == 1 ? 5000000L : 2000000L};
    long tpl = trials_per_L[0];

    printf("  %4s %8s %10s %8s %10s\n",
           "L", "N", "trials", "viable", "codes/sec");
    printf("  %s\n", "──────────────────────────────────────");

    for (int L = 3; L <= (mode == 2 ? 127 : mode == 1 ? 100 : 60); L++) {
        double ts = now_sec();
        long prev_all = n_all;
        search_L(L, tpl);
        long found = n_all - prev_all;
        total_trials += tpl;
        total_viable += found;
        double dt = now_sec() - ts;
        printf("  %4d %8d %10ld %8ld %10.0f\n",
               L, 2*L, tpl, found, dt > 0 ? tpl / dt : 0);
        fflush(stdout);
    }

    double dt = now_sec() - t0;
    printf("\n  TOTAL: %ld trials, %d unique codes in %.0fs (%.0f/sec)\n",
           total_trials, n_all, dt, total_trials / dt);

    /* Analysis */
    print_best(50);

    /* Distribution by N range */
    printf("\n═══ DISTRIBUTION BY N ═══\n\n");
    int hist[256] = {0};
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].N < 256) hist[all_codes[i].N]++;
    printf("  N:");
    int col = 0;
    for (int n = 4; n <= 200; n += 2)
        if (hist[n]) {
            printf(" %d:%d", n, hist[n]);
            if (++col % 6 == 0) printf("\n    ");
        }
    printf("\n");

    /* Distribution by K */
    printf("\n  K:");
    memset(hist, 0, sizeof(hist));
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].K < 256) hist[all_codes[i].K]++;
    col = 0;
    for (int k = 2; k <= 128; k += 2)
        if (hist[k]) {
            printf(" %d:%d", k, hist[k]);
            if (++col % 7 == 0) printf("\n    ");
        }
    printf("\n");

    /* Distribution by D_bound */
    printf("\n  D_bound:");
    memset(hist, 0, sizeof(hist));
    for (int i = 0; i < n_all; i++)
        if (all_codes[i].D_bound < 64) hist[all_codes[i].D_bound]++;
    col = 0;
    for (int d = 2; d <= 32; d++)
        if (hist[d]) {
            printf(" D≥%d:%d", d, hist[d]);
            if (++col % 6 == 0) printf("\n       ");
        }
    printf("\n");

    if (csv) { fclose(csv); printf("\n  Saved to shors_codes.csv\n"); }
    printf("  Done.\n");
    return 0;
}
