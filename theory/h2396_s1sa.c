// h2396_s1sa.c -- stochastic local search for stage 1 of the Paley plug-in Legendre-quadruple slice on Z_n x F_q
// (report/prompt_599_plugin.md, Part A).  Stage 1 is the F_q-compression: four sequences u_i = x0_i + (q-1) e_i on Z_n,
// x0 in {+-1}^n, e in {0,+-1}^n (the defects), with
//     sum_i PAF_{u_i}(d) = 4m + 4 - 4q  (d = 0),   -4q  (d != 0),       m = q n.
// The number of defects E and their types (e_plus: x0 = e, letter +-q; e_minus: x0 = -e, letter +-(q-2)) come from the
// defect identity (q+1) e_plus + (q-3) e_minus = 4(n-1); E is kept fixed by the moves, and the lag-0 term fixes the type
// split.  Cost = sum over all n lags of (sum_i PAF_{u_i}(d) - T(d))^2.  Moves: flip x0 at a position; move a defect
// (with a random sign) to a non-defect position of any block; flip the sign of a defect.  Every reported solution is
// re-checked from scratch.  This is a sampler (not uniform, not exhaustive): it only produces instances.
// usage: h2396_s1sa n q e_plus e_minus nsol seed maxsteps_per_restart [timelimit_s] [T0] [Tfloor] [Tdecay]  (S1SA_VERBOSE=1: costs)
// output: one line per solution: "SOL u_0(0..n-1) | u_1 | u_2 | u_3"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#define MAXN 64
static int n, q, E, ep, em, T[MAXN];
static int x0[4][MAXN], ee[4][MAXN], u[4][MAXN];
static long P[MAXN];                 // P[d] = sum_i PAF_{u_i}(d), d = 0..n-1
static uint64_t rs = 88172645463325252ULL;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static double T0init = 2000.0, Tfloor = 20.0, Tdecay = 0.995; static int verbose = 0;
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static void full_paf(void) {
    for (int d = 0; d < n; d++) {
        long s = 0;
        for (int i = 0; i < 4; i++) for (int k = 0; k < n; k++) s += (long)u[i][k] * u[i][(k + d) % n];
        P[d] = s;
    }
}
static long cost_of(const long *PP) { long c = 0; for (int d = 0; d < n; d++) { long z = PP[d] - T[d]; c += z * z; } return c; }

// change u[i][p] to newv: delta of P over all lags (d and n-d both updated)
static void delta_P(int i, int p, int newv, long *dP) {
    int old = u[i][p], D = newv - old;
    for (int d = 0; d < n; d++) dP[d] = 0;
    if (!D) return;
    dP[0] = (long)newv * newv - (long)old * old;
    for (int d = 1; d < n; d++) {
        int a = u[i][(p + d) % n], b = u[i][(p - d + n) % n];
        // PAF(d) = sum_k u(k) u(k+d): position p enters as k = p (partner p+d) and k = p-d (partner p)
        if ((p + d) % n == p) continue;
        long c = (long)D * (a + b);
        if ((p + d) % n == (p - d + n) % n) c = (long)D * (a + b);      // d = n/2: partner counted twice, a == b
        dP[d] = c;
    }
}

static void set_u(int i, int p, int v, long *dP) { for (int d = 0; d < n; d++) P[d] += dP[d]; u[i][p] = v; }

static void random_start(void) {
    for (int i = 0; i < 4; i++) for (int k = 0; k < n; k++) { x0[i][k] = (rnd() & 1) ? 1 : -1; ee[i][k] = 0; }
    int placed = 0;
    while (placed < E) {
        int i = rnd() % 4, k = rnd() % n;
        if (ee[i][k]) continue;
        int s = (rnd() & 1) ? 1 : -1; ee[i][k] = s;
        x0[i][k] = (placed < ep) ? s : -s;          // start with the right type split
        placed++;
    }
    for (int i = 0; i < 4; i++) for (int k = 0; k < n; k++) u[i][k] = x0[i][k] + (q - 1) * ee[i][k];
    full_paf();
}

static int check_solution(void) {
    long PP[MAXN];
    for (int d = 0; d < n; d++) { long s = 0; for (int i = 0; i < 4; i++) for (int k = 0; k < n; k++) s += (long)u[i][k] * u[i][(k + d) % n]; PP[d] = s; }
    for (int d = 0; d < n; d++) if (PP[d] != T[d]) return 0;
    int cnt = 0, cp = 0;
    for (int i = 0; i < 4; i++) for (int k = 0; k < n; k++) if (ee[i][k]) { cnt++; if (x0[i][k] == ee[i][k]) cp++; }
    return cnt == E && cp == ep;
}

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: h2396_s1sa n q e_plus e_minus nsol seed maxsteps [timelimit]\n"); return 1; }
    n = atoi(argv[1]); q = atoi(argv[2]); ep = atoi(argv[3]); em = atoi(argv[4]); E = ep + em;
    int nsol = atoi(argv[5]); rs ^= (uint64_t)atoll(argv[6]) * 0x9E3779B97F4A7C15ULL; for (int z = 0; z < 20; z++) rnd();
    long maxsteps = atol(argv[7]); double tlim = argc > 8 ? atof(argv[8]) : 1e18;
    if (argc > 9) T0init = atof(argv[9]);
    if (argc > 10) Tfloor = atof(argv[10]);
    if (argc > 11) Tdecay = atof(argv[11]);
    if (getenv("S1SA_VERBOSE")) verbose = 1;
    if ((q + 1) * ep + (q - 3) * em != 4 * (n - 1)) { fprintf(stderr, "allocation violates the defect identity\n"); return 2; }
    long m = (long)q * n;
    for (int d = 0; d < n; d++) T[d] = d == 0 ? (int)(4 * m + 4 - 4 * q) : -4 * q;
    if (getenv("S1SA_PLANT")) {                    // planted control: target = PAF sum of a random configuration
        uint64_t keep = rs; rs ^= (uint64_t)atoll(getenv("S1SA_PLANT")) * 0xD1B54A32D192ED03ULL; for (int z = 0; z < 20; z++) rnd();
        random_start(); for (int d = 0; d < n; d++) T[d] = (int)P[d]; rs = keep;
        fprintf(stderr, "planted target:"); for (int d = 0; d <= n / 2; d++) fprintf(stderr, " %d", T[d]); fprintf(stderr, "\n");
    }
    double t0 = now(); int found = 0; long restarts = 0, totsteps = 0;
    long dP[MAXN], dP2[MAXN], PP[MAXN];
    while (found < nsol && now() - t0 < tlim) {
        random_start(); restarts++;
        long c = cost_of(P), cmin = c;
        double Temp = T0init;
        for (long step = 0; step < maxsteps && c > 0; step++) {
            totsteps++;
            if ((step & 1023) == 0) { Temp *= Tdecay; if (Temp < Tfloor) Temp = Tfloor; }
            int mv = rnd() % 10;
            if (mv < 7) {                              // flip x0 at (i, p)
                int i = rnd() % 4, p = rnd() % n;
                int nx = -x0[i][p], nv = nx + (q - 1) * ee[i][p];
                // keep the type split: flipping x0 at a defect changes its type; allow, lag 0 penalises
                delta_P(i, p, nv, dP);
                for (int d = 0; d < n; d++) PP[d] = P[d] + dP[d];
                long c2 = cost_of(PP);
                if (c2 <= c || urand() < exp((double)(c - c2) / Temp)) { set_u(i, p, nv, dP); x0[i][p] = nx; c = c2; }
                if (c < cmin) cmin = c;
            } else if (mv < 9) {                       // move a defect
                int i, p; do { i = rnd() % 4; p = rnd() % n; } while (!ee[i][p]);
                int i2, p2; do { i2 = rnd() % 4; p2 = rnd() % n; } while (ee[i2][p2]);
                int typ_plus = (x0[i][p] == ee[i][p]);
                int s2 = (rnd() & 1) ? 1 : -1;
                int nx2 = typ_plus ? s2 : -s2;             // keep the type of the moved defect
                int v1 = x0[i][p], v2 = nx2 + (q - 1) * s2;
                long savedP[MAXN]; memcpy(savedP, P, sizeof(long) * n);
                int ou1 = u[i][p], ou2 = u[i2][p2];
                delta_P(i, p, v1, dP); set_u(i, p, v1, dP);
                delta_P(i2, p2, v2, dP2); set_u(i2, p2, v2, dP2);
                long c2 = cost_of(P);
                if (c2 <= c || urand() < exp((double)(c - c2) / Temp)) {
                    ee[i][p] = 0; ee[i2][p2] = s2; x0[i2][p2] = nx2; c = c2;
                } else { memcpy(P, savedP, sizeof(long) * n); u[i][p] = ou1; u[i2][p2] = ou2; }
            } else {                                   // flip the sign of a defect, keeping its type
                int i, p; do { i = rnd() % 4; p = rnd() % n; } while (!ee[i][p]);
                int nv = -u[i][p];
                delta_P(i, p, nv, dP);
                for (int d = 0; d < n; d++) PP[d] = P[d] + dP[d];
                long c2 = cost_of(PP);
                if (c2 <= c || urand() < exp((double)(c - c2) / Temp)) { set_u(i, p, nv, dP); x0[i][p] = -x0[i][p]; ee[i][p] = -ee[i][p]; c = c2; }
            }
        }
        if (verbose) fprintf(stderr, "restart %ld: final cost %ld, min %ld\n", restarts, c, cmin);
        if (c == 0 && check_solution()) {
            found++;
            printf("SOL");
            for (int i = 0; i < 4; i++) { if (i) printf(" |"); for (int k = 0; k < n; k++) printf(" %d", u[i][k]); }
            printf("\n"); fflush(stdout);
        }
    }
    fprintf(stderr, "n=%d q=%d alloc=(%d,%d) found=%d restarts=%ld steps=%ld secs=%.1f\n", n, q, ep, em, found, restarts, totsteps, now() - t0);
    return 0;
}
