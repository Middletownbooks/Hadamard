/*
 * Low-energy list generator, fast version of theory/skew_lg.c (same input/output and command line).
 * Written for this repository; the energy is the one of K. Kumar's search/cp_lowe.py and the Metropolis list
 * generation follows the idea of his search/listgen (Apache-2.0); no code is copied.
 *
 * State: orbit signs e, the whitened vector z = W (P - c) and, for every orbit o, the whitened field
 * h[o] = W g(o), g_s(o) = sum_o' e_o' M_s(o, o'), M = N + N^T. Flipping orbits a_1, ..., a_m in turn changes z by
 *   dz = sum_a -2 e_a (h[a] - e_a WM[a][a] - sum_{b < a} 2 e_b WM[a][b]),   WM[a][b] = W M(a, b)  (a C-vector),
 * so a proposal costs O(C) and a committed move O(r C). The exact integer P is kept alongside (integer field g).
 * usage: skew_lg2 in.bin out.bin NKEEP TAU T MAXMOVES SEED THREADS RESTART_EVERY [MAXSECS]
 *   cc -O3 -march=native -fopenmp skew_lg2.c -o skew_lg2 -lm
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int r, C, U, skew, nstart;
static int32_t *N, *M, *uo, *su, *wt; static double *W, *cc; static int8_t *starts;
static float *WM;              /* [r][r][C] */
static int ord[512][2], nord[512];
typedef struct { uint64_t a, b; } H2;
static H2 *table; static uint64_t tmask; static long nfound = 0, nkeep;
static int16_t *okeys; static H2 *omask; static float *oE;

static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static int insert(H2 h) {
    uint64_t i = mix(h.a ^ mix(h.b)) & tmask;
    while (1) {
        if (table[i].a == 0 && table[i].b == 0) { table[i] = h; return 1; }
        if (table[i].a == h.a && table[i].b == h.b) return 0;
        i = (i + 1) & tmask;
    }
}
static inline double rnd(uint64_t *s) { *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return (*s >> 11) * (1.0 / 9007199254740992.0); }

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage\n"); return 1; }
    FILE *f = fopen(argv[1], "rb"); int32_t hd[5];
    if (fread(hd, 4, 5, f) != 5) return 2;
    r = hd[0]; C = hd[1]; U = hd[2]; skew = hd[3]; nstart = hd[4];
    N = malloc(4L * C * r * r); W = malloc(8L * C * C); cc = malloc(8L * C);
    uo = malloc(4 * r); su = malloc(4 * r); wt = malloc(4 * U); starts = malloc((size_t)nstart * U);
    if (fread(N, 4, (size_t)C * r * r, f) != (size_t)C * r * r || fread(W, 8, (size_t)C * C, f) != (size_t)C * C ||
        fread(cc, 8, C, f) != (size_t)C || fread(uo, 4, r, f) != (size_t)r || fread(su, 4, r, f) != (size_t)r ||
        fread(wt, 4, U, f) != (size_t)U || fread(starts, 1, (size_t)nstart * U, f) != (size_t)nstart * U) return 3;
    fclose(f);
    nkeep = atol(argv[3]); double tau = atof(argv[4]), T = atof(argv[5]); long maxmoves = atol(argv[6]);
    uint64_t seed = strtoull(argv[7], 0, 10); int threads = atoi(argv[8]); long restart = atol(argv[9]);
    double maxsecs = argc > 10 ? atof(argv[10]) : 1e9;
    M = malloc(4L * C * r * r);
    for (int s = 0; s < C; s++) for (int i = 0; i < r; i++) for (int j = 0; j < r; j++)
        M[((long)s * r + i) * r + j] = N[((long)s * r + i) * r + j] + N[((long)s * r + j) * r + i];
    WM = malloc(4L * r * r * C);
    for (int a = 0; a < r; a++) for (int b = 0; b < r; b++) for (int i = 0; i < C; i++) {
        double z = 0; for (int s = 0; s < C; s++) z += W[i * C + s] * M[((long)s * r + a) * r + b];
        WM[((long)a * r + b) * C + i] = (float)z;
    }
    for (int u = 0; u < U; u++) nord[u] = 0;
    for (int o = 0; o < r; o++) if (uo[o] >= 0) ord[uo[o]][nord[uo[o]]++] = o;
    long tsize = 1; while (tsize < 3 * nkeep) tsize <<= 1; tmask = tsize - 1;
    table = calloc(tsize, sizeof(H2));
    okeys = malloc(2L * C * nkeep); omask = malloc(sizeof(H2) * nkeep); oE = malloc(4L * nkeep);
    volatile int done = 0; long moves_total = 0; double t0 = omp_get_wtime(); long lastn = 0, lastm = 0;
    omp_set_num_threads(threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint64_t rs = mix(seed * 1000003ULL + tid + 1);
        int8_t *e = malloc(r), *us = malloc(U);
        int32_t *g = malloc(4L * C * r);          /* integer field g_s(o), layout [o][s] */
        float *h = malloc(4L * r * C);            /* whitened field, layout [o][i] */
        double P[512]; float z[512], dz[512];
        long it = 0; int si = tid % nstart;
        while (!done) {
            memcpy(us, starts + (long)si * U, U); si = (si + threads) % nstart;
            for (int o = 0; o < r; o++) e[o] = uo[o] < 0 ? (int8_t)su[o] : (int8_t)(us[uo[o]] * su[o]);
            for (int s = 0; s < C; s++) {
                long ps = 0;
                for (int o = 0; o < r; o++) {
                    long gg = 0; const int32_t *Mr = M + ((long)s * r + o) * r;
                    for (int o2 = 0; o2 < r; o2++) gg += e[o2] * Mr[o2];
                    g[(long)o * C + s] = (int32_t)gg;
                    const int32_t *Nr = N + ((long)s * r + o) * r;
                    for (int o2 = 0; o2 < r; o2++) ps += e[o] * e[o2] * Nr[o2];
                }
                P[s] = (double)ps;
            }
            for (int o = 0; o < r; o++) for (int i = 0; i < C; i++) {
                double t = 0; for (int s = 0; s < C; s++) t += W[i * C + s] * g[(long)o * C + s];
                h[(long)o * C + i] = (float)t;
            }
            double E = 0;
            for (int i = 0; i < C; i++) { double t = 0; for (int s = 0; s < C; s++) t += W[i * C + s] * (P[s] - cc[s]); z[i] = (float)t; E += t * t; }
            for (long k = 0; k < restart && !done; k++) {
                int u1 = (int)(rnd(&rs) * U), u2 = -1;
                if (!skew) {
                    u2 = (int)(rnd(&rs) * U);
                    if (u1 == u2 || wt[u1] != wt[u2] || us[u1] == us[u2]) continue;
                }
                int fl[4], nf = 0;
                for (int i = 0; i < nord[u1]; i++) fl[nf++] = ord[u1][i];
                if (u2 >= 0) for (int i = 0; i < nord[u2]; i++) fl[nf++] = ord[u2][i];
                for (int i = 0; i < C; i++) dz[i] = 0;
                for (int a = 0; a < nf; a++) {
                    int oa = fl[a]; float ea = e[oa];
                    const float *ha = h + (long)oa * C, *wd = WM + ((long)oa * r + oa) * C;
                    for (int i = 0; i < C; i++) {
                        float t = ha[i] - ea * wd[i];
                        for (int b = 0; b < a; b++) t -= 2.0f * e[fl[b]] * WM[((long)oa * r + fl[b]) * C + i];
                        dz[i] += -2.0f * ea * t;
                    }
                }
                double En = 0; for (int i = 0; i < C; i++) { double t = z[i] + dz[i]; En += t * t; }
                double dE = En - E;
                if (dE <= 0 || rnd(&rs) < exp(-dE / T)) {
                    for (int a = 0; a < nf; a++) {
                        int oa = fl[a]; int ea = e[oa];
                        /* exact integer P: dP_s = -2 e_a (g_s(a) - e_a M_s(a, a)) */
                        for (int s = 0; s < C; s++) P[s] += -2.0 * ea * (g[(long)oa * C + s] - ea * M[((long)s * r + oa) * r + oa]);
                        for (int o = 0; o < r; o++) {
                            const float *wv = WM + ((long)o * r + oa) * C; float *ho = h + (long)o * C;
                            for (int i = 0; i < C; i++) ho[i] -= 2.0f * ea * wv[i];
                            int32_t *go = g + (long)o * C;
                            for (int s = 0; s < C; s++) go[s] -= 2 * ea * M[((long)s * r + o) * r + oa];
                        }
                        e[oa] = -e[oa];
                    }
                    us[u1] = -us[u1]; if (u2 >= 0) us[u2] = -us[u2];
                    for (int i = 0; i < C; i++) z[i] += dz[i];
                    E = En;
                    if (E <= tau) {
                        /* recompute E exactly from P to avoid float drift */
                        double Ex = 0;
                        for (int i = 0; i < C; i++) { double t = 0; for (int s = 0; s < C; s++) t += W[i * C + s] * (P[s] - cc[s]); z[i] = (float)t; Ex += t * t; }
                        E = Ex;
                        if (E <= tau) {
                            H2 hh = {0, 0};
                            for (int u = 0; u < U; u++) if (us[u] < 0) { if (u < 64) hh.a |= 1ULL << u; else hh.b |= 1ULL << (u - 64); }
                            hh.b |= 1ULL << 63;
                            #pragma omp critical
                            {
                                if (!done && insert(hh)) {
                                    long i = nfound++;
                                    for (int s = 0; s < C; s++) okeys[i * C + s] = (int16_t)lrint(P[s]);
                                    omask[i] = hh; oE[i] = (float)E;
                                    if (nfound >= nkeep) done = 1;
                                }
                            }
                        }
                    }
                }
                if (++it % 1000000 == 0) {
                    long mt;
                    #pragma omp atomic capture
                    mt = moves_total += 1000000;
                    if (mt >= maxmoves) done = 1;
                    if (tid == 0) {
                        if (omp_get_wtime() - t0 > maxsecs) done = 1;
                        if (mt - lastm >= 400000000L) {
                            if (nfound - lastn < nfound / 500 + 1) done = 1;
                            lastn = nfound; lastm = mt;
                            fprintf(stderr, "moves %ld kept %ld (%.0fs)\n", mt, nfound, omp_get_wtime() - t0);
                        }
                    }
                }
            }
        }
    }
    FILE *o = fopen(argv[2], "wb");
    int64_t n = nfound; int32_t Cc = C;
    fwrite(&n, 8, 1, o); fwrite(&Cc, 4, 1, o);
    fwrite(okeys, 2, (size_t)n * C, o); fwrite(omask, sizeof(H2), n, o); fwrite(oE, 4, n, o);
    fclose(o);
    fprintf(stderr, "done: kept %ld moves %ld in %.0fs\n", nfound, moves_total, omp_get_wtime() - t0);
    return 0;
}
