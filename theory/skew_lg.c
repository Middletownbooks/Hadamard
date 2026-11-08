/*
 * Low-energy list generator for one block of a skew Goethals-Seidel multiplier slice (any abelian group).
 * Written for this repository; the energy (whitened distance of the PAF class vector to the conditional mean of a
 * solution) is the one of K. Kumar's search/cp_lowe.py and search/listgen (Apache-2.0); no code is copied.
 *
 * Variables: U units; unit u flips the signs of its orbits (free block: one orbit per unit; skew block: the pair
 * {O, -O} with opposite signs). PAF at class representative s: P_s = sum_{O,O'} e_O e_O' N_s(O, O').
 * Moves: free block: swap two units of equal weight and opposite sign (keeps the row sum); skew block: flip one unit.
 * Metropolis at temperature T on E = |W (P - c)|^2; every visited state with E <= tau is recorded once (64+64-bit hash
 * of the sign vector) until NKEEP distinct states, or MAXMOVES moves. Restarts from the supplied starting states.
 *
 * input (little endian): int32 r, C, U, skew, nstart; int32 N[C][r][r]; double W[C][C]; double c[C];
 *   int32 unit_of_orbit[r]; int32 sign_in_unit[r] (+1/-1: orbit sign = unit sign * this); int32 weight[U];
 *   int8 starts[nstart][U] (unit signs)
 * usage: skew_lg in.bin out.bin NKEEP TAU T MAXMOVES SEED THREADS RESTART_EVERY [MAXSECS]
 * Stops also after MAXSECS, or when 10^8 consecutive moves add fewer than 0.2% new states (saturation).
 * output: int64 n; int32 C; int16 keys[n][C]; uint64 masks[n][2] (unit sign bits, 1 = negative), then float E[n]
 *   cc -O3 -march=native -fopenmp skew_lg.c -o skew_lg -lm
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int r, C, U, skew, nstart;
static int32_t *N; static double *W, *cc; static int32_t *uo, *su, *wt; static int8_t *starts;
static int32_t *M;   /* M[s][O][O'] = N + N^T */
static int ord[256][256]; /* orbits of a unit */
static int nord[256];

typedef struct { uint64_t a, b; } H2;
static H2 *table; static uint64_t tmask; static long nfound = 0, nkeep;
static int16_t *okeys; static H2 *omask; static float *oE;

static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }

static int insert(H2 h) {           /* 1 if new; thread safe by critical section at caller */
    uint64_t i = mix(h.a ^ mix(h.b)) & tmask;
    while (1) {
        if (table[i].a == 0 && table[i].b == 0) { table[i] = h; return 1; }
        if (table[i].a == h.a && table[i].b == h.b) return 0;
        i = (i + 1) & tmask;
    }
}

static inline double rnd(uint64_t *s) { *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return (*s >> 11) * (1.0 / 9007199254740992.0); }

static double energy(const double *P) {
    double E = 0, d[512];
    for (int i = 0; i < C; i++) d[i] = P[i] - cc[i];
    for (int i = 0; i < C; i++) { double z = 0; for (int j = 0; j < C; j++) z += W[i * C + j] * d[j]; E += z * z; }
    return E;
}

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int32_t hd[5]; if (fread(hd, 4, 5, f) != 5) return 2;
    r = hd[0]; C = hd[1]; U = hd[2]; skew = hd[3]; nstart = hd[4];
    N = malloc(4L * C * r * r); W = malloc(8L * C * C); cc = malloc(8L * C);
    uo = malloc(4 * r); su = malloc(4 * r); wt = malloc(4 * U); starts = malloc((size_t)nstart * U);
    if (fread(N, 4, (size_t)C * r * r, f) != (size_t)C * r * r) return 3;
    if (fread(W, 8, (size_t)C * C, f) != (size_t)C * C) return 3;
    if (fread(cc, 8, C, f) != (size_t)C) return 3;
    if (fread(uo, 4, r, f) != (size_t)r || fread(su, 4, r, f) != (size_t)r || fread(wt, 4, U, f) != (size_t)U) return 3;
    if (fread(starts, 1, (size_t)nstart * U, f) != (size_t)nstart * U) return 3;
    fclose(f);
    nkeep = atol(argv[3]); double tau = atof(argv[4]), T = atof(argv[5]); long maxmoves = atol(argv[6]);
    uint64_t seed = strtoull(argv[7], 0, 10); int threads = atoi(argv[8]); long restart = atol(argv[9]);
    double maxsecs = argc > 10 ? atof(argv[10]) : 1e9; double t0 = omp_get_wtime(); long lastn = 0, lastm = 0;
    M = malloc(4L * C * r * r);
    for (int s = 0; s < C; s++) for (int i = 0; i < r; i++) for (int j = 0; j < r; j++)
        M[((long)s * r + i) * r + j] = N[((long)s * r + i) * r + j] + N[((long)s * r + j) * r + i];
    for (int u = 0; u < U; u++) nord[u] = 0;
    for (int o = 0; o < r; o++) if (uo[o] >= 0) ord[uo[o]][nord[uo[o]]++] = o;
    long tsize = 1; while (tsize < 3 * nkeep) tsize <<= 1; tmask = tsize - 1;
    table = calloc(tsize, sizeof(H2));
    okeys = malloc(2L * C * nkeep); omask = malloc(sizeof(H2) * nkeep); oE = malloc(4L * nkeep);
    long moves_total = 0; int done = 0;
    omp_set_num_threads(threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint64_t rs = mix(seed * 1000003ULL + tid + 1);
        int8_t *e = malloc(r), *us = malloc(U);
        double *g = malloc(8L * C * r), P[512], Pn[512];
        long it = 0; int si = tid % (nstart > 0 ? nstart : 1);
        while (!done) {
            /* (re)start */
            memcpy(us, starts + (long)si * U, U); si = (si + threads) % nstart;
            for (int o = 0; o < r; o++) e[o] = uo[o] < 0 ? (int8_t)su[o] : (int8_t)(us[uo[o]] * su[o]);
            for (int s = 0; s < C; s++) {
                double ps = 0;
                for (int o = 0; o < r; o++) {
                    double gg = 0; const int32_t *Mr = M + ((long)s * r + o) * r;
                    for (int o2 = 0; o2 < r; o2++) gg += e[o2] * Mr[o2];
                    g[(long)s * r + o] = gg;
                    const int32_t *Nr = N + ((long)s * r + o) * r;
                    for (int o2 = 0; o2 < r; o2++) ps += e[o] * e[o2] * Nr[o2];
                }
                P[s] = ps;
            }
            double E = energy(P);
            for (long k = 0; k < restart && !done; k++) {
                int u1 = (int)(rnd(&rs) * U), u2 = -1;
                if (!skew) {
                    u2 = (int)(rnd(&rs) * U);
                    if (u1 == u2 || wt[u1] != wt[u2] || us[u1] == us[u2]) continue;
                }
                /* delta of P for flipping the orbits of u1 (and u2) sequentially, without committing */
                int fl[8], nf = 0;
                for (int i = 0; i < nord[u1]; i++) fl[nf++] = ord[u1][i];
                if (u2 >= 0) for (int i = 0; i < nord[u2]; i++) fl[nf++] = ord[u2][i];
                for (int s = 0; s < C; s++) {
                    double d = 0;
                    for (int a = 0; a < nf; a++) {
                        int o = fl[a]; const int32_t *Mr = M + ((long)s * r + o) * r;
                        double ga = g[(long)s * r + o] - e[o] * Mr[o];
                        /* corrections for orbits already flipped in this move */
                        for (int b = 0; b < a; b++) ga -= 2.0 * e[fl[b]] * Mr[fl[b]];
                        d += -2.0 * e[o] * ga;
                    }
                    Pn[s] = P[s] + d;
                }
                double En = energy(Pn), dE = En - E;
                if (dE <= 0 || rnd(&rs) < exp(-dE / T)) {
                    for (int a = 0; a < nf; a++) {
                        int o = fl[a];
                        for (int s = 0; s < C; s++) {
                            const int32_t *Mc = M + (long)s * r * r;
                            double *gs = g + (long)s * r;
                            for (int o2 = 0; o2 < r; o2++) gs[o2] -= 2.0 * e[o] * Mc[(long)o2 * r + o];
                        }
                        e[o] = -e[o];
                    }
                    us[u1] = -us[u1]; if (u2 >= 0) us[u2] = -us[u2];
                    memcpy(P, Pn, 8L * C); E = En;
                    if (E <= tau) {
                        H2 h = {0, 0};
                        for (int u = 0; u < U; u++) if (us[u] < 0) { if (u < 64) h.a |= 1ULL << u; else h.b |= 1ULL << (u - 64); }
                        h.b |= 1ULL << 63;  /* never all-zero */
                        #pragma omp critical
                        {
                            if (!done && insert(h)) {
                                long i = nfound++;
                                for (int s = 0; s < C; s++) okeys[i * C + s] = (int16_t)lrint(P[s]);
                                omask[i] = h; oE[i] = (float)E;
                                if (nfound >= nkeep) done = 1;
                            }
                        }
                    }
                }
                if (++it % 1000000 == 0) {
                    #pragma omp atomic
                    moves_total += 1000000;
                    if (moves_total >= maxmoves) done = 1;
                    if (tid == 0) {
                        if (omp_get_wtime() - t0 > maxsecs) done = 1;
                        if (moves_total - lastm >= 100000000L) {
                            if (nfound - lastn < nfound / 500 + 1) done = 1;
                            lastn = nfound; lastm = moves_total;
                        }
                    }
                    if (tid == 0 && it % 50000000 == 0) { fprintf(stderr, "moves %ld kept %ld\n", moves_total, nfound); }
                }
            }
        }
    }
    FILE *o = fopen(argv[2], "wb");
    int64_t n = nfound; int32_t Cc = C;
    fwrite(&n, 8, 1, o); fwrite(&Cc, 4, 1, o);
    fwrite(okeys, 2, (size_t)n * C, o); fwrite(omask, sizeof(H2), n, o); fwrite(oE, 4, n, o);
    fclose(o);
    fprintf(stderr, "done: kept %ld moves %ld\n", nfound, moves_total);
    return 0;
}
