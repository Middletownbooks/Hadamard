/*
 * Exhaustive low-energy enumeration of one block of a skew Goethals-Seidel multiplier slice (any abelian group).
 * Written for this repository; the filter is the energy of K. Kumar's search/cp_enum.c / cp_lowe.py (Apache-2.0),
 * |W (P - c)|^2 <= tau, but no code is copied: here the units (orbits, or skew pairs {O, -O}) are cut into G groups,
 * the PAF class vector is  P = K0 + sum_g A_g[p_g] + sum_{g<h} B_gh[p_g, p_h]  (exact integer tables), and the
 * whitened vector W P is additive as well (float tables), so every member costs G + G(G-1)/2 vector additions.
 * Members of the last group are bucketed by row sum, so only members with the target row sum are visited.
 *
 * input: int32 r, C, U, G, skew, K; int32 N[C][r][r]; double W[C][C]; double c[C];
 *        int32 unit_of_orbit[r] (-1 = fixed); int32 sign_in_unit[r] (fixed orbits: their sign);
 *        int32 lin[r][K] (linear constraint coefficients of each orbit); int32 target[K]; int32 group_of_unit[U]
 *   The members kept satisfy sum_o e_o lin[o][.] = target exactly (row sum: K = 1, lin = orbit sizes; a stage-1 class:
 *   lin = the contribution of each orbit to the compressed coordinates).
 * usage: skew_enum_lowe in.bin out.bin TAU MAXOUT THREADS
 * output: int64 n; int64 visited; int32 C; int16 keys[n][C]; uint64 masks[n][2]; float E[n]
 *   (stops writing at MAXOUT and reports overflow by visited < 0)
 *   cc -O3 -march=native -fopenmp skew_enum_lowe.c -o skew_enum_lowe -lm
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CP 32   /* padded class dimension (C <= 32) */
static int r, C, U, G, skew, K;
#define KMAX 24
static int32_t *lin, tgtv[KMAX];
static int32_t *N, *uo, *su, *gu; static double *W, *cc;
static int gsz[16], gunits[16][24];
static int npat[16];
static int32_t *Ai[16];  /* [npat][CP] int */
static float *Af[16];
static int32_t *Bi[16][16]; static float *Bf[16][16];   /* [npat_g][npat_h][CP] */
static int32_t *rsg[16];  /* [npat][K] */
static uint64_t *rsh[16]; static uint64_t Rk[KMAX];
static int32_t K0i[CP]; static float K0f[CP];

static uint64_t *cmp_h; static int cmpi(const void *a, const void *b) { uint64_t x = cmp_h[*(const int *)a], y = cmp_h[*(const int *)b]; return x < y ? -1 : x > y; }
static int orb_sign(int o, int g, int p) {   /* sign of orbit o under pattern p of group g (o's unit in g) */
    int u = uo[o]; if (u < 0) return su[o];
    int k = -1; for (int i = 0; i < gsz[g]; i++) if (gunits[g][i] == u) k = i;
    if (k < 0) return 0;
    return ((p >> k) & 1) ? -su[o] : su[o];
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); int32_t hd[6];
    if (fread(hd, 4, 6, f) != 6) return 2;
    r = hd[0]; C = hd[1]; U = hd[2]; G = hd[3]; skew = hd[4]; K = hd[5];
    if (K > KMAX) return 2;
    if (C > CP || G > 16) { fprintf(stderr, "C or G too large\n"); return 2; }
    N = malloc(4L * C * r * r); W = malloc(8L * C * C); cc = malloc(8L * C);
    uo = malloc(4 * r); su = malloc(4 * r); gu = malloc(4 * U); lin = malloc(4L * r * K);
    if (fread(N, 4, (size_t)C * r * r, f) != (size_t)C * r * r || fread(W, 8, (size_t)C * C, f) != (size_t)C * C ||
        fread(cc, 8, C, f) != (size_t)C || fread(uo, 4, r, f) != (size_t)r || fread(su, 4, r, f) != (size_t)r ||
        fread(lin, 4, (size_t)r * K, f) != (size_t)r * K || fread(tgtv, 4, K, f) != (size_t)K || fread(gu, 4, U, f) != (size_t)U) return 3;
    for (int j = 0; j < K; j++) { uint64_t x = 0x9E3779B97F4A7C15ULL * (j + 1); x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29; Rk[j] = x | 1; }
    fclose(f);
    double tau = atof(argv[3]); long maxout = atol(argv[4]); int threads = atoi(argv[5]);
    for (int g = 0; g < G; g++) gsz[g] = 0;
    for (int u = 0; u < U; u++) { int g = gu[u]; gunits[g][gsz[g]++] = u; }
    for (int g = 0; g < G; g++) { if (gsz[g] > 20) { fprintf(stderr, "group too large\n"); return 2; } npat[g] = 1 << gsz[g]; }
    /* orbits of each group */
    int *og = malloc(4 * r); for (int o = 0; o < r; o++) og[o] = uo[o] < 0 ? -1 : gu[uo[o]];
    /* constant part: fixed orbits */
    for (int s = 0; s < C; s++) {
        long k = 0;
        for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) if (og[i] < 0 && og[j] < 0) k += (long)su[i] * su[j] * N[((long)s * r + i) * r + j];
        K0i[s] = (int32_t)k;
    }
    int32_t rs0[KMAX] = {0}; for (int o = 0; o < r; o++) if (og[o] < 0) for (int j = 0; j < K; j++) rs0[j] += su[o] * lin[o * K + j];
    /* self tables (including cross terms with the fixed orbits) */
    for (int g = 0; g < G; g++) {
        Ai[g] = calloc((size_t)npat[g] * CP, 4); rsg[g] = calloc((size_t)npat[g] * K, 4); rsh[g] = calloc(npat[g], 8);
        for (int p = 0; p < npat[g]; p++) {
            for (int o = 0; o < r; o++) if (og[o] == g) { int e_ = orb_sign(o, g, p); for (int j = 0; j < K; j++) rsg[g][(long)p * K + j] += e_ * lin[o * K + j]; }
            uint64_t hh = 0; for (int j = 0; j < K; j++) hh += Rk[j] * (uint64_t)(int64_t)rsg[g][(long)p * K + j]; rsh[g][p] = hh;
            for (int s = 0; s < C; s++) {
                long k = 0;
                for (int i = 0; i < r; i++) {
                    if (og[i] != g) continue; int ei = orb_sign(i, g, p);
                    for (int j = 0; j < r; j++) {
                        long nij = N[((long)s * r + i) * r + j], nji = N[((long)s * r + j) * r + i];
                        if (og[j] == g) k += (long)ei * orb_sign(j, g, p) * nij;
                        else if (og[j] < 0) k += (long)ei * su[j] * (nij + nji);
                    }
                }
                Ai[g][(long)p * CP + s] = (int32_t)k;
            }
        }
    }
    for (int g = 0; g < G; g++) for (int h = g + 1; h < G; h++) {
        Bi[g][h] = calloc((size_t)npat[g] * npat[h] * CP, 4);
        /* bilinear: B[p][q][s] = sum_{i in g, j in h} e_i(p) e_j(q) (N_ij + N_ji); compute via per-orbit sign vectors */
        int ng = 0, nh = 0, oi[64], oj[64];
        for (int o = 0; o < r; o++) { if (og[o] == g) oi[ng++] = o; if (og[o] == h) oj[nh++] = o; }
        int8_t *Eg = malloc((size_t)npat[g] * ng), *Eh = malloc((size_t)npat[h] * nh);
        for (int p = 0; p < npat[g]; p++) for (int a = 0; a < ng; a++) Eg[(long)p * ng + a] = orb_sign(oi[a], g, p);
        for (int q = 0; q < npat[h]; q++) for (int b = 0; b < nh; b++) Eh[(long)q * nh + b] = orb_sign(oj[b], h, q);
        int32_t *Mx = malloc(4L * C * ng * nh);
        for (int s = 0; s < C; s++) for (int a = 0; a < ng; a++) for (int b = 0; b < nh; b++)
            Mx[((long)s * ng + a) * nh + b] = N[((long)s * r + oi[a]) * r + oj[b]] + N[((long)s * r + oj[b]) * r + oi[a]];
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for (int p = 0; p < npat[g]; p++) {
            int32_t tmp[CP * 64];
            for (int s = 0; s < C; s++) for (int b = 0; b < nh; b++) {
                int32_t t = 0; for (int a = 0; a < ng; a++) t += Eg[(long)p * ng + a] * Mx[((long)s * ng + a) * nh + b];
                tmp[s * 64 + b] = t;
            }
            for (int q = 0; q < npat[h]; q++) for (int s = 0; s < C; s++) {
                int32_t t = 0; for (int b = 0; b < nh; b++) t += tmp[s * 64 + b] * Eh[(long)q * nh + b];
                Bi[g][h][((long)p * npat[h] + q) * CP + s] = t;
            }
        }
        free(Eg); free(Eh); free(Mx);
    }
    /* whitened float tables: W (x) for each table; constant gets -W c */
    float Wf[CP * CP]; for (int i = 0; i < C; i++) for (int j = 0; j < C; j++) Wf[i * CP + j] = (float)W[i * C + j];
    #define WMUL(dst, src) do { for (int i_ = 0; i_ < CP; i_++) { float z_ = 0; if (i_ < C) for (int j_ = 0; j_ < C; j_++) z_ += Wf[i_ * CP + j_] * (float)(src)[j_]; (dst)[i_] = z_; } } while (0)
    { double t[CP]; for (int s = 0; s < C; s++) t[s] = K0i[s] - cc[s];
      for (int i = 0; i < CP; i++) { double z = 0; if (i < C) for (int j = 0; j < C; j++) z += W[i * C + j] * t[j]; K0f[i] = (float)z; } }
    for (int g = 0; g < G; g++) { Af[g] = malloc(4L * npat[g] * CP); for (long p = 0; p < npat[g]; p++) WMUL(Af[g] + p * CP, Ai[g] + p * CP); }
    for (int g = 0; g < G; g++) for (int h = g + 1; h < G; h++) {
        long n2 = (long)npat[g] * npat[h]; Bf[g][h] = malloc(4L * n2 * CP);
        #pragma omp parallel for num_threads(threads)
        for (long t = 0; t < n2; t++) WMUL(Bf[g][h] + t * CP, Bi[g][h] + t * CP);
    }
    /* per-coordinate suffix bounds, and the last group's patterns sorted by the hash of their constraint vector */
    int L = G - 1;
    static int32_t mn[17][KMAX], mx[17][KMAX], sufmin[18][KMAX], sufmax[18][KMAX];
    for (int g = 0; g < G; g++) for (int j = 0; j < K; j++) {
        mn[g][j] = 1 << 30; mx[g][j] = -(1 << 30);
        for (int p = 0; p < npat[g]; p++) { int x = rsg[g][(long)p * K + j]; if (x < mn[g][j]) mn[g][j] = x; if (x > mx[g][j]) mx[g][j] = x; }
    }
    for (int j = 0; j < K; j++) { sufmin[G][j] = sufmax[G][j] = 0; }
    for (int g = G - 1; g >= 0; g--) for (int j = 0; j < K; j++) { sufmin[g][j] = sufmin[g + 1][j] + mn[g][j]; sufmax[g][j] = sufmax[g + 1][j] + mx[g][j]; }
    int *bidx = malloc(4 * npat[L]); for (int p = 0; p < npat[L]; p++) bidx[p] = p;
    cmp_h = rsh[L]; qsort(bidx, npat[L], sizeof(int), cmpi);
    uint64_t *bh = malloc(8 * npat[L]); for (int p = 0; p < npat[L]; p++) bh[p] = rsh[L][bidx[p]];
    long outn = 0, visited = 0; int overflow = 0;
    int16_t *okeys = malloc(2L * C * maxout); uint64_t *omask = malloc(16L * maxout); float *oE = malloc(4L * maxout);
    long total_outer = npat[0];
    #pragma omp parallel num_threads(threads) reduction(+:visited)
    {
        int p[16]; float base[16][CP];
        #pragma omp for schedule(dynamic, 1)
        for (long p0 = 0; p0 < total_outer; p0++) {
            if (overflow) continue;
            p[0] = (int)p0;
            /* iterative nested loops over groups 1..L-1 with partial whitened sums base[g] */
            for (int i = 0; i < CP; i++) base[0][i] = K0f[i] + Af[0][(long)p0 * CP + i];
            int depth = 1; p[1] = -1;
            int32_t rsp[17][KMAX]; for (int j = 0; j < K; j++) rsp[0][j] = rs0[j] + rsg[0][(long)p0 * K + j];
            { int ok = 1; for (int j = 0; j < K; j++) { int lo = rsp[0][j] + sufmin[1][j], hi = rsp[0][j] + sufmax[1][j]; if (tgtv[j] < lo || tgtv[j] > hi) ok = 0; } if (!ok) continue; }
            while (depth >= 1) {
                if (depth == L) {
                    int32_t need[KMAX]; uint64_t hn = 0;
                    for (int j = 0; j < K; j++) { need[j] = tgtv[j] - rsp[L - 1][j]; hn += Rk[j] * (uint64_t)(int64_t)need[j]; }
                    long lo = 0, hi = npat[L];
                    while (lo < hi) { long mid = (lo + hi) >> 1; if (bh[mid] < hn) lo = mid + 1; else hi = mid; }
                    int b0 = (int)lo, b1 = b0; while (b1 < npat[L] && bh[b1] == hn) b1++;
                    for (int bk = b0; bk < b1; bk++) {
                        int q = bidx[bk];
                        { int ok = 1; for (int j = 0; j < K; j++) if (rsg[L][(long)q * K + j] != need[j]) ok = 0; if (!ok) continue; }
                        /* energy in chunks of 8 coordinates with early exit (the partial sum is a lower bound) */
                        float E = 0; int alive = 1;
                        for (int c0 = 0; c0 < CP && alive; c0 += 8) {
                            float t[8];
                            for (int i = 0; i < 8; i++) t[i] = base[L - 1][c0 + i] + Af[L][(long)q * CP + c0 + i];
                            for (int g = 0; g < L; g++) { const float *bb = Bf[g][L] + ((long)p[g] * npat[L] + q) * CP + c0; for (int i = 0; i < 8; i++) t[i] += bb[i]; }
                            for (int i = 0; i < 8; i++) E += t[i] * t[i];
                            if (E > tau) alive = 0;
                        }
                        visited++;
                        if (E <= tau) {
                            int32_t key[CP]; for (int s = 0; s < C; s++) key[s] = K0i[s];
                            p[L] = q;
                            for (int g = 0; g <= L; g++) for (int s = 0; s < C; s++) key[s] += Ai[g][(long)p[g] * CP + s];
                            for (int g = 0; g <= L; g++) for (int h = g + 1; h <= L; h++) for (int s = 0; s < C; s++) key[s] += Bi[g][h][((long)p[g] * npat[h] + p[h]) * CP + s];
                            uint64_t m0 = 0, m1 = 0;
                            for (int g = 0; g <= L; g++) for (int k = 0; k < gsz[g]; k++) if ((p[g] >> k) & 1) { int u = gunits[g][k]; if (u < 64) m0 |= 1ULL << u; else m1 |= 1ULL << (u - 64); }
                            long i;
                            #pragma omp atomic capture
                            i = outn++;
                            if (i < maxout) { for (int s = 0; s < C; s++) okeys[i * C + s] = (int16_t)key[s]; omask[2 * i] = m0; omask[2 * i + 1] = m1; oE[i] = E; }
                            else overflow = 1;
                        }
                    }
                    depth--; continue;
                }
                p[depth]++;
                if (p[depth] >= npat[depth]) { p[depth] = -1; depth--; continue; }
                { int ok = 1;
                  for (int j = 0; j < K; j++) { int x = rsp[depth - 1][j] + rsg[depth][(long)p[depth] * K + j]; rsp[depth][j] = x;
                      if (tgtv[j] < x + sufmin[depth + 1][j] || tgtv[j] > x + sufmax[depth + 1][j]) ok = 0; }
                  if (!ok) continue; }
                for (int i = 0; i < CP; i++) {
                    float z = base[depth - 1][i] + Af[depth][(long)p[depth] * CP + i];
                    for (int g = 0; g < depth; g++) z += Bf[g][depth][((long)p[g] * npat[depth] + p[depth]) * CP + i];
                    base[depth][i] = z;
                }
                depth++; if (depth < L) p[depth] = -1;
            }
        }
    }
    long n = outn < maxout ? outn : maxout;
    FILE *o = fopen(argv[2], "wb"); int64_t nn = n, vv = overflow ? -visited : visited; int32_t Cc = C;
    fwrite(&nn, 8, 1, o); fwrite(&vv, 8, 1, o); fwrite(&Cc, 4, 1, o);
    fwrite(okeys, 2, (size_t)n * C, o); fwrite(omask, 16, n, o); fwrite(oE, 4, n, o); fclose(o);
    fprintf(stderr, "enum: visited %ld kept %ld%s\n", visited, n, overflow ? " (OVERFLOW)" : "");
    return 0;
}
