// h2396_stage2.c -- exact stage-2 count of the Paley plug-in Legendre-quadruple slice on Z_n x F_q, n even
// (report/prompt_599_plugin.md, Part A, gate A2).
// Given a stage-1 solution u_0..u_3 on Z_n (u = x0 + (q-1) e; letters +-1 clean, +-q (x0 = e), +-(q-2) (x0 = -e)),
// s = x0 - e and the clean sets are determined.  A completion is t_i in {+-1} on the clean positions of block i with
//   (I')  sum_i PAF_{t_i}(h) = T2(h) := -(1/q) sum_i PAF_{s_i}(h)      for h != 0,
//   (Q)   sum_i [CCF(s_i, t_i) + CCF(t_i, s_i)](h) = 0                 for all h,
// which together with stage 1 is the whole system (I), (J), (Q) [P: theory/h2396_plugin_system.py].
// Method (exact): split along L = {0, k}, k = n/2.  S_i(j) = t_i(j) + t_i(j + k) on Z_k; the S-level is the image of
// (I'), (Q) under the compression Z[Z_n] -> Z[Z_k] (a ring map): sum PAF_{S} = T_S, sum [CCF(sig, S) + CCF(S, sig)] = 0
// with sig(j) = s(j) + s(j + k).  S-level: per-block lists (PSD-filtered: |S^(phi)|^2 <= T_S^(phi)), 2+2 join by a hash
// of the key.  For every S-solution, the D-check enumerates the remaining signs (pairs with S = 0, both clean) of each
// block and joins 2+2 on the full K-level key (PAF_t(h), h = 1..n/2; Q_t(h), h = 0..n/2).  Every completion found is
// re-checked against (I') and (Q) from scratch.
// usage: h2396_stage2 n q < instance        (instance: 4 lines of n integers u_i(0..n-1))
//        prints "COUNT N  sols_S M  ..." ; env ST2_PRINT=K prints up to K completions as "T t_0 | t_1 | t_2 | t_3"
//        planted control: after the 4 lines, "PLANT" then T2(1..n-1) and the Q target (0..n-1); u then only fixes s and
//        the clean sets (letters as above).  (Q) carries eps = chi(-1): CCF(s,t) + eps CCF(t,s).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAXN 64
static int n, q, k, hS, hK, eps = 1, planted = 0;
static long QT[MAXN];               // Q target per lag (0 unless planted)
static int u[4][MAXN], s[4][MAXN], cl[4][MAXN];
static long T2[MAXN];
static int sig[4][MAXN];
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

// ---------------- S-level lists
typedef struct { int nk; long n, cap; int8_t *S; int32_t *key; } List;   // S: n x k values; key: n x nk
static List LS[4];
static int nkS;                       // key length at S-level: PAF_S(0..hS) and Q_S(0..hS)
static long TS[MAXN];                 // S-level PAF target on Z_k
static double TShat[MAXN];

static void list_push(List *L, const int8_t *S, const int32_t *key) {
    if (L->n == L->cap) {
        L->cap = L->cap ? 2 * L->cap : 1024;
        L->S = realloc(L->S, (size_t)L->cap * k); L->key = realloc(L->key, sizeof(int32_t) * (size_t)L->cap * L->nk);
        if (!L->S || !L->key) { fprintf(stderr, "oom list\n"); exit(3); }
    }
    memcpy(L->S + (size_t)L->n * k, S, k); memcpy(L->key + (size_t)L->n * L->nk, key, sizeof(int32_t) * L->nk); L->n++;
}

static int curS[MAXN];
static void s_key(int b, const int *S, int32_t *key) {
    for (int d = 0; d <= hS; d++) {
        long p = 0, qq = 0;
        for (int j = 0; j < k; j++) { p += S[j] * S[(j + d) % k]; qq += (long)sig[b][j] * S[(j + d) % k] + (long)eps * S[j] * sig[b][(j + d) % k]; }
        key[d] = (int32_t)p; key[hS + 1 + d] = (int32_t)qq;
    }
}
static int psd_ok(const int *S) {
    for (int l = 1; l <= k / 2; l++) {
        double re = 0, im = 0;
        for (int j = 0; j < k; j++) { re += S[j] * cos(2 * M_PI * l * j / k); im += S[j] * sin(2 * M_PI * l * j / k); }
        if (re * re + im * im > TShat[l] + 1e-6) return 0;
    }
    long r = 0; for (int j = 0; j < k; j++) r += S[j];
    if ((double)r * r > TShat[0] + 1e-6) return 0;
    return 1;
}
static void gen_S(int b, int j) {
    if (j == k) {
        if (!psd_ok(curS)) return;
        int8_t S8[MAXN]; int32_t key[2 * MAXN];
        for (int z = 0; z < k; z++) S8[z] = (int8_t)curS[z];
        s_key(b, curS, key); list_push(&LS[b], S8, key); return;
    }
    int c0 = cl[b][j], c1 = cl[b][j + k];
    if (c0 && c1) { int opt[3] = {-2, 0, 2}; for (int o = 0; o < 3; o++) { curS[j] = opt[o]; gen_S(b, j + 1); } }
    else if (c0 || c1) { curS[j] = 1; gen_S(b, j + 1); curS[j] = -1; gen_S(b, j + 1); }
    else { curS[j] = 0; gen_S(b, j + 1); }
}

// ---------------- hash of S-level pair keys (open addressing on 64-bit fingerprints)
static uint64_t fp_key(const int32_t *a, const int32_t *b, int nk, int sign, const long *T) {
    // fingerprint of (a + b) if sign = +1, of (T - a - b) if sign = -1
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < nk; i++) {
        long v = sign > 0 ? (long)a[i] + b[i] : T[i] - a[i] - b[i];
        h ^= (uint64_t)(v + 0x40000000L); h *= 1099511628211ULL; h ^= h >> 29;
    }
    return h | 1;
}
typedef struct { uint64_t fp; uint32_t a, b; } Ent;

// ---------------- D-check: full K-level key per block, given S
static int nkK;                       // PAF_t(1..hK), Q_t(0..hK)
static long TK[2 * MAXN];
static long count_total = 0, sols_S = 0, dcands = 0;
static int print_left = 0;

static void t_key(int b, const int *t, long *key) {
    for (int h = 1; h <= hK; h++) { long p = 0; for (int x = 0; x < n; x++) p += t[x] * t[(x + h) % n]; key[h - 1] = p; }
    for (int h = 0; h <= hK; h++) {
        long qq = 0; for (int x = 0; x < n; x++) qq += (long)s[b][x] * t[(x + h) % n] + (long)eps * t[x] * s[b][(x + h) % n];
        key[hK + h] = qq;
    }
}
typedef struct { long n; long *key; int8_t *t; } DList;
static void build_D(int b, const int8_t *S, DList *D) {
    int zpos[MAXN], z = 0, t[MAXN];
    for (int j = 0; j < k; j++) {
        int c0 = cl[b][j], c1 = cl[b][j + k];
        if (c0 && c1) {
            if (S[j] == 2) { t[j] = 1; t[j + k] = 1; }
            else if (S[j] == -2) { t[j] = -1; t[j + k] = -1; }
            else { zpos[z++] = j; t[j] = 1; t[j + k] = -1; }
        } else if (c0) { t[j] = S[j]; t[j + k] = 0; }
        else if (c1) { t[j] = 0; t[j + k] = S[j]; }
        else { t[j] = 0; t[j + k] = 0; }
    }
    long N = 1L << z; D->n = N; D->key = malloc(sizeof(long) * N * nkK); D->t = malloc((size_t)N * n);
    for (long m = 0; m < N; m++) {
        for (int i = 0; i < z; i++) { int j = zpos[i]; int sg = (m >> i) & 1 ? -1 : 1; t[j] = sg; t[j + k] = -sg; }
        t_key(b, t, D->key + m * nkK);
        for (int x = 0; x < n; x++) D->t[m * n + x] = (int8_t)t[x];
    }
}
static int verify_t(int8_t *tt[4]) {
    for (int h = 1; h < n; h++) { long p = 0; for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) p += tt[b][x] * tt[b][(x + h) % n]; if (p != T2[h]) return 0; }
    for (int h = 0; h < n; h++) { long qq = 0; for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) qq += (long)s[b][x] * tt[b][(x + h) % n] + (long)eps * tt[b][x] * s[b][(x + h) % n]; if (qq != QT[h]) return 0; }
    for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) if ((tt[b][x] != 0) != cl[b][x]) return 0;
    return 1;
}
static void d_check(const int8_t *SA, const int8_t *SB, const int8_t *SC, const int8_t *SD) {
    DList D[4]; const int8_t *SS[4] = {SA, SB, SC, SD};
    for (int b = 0; b < 4; b++) build_D(b, SS[b], &D[b]);
    // 2+2 join: hash (A,B) sums; probe (C,D)
    long NAB = D[0].n * D[1].n; long cap = 1; while (cap < 2 * NAB) cap <<= 1;
    Ent *H = calloc(cap, sizeof(Ent));
    long Tz[2 * MAXN]; for (int i = 0; i < nkK; i++) Tz[i] = TK[i];
    for (long a = 0; a < D[0].n; a++) for (long bb = 0; bb < D[1].n; bb++) {
        uint64_t h = 1469598103934665603ULL;
        for (int i = 0; i < nkK; i++) { long v = D[0].key[a * nkK + i] + D[1].key[bb * nkK + i]; h ^= (uint64_t)(v + 0x40000000L); h *= 1099511628211ULL; h ^= h >> 29; }
        h |= 1; long pos = h & (cap - 1);
        while (H[pos].fp) pos = (pos + 1) & (cap - 1);
        H[pos].fp = h; H[pos].a = (uint32_t)a; H[pos].b = (uint32_t)bb;
    }
    for (long c = 0; c < D[2].n; c++) for (long d = 0; d < D[3].n; d++) {
        uint64_t h = 1469598103934665603ULL;
        for (int i = 0; i < nkK; i++) { long v = Tz[i] - D[2].key[c * nkK + i] - D[3].key[d * nkK + i]; h ^= (uint64_t)(v + 0x40000000L); h *= 1099511628211ULL; h ^= h >> 29; }
        h |= 1; long pos = h & (cap - 1);
        while (H[pos].fp) {
            if (H[pos].fp == h) {
                long a = H[pos].a, bb = H[pos].b; int ok = 1;
                for (int i = 0; i < nkK && ok; i++) if (D[0].key[a * nkK + i] + D[1].key[bb * nkK + i] + D[2].key[c * nkK + i] + D[3].key[d * nkK + i] != Tz[i]) ok = 0;
                dcands++;
                if (ok) {
                    int8_t *tt[4] = {D[0].t + a * n, D[1].t + bb * n, D[2].t + c * n, D[3].t + d * n};
                    if (!verify_t(tt)) { fprintf(stderr, "VERIFY FAILED\n"); exit(5); }
                    count_total++;
                    if (print_left > 0) {
                        print_left--; printf("T");
                        for (int b = 0; b < 4; b++) { if (b) printf(" |"); for (int x = 0; x < n; x++) printf(" %d", tt[b][x]); }
                        printf("\n");
                    }
                }
            }
            pos = (pos + 1) & (cap - 1);
        }
    }
    free(H); for (int b = 0; b < 4; b++) { free(D[b].key); free(D[b].t); }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: h2396_stage2 n q < instance\n"); return 1; }
    n = atoi(argv[1]); q = atoi(argv[2]); if (n % 2) { fprintf(stderr, "n must be even\n"); return 1; }
    k = n / 2; hS = k / 2; hK = n / 2; nkS = 2 * (hS + 1); nkK = hK + (hK + 1);
    if (getenv("ST2_PRINT")) print_left = atoi(getenv("ST2_PRINT"));
    for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) if (scanf("%d", &u[b][x]) != 1) { fprintf(stderr, "bad input\n"); return 1; }
    // decode
    for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) {
        int v = u[b][x];
        if (v == 1 || v == -1) { s[b][x] = v; cl[b][x] = 1; }
        else if (v == q || v == -q) { s[b][x] = 0; cl[b][x] = 0; }             // x0 = e
        else if (v == q - 2 || v == -(q - 2)) { s[b][x] = -2 * (v / (q - 2)); cl[b][x] = 0; }   // x0 = -e: s = -2e
        else { fprintf(stderr, "bad letter %d\n", v); return 1; }
    }
    eps = ((q - 1) / 2) % 2 ? -1 : 1;
    char tag[16] = "";
    if (scanf("%15s", tag) == 1 && !strcmp(tag, "PLANT")) {
        planted = 1;
        for (int h = 1; h < n; h++) if (scanf("%ld", &T2[h]) != 1) return 1;
        for (int h = 0; h < n; h++) if (scanf("%ld", &QT[h]) != 1) return 1;
    }
    // stage-1 check
    long m = (long)q * n;
    for (int h = 0; h < n && !planted; h++) {
        long p = 0; for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) p += (long)u[b][x] * u[b][(x + h) % n];
        long T = h == 0 ? 4 * m + 4 - 4 * q : -4 * q;
        if (p != T) { fprintf(stderr, "not a stage-1 solution (lag %d: %ld vs %ld)\n", h, p, T); return 2; }
    }
    // targets
    for (int h = 0; h < n; h++) {
        long p = 0; for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) p += (long)s[b][x] * s[b][(x + h) % n];
        if (h == 0) { long c = 0; for (int b = 0; b < 4; b++) for (int x = 0; x < n; x++) c += cl[b][x]; T2[0] = c;
                      if (!planted && p + (long)q * c != 4 * m + 4) { fprintf(stderr, "lag-0 identity fails\n"); return 2; } }
        else if (!planted) { if (p % q) { printf("COUNT 0  reason T2-not-integral lag %d\n", h); return 0; } T2[h] = -p / q; }
    }
    for (int d = 0; d < k; d++) TS[d] = T2[d] + T2[d + k];
    for (int l = 0; l <= k / 2; l++) { double a = 0; for (int d = 0; d < k; d++) a += TS[d] * cos(2 * M_PI * l * d / k); TShat[l] = a; }
    for (int l = 0; l <= k / 2; l++) if (TShat[l] < -1e-6) { printf("COUNT 0  reason S-target-not-PSD char %d\n", l); return 0; }
    for (int b = 0; b < 4; b++) for (int j = 0; j < k; j++) sig[b][j] = s[b][j] + s[b][j + k];
    long TSk[2 * MAXN]; for (int d = 0; d <= hS; d++) { TSk[d] = TS[d]; TSk[hS + 1 + d] = QT[d] + QT[d + k]; }
    for (int h = 1; h <= hK; h++) TK[h - 1] = T2[h];
    for (int h = 0; h <= hK; h++) TK[hK + h] = QT[h];
    double t0 = now();
    for (int b = 0; b < 4; b++) { LS[b].nk = nkS; gen_S(b, 0); }
    fprintf(stderr, "S-lists %ld %ld %ld %ld (%.1fs)\n", LS[0].n, LS[1].n, LS[2].n, LS[3].n, now() - t0);
    // S-level 2+2 join
    long NAB = LS[0].n * LS[1].n; long cap = 1; while (cap < 2 * NAB) cap <<= 1;
    Ent *H = calloc(cap, sizeof(Ent)); if (!H) { fprintf(stderr, "oom S hash (%ld)\n", cap); return 3; }
    for (long a = 0; a < LS[0].n; a++) for (long b = 0; b < LS[1].n; b++) {
        uint64_t h = fp_key(LS[0].key + a * nkS, LS[1].key + b * nkS, nkS, 1, NULL);
        long pos = h & (cap - 1); while (H[pos].fp) pos = (pos + 1) & (cap - 1);
        H[pos].fp = h; H[pos].a = (uint32_t)a; H[pos].b = (uint32_t)b;
    }
    for (long c = 0; c < LS[2].n; c++) for (long d = 0; d < LS[3].n; d++) {
        uint64_t h = fp_key(LS[2].key + c * nkS, LS[3].key + d * nkS, nkS, -1, TSk);
        long pos = h & (cap - 1);
        while (H[pos].fp) {
            if (H[pos].fp == h) {
                long a = H[pos].a, b = H[pos].b; int ok = 1;
                for (int i = 0; i < nkS && ok; i++) if ((long)LS[0].key[a * nkS + i] + LS[1].key[b * nkS + i] + LS[2].key[c * nkS + i] + LS[3].key[d * nkS + i] != TSk[i]) ok = 0;
                if (ok) { sols_S++; d_check(LS[0].S + a * k, LS[1].S + b * k, LS[2].S + c * k, LS[3].S + d * k); }
            }
            pos = (pos + 1) & (cap - 1);
        }
    }
    printf("COUNT %ld  sols_S %ld  dcands %ld  lists %ld %ld %ld %ld  secs %.2f\n", count_total, sols_S, dcands, LS[0].n, LS[1].n, LS[2].n, LS[3].n, now() - t0);
    return 0;
}
