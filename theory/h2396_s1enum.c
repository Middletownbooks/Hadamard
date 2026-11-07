// h2396_s1enum.c -- exact enumeration of stage 1 of the Paley plug-in Legendre-quadruple slice on Z_n x F_q, n even
// (report/prompt_599_plugin.md, gate A2).  Shares no code with theory/open_screen_stage1_exact.py or ab_lib.
// Stage 1: u_i = x0_i + (q-1) e_i on Z_n (x0 in {+-1}^n, e in {0,+-1}^n the defects), sum_i PAF_{u_i} = T,
// T(0) = 4m + 4 - 4q, T(d != 0) = -4q, m = qn; row sums (+-2, 0, 0, 0); #defects = e_plus + e_minus and
// #{x0 = e on a defect} = e_plus (both forced by T(0) once E is fixed).
// Method.  Per block, the defect pattern e is taken up to the group {k -> +-k + a} x {(x0, e) -> -(x0, e)} (order 4n),
// which maps (x0, e) block solutions to block solutions; canonical patterns carry their orbit sizes.  For every ordered
// 4-tuple of canonical patterns with sum D_i = E, every placement of the row sum (+-2 on one block) and every block:
// the x0 list is every x0 with sum x0 = r_i - (q-1) e^_i whose u passes the PSD bound |u^(phi)|^2 <= 4m + 4 for all
// phi != 1; then a 2+2 join on (PAF_u(d))_{d=0..n/2}.  Ordered stage-1 count = sum over tuples of prod_i |orbit(e_i)|
// x #(x0 solutions) [P: the count is constant on orbits].  Solutions are written (canonical e-tuple, x0, weight) for
// stage 2.  Each solution is re-checked by a direct PAF.
// usage: h2396_s1enum n q e_plus e_minus [outfile] [--zqpre]
// --zqpre: only (tuple, row placement) whose profile (D_i, e^_i, r_i) passes the Z_q level (theory/h2396_zq_profiles.py)
//          are joined; the count is then of the stage-1 solutions that can still complete [P: the Z_q level is necessary].
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAXN 32
static int n, q, ep, em, E, hK;
static long T[MAXN];
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

typedef struct { int D, esum; int8_t e[MAXN]; long orb; } Pat;
static Pat *pats[MAXN + 1]; static long npats[MAXN + 1];

// canonical form of a signed pattern under dihedral x sign: lexicographically smallest image (as int8 vector)
static void apply(const int8_t *e, int a, int refl, int sg, int8_t *out) {
    for (int x = 0; x < n; x++) { int y = refl ? ((a - x) % n + n) % n : (x + a) % n; out[y] = (int8_t)(sg * e[x]); }
}
static int cmp8(const int8_t *a, const int8_t *b) { for (int i = 0; i < n; i++) { if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1; } return 0; }
static long orbit_size(const int8_t *e) {
    // number of distinct images under the 4n group
    int8_t img[4 * MAXN * 2][MAXN]; int ni = 0;
    for (int a = 0; a < n; a++) for (int r = 0; r < 2; r++) for (int s = 0; s < 2; s++) {
        int8_t o[MAXN]; apply(e, a, r, s ? -1 : 1, o);
        int dup = 0; for (int t = 0; t < ni; t++) if (!cmp8(img[t], o)) { dup = 1; break; }
        if (!dup) { memcpy(img[ni], o, n); ni++; }
    }
    return ni;
}
static int is_canonical(const int8_t *e) {
    for (int a = 0; a < n; a++) for (int r = 0; r < 2; r++) for (int s = 0; s < 2; s++) {
        int8_t o[MAXN]; apply(e, a, r, s ? -1 : 1, o);
        if (cmp8(o, e) < 0) return 0;
    }
    return 1;
}
static void gen_pats(int D) {
    long cap = 1024; pats[D] = malloc(sizeof(Pat) * cap); npats[D] = 0;
    // choose D positions (combination) and signs
    int pos[MAXN];
    for (int i = 0; i < D; i++) pos[i] = i;
    while (1) {
        for (int sm = 0; sm < (1 << D); sm++) {
            int8_t e[MAXN]; memset(e, 0, sizeof(e)); int es = 0;
            for (int i = 0; i < D; i++) { e[pos[i]] = (sm >> i) & 1 ? -1 : 1; es += e[pos[i]]; }
            if (D && !is_canonical(e)) continue;
            if (npats[D] == cap) { cap *= 2; pats[D] = realloc(pats[D], sizeof(Pat) * cap); }
            Pat *P = &pats[D][npats[D]++]; P->D = D; P->esum = es; memcpy(P->e, e, n); P->orb = D ? orbit_size(e) : 2;
        }
        if (D == 0) break;
        int i = D - 1; while (i >= 0 && pos[i] == n - D + i) i--;
        if (i < 0) break;
        pos[i]++; for (int j = i + 1; j < D; j++) pos[j] = pos[j - 1] + 1;
    }
}

// x0 lists for (pattern, row sum)
typedef struct { long n; int32_t *key; uint32_t *x; } XL;      // x: bitmask of x0 = -1
static double cs[MAXN][MAXN], sn[MAXN][MAXN];
static void build_xl(const Pat *P, int r, XL *L) {
    long cap = 4096; L->n = 0; L->key = malloc(sizeof(int32_t) * cap * (hK + 1)); L->x = malloc(sizeof(uint32_t) * cap);
    int target_sum = r - (q - 1) * P->esum;          // sum x0
    if ((target_sum - n) % 2 || abs(target_sum) > n) return;
    int nneg = (n - target_sum) / 2;
    // enumerate subsets of size nneg (positions with x0 = -1)
    int pos[MAXN]; for (int i = 0; i < nneg; i++) pos[i] = i;
    double lim = 4.0 * q * n + 4 + 1e-6;
    while (1) {
        uint32_t mask = 0; for (int i = 0; i < nneg; i++) mask |= 1u << pos[i];
        int u[MAXN];
        for (int x = 0; x < n; x++) u[x] = ((mask >> x) & 1 ? -1 : 1) + (q - 1) * P->e[x];
        int ok = 1;
        for (int l = 1; l <= n / 2 && ok; l++) {
            double re = 0, im = 0; for (int x = 0; x < n; x++) { re += u[x] * cs[l][x]; im += u[x] * sn[l][x]; }
            if (re * re + im * im > lim) ok = 0;
        }
        if (ok) {
            if (L->n == cap) { cap *= 2; L->key = realloc(L->key, sizeof(int32_t) * cap * (hK + 1)); L->x = realloc(L->x, sizeof(uint32_t) * cap); }
            for (int d = 0; d <= hK; d++) { long p = 0; for (int x = 0; x < n; x++) p += u[x] * u[(x + d) % n]; L->key[L->n * (hK + 1) + d] = (int32_t)p; }
            L->x[L->n++] = mask;
        }
        if (nneg == 0) break;
        int i = nneg - 1; while (i >= 0 && pos[i] == n - nneg + i) i--;
        if (i < 0) break;
        pos[i]++; for (int j = i + 1; j < nneg; j++) pos[j] = pos[j - 1] + 1;
    }
}

typedef struct { uint64_t fp; uint32_t a, b; } Ent;
static uint64_t fpk(const int32_t *a, const int32_t *b, int sign) {
    uint64_t h = 1469598103934665603ULL;
    for (int d = 0; d <= hK; d++) { long v = sign > 0 ? (long)a[d] + b[d] : T[d] - a[d] - b[d]; h ^= (uint64_t)(v + 0x40000000L); h *= 1099511628211ULL; h ^= h >> 29; }
    return h | 1;
}

static int zqpre = 0;                 // --zqpre: skip (tuple, rows) whose (D_i, e^_i, r_i) profile fails the Z_q level
static int zq_profile_ok(Pat *Pt[4], int rr[4]) {
    // exists t^_i, |t^_i| <= n - D_i, t^_i = n - D_i (mod 2): sum (s^_i^2 + q t^_i^2) = 4qn + 4, sum s^_i t^_i = 0
    long sh[4], rest = 4L * q * n + 4;
    for (int b = 0; b < 4; b++) { sh[b] = rr[b] - (long)q * Pt[b]->esum; rest -= sh[b] * sh[b]; }
    if (rest < 0 || rest % q) return 0;
    long R2 = rest / q; int c[4]; for (int b = 0; b < 4; b++) c[b] = n - Pt[b]->D;
    for (int t0 = -c[0]; t0 <= c[0]; t0 += 2) { long q0 = R2 - (long)t0 * t0; if (q0 < 0) continue;
    for (int t1 = -c[1]; t1 <= c[1]; t1 += 2) { long q1 = q0 - (long)t1 * t1; if (q1 < 0) continue;
    for (int t2 = -c[2]; t2 <= c[2]; t2 += 2) { long q2 = q1 - (long)t2 * t2; if (q2 < 0) continue;
    for (int t3 = -c[3]; t3 <= c[3]; t3 += 2) { if ((long)t3 * t3 != q2) continue;
        if (((q - 1) / 2) % 2 || sh[0] * t0 + sh[1] * t1 + sh[2] * t2 + sh[3] * t3 == 0) return 1; } } } }   // sum s^ t^ = 0 only for eps = +1
    return 0;
}
static long double total_count = 0; static long nsol_written = 0, ntuples = 0, njoins = 0, nzqskip = 0;
static FILE *fo = NULL;

#define HLIM (1L << 26)                 // max hashed pairs per chunk (memory: 2 x HLIM x 16 bytes)
static void join_tuple(Pat *Pt0[4], int rr0[4]) {
    // order the blocks so that the two smallest x0 lists are hashed; the sum over blocks is symmetric
    XL L0[4]; int perm[4] = {0, 1, 2, 3};
    for (int b = 0; b < 4; b++) { build_xl(Pt0[b], rr0[b], &L0[b]); if (!L0[b].n) { for (int c = 0; c <= b; c++) { free(L0[c].key); free(L0[c].x); } return; } }
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (L0[perm[j]].n < L0[perm[i]].n) { int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    XL L[4]; Pat *Pt[4];
    for (int b = 0; b < 4; b++) { L[b] = L0[perm[b]]; Pt[b] = Pt0[perm[b]]; }
    njoins++;
    long chunk = L[0].n; while (chunk > 1 && chunk * L[1].n > HLIM) chunk = (chunk + 1) / 2;
    long nsol = 0;
    for (long a0 = 0; a0 < L[0].n; a0 += chunk) {
        long a1 = a0 + chunk < L[0].n ? a0 + chunk : L[0].n;
        long NAB = (a1 - a0) * L[1].n; long cap = 1; while (cap < 2 * NAB) cap <<= 1;
        Ent *H = calloc(cap, sizeof(Ent));
        if (!H) { fprintf(stderr, "oom join %ld\n", cap); exit(3); }
        for (long a = a0; a < a1; a++) for (long b = 0; b < L[1].n; b++) {
            uint64_t h = fpk(L[0].key + a * (hK + 1), L[1].key + b * (hK + 1), 1); long p = h & (cap - 1);
            while (H[p].fp) p = (p + 1) & (cap - 1);
            H[p].fp = h; H[p].a = (uint32_t)a; H[p].b = (uint32_t)b;
        }
        for (long c = 0; c < L[2].n; c++) for (long d = 0; d < L[3].n; d++) {
            uint64_t h = fpk(L[2].key + c * (hK + 1), L[3].key + d * (hK + 1), -1); long p = h & (cap - 1);
            while (H[p].fp) {
                if (H[p].fp == h) {
                    long a = H[p].a, b = H[p].b; int ok = 1;
                    for (int z = 0; z <= hK && ok; z++) if ((long)L[0].key[a * (hK + 1) + z] + L[1].key[b * (hK + 1) + z] + L[2].key[c * (hK + 1) + z] + L[3].key[d * (hK + 1) + z] != T[z]) ok = 0;
                    if (ok) {
                        uint32_t xs[4] = {L[0].x[a], L[1].x[b], L[2].x[c], L[3].x[d]};
                        int u[4][MAXN];
                        for (int i = 0; i < 4; i++) for (int x = 0; x < n; x++) u[i][x] = ((xs[i] >> x) & 1 ? -1 : 1) + (q - 1) * Pt[i]->e[x];
                        int good = 1, cp = 0, cd = 0;
                        for (int h2 = 0; h2 < n && good; h2++) { long s = 0; for (int i = 0; i < 4; i++) for (int x = 0; x < n; x++) s += (long)u[i][x] * u[i][(x + h2) % n]; if (s != (h2 ? -4 * q : 4L * q * n + 4 - 4 * q)) good = 0; }
                        for (int i = 0; i < 4; i++) for (int x = 0; x < n; x++) if (Pt[i]->e[x]) { cd++; if (((xs[i] >> x) & 1 ? -1 : 1) == Pt[i]->e[x]) cp++; }
                        if (!good || cd != E || cp != ep) { fprintf(stderr, "RECHECK FAILED\n"); exit(5); }
                        nsol++;
                        if (fo) {
                            long double w = 1; for (int i = 0; i < 4; i++) w *= Pt[i]->orb;
                            fprintf(fo, "%.0Lf", w);
                            // write the blocks in their original order (block 0 carries the row sum as placed)
                            int inv[4]; for (int b2 = 0; b2 < 4; b2++) inv[perm[b2]] = b2;
                            for (int i = 0; i < 4; i++) { fprintf(fo, " |"); for (int x = 0; x < n; x++) fprintf(fo, " %d", u[inv[i]][x]); }
                            fprintf(fo, "\n"); nsol_written++;
                        }
                    }
                }
                p = (p + 1) & (cap - 1);
            }
        }
        free(H);
    }
    if (nsol) { long double w = 1; for (int i = 0; i < 4; i++) w *= Pt[i]->orb; total_count += w * nsol; }
    for (int b = 0; b < 4; b++) { free(L0[b].key); free(L0[b].x); }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: h2396_s1enum n q e_plus e_minus [outfile]\n"); return 1; }
    n = atoi(argv[1]); q = atoi(argv[2]); ep = atoi(argv[3]); em = atoi(argv[4]); E = ep + em; hK = n / 2;
    if (n % 2 || n > 32) { fprintf(stderr, "n even <= 32\n"); return 1; }
    if ((q + 1) * ep + (q - 3) * em != 4 * (n - 1)) { fprintf(stderr, "allocation violates the defect identity\n"); return 2; }
    if (argc > 5 && strcmp(argv[5], "--zqpre")) fo = fopen(argv[5], "w");
    for (int a = 5; a < argc; a++) if (!strcmp(argv[a], "--zqpre")) zqpre = 1;
    for (int d = 0; d < n; d++) T[d] = d ? -4 * q : 4L * q * n + 4 - 4 * q;
    for (int l = 0; l < n; l++) for (int x = 0; x < n; x++) { cs[l][x] = cos(2 * M_PI * l * x / n); sn[l][x] = sin(2 * M_PI * l * x / n); }
    double t0 = now();
    for (int D = 0; D <= E; D++) gen_pats(D);
    long double sumorb[MAXN + 1];
    for (int D = 0; D <= E; D++) { sumorb[D] = 0; for (long i = 0; i < npats[D]; i++) sumorb[D] += pats[D][i].orb; }
    fprintf(stderr, "patterns:"); for (int D = 0; D <= E; D++) fprintf(stderr, " D=%d:%ld", D, npats[D]); fprintf(stderr, " (%.1fs)\n", now() - t0);
    // D = 0: orbit size 2 was stored for the empty pattern only to keep a weight; the empty pattern is fixed by
    // the group, so its orbit is 1.
    for (long i = 0; i < npats[0]; i++) pats[0][i].orb = 1;
    int Ds[4];
    for (Ds[0] = 0; Ds[0] <= E; Ds[0]++) for (Ds[1] = 0; Ds[0] + Ds[1] <= E; Ds[1]++) for (Ds[2] = 0; Ds[0] + Ds[1] + Ds[2] <= E; Ds[2]++) {
        Ds[3] = E - Ds[0] - Ds[1] - Ds[2];
        long i0, i1, i2, i3;
        for (i0 = 0; i0 < npats[Ds[0]]; i0++) for (i1 = 0; i1 < npats[Ds[1]]; i1++) for (i2 = 0; i2 < npats[Ds[2]]; i2++) for (i3 = 0; i3 < npats[Ds[3]]; i3++) {
            Pat *Pt[4] = {&pats[Ds[0]][i0], &pats[Ds[1]][i1], &pats[Ds[2]][i2], &pats[Ds[3]][i3]};
            // quick row-sum feasibility per block for r in {0, +-2}
            int okr[4][3];
            for (int b = 0; b < 4; b++) for (int ri = 0; ri < 3; ri++) {
                int r = (ri - 1) * 2, ts = r - (q - 1) * Pt[b]->esum;
                okr[b][ri] = (abs(ts) <= n && (ts - n) % 2 == 0);
            }
            for (int big = 0; big < 4; big++) for (int sg = -1; sg <= 1; sg += 2) {
                int rr[4]; int ok = 1;
                for (int b = 0; b < 4; b++) { rr[b] = b == big ? 2 * sg : 0; if (!okr[b][rr[b] / 2 + 1]) ok = 0; }
                if (!ok) continue;
                if (zqpre && !zq_profile_ok(Pt, rr)) { nzqskip++; continue; }
                ntuples++;
                join_tuple(Pt, rr);
            }
        }
    }
    printf("n=%d q=%d alloc=(%d,%d) %s=%.0Lf tuples=%ld joins=%ld written=%ld zq_skipped=%ld secs=%.1f\n", n, q, ep, em, zqpre ? "stage1_ordered_count_in_Zq_alive_profiles" : "stage1_ordered_count", total_count, ntuples, njoins, nsol_written, nzqskip, now() - t0);
    if (fo) fclose(fo);
    return 0;
}
