r"""[P, M] The negation split at v = 863 (S7 host, Codex certificate y): S = x + x(-.), D = x - x(-.).
[P] S is symmetric and D antisymmetric, so CC(S, D) = -CC(D, S) and PAF_x = (PAF_S + PAF_D) / 4 for every block: every
class equation reads (S-part) + (D-part) = 4 F; no equation involves S alone, so there is no S-level to solve first and
all of stage 2 (275.6 bits given y) falls on the D-level.  This script computes h_D = H(D | S, y) exactly:
cells (q, 0) pair their own classes (C_j, -C_j); cells (q, u), (q, -u) with u != 0 pair (h, C_j) with (-h, -C_j);
singletons are fixed by y.  D_def = 275.6 - h_D (all equations on the D side).
usage: python3 theory/h863_negsplit.py"""
import sys, os, json, math, itertools as it
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L
Y, E, S = L.load_cert()
rows = {s: [r for r in it.product((-1, 1), repeat=6) if sum(r) == s] for s in (-2, 0, 2)}
NEG = L.NEG                                      # class index of -C_j


def H(cnt, tot):
    return -sum(v / tot * math.log2(v / tot) for v in cnt.values())


hD = 0.0
for i in range(4):
    for q in range(4):
        # u = 0: one cell, pairs (j, NEG[j]) inside the row
        R = rows[int(S[i, q * 5])]
        cnt = {}
        for r in R:
            k = tuple(r[a] + r[b] for a, b in L.PAIRS); cnt[k] = cnt.get(k, 0) + 1
        hD += math.log2(len(R)) - H(cnt, len(R))
        # u = 1, 4 and u = 2, 3: two cells, pairs (h, C_j) <-> (-h, C_NEG[j])
        for u, v in ((1, 4), (2, 3)):
            R0, R1 = rows[int(S[i, q * 5 + u])], rows[int(S[i, q * 5 + v])]
            cnt = {}
            for r in R0:
                for t in R1:
                    k = tuple(r[j] + t[NEG[j]] for j in range(6)); cnt[k] = cnt.get(k, 0) + 1
            hD += math.log2(len(R0) * len(R1)) - H(cnt, len(R0) * len(R1))
res = dict(h_D=round(hD, 2), C_total_given_y=275.6, D_def=round(275.6 - hD, 2), log2_fibre_y=339.11)
print(json.dumps(res))
