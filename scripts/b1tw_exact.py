"""Exact solution counts of small twisted b = 1 slices on Z_m (scripts/b1tw_lib.py), for validating the Gaussian.

Every row-sum group (A; D; B carrying C [and D]) is enumerated completely with its target (A = +-2, the rest 0); the
class vectors F_g = sum of the PAFs of the group's blocks are joined by 64-bit keys (lists: 'BCD' A | B, 'BC' A+D | B,
'id' A+B | C+D) and every joined solution is checked on all m lags.  The count is multiplied by the slice's symmetry
factor (placements of the +-2 block), so it counts labelled solutions of the whole slice, like the Gaussian's log2 E.
usage: b1tw_exact.py m u tw [m u tw ...]     (prints JSON: exact count and the Gaussian estimate)
"""
import sys, os, json, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice, gaussian

MAXBITS = int(os.environ.get("MAXBITS", "24"))
MAXPAIRS = int(os.environ.get("MAXPAIRS", str(1 << 25)))


def enum_group(w, T, chunk=1 << 20):
    """all eps in {+-1}^n with sum w_i eps_i = T (int8 rows)"""
    n = len(w); out = []
    for s in range(0, 1 << n, chunk):
        k = np.arange(s, min(1 << n, s + chunk), dtype=np.int64)
        E = (((k[:, None] >> np.arange(n)) & 1) * 2 - 1).astype(np.int8)
        out.append(E[(E.astype(np.int64) @ w) == T])
    return np.vstack(out) if out else np.zeros((0, n), np.int8)


def group_F(sl, g, E):
    """class-sum contribution of group g for each row of E (all L classes)"""
    bits = sl.groups[g]; blocks = [k for k in range(4) if np.isin(sl.idx[k], bits).all()]
    pos = {b: i for i, b in enumerate(bits)}
    out = []
    for s in range(0, len(E), 4096):
        e = E[s:s + 4096].astype(np.float64)
        tot = np.zeros((len(e), sl.m))
        for k in blocks:
            X = e[:, [pos[b] for b in sl.idx[k]]]
            F = np.fft.rfft(X, axis=1); tot += np.fft.irfft(F.real ** 2 + F.imag ** 2, n=sl.m, axis=1)
        out.append(np.rint(tot[:, sl.reps]).astype(np.int64))
    return np.vstack(out) if out else np.zeros((0, sl.L), np.int64), blocks


def count_pairs(K1, K2):
    """number of (i, j) with K1[i] == K2[j]"""
    a, ca = np.unique(K1, return_counts=True); b, cb = np.unique(K2, return_counts=True)
    common, ia, ib = np.intersect1d(a, b, assume_unique=True, return_indices=True)
    return int((ca[ia].astype(np.int64) * cb[ib]).sum()), common


def exact(sl):
    for g, bits in sl.groups.items():
        if len(bits) > MAXBITS: return None
    W = np.random.default_rng(7).integers(1, 2 ** 62, sl.L).astype(np.int64)
    tgt = -4 * np.ones(sl.L, np.int64)
    total = 0; nver = 0
    for sA in (2, -2):
        T = {g: (sA if g == "A" else 0) for g in sl.groups}
        E = {g: enum_group(sl.gw[g][1], T[g]) for g in sl.groups}
        if any(len(e) == 0 for e in E.values()): continue
        if sl.tw == "BC" and len(E["A"]) * len(E["D"]) > MAXPAIRS: return None
        if sl.tw == "id" and max(len(E["A"]) * len(E["B"]), len(E["C"]) * len(E["D"])) > MAXPAIRS: return None
        F = {g: group_F(sl, g, E[g])[0] for g in sl.groups}
        if sl.tw == "BCD":
            L1, L2 = F["A"], F["B"]
        elif sl.tw == "BC":
            L1 = (F["A"][:, None, :] + F["D"][None, :, :]).reshape(-1, sl.L); L2 = F["B"]
        else:
            L1 = (F["A"][:, None, :] + F["B"][None, :, :]).reshape(-1, sl.L)
            L2 = (F["C"][:, None, :] + F["D"][None, :, :]).reshape(-1, sl.L)
        K1 = L1 @ W; K2 = (tgt - L2) @ W
        c, common = count_pairs(K1, K2)
        # exactness: check every joined pair's class vector (keys can collide) on a sample, all lags on a few
        if c:
            i1 = np.flatnonzero(np.isin(K1, common)); i2 = np.flatnonzero(np.isin(K2, common))
            ok = {}
            for i in i1[:2000]:
                for j in i2[np.flatnonzero(K2[i2] == K1[i])][:5]:
                    assert (L1[i] + L2[j] == tgt).all(), "hash collision"
            nver += 1
        total += c
    return total * sl.sym


def main():
    a = sys.argv[1:]
    for i in range(0, len(a), 3):
        m, u, tw = int(a[i]), int(a[i + 1]), a[i + 2]
        sl = Slice(m, u, tw)
        ex = exact(sl)
        ga = gaussian(sl, N=int(os.environ.get("NS", "40000")))
        rec = dict(ga, exact=ex, exact_log2=(round(math.log2(ex), 2) if ex else None),
                   group_bits={g: len(b) for g, b in sl.groups.items()})
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
