"""Exact quotient (lag-0) screen for twisted b = 1 slices (scripts/b1tw_lib.py) [the argument is P].

For every subgroup H of Z_m (order K, 1 < K < m) compress each block onto Z_F = Z_m / H (F = m / K):
y_k(b) = sum_{x = b mod F} X_k(x).  Summing sum_k PAF_k over H gives the lag-0 identity
    sum_k sum_b y_k(b)^2 = 4m - 4(K - 1),            with row sums sum_b y_k(b) = s_k.
In a twisted slice every row-sum group is a representative block R (A, D, or B) times its multiplicity (B, C [, D] have
equal compressed energies because C = B o u^-1, D = B o u^-2 and u permutes the fibres).  R is invariant under a
subgroup V (U for A and D, U^2 for B under 'BC', U^3 under 'BCD'); its compression is constant on the V-orbits of Z_F
("cells"), and each free bit lies in exactly one cell (the image of a V-orbit of Z_m is one V-orbit of Z_F).  So the
energy is separable: cell value y_c = sum_{bits in c} eps_i n_i (n_i = positions of bit i per fibre element), energy
w_c y_c^2, row sum w_c y_c.  A DP over the cells decides exactly whether energy T0 = 4m - 4(K - 1) with row sums
(A, D, B) = (+-2, 0, 0) (and the symmetric placements) is attainable.
usage: b1tw_lag0.py m u tw [m u tw ...]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice


def divisors(n):
    return [k for k in range(1, n + 1) if n % k == 0]


def cells_of_group(sl, g, F):
    """for group g: representative block k, multiplicity, and per cell (weight w_c, list of bit counts n_i)"""
    bits = sl.groups[g]
    rep = {"A": 0, "D": 3, "B": 1, "C": 2}[g]
    mult = {"A": 1, "D": 1, "C": 1, "B": {"id": 1, "BC": 2, "BCD": 3}[sl.tw]}[g]
    vpow = {"A": 1, "D": 1, "C": 1, "B": {"id": 1, "BC": 2, "BCD": 3}[sl.tw]}[g]
    m, u = sl.m, sl.u; v = pow(u, vpow, m)
    # V-orbits of Z_F (cells)
    lab = -np.ones(F, np.int64); cells = []
    for b in range(F):
        if lab[b] >= 0: continue
        o, y = [b], b * v % F
        while y != b:
            o.append(y); y = y * v % F
        for y in o: lab[y] = len(cells)
        cells.append(o)
    idx = sl.idx[rep]
    out = {}
    for i in bits:
        pos = np.flatnonzero(idx == i)
        cs = set(lab[pos % F].tolist()); assert len(cs) == 1, "bit spans two cells"
        c = cs.pop(); w = len(cells[c]); assert len(pos) % w == 0
        out.setdefault(c, []).append(len(pos) // w)
    return rep, mult, [(len(cells[c]), ns) for c, ns in sorted(out.items())]


def letters(ns):
    S = {0}
    for n in ns: S = {s + n for s in S} | {s - n for s in S}
    return sorted(S)


def group_table(cells, Emax, mult):
    """boolean table T[E, s + S] : energy E (already times mult) and row sum s attainable by the representative block"""
    S = sum(w * sum(ns) for w, ns in cells)
    T = np.zeros((Emax + 1, 2 * S + 1), bool); T[0, S] = True
    for w, ns in cells:
        N = np.zeros_like(T)
        for y in letters(ns):
            e = mult * w * y * y; s = w * y
            if e > Emax: continue
            src = T[:Emax + 1 - e]
            if s >= 0: N[e:, s:] |= src[:, :2 * S + 1 - s]
            else: N[e:, :2 * S + 1 + s] |= src[:, -s:]
        T = N
    return T, S


def screen(sl, Ks=None):
    m = sl.m; res = []
    for K in (Ks or [k for k in divisors(m) if 1 < k < m]):
        F = m // K; T0 = 4 * m - 4 * (K - 1)
        tabs = {}
        for g in sl.groups:
            rep, mult, cells = cells_of_group(sl, g, F)
            tabs[g] = group_table(cells, T0, mult)
        # combine: placements of the +-2 block among the pi-fixed groups
        fixed = {"id": ["A", "B", "C", "D"], "BC": ["A", "D"], "BCD": ["A"]}[sl.tw]
        ok = False
        for g2 in fixed:
            for s2 in (2, -2):
                # energy-reachability with row sums fixed per group: E-vectors
                acc = np.zeros(T0 + 1, bool); acc[0] = True
                for g, (T, S) in tabs.items():
                    want = s2 if g == g2 else 0
                    if abs(want) > S: acc[:] = False; break
                    col = T[:, want + S]
                    nxt = np.zeros_like(acc)
                    for e in np.flatnonzero(col):
                        nxt[e:] |= acc[:T0 + 1 - e]
                    acc = nxt
                if acc[T0]: ok = True; break
            if ok: break
        res.append(dict(K=K, F=F, T0=T0, attainable=bool(ok)))
        if not ok: break
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    for i in range(0, len(a), 3):
        sl = Slice(int(a[i]), int(a[i + 1]), a[i + 2])
        r = screen(sl)
        dead = [x for x in r if not x["attainable"]]
        print(json.dumps(dict(m=sl.m, u=sl.u, d=sl.d, tw=sl.tw, alive=not dead, killed_by=dead[0] if dead else None,
                              checked=[x["K"] for x in r])), flush=True)


def quotient_screen(sl, K, maxvec=2_000_000):
    """exact all-lags quotient test for one subgroup H of order K: the compressed blocks y_k on Z_F (F = m/K) must
    satisfy sum_k PAF_{y_k}(j) = -4K + [j = 0](4m + 4) for every j, with row sums (+-2, 0, 0, 0) placed as in the slice.
    Enumerates every attainable compressed representative block per group (cells x letters).  Returns True/False, or
    None if a group has more than maxvec compressed vectors."""
    import itertools
    m, F = sl.m, sl.m // K
    tgt = np.full(F, -4 * K, np.int64); tgt[0] += 4 * m + 4
    vpow = {"id": 1, "BC": 2, "BCD": 3}[sl.tw]
    contrib = {}
    for g in sl.groups:
        rep, mult, cells = cells_of_group(sl, g, F)
        # cell membership on Z_F for the representative block
        v = pow(sl.u, vpow if g == "B" else 1, m)
        lab = -np.ones(F, np.int64); co = []
        for b in range(F):
            if lab[b] >= 0: continue
            o, y = [b], b * v % F
            while y != b:
                o.append(y); y = y * v % F
            for y in o: lab[y] = len(co)
            co.append(o)
        # cells_of_group returns cells in sorted cell-index order with bits; rebuild the letter sets per cell index
        bits = sl.groups[g]; idx = sl.idx[rep]; per = {}
        for i in bits:
            pos = np.flatnonzero(idx == i); c = int(lab[pos[0] % F]); per.setdefault(c, []).append(len(pos) // len(co[c]))
        cids = sorted(per); L = [letters(per[c]) for c in cids]
        n = 1
        for x in L: n *= len(x)
        if n > maxvec: return None
        Y = np.zeros((n, F), np.int64)
        grid = np.array(list(itertools.product(*L)), np.int64) if L else np.zeros((1, 0), np.int64)
        for j, c in enumerate(cids): Y[:, co[c]] = grid[:, [j]]
        # blocks of the group on Z_F: representative, then images under u^-1 (C) and u^-2 (D)
        ui = pow(sl.u, -1, m) % F if math.gcd(sl.u, F) == 1 else None
        blocks = [Y]
        if g == "B" and sl.tw in ("BC", "BCD"):
            perm1 = np.array([(ui * b) % F for b in range(F)])
            blocks.append(Y[:, perm1])
            if sl.tw == "BCD": blocks.append(Y[:, (perm1[perm1])])
        P = np.zeros((n, F), np.int64)
        for Z in blocks:
            Fz = np.fft.fft(Z.astype(float), axis=1); P += np.rint(np.fft.ifft(np.abs(Fz) ** 2, axis=1).real).astype(np.int64)
        contrib[g] = (Y.sum(1), P)
    fixed = {"id": ["A", "B", "C", "D"], "BC": ["A", "D"], "BCD": ["A"]}[sl.tw]
    W = np.random.default_rng(3).integers(1, 2 ** 62, F).astype(np.int64)
    gs = list(contrib)
    for g2 in fixed:
        for s2 in (2, -2):
            sel = {g: contrib[g][1][contrib[g][0] == (s2 if g == g2 else 0)] for g in gs}
            if any(len(x) == 0 for x in sel.values()): continue
            # join: last group against the sum of the others
            acc = np.unique(sel[gs[0]], axis=0)
            for g in gs[1:-1]:
                b = np.unique(sel[g], axis=0)
                acc = np.unique((acc[:, None, :] + b[None, :, :]).reshape(-1, F), axis=0)
                if len(acc) > 5_000_000: return None
            last = np.unique(sel[gs[-1]], axis=0)
            k1 = set((acc @ W).tolist()); k2 = ((tgt - last) @ W).tolist()
            if any(k in k1 for k in k2): return True
    return False
