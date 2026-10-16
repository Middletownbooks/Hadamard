"""Twisted b = 1 slices on a finite abelian group G = Z_{n1} x ... x Z_{nk} (noncyclic allowed), for the open orders
with square factors in v - 1 (report/prompt_b1_twist_design.md §3, stage B).

An automorphism u of G is given by its index permutation (theory/abgroup.py AbGroup / aut_perm).  Twist types and bit
layout are those of scripts/b1tw_lib.py ('id', 'BC', 'BCD'), with orbits of u in place of multiplier orbits.
PAF over G by an n-dimensional FFT; equations = classes of <u, -1> on G minus {0}.
Quotient screens [P] for every u-stable subgroup H given by generators: compress each block onto the cosets of H, the
lag-0 and all-lags conditions sum_k PAF_{y_k}(j) = -4|H| + [j = 0](4|G| + 4) on G/H, with the twisted structure
(C = B o u^-1, D = B o u^-2 on G/H).
"""
import os, sys, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "theory"))
from b1tw_lib import count_signed, SignedSampler
from abgroup import AbGroup


def perm_order(P):
    k, Q = 1, P.copy()
    while not np.array_equal(Q, np.arange(len(P))):
        Q = P[Q]; k += 1
    return k


class AbSlice:
    def __init__(self, G, uperm, tw, label=""):
        self.G, self.u, self.tw, self.label = G, np.asarray(uperm), tw, label
        n = G.n; self.n = n; self.shape = G.mods
        self.d = perm_order(self.u)
        if tw == "BC": assert self.d % 2 == 0
        if tw == "BCD": assert self.d % 3 == 0
        seen = np.zeros(n, bool); self.orbits = []
        for x in range(n):
            if seen[x]: continue
            o = [x]; y = int(self.u[x])
            while y != x:
                o.append(y); y = int(self.u[y])
            seen[o] = True; self.orbits.append(o)
        neg = G.neg_perm(); self.neg = neg
        seen = np.zeros(n, bool); seen[0] = True; reps, sizes = [], []
        for x in range(1, n):
            if seen[x]: continue
            c, stack = {x}, [x]
            while stack:
                y = stack.pop()
                for z in (int(self.u[y]), int(neg[y])):
                    if z not in c: c.add(z); stack.append(z)
            seen[list(c)] = True; reps.append(x); sizes.append(len(c))
        self.reps, self.csize = np.array(reps), np.array(sizes)
        idx = -np.ones((4, n), np.int64); groups = {k: [] for k in "ABCD"}; nb = 0
        def new(g):
            nonlocal nb
            groups[g].append(nb); nb += 1; return nb - 1
        for o in self.orbits:
            s = len(o); a = new("A")
            for y in o: idx[0, y] = a
            if tw == "id":
                for k, g in ((1, "B"), (2, "C"), (3, "D")):
                    b = new(g)
                    for y in o: idx[k, y] = b
            elif tw == "BC":
                dd = new("D")
                for y in o: idx[3, y] = dd
                if s % 2 == 0:
                    b, c = new("B"), new("B")
                    for j, y in enumerate(o): idx[1, y], idx[2, y] = (b, c) if j % 2 == 0 else (c, b)
                else:
                    b = new("B")
                    for y in o: idx[1, y] = idx[2, y] = b
            else:
                if s % 3 == 0:
                    b, c, dd = new("B"), new("B"), new("B")
                    for j, y in enumerate(o):
                        idx[1, y] = (b, dd, c)[j % 3]; idx[2, y] = (c, b, dd)[j % 3]; idx[3, y] = (dd, c, b)[j % 3]
                else:
                    b = new("B")
                    for y in o: idx[1, y] = idx[2, y] = idx[3, y] = b
        assert (idx >= 0).all()
        self.idx, self.nbits = idx, nb
        self.groups = {g: np.array(v, np.int64) for g, v in groups.items() if v}
        blk = {"A": 0, "B": 1, "C": 2, "D": 3}
        self.gw = {g: (blk[g], np.bincount(idx[blk[g]], minlength=nb)[bits].astype(np.int64)) for g, bits in self.groups.items()}
        self.sym = {"id": 4, "BC": 2, "BCD": 1}[tw]
        self.r, self.L = len(self.orbits), len(self.reps)
        uu = np.arange(n); pows = []
        for _ in range(self.d): pows.append(uu.copy()); uu = self.u[uu]
        self.minus1_in_U = any(np.array_equal(p, neg) for p in pows)
        # classes of size 1: elements with x = -x and u x = x (2-torsion fixed by u); their values are not all implied
        self.small_classes = [c for c in range(self.L) if self.csize[c] == 1]

    def paf_classes(self, beta, batch=2000):
        out = []
        for s in range(0, len(beta), batch):
            X = beta[s:s + batch][:, self.idx].astype(np.float64)          # (b, 4, n)
            X = X.reshape((len(X), 4) + tuple(self.shape))
            ax = tuple(range(2, 2 + len(self.shape)))
            Fx = np.fft.fftn(X, axes=ax)
            P = np.fft.ifftn(np.abs(Fx) ** 2, axes=ax).real.reshape(len(X), 4, self.n).sum(1)
            out.append(np.rint(P[:, self.reps]).astype(np.int64))
        return np.vstack(out)


def log2_rowsum_configs(sl):
    tot = 0
    for sA in (2, -2):
        c = 1
        for g, (k, w) in sl.gw.items(): c *= count_signed(w, sA if g == "A" else 0)
        tot += c
    return math.log2(tot) if tot else float("-inf")


def sample(sl, N, rng):
    beta = np.empty((N, sl.nbits), np.int8); k2 = int(rng.binomial(N, 0.5)); row = 0
    for sA, k in ((2, k2), (-2, N - k2)):
        if k == 0: continue
        for g, bits in sl.groups.items():
            beta[row:row + k, bits] = SignedSampler(sl.gw[g][1], sA if g == "A" else 0).sample(k, rng)
        row += k
    return beta


def gaussian(sl, N=20000, seed=1):
    rec = dict(G=repr(sl.G), label=sl.label, d=sl.d, tw=sl.tw, orbits=sl.r, classes=sl.L, bits=sl.nbits,
               minus1_in_U=sl.minus1_in_U, group_bits={g: len(b) for g, b in sl.groups.items()})
    lrs = log2_rowsum_configs(sl); rec["log2_rowsum_configs"] = round(lrs, 2)
    if lrs == float("-inf"):
        rec.update(feasible=False, why="row sums"); return rec
    F = sl.paf_classes(sample(sl, N, np.random.default_rng(seed)))
    # drop one size-1 class (implied by the row sums) if there is one; otherwise drop the smallest class
    drop = sl.small_classes[0] if sl.small_classes else int(np.argmin(sl.csize))
    keep = [c for c in range(sl.L) if c != drop]
    F = F[:, keep]; tgt = -4 * np.ones(len(keep))
    steps = np.array([max(1, int(np.gcd.reduce(np.abs(F[:, c] - F[0, c])))) for c in range(len(keep))])
    if not all(((tgt[c] - F[0, c]) % steps[c]) == 0 for c in range(len(keep))):
        rec.update(feasible=False, why="target off the per-class lattice"); return rec
    mu = F.mean(0); S = np.atleast_2d(np.cov(F.T)); sign, logdet = np.linalg.slogdet(S)
    if sign <= 0:
        rec.update(feasible=None, why="degenerate covariance"); return rec
    dl = tgt - mu; q = dl @ np.linalg.solve(S, dl)
    log2P = (np.log(steps).sum() - 0.5 * len(keep) * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
    rec.update(feasible=True, n_cond=round(-log2P, 1), log2E=round(math.log2(sl.sym) + lrs + log2P, 1),
               steps=sorted(set(steps.tolist())), mahal=round(float(q), 2), dropped_class_size=int(sl.csize[drop]))
    return rec


# ------------------------------------------------------------------ quotient screens
def cosets(G, gens):
    """coset label of every element for the subgroup generated by gens"""
    n = G.n; H = {0}; frontier = [0]
    add = [G.add_perm(g) for g in gens]
    while frontier:
        nf = []
        for h in frontier:
            for P in add:
                y = int(P[h])
                if y not in H: H.add(y); nf.append(y)
        frontier = nf
    H = sorted(H); lab = -np.ones(n, np.int64); reps = []
    for x in range(n):
        if lab[x] >= 0: continue
        cos = G.add_perm(x)[H]; lab[cos] = len(reps); reps.append(x)
    return lab, reps, len(H)


def letters(ns):
    S = {0}
    for k in ns: S = {s + k for s in S} | {s - k for s in S}
    return sorted(S)


def quotient_screen(sl, gens, maxvec=300_000, energy_only_fallback=True):
    """exact test on G/H, H = <gens> (must be u-stable).  Returns (verdict, test) with verdict True/False/None."""
    G = sl.G; lab, creps, K = cosets(G, gens); Fn = len(creps); n = sl.n
    # u acts on cosets
    ucos = np.array([lab[sl.u[creps[b]]] for b in range(Fn)])
    assert all(lab[sl.u[x]] == ucos[lab[x]] for x in range(n)), "H not u-stable"
    # coset addition table and negation
    addc = np.array([[lab[G.idx(G.coords[creps[a]] + G.coords[creps[b]])] for b in range(Fn)] for a in range(Fn)])
    tgt = np.full(Fn, -4 * K, np.int64); tgt[lab[0]] += 4 * n + 4
    vpow = {"id": 1, "BC": 2, "BCD": 3}[sl.tw]
    contrib = {}; E_tables = {}
    for g in sl.groups:
        rep = {"A": 0, "B": 1, "C": 2, "D": 3}[g]; vp = vpow if g == "B" else 1
        mult = vpow if (g == "B" and sl.tw != "id") else 1
        vc = np.arange(Fn)
        for _ in range(vp): vc = ucos[vc]
        cl = -np.ones(Fn, np.int64); cells = []
        for b in range(Fn):
            if cl[b] >= 0: continue
            o, y = [b], int(vc[b])
            while y != b:
                o.append(y); y = int(vc[y])
            cl[o] = len(cells); cells.append(o)
        per = {}
        for i in sl.groups[g]:
            pos = np.flatnonzero(sl.idx[rep] == i); cs = set(cl[lab[pos]].tolist()); assert len(cs) == 1
            c = cs.pop(); w = len(cells[c]); assert len(pos) % w == 0
            per.setdefault(c, []).append(len(pos) // w)
        cids = sorted(per); Ls = [letters(per[c]) for c in cids]
        E_tables[g] = (mult, [(len(cells[c]), Ls[j]) for j, c in enumerate(cids)])
        tot = 1
        for x in Ls: tot *= len(x)
        if tot > maxvec: contrib = None; continue
        if contrib is None: continue
        Y = np.zeros((tot, Fn), np.int64)
        grid = np.array(list(itertools.product(*Ls)), np.int64)
        for j, c in enumerate(cids): Y[:, cells[c]] = grid[:, [j]]
        blocks = [Y]
        if g == "B" and sl.tw in ("BC", "BCD"):
            uinv = np.argsort(ucos)                    # coset map of u^-1
            blocks.append(Y[:, uinv])                  # y_C(b) = y_B(u^-1 b)
            if sl.tw == "BCD": blocks.append(Y[:, uinv[uinv]])
        P = np.zeros((tot, Fn), np.int64)
        for Z in blocks:
            for j in range(Fn): P[:, j] += (Z * Z[:, addc[:, j]]).sum(1)
        contrib[g] = (Y.sum(1), P)
    fixed = {"id": ["A", "B", "C", "D"], "BC": ["A", "D"], "BCD": ["A"]}[sl.tw]
    if contrib is not None:
        W = np.random.default_rng(3).integers(1, 2 ** 62, Fn).astype(np.int64); gs = list(contrib)
        unknown = False
        for g2 in fixed:
            for s2 in (2, -2):
                sel = {g: contrib[g][1][contrib[g][0] == (s2 if g == g2 else 0)] for g in gs}
                if any(len(x) == 0 for x in sel.values()): continue
                acc = np.unique(sel[gs[0]] @ W)
                for g in gs[1:-1]:
                    b = np.unique(sel[g] @ W)
                    if len(acc) * len(b) > 50_000_000: unknown = True; break
                    acc = np.unique((acc[:, None] + b[None, :]).ravel())
                if unknown: break
                if np.isin(np.unique((tgt - sel[gs[-1]]) @ W), acc).any(): return True, "all-lags"
            if unknown: break
        if not unknown: return False, "all-lags"
    if not energy_only_fallback: return None, "too large"
    T0 = 4 * n - 4 * (K - 1)
    acc = np.zeros(T0 + 1, bool); acc[0] = True
    for g, (mult, cells) in E_tables.items():
        for w, Ls in cells:
            nxt = np.zeros_like(acc)
            for e in sorted({mult * w * y * y for y in Ls}):
                if e <= T0: nxt[e:] |= acc[:T0 + 1 - e]
            acc = nxt
    return bool(acc[T0]), "lag0-energy"
