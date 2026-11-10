#!/usr/bin/env python3
"""Hashed meet-in-the-middle version of theory/fable_tmult.py for t = 5 (and 3): exhaustive over all assignments
of signed letters to the 4 x t slots, for a letter-relation type given as in fable_tmult_types.py
(e.g. "f s s s", "f f s =2").  Blocks are enumerated up to sign and cyclic shift; pair sums of 64-bit hashes are
sorted and matched; every hash match is re-checked exactly.  [M]
Usage: fable_tmult_h.py t "f s s s" ["f f s =2" ...]
"""
import sys, itertools, random, time
import numpy as np
sys.path.insert(0, 'theory')
import fable_tmult_types as TT

def run(t, ty, show=2):
    sym, conj, letters = TT.make(ty)
    nsym = len(conj)
    mono = {}
    for a in range(nsym):
        for b in range(a, nsym):
            mono[(a, b)] = len(mono)
    nm = len(mono); half = (t - 1) // 2
    cls = {}
    for s in range(nsym):
        cls.setdefault(min(s, conj[s]), len(cls))
    ncls = len(cls)
    target = np.zeros(ncls, dtype=np.int64)
    for k in range(4):
        s = sym[(k, 0)]
        target[cls[min(s, conj[s])]] += t
    L = [(e, l) for e in (1, -1) for l in range(len(letters))]
    nL = len(L)
    lsym = np.array([sym[letters[l]] for e, l in L]); lsgn = np.array([e for e, l in L])
    lcon = np.array([conj[s] for s in lsym])
    lcls = np.array([cls[min(s, conj[s])] for s in lsym])
    # all configs as array (nL^t, t)
    cfg = np.array(list(itertools.product(range(nL), repeat=t)), dtype=np.int16)
    # canonical: keep cfg that is lexicographically minimal among shifts and global sign flips
    flip = np.array([(i + nL // 2) % nL for i in range(nL)])     # sign flip of a letter index
    keep = np.ones(len(cfg), dtype=bool)
    base = (cfg.astype(np.int64) * (nL ** np.arange(t - 1, -1, -1))).sum(1)
    for sh in range(t):
        for fl in (0, 1):
            if sh == 0 and fl == 0: continue
            c2 = np.roll(cfg, -sh, axis=1)
            if fl: c2 = flip[c2]
            code = (c2.astype(np.int64) * (nL ** np.arange(t - 1, -1, -1))).sum(1)
            keep &= base <= code
    cfg = cfg[keep]
    n = len(cfg)
    # exact vectors
    V = np.zeros((n, half * nm + ncls), dtype=np.int16)
    for d in range(1, half + 1):
        for tau in range(t):
            i1 = cfg[:, (tau + d) % t]; i0 = cfg[:, tau]
            a = lsym[i1]; b = lcon[i0]
            lo = np.minimum(a, b); hi = np.maximum(a, b)
            # monomial index for (lo, hi)
            midx = np.array([[mono[(x, y)] if x <= y else -1 for y in range(nsym)] for x in range(nsym)])[lo, hi]
            np.add.at(V, (np.arange(n), (d - 1) * nm + midx), (lsgn[i1] * lsgn[i0]).astype(np.int16))
    for tau in range(t):
        np.add.at(V, (np.arange(n), half * nm + lcls[cfg[:, tau]]), 1)
    rng = np.random.default_rng(12345)
    w = rng.integers(-(1 << 61), 1 << 61, size=V.shape[1], dtype=np.int64)
    np.seterr(over='ignore')
    h = (V.astype(np.int64) * w).sum(1)
    tgt = np.concatenate([np.zeros(half * nm, dtype=np.int64), target])
    ht = int((tgt * w).sum())
    # pair sums i <= j
    t0 = time.time()
    iu, ju = np.triu_indices(n)
    A = h[iu] + h[ju]
    order = np.argsort(A, kind='stable'); As = A[order]
    B = np.int64(ht) - As                       # need another pair with sum B
    pos = np.searchsorted(As, B)
    ok = (pos < len(As)) & (As[np.minimum(pos, len(As) - 1)] == B)
    cand = np.nonzero(ok)[0]
    nsol = 0; ex = []
    names = 'ABCD'
    for c in cand:
        p = order[c]
        lo = np.searchsorted(As, B[c], 'left'); hi_ = np.searchsorted(As, B[c], 'right')
        for q in order[lo:hi_]:
            if q < p: continue
            tot = V[iu[p]].astype(np.int64) + V[ju[p]] + V[iu[q]] + V[ju[q]]
            if np.array_equal(tot, tgt):
                nsol += 1
                if len(ex) < show: ex.append((iu[p], ju[p], iu[q], ju[q]))
    dim = sum(2 if e == 'f' else 1 if e == 's' else 0 for e in ty)
    print(f"t={t} type={' '.join(ty)} dim={dim}: {nL} signed letters, {n} block classes, {len(A)} pairs, solutions={nsol} ({time.time()-t0:.0f}s)", flush=True)
    for e4 in ex:
        for i in e4:
            print('    ', ' '.join(('+' if L[x][0] > 0 else '-') + names[letters[L[x][1]][0]] + ('*' if letters[L[x][1]][1] else ' ') for x in cfg[i]))
        print('     --')
    return nsol

if __name__ == '__main__':
    t = int(sys.argv[1])
    for s in sys.argv[2:]:
        run(t, s.split())
