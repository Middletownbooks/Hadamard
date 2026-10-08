"""Exact quotient (compression) screen for symmetric b1 / b0 / lp slices over an abelian group G.

Choose a subgroup H that is a union of whole components of G (so every symmetry preserves it); Q = G/H.
For S-invariant blocks X_k, the compressions Y_k(q) = sum_{x in q} X_k(x) satisfy
        sum_k PAF_{Y_k}(q) = sum_{h in H} T(q + h)      (q in Q)
where T is the target on G (b1: 4m at 0, -4 elsewhere; b0: 4m at 0, 0 elsewhere; lp: 2m at 0, -2 elsewhere),
plus the row-sum condition.  The achievable Y_k form a product over the S-bar orbits R of Q of the sets
V_R (signed sums of the S-orbits lying over R), computed exactly by DP.  Each Y_k is filtered by the PSD
bound over Q, then an exact multiset join is run.  'False' means the slice has no solutions [exact].
usage: python3 compress_screen.py "comps" "gens(json)" kind H_components
"""
import itertools
import json
import sys

import numpy as np

from abgroup import AbGroup, AutGroup, slice_structure


def quotient_data(G, auts, Hcomp):
    """Q = product of the components not in Hcomp; projection of G -> Q as index array"""
    keep = [i for i in range(len(G.comps)) if i not in Hcomp]
    Qcomps = [G.comps[i] for i in keep]
    Q = AbGroup(Qcomps)
    # coordinate columns of G belonging to kept components
    cols = []
    c0 = 0
    for i, c in enumerate(G.comps):
        k = c[2] if c[0] == 'E' else 1
        if i in keep:
            cols += list(range(c0, c0 + k))
        c0 += k
    proj = Q.idx(G.coords[:, cols])
    Qauts = [[a[i] for i in keep] for a in auts]
    return Q, proj, Qauts


def screen(G, auts, kind, Hcomp, max_vectors=1_200_000, verbose=True):
    orb, lag, frq = slice_structure(G, auts)
    r = int(orb.max()) + 1
    Q, proj, Qauts = quotient_data(G, auts, Hcomp)
    qorb, _, _ = slice_structure(Q, Qauts)
    nQ = Q.n
    m = G.n
    H = m // nQ
    # target on Q
    if kind == 'b1':
        T0, Toff, bound = 4 * m, -4, 4 * (m + 1)
    elif kind == 'b0':
        T0, Toff, bound = 4 * m, 0, 4 * m
    else:
        T0, Toff, bound = 2 * m, -2, 2 * m + 2
    tq = np.zeros(nQ, dtype=np.int64)
    tq[0] = T0 + Toff * (H - 1)
    tq[1:] = Toff * H
    # per S-orbit of G: count vector over Q
    cnt = np.zeros((r, nQ), dtype=np.int64)
    np.add.at(cnt, (orb, proj), 1)
    # group S-orbits by the Q-orbit they lie over (an S-orbit maps into one S-bar orbit)
    Rlab = qorb
    groups = {}
    for o in range(r):
        qs = np.flatnonzero(cnt[o])
        R = Rlab[qs[0]]
        assert all(Rlab[q] == R for q in qs)
        groups.setdefault(R, []).append(o)
    # V_R: achievable count vectors restricted to the points of R
    Rpts = {R: np.flatnonzero(Rlab == R) for R in groups}
    VR = {}
    for R, os_ in groups.items():
        pts = Rpts[R]
        S = {tuple(np.zeros(len(pts), dtype=np.int64))}
        for o in os_:
            c = cnt[o, pts]
            S = {tuple(np.array(s) + c) for s in S} | {tuple(np.array(s) - c) for s in S}
        VR[R] = sorted(S)
    nvec = int(np.prod([len(v) for v in VR.values()]))
    if verbose:
        print(f'  screen: G={G} Q={Q} |H|={H} r={r} Q-orbits={len(groups)} vectors/block(unpruned)={nvec}', flush=True)
    budget = int(tq[0])  # sum over the 4 (or 2) blocks of sum_q Y(q)^2 equals tq[0] exactly
    # budgeted vectorised enumeration over the S-bar orbits (largest value sets first)
    Rkeys = sorted(VR, key=lambda R: -len(VR[R]))
    Y = np.zeros((1, nQ), dtype=np.int16)
    ss = np.zeros(1, dtype=np.int64)
    for R in Rkeys:
        pts = Rpts[R]
        vals = np.array(VR[R], dtype=np.int16)          # (nv, |R|)
        vs = (vals.astype(np.int64) ** 2).sum(1)
        newY, newS = [], []
        for t in range(len(vals)):
            keep = ss + vs[t] <= budget
            if not keep.any():
                continue
            Z = Y[keep].copy()
            Z[:, pts] = vals[t]
            newY.append(Z)
            newS.append(ss[keep] + vs[t])
        if not newY:
            return False, dict(blocks=0)
        Y = np.concatenate(newY)
        ss = np.concatenate(newS)
        if len(Y) > max_vectors:
            return None, dict(reason=f'too many compressed vectors after pruning (>{max_vectors})')
    rs = Y.sum(1).astype(np.int64)
    if kind == 'b1':
        okr = np.isin(rs, [0, 2, -2]) if m % 2 == 0 else np.isin(rs, [1, -1])
    elif kind == 'lp':
        okr = np.isin(rs, [1, -1])
    else:
        okr = rs * rs <= 4 * m
    Y = Y[okr]
    rs = rs[okr]
    keepP, pafs = [], []
    ax = tuple(range(1, len(Q.mods) + 1))
    for c0 in range(0, len(Y), 200000):
        Yc = Y[c0:c0 + 200000].astype(float).reshape((-1,) + Q.mods)
        Fc = np.fft.fftn(Yc, axes=ax)
        Pc = (np.abs(Fc) ** 2).reshape(len(Yc), -1)
        ok = (Pc[:, 1:] <= bound + 1e-6).all(1)
        keepP.append(np.flatnonzero(ok) + c0)
        if ok.any():
            pafs.append(np.rint(np.fft.ifftn(Pc[ok].reshape((-1,) + Q.mods), axes=ax).real.reshape(int(ok.sum()), -1)).astype(np.int64))
    keepP = np.concatenate(keepP) if keepP else np.zeros(0, dtype=int)
    Y, rs = Y[keepP], rs[keepP]
    paf = np.concatenate(pafs) if pafs else np.zeros((0, nQ), dtype=np.int64)
    if verbose:
        print(f'  screen: compressed blocks passing row-sum+PSD: {len(Y)}', flush=True)
    # exact join on PAF over Q (all q including 0)
    rng = np.random.default_rng(7)
    Rh = rng.integers(1, 2 ** 62, size=nQ, dtype=np.int64).astype(np.uint64)
    h = (paf.astype(np.uint64) * Rh[None, :]).sum(1)
    ht = (tq.astype(np.uint64) * Rh).sum()
    idx = np.arange(len(Y))
    if kind == 'lp':
        need = ht - h
        hs = np.sort(h)
        pos = np.searchsorted(hs, need)
        found = bool(((pos < len(hs)) & (hs[np.minimum(pos, len(hs) - 1)] == need)).any())
        return found, dict(blocks=len(Y))
    if kind == 'b1' and m % 2 == 0:
        A_idx, rest = idx[np.abs(rs) == 2], idx[rs == 0]
        lists = (A_idx, rest, rest, rest)
    elif kind == 'b1':
        lists = (idx, idx, idx, idx)
    else:
        lists = None
    if lists is not None:
        IA, IB, IC, ID = lists
        if min(map(len, lists)) == 0:
            return False, dict(blocks=len(Y))
        keysCD = (h[IC][:, None] + h[ID][None, :]).ravel() if len(IC) * len(ID) < 4e8 else None
        if keysCD is None:
            return None, dict(reason='join too large')
        keysCD = np.unique(keysCD)
        for s in range(0, len(IA), 256):
            A = IA[s:s + 256]
            need = (ht - h[A][:, None] - h[IB][None, :]).ravel()
            pos = np.searchsorted(keysCD, need)
            if ((pos < len(keysCD)) & (keysCD[np.minimum(pos, len(keysCD) - 1)] == need)).any():
                return True, dict(blocks=len(Y))
        return False, dict(blocks=len(Y))
    # b0: any four blocks with sum a^2 = 4m
    vals = sorted(set(np.abs(rs).tolist()))
    for a, b, c in itertools.combinations_with_replacement(vals, 3):
        d2 = 4 * m - a * a - b * b - c * c
        d = int(round(np.sqrt(max(d2, 0))))
        if d2 < 0 or d * d != d2 or d not in vals or d < c:
            continue
        IA, IB, IC, ID = [idx[np.abs(rs) == x] for x in (a, b, c, d)]
        keysCD = np.unique((h[IC][:, None] + h[ID][None, :]).ravel())
        need = (ht - h[IA][:, None] - h[IB][None, :]).ravel()
        pos = np.searchsorted(keysCD, need)
        if ((pos < len(keysCD)) & (keysCD[np.minimum(pos, len(keysCD) - 1)] == need)).any():
            return True, dict(blocks=len(Y))
    return False, dict(blocks=len(Y))


if __name__ == '__main__':
    comps = eval(sys.argv[1])
    gens = json.loads(sys.argv[2])
    kind = sys.argv[3]
    Hcomp = eval(sys.argv[4])
    G = AbGroup(comps)
    AG = AutGroup(G)
    auts = [AG.to_aut(tuple(g)) for g in gens]
    print(screen(G, auts, kind, Hcomp))
