r"""Exact stage 1 by direct enumeration, for slices with few G'-orbits and a tight lag-0 energy [M].

Per block: every letter vector y (one letter per G'-orbit) with weighted energy sum_b w_b y_b^2 <= T1(0) - 3 Emin,
built orbit by orbit with energy pruning; one representative of {y, -y}.  Then the energy-split hashed 4-list join
(theory/ab_s1sd.ejoin4) of the PAF over G' at its lag classes against T1 (the b >= 3 or b = 0 stage-1 target), every
match re-checked exactly.  Output: all stage-1 solutions (4 x |G'| arrays), canonical up to block order and signs.
usage: python3 theory/ab_s1_direct.py '<cfg json>' --out FILE [--max-out N]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD


def block_vectors(sl, Emax):
    L = [np.array(sorted(l), np.int64) for l in sl.letters]; w = sl.bsz.astype(np.int64); nb = len(L)
    minE = [int(w[j] * (L[j] ** 2).min()) for j in range(nb)]
    rest = np.concatenate([np.cumsum(minE[::-1])[::-1][1:], [0]])           # min energy of orbits j+1..
    cur = np.zeros((1, 0), np.int64); E = np.zeros(1, np.int64)
    for j in range(nb):
        e_j = w[j] * L[j] ** 2
        En = (E[:, None] + e_j[None, :]).ravel()
        ok = En + rest[j] <= Emax
        idx = np.nonzero(ok)[0]
        cur = np.concatenate([np.repeat(cur, len(L[j]), 0)[idx], np.tile(L[j], len(E))[idx][:, None]], 1)
        E = En[idx]
    # one of y, -y: first nonzero letter positive (letters are odd, so the first entry decides)
    keep = cur[:, 0] > 0
    return cur[keep], E[keep]


def ejoin4_sym(K, T, chunk=1 << 22, max_out=10 ** 7):
    """unordered quadruples {i, j, k, l} (with repetition) of rows of K summing to T; coordinate 0 = energy.
    Left pairs i <= j with energy EL <= T0/2, right pairs k <= l with energy T0 - EL; per energy split the smaller
    side is sorted and the larger streamed in chunks; every match re-checked exactly"""
    n, d = K.shape
    wv = np.random.default_rng(99).integers(1, 2 ** 62, size=d, dtype=np.int64)
    H = K @ wv; E = K[:, 0]; T0 = int(T[0]); th = np.int64(int(T.astype(np.int64) @ wv))
    order = np.argsort(E, kind='stable'); Es = E[order]
    ue, st = np.unique(Es, return_index=True); en = np.r_[st[1:], len(Es)]
    grp = {int(a): order[b:c] for a, b, c in zip(ue, st, en)}
    es = sorted(grp)
    def pairs_with_sum(S):
        out = []
        for ea in es:
            eb = S - ea
            if eb < ea: break
            if eb in grp: out.append((ea, eb))
        return out
    def pair_arrays(plist):
        """yield (hash, i, j) chunks for all pairs i <= j (i in group ea, j in group eb)"""
        for ea, eb in plist:
            A = grp[ea]; B = grp[eb]
            for c0 in range(0, len(A), max(1, chunk // max(1, len(B)))):
                a = A[c0:c0 + max(1, chunk // max(1, len(B)))]
                I = np.repeat(a, len(B)); J = np.tile(B, len(a))
                if ea == eb:
                    keep = I <= J; I, J = I[keep], J[keep]
                yield H[I] + H[J], I, J
    def size(plist):
        s = 0
        for ea, eb in plist:
            na, nb = len(grp[ea]), len(grp[eb]); s += na * (na + 1) // 2 if ea == eb else na * nb
        return s
    out = []; sums = sorted({ea + eb for ea in es for eb in es if ea <= eb})
    for EL in sums:
        if EL > T0 - EL: break
        L = pairs_with_sum(EL); R = pairs_with_sum(T0 - EL)
        if not L or not R: continue
        if size(L) > size(R): L, R, flip = R, L, True
        else: flip = False
        parts = list(pair_arrays(L))
        LH = np.concatenate([p[0] for p in parts]); LI = np.concatenate([p[1] for p in parts]); LJ = np.concatenate([p[2] for p in parts])
        o = np.argsort(LH); LS = LH[o]
        for RH0, RI, RJ in pair_arrays(R):
            RH = th - RH0
            lo = np.searchsorted(LS, RH, 'left'); hi = np.searchsorted(LS, RH, 'right')
            for r in np.nonzero(hi > lo)[0]:
                for p in o[lo[r]:hi[r]]:
                    q = (int(LI[p]), int(LJ[p]), int(RI[r]), int(RJ[r]))
                    if ((K[q[0]] + K[q[1]] + K[q[2]] + K[q[3]]) == T).all():
                        out.append(q)
                        if len(out) > max_out: raise RuntimeError("output cap")
    return out


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(a[0]); out = a[a.index('--out') + 1]
    max_out = int(a[a.index('--max-out') + 1]) if '--max-out' in a else 10 ** 7
    sl = W.build(cfg); t0 = time.time()
    T1 = sl.T1; Gp = sl.Gp
    Up = [Gp.act(u) for u in sl.Sp]
    reps, _ = AGs.lag_reps(Gp, Up)
    w = sl.bsz.astype(np.int64)
    Emin = int(sum(w[j] * min(y * y for y in sl.letters[j]) for j in range(len(sl.breps))))
    Emax = int(T1[0]) - 3 * Emin
    Y, E = block_vectors(sl, Emax)
    full = Y[:, sl.borb]
    P = np.zeros((len(Y), len(reps)), np.int64)
    for c0 in range(0, len(Y), 1 << 16):
        P[c0:c0 + (1 << 16)] = Gp.paf(full[c0:c0 + (1 << 16)])[:, reps]
    assert (P[:, 0] == E).all()
    print(json.dumps(dict(stage="blocks", v=cfg['v'], orbits=len(sl.breps), Emax=Emax, vectors=len(Y), secs=round(time.time() - t0, 1))), flush=True)
    # dedupe identical PAF signatures for the join, then expand
    keys, inv = np.unique(P, axis=0, return_inverse=True)
    inv = inv.ravel()
    groups = [[] for _ in range(len(keys))]
    for i, g in enumerate(inv.tolist()): groups[g].append(i)
    quads = ejoin4_sym(keys, T1[reps], max_out=max_out)
    sols = set()
    for (i0, i1, i2, i3) in quads:
        for a0 in groups[i0]:
            for a1 in groups[i1]:
                for a2 in groups[i2]:
                    for a3 in groups[i3]:
                        q = tuple(sorted((a0, a1, a2, a3)))
                        sols.add(q)
    res = []
    for q in sorted(sols):
        Yq = full[list(q)]
        if (Gp.paf(Yq).sum(0) == T1).all(): res.append(Yq.tolist())
    print(json.dumps(dict(stage="done", signature_quads=len(quads), stage1_solutions=len(res), all_verified=True,
                          secs=round(time.time() - t0, 1))), flush=True)
    json.dump(dict(cfg=json.loads(a[0]), examples=res), open(out, 'w'))


if __name__ == "__main__":
    main()
