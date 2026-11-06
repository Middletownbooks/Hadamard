r"""Second, independent stage-1 enumerator for any b1 / k4 slice key of the open-order screen [M].

Shares no code with theory/open_screen_stage1_exact.py, ab_lib, ab_sweep, ab_b0_stage1 (only the key format).
Key ["fam", v, Q, H1, q, s, U]: G = Q x H1 x Z_q, G' = Q x H1 (axes in that order, row-major index), M of index
s in F_q^*, f = (q-1)/s.  U entries: [axis, order] (the subgroup of that order of the units of the axis -- the
first generator found, as the screen's 2-entry keys; only used on axes where it is unique, i.e. cyclic units),
[axis, order, gen] (generator given), or ["neg", [axes]] (x -> -x on those axes).
Stage-1 target T1(h) = sum_z T(h, z), from the family's target on G (|G| = m):
  b1: T = 4m at 0, -4 elsewhere;   k4: T = 4m at 0, -12 on K \ 0 (K = H1 x Z_q, i.e. Q-part 0), +4 off K.
Letters alpha + f sigma, alpha = +-1, sigma a sum of s signs.  y is constant on U-orbits of G'.
Enumeration: level-wise over the U-orbits (largest first), keeping partial vectors whose energy plus the least
possible rest is <= T1(0) - 3 E_min and whose |row| minus the largest possible rest is <= isqrt(sum T1); at the end
r^2 <= sum T1.  [P: each block of a solution has energy T1(0) minus the other three energies, each >= E_min, and
r_k^2 <= sum_k r_k^2 = sum_h T1(h).]  Block sign normalised (y[0] > 0).
PAF on one lag per U-orbit (PAF is constant on U-orbits of lags): sum_{i,j} v_i v_j N_ij(h), N counted on G'.
Join: for every row multiset {|r_a|, |r_b|, |r_c|, |r_d|} with sum of squares = sum T1 and every split of the lag-0
energy T1(0) = e_left + e_right (bucket by bucket), unordered pairs of distinct block-PAF vectors, a 64-bit linear hash (exact mod 2^64) sorted on the right, the left searched, hits re-checked on
the exact vectors; solutions expanded to blocks; every solution re-checked on all lags of G' by the definition.
Caps: --cap blocks (default 3e6) and --maxpairs per energy bucket (default 4e8); hitting one reports UNDECIDED.
Planted: --planted N quadruples of enumerated blocks with sum of rows^2 = sum T1 and total energy <= T1(0), so the
list is complete for the planted target; recovery of the planted canonical quadruple is required.
usage: python3 theory/h4892_stage1_any.py 'KEY' [--planted N] [--out FILE] [--full FILE] [--cap C] [--maxpairs P]"""
import sys, json, time, itertools, math
from collections import defaultdict
import numpy as np


def order_mod(u, n):
    x, k = u % n, 1
    while x != 1: x = x * u % n; k += 1
    return k


def setup(key):
    fam, v, Q, H1, q, s, U = key
    axes = list(Q) + list(H1); nQ = len(Q)
    els = list(itertools.product(*[range(a) for a in axes])) if axes else [()]
    n = len(els); idx = {e: i for i, e in enumerate(els)}
    mults = []
    for g in U:
        mu = [1] * len(axes)
        if g[0] == 'neg':
            for a in g[1]: mu[a] = -1
        elif len(g) == 3: mu[g[0]] = g[2]; assert order_mod(g[2], axes[g[0]]) == g[1]
        else:
            a, o = g; na = axes[a]
            mu[a] = next(u for u in range(2, na) if math.gcd(u, na) == 1 and order_mod(u, na) == o)
        mults.append(mu)
    act = lambda mu, e: tuple((x * c) % na for x, c, na in zip(e, mu, axes))
    oid = [-1] * n; orbs = []
    for i, e in enumerate(els):
        if oid[i] >= 0: continue
        S = {e}; st = [e]
        while st:
            x = st.pop()
            for mu in mults:
                y = act(mu, x)
                if y not in S: S.add(y); st.append(y)
        for x in S: oid[idx[x]] = len(orbs)
        orbs.append(sorted(idx[x] for x in S))
    f = (q - 1) // s
    let = sorted({al + f * sg for al in (1, -1) for sg in range(-s, s + 1, 2)})
    m = n * q
    T1 = []
    for e in els:
        t = 0
        for z in range(q):
            if all(x == 0 for x in e) and z == 0: t += 4 * m
            elif fam == 'b1': t += -4
            elif fam == 'k4': t += -12 if all(x == 0 for x in e[:nQ]) else 4
            else: raise ValueError(fam)
        T1.append(t)
    addi = [[idx[tuple((x + y) % na for x, y, na in zip(e, h, axes))] for e in els] for h in els]
    k = len(orbs); lag = [o[0] for o in orbs]
    N = np.zeros((k, k, k), np.int64)
    for l, h in enumerate(lag):
        for g in range(n): N[l, oid[g], oid[addi[h][g]]] += 1
    return dict(n=n, orbs=orbs, oid=oid, let=let, T1=T1, addi=addi, N=N, W=np.array([len(o) for o in orbs]), k=k,
                lagT=np.array([T1[h] for h in lag], np.int64))


def blocks(st, T0, R, cap):
    W = st['W']; L = np.array(st['let'], np.int64); k = st['k']
    m0 = int(np.min(L * L)); big = int(np.max(np.abs(L)))
    Emax = T0 - 3 * m0 * st['n']; Rb = math.isqrt(R)
    order = np.argsort(-W, kind='stable')
    rest = [int(W[order[t + 1:]].sum()) for t in range(k)]
    V = np.zeros((1, 0), np.int64); E = np.zeros(1, np.int64); r = np.zeros(1, np.int64)
    for t in range(k):
        w = int(W[order[t]]); nV, nE, nr = [], [], []; tot = 0
        for c0 in range(0, len(V), 1 << 16):                  # chunked expansion (memory)
            v_, e_, r_ = V[c0:c0 + (1 << 16)], E[c0:c0 + (1 << 16)], r[c0:c0 + (1 << 16)]
            v2 = np.concatenate([np.repeat(v_, len(L), 0), np.tile(L, len(v_))[:, None]], 1)
            e2 = np.repeat(e_, len(L)) + w * np.tile(L * L, len(e_)); r2 = np.repeat(r_, len(L)) + w * np.tile(L, len(r_))
            keep = (e2 + m0 * rest[t] <= Emax) & (np.abs(r2) - big * rest[t] <= Rb)
            nV.append(v2[keep]); nE.append(e2[keep]); nr.append(r2[keep]); tot += int(keep.sum())
            if tot > cap: raise OverflowError(f'more than {cap} partial blocks at level {t}')
        V, E, r = np.concatenate(nV), np.concatenate(nE), np.concatenate(nr)
    keep = r * r <= R
    V = V[keep]; r = r[keep]
    out = np.zeros_like(V); out[:, order] = V
    out = out[out[:, st['oid'][0]] > 0]                                  # block sign: y[0] > 0
    return out


def solve(st, lagT, R, cap, maxpairs):
    V = blocks(st, int(lagT[0]), R, cap)
    P = np.einsum('bi,bj,lij->bl', V, V, st['N'])
    rows = np.abs(V @ st['W'])
    keys, inv = np.unique(P, axis=0, return_inverse=True); inv = inv.ravel()
    krow = np.zeros(len(keys), np.int64); krow[inv] = rows
    members = defaultdict(list)
    for i, kk in enumerate(inv): members[int(kk)].append(i)
    rng = np.random.default_rng(7); hv = rng.integers(-2 ** 62, 2 ** 62, size=st['k'], dtype=np.int64)
    H = (keys * hv).sum(1); tH = np.int64((lagT * hv).sum())
    rset = sorted(set(krow.tolist()))
    combos = {tuple(sorted(c)) for c in itertools.product(rset, repeat=4) if sum(x * x for x in c) == R}
    sols = set(); stats = []
    ken = keys[:, 0]                                           # lag-0 value = block energy
    grp = defaultdict(list)                                    # (row, energy) -> key indices
    for i in range(len(keys)): grp[(int(krow[i]), int(ken[i]))].append(i)
    grp = {g: np.array(v, np.int64) for g, v in grp.items()}
    T0 = int(lagT[0])
    def pairs_e(ra, rb, e):
        """all unordered key pairs with rows (ra, rb) and energy sum e (memory: one energy bucket at a time)"""
        I, J = [], []
        for (r1, e1), A in grp.items():
            if r1 != ra or (rb, e - e1) not in grp: continue
            B = grp[(rb, e - e1)]
            X, Y = np.meshgrid(A, B, indexing='ij'); X = X.ravel(); Y = Y.ravel()
            if ra == rb: mk = X <= Y; X, Y = X[mk], Y[mk]
            I.append(X); J.append(Y)
        if not I: return np.zeros(0, np.int64), np.zeros(0, np.int64)
        return np.concatenate(I), np.concatenate(J)
    for c in sorted(combos):
        for (ra, rb), (rc, rd) in {((c[0], c[1]), (c[2], c[3])), ((c[0], c[2]), (c[1], c[3])), ((c[0], c[3]), (c[1], c[2]))}:
            eL = sorted({e1 + e2 for (r1, e1) in grp if r1 == ra for (r2, e2) in grp if r2 == rb})
            eR = {e1 + e2 for (r1, e1) in grp if r1 == rc for (r2, e2) in grp if r2 == rd}
            nL = nR = 0
            for e in eL:
                if T0 - e not in eR: continue
                I, J = pairs_e(rc, rd, T0 - e)
                if len(I) > maxpairs: raise OverflowError(f'{len(I)} pairs in one energy bucket')
                rh = H[I] + H[J]; o = np.argsort(rh); rs = rh[o]; nR += len(rs)
                A_, B_ = pairs_e(ra, rb, e); nL += len(A_)
                if len(rs) == 0 or len(A_) == 0: continue
                for c0 in range(0, len(A_), 1 << 22):
                    a_, b_ = A_[c0:c0 + (1 << 22)], B_[c0:c0 + (1 << 22)]
                    need = tH - (H[a_] + H[b_])
                    pos = np.minimum(np.searchsorted(rs, need), len(rs) - 1)
                    for h in np.nonzero(rs[pos] == need)[0]:
                        lo, hi = np.searchsorted(rs, need[h], 'left'), np.searchsorted(rs, need[h], 'right')
                        for qq in o[lo:hi]:
                            ks = (a_[h], b_[h], I[qq], J[qq])
                            if not (sum(keys[x] for x in ks) == lagT).all(): continue
                            for t4 in itertools.product(*[members[int(x)] for x in ks]):
                                sols.add(tuple(sorted(tuple(int(z) for z in V[x]) for x in t4)))
            stats.append(dict(rows=[ra, rb, rc, rd], left=int(nL), right=int(nR)))
    return sols, V, stats


def expand(st, v):
    y = [0] * st['n']
    for j, o in enumerate(st['orbs']):
        for g in o: y[g] = int(v[j])
    return y


def main():
    a = sys.argv[1:]; key = json.loads(a[0])
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 2
    cap = int(float(a[a.index('--cap') + 1])) if '--cap' in a else 3 * 10 ** 6
    mp = int(float(a[a.index('--maxpairs') + 1])) if '--maxpairs' in a else 4 * 10 ** 8
    outp = a[a.index('--out') + 1] if '--out' in a else None
    t0 = time.time(); st = setup(key); R = int(sum(st['T1']))
    rec = dict(key=key, script='theory/h4892_stage1_any.py', U_orbits=st['k'], letters=len(st['let']), sumT1=R)
    try:
        sols, V, stats = solve(st, st['lagT'], R, cap, mp)
        for Y in sols:
            ys = [expand(st, y) for y in Y]
            tot = [sum(y[g] * y[st['addi'][h][g]] for y in ys for g in range(st['n'])) for h in range(st['n'])]
            assert tot == st['T1']
        rec.update(blocks=int(len(V)), solutions=len(sols), join=stats, secs=round(time.time() - t0, 1),
                   canonical=[[expand(st, y) for y in Y] for Y in sorted(sols)][:50])
        if '--full' in a:                                      # every canonical solution, as U-orbit value vectors
            json.dump(dict(key=key, orbits=st['orbs'], solutions=[list(map(list, Y)) for Y in sorted(sols)]),
                      open(a[a.index('--full') + 1], 'w'))
    except (OverflowError, MemoryError) as e:
        rec.update(undecided=repr(e)[:150], secs=round(time.time() - t0, 1)); V = None
    print(json.dumps({x: y for x, y in rec.items() if x != 'canonical'}), flush=True)
    pl = []
    if V is not None and len(V):
        rng = np.random.default_rng(99); rows = V @ st['W']; EE = (V * V) @ st['W']
        for t in range(npl):
            got = None
            for _ in range(200000):
                ix = rng.integers(len(V), size=4)
                if int((rows[ix] ** 2).sum()) == R and int(EE[ix].sum()) <= st['lagT'][0]: got = ix; break
            if got is None: pl.append(dict(skipped='no plant with sum r^2 = sum T1 and energy <= T1(0) found')); continue
            Y = [V[i] * (1 if rng.random() < .5 else -1) for i in got]
            Pp = sum(np.einsum('i,j,lij->l', y, y, st['N']) for y in Y)
            try:
                sp, _, _ = solve(st, Pp, R, cap, mp)
                canon = tuple(sorted(tuple(int(z) for z in (y if y[st['oid'][0]] > 0 else -y)) for y in Y))
                pl.append(dict(recovered=canon in sp, solutions=len(sp)))
            except (OverflowError, MemoryError) as e: pl.append(dict(undecided=repr(e)[:100]))
    rec['planted'] = pl
    print(json.dumps(dict(planted=pl)), flush=True)
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
