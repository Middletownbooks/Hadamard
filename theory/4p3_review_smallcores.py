#!/usr/bin/env python3
r"""[M] Small random k4 cores for the review of report/openai_4p3_answer.md (Claim 2, the universal border, by the OpenAI
model): cores on groups whose odd part K is NOT cyclic (Z3 x Z3, Z3 x Z3 x Z5 ...), found by simulated annealing, so
that the border theorem is also tested end to end where K is non-cyclic.  New code (Claude Code review); imports no
repository code.  A found core is checked here by direct integer PAF at every lag; its border and matrix are then
handled by theory/4p3_review_border.py --elliptic FILE and theory/4p3_review_hverify.py.

Objective: sum over g != 0 of (sum_i PAF_i(g) - F(g))^2, F = -12 on K \ 0 and +4 off K; single-bit flips, Metropolis.
usage: python3 theory/4p3_review_smallcores.py --mods 2,2,3,3 --q v4 --out FILE [--seconds 300] [--seed 1]
"""
import argparse, itertools, json, math, sys, time
import numpy as np


def group(mods):
    n = int(np.prod(mods))
    co = np.array(np.unravel_index(np.arange(n), mods)).T
    M = np.array(mods)
    idx = lambda c: np.ravel_multi_index(tuple(np.moveaxis(np.mod(c, M), -1, 0)), mods)
    add = idx(co[:, None, :] + co[None, :, :])
    neg = idx(-co)
    return n, co, add, neg


def slice_search(mods, q, M, seed=1):
    """exhaustive search inside the slice of <M>-invariant sequences (M an automorphism of G given as an integer
    matrix on the coordinates): all 2^r orbit-sign vectors per block, distinct PAF vectors, 2 + 2 join on a linear
    64-bit hash (exact under wraparound), exact check.  Returns (r, number of distinct PAF vectors, a core or None)."""
    n, co, add, neg = group(mods)
    Mx = np.array(M, np.int64)
    Marr = np.array(mods)
    img = np.ravel_multi_index(tuple(np.moveaxis(np.mod(co @ Mx.T, Marr), -1, 0)), mods)
    assert len(set(img.tolist())) == n
    orb = -np.ones(n, np.int64); r = 0
    for g in range(n):
        if orb[g] < 0:
            h = g
            while orb[h] < 0:
                orb[h] = r; h = img[h]
            r += 1
    qidx = (2 * co[:, 0] + co[:, 1]) if q == 'v4' else co[:, 0] % 4
    F = np.where(qidx == 0, -12, 4).astype(np.int64); F[0] = 4 * n
    E = np.array(list(itertools.product((1, -1), repeat=r)), np.int64)
    X = E[:, orb]                                                   # all invariant sequences
    P = np.stack([(X * X[:, add[:, g]]).sum(1) for g in range(n)], 1)
    U, first = np.unique(P, axis=0, return_index=True)
    rng = np.random.default_rng(seed)
    coef = rng.integers(1, 2 ** 62, size=n, dtype=np.int64)
    hU = U @ coef
    I, J = np.triu_indices(len(U))
    hp = hU[I] + hU[J]
    o = np.argsort(hp, kind='stable'); hs = hp[o]
    need = np.int64(int(F @ coef)) - hp
    lo = np.searchsorted(hs, need, 'left'); hi = np.searchsorted(hs, need, 'right')
    for x in np.flatnonzero(hi > lo):
        for y in o[lo[x]:hi[x]]:
            a, b, c, d = int(I[x]), int(J[x]), int(I[y]), int(J[y])
            if np.array_equal(U[a] + U[b] + U[c] + U[d], F):
                xs = [X[first[k]] for k in (a, b, c, d)]
                S = sum((v[:, None] * v[add]).sum(0) for v in xs)
                assert np.array_equal(S, F)
                return r, len(U), xs
    return r, len(U), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slice', default='', help='exhaustive slice search: automorphism matrix rows, e.g. '
                    '"1,0,0,0;0,1,0,0;0,0,0,2;0,0,1,0"')
    ap.add_argument('--mods', required=True)
    ap.add_argument('--q', choices=('v4', 'z4'), required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seconds', type=float, default=300)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument("--want", type=int, default=1)
    ap.add_argument("--tmin", type=float, default=6.0)
    ap.add_argument("--iters", type=int, default=2000000)
    a = ap.parse_args()
    mods = tuple(int(x) for x in a.mods.split(','))
    if a.slice:
        t0 = time.time()
        M = [[int(v) for v in row.split(',')] for row in a.slice.split(';')]
        r, nu, xs = slice_search(mods, a.q, M, seed=a.seed)
        res = dict(script='theory/4p3_review_smallcores.py --slice', command=' '.join(sys.argv), mods=list(mods), q=a.q,
                   automorphism=M, orbits_per_block=r, distinct_paf_vectors=nu, found=int(xs is not None),
                   secs=round(time.time() - t0, 1),
                   cores=[] if xs is None else [dict(name='slice core %s %s' % (a.q, 'x'.join(map(str, mods))),
                                                     mods=list(mods), q=a.q, sequences=[v.tolist() for v in xs],
                                                     source='theory/4p3_review_smallcores.py --slice %s' % a.slice)])
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
        print(json.dumps({k: v for k, v in res.items() if k != 'cores'}))
        return
    n, co, add, neg = group(mods)
    qidx = (2 * co[:, 0] + co[:, 1]) if a.q == 'v4' else co[:, 0] % 4
    assert (a.q == 'v4' and mods[:2] == (2, 2)) or (a.q == 'z4' and mods[0] == 4)
    F = np.where(qidx == 0, -12, 4).astype(np.int64)
    F[0] = 4 * n
    rng = np.random.default_rng(a.seed)
    t0 = time.time()
    found = []
    restarts = 0
    while time.time() - t0 < a.seconds and len(found) < a.want:
        restarts += 1
        X = rng.choice(np.array([-1, 1], np.int64), size=(4, n))
        S = sum((X[i][:, None] * X[i][add]).sum(0) for i in range(4))
        D = S - F
        D[0] = 0
        cost = int((D * D).sum())
        T = 40.0
        for it in range(a.iters):
            i = int(rng.integers(4)); h = int(rng.integers(n))
            x = X[i]
            # dS[g] for g != 0: -2 x(h) (x(h + g) + x(h - g)), with h + g = add[h, g] and h - g = add[h, neg[g]]
            dS = -2 * x[h] * (x[add[h]] + x[add[h][neg]])
            dS[0] = 0
            nd = D + dS
            nc = int((nd * nd).sum())
            if nc <= cost or rng.random() < math.exp((cost - nc) / T):
                x[h] = -x[h]
                D = nd; cost = nc
                if cost == 0:
                    break
            T = max(a.tmin, T * 0.99997)
        if cost == 0:
            xs = [X[i].copy() for i in range(4)]
            S = sum((x[:, None] * x[add]).sum(0) for x in xs)
            assert np.array_equal(S, F), 'annealer returned a non-core'
            found.append([x.tolist() for x in xs])
    res = dict(script='theory/4p3_review_smallcores.py', command=' '.join(sys.argv), seed=a.seed, mods=list(mods),
               q=a.q, restarts=restarts, found=len(found), secs=round(time.time() - t0, 1),
               cores=[dict(name='random core %s %s #%d' % (a.q, 'x'.join(map(str, mods)), k), mods=list(mods), q=a.q,
                           sequences=s, source='theory/4p3_review_smallcores.py seed %d' % a.seed)
                      for k, s in enumerate(found)])
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'cores'}))


if __name__ == '__main__':
    main()
