#!/usr/bin/env python3
"""The slice 599 = 16*37 + 7 over G = Z_2^4 x Z_37 with H = <sigma> x <10> (sigma of order 3 in GL(4,2) without
nonzero fixed points, 10 of order 3 mod 37; 78 orbits per block) and a non-classical character set T:
all Q-level solutions, written out and re-checked from scratch.  [M]
Q-level: cell sums y_k(q) (constant on the sigma-orbits of Q, values = +-1 + 3m, |m| <= 12, m even) with
sum_k yhat_k(rho)^2 = 28 on T and 2396 off T (16-point Walsh transforms).
Usage: fable_q599.py [T as comma list]   -> logs/fable_fresh/q599_solutions.json
"""
import sys, json, itertools, math, collections
import numpy as np
sys.path.insert(0, 'theory')
import fable_tcensus2 as C2
from fable_tclasses import kind

def main():
    p = 599; t = 7; n = 4 * p; K = 16
    T = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0, 1, 2, 3, 4, 8, 12]
    sig = [x for x in C2.gl_odd_classes() if (x[0], x[1]) == (3, 1)][0]
    o, fix, img, dual = sig
    qorb = C2.perm_orbits(img); corb = C2.perm_orbits(dual)
    print('sigma on cells:', qorb); print('sigma on characters:', corb); print('T =', T, kind(T, K))
    Tset = set(T)
    assert all((set(c) <= Tset) or not (set(c) & Tset) for c in corb), 'T not sigma-invariant'
    Y = [e + 3 * m for e in (1, -1) for m in range(-12, 13, 2)]
    Y = sorted(set(Y))
    N = len(qorb)
    walsh = lambda y: [sum((-1) ** bin(r & q).count('1') * y[q] for q in range(K)) for r in range(K)]
    # enumerate blocks: y constant on sigma-orbits
    lowb = 5
    blocks = []
    Yv = np.array(Y)
    # N = 6, |Y| = 26: 3.1e8 -> enumerate in chunks over the first coordinate
    chars = np.array([[(-1) ** bin(r & q).count('1') for q in range(K)] for r in range(K)], dtype=np.int64)
    cellidx = np.zeros(K, dtype=int)
    for j, ob in enumerate(qorb):
        for q in ob: cellidx[q] = j
    P = np.zeros((K, N), dtype=np.int64)
    for r in range(K):
        for q in range(K): P[r, cellidx[q]] += chars[r, q]
    Tm = np.array([r in Tset for r in range(K)])
    rest = np.array(list(itertools.product(range(len(Y)), repeat=N - 2)), dtype=np.int16)
    Vrest = Yv[rest].astype(np.int32)
    Wrest = Vrest @ P[:, 2:].T.astype(np.int32)
    for y0 in Y:
        for y1 in Y:
            W = Wrest + (y0 * P[:, 0] + y1 * P[:, 1]).astype(np.int32)[None, :]
            ok = np.all(np.abs(W[:, Tm]) <= lowb, axis=1)
            idxs = np.nonzero(ok)[0]
            if len(idxs) == 0: continue
            W2 = W[idxs]
            ok2 = np.all(W2[:, ~Tm].astype(np.int64) ** 2 <= n, axis=1)
            for idx, w in zip(idxs[ok2], W2[ok2]):
                blocks.append((tuple([y0, y1] + Vrest[idx].tolist()), tuple(int(x) for x in w)))
    print(len(blocks), 'admissible block vectors')
    E = np.array([[w * w for w in b[1]] for b in blocks], dtype=np.int64)
    tgt = np.array([28 if r in Tset else n for r in range(K)], dtype=np.int64)
    # 2 + 2 join (exact, dictionary on pair energy vectors)
    d = collections.defaultdict(list)
    for i in range(len(blocks)):
        S2 = E[i][None, :] + E[i:]
        good = np.nonzero(np.all(S2 <= tgt[None, :], axis=1))[0]
        for j in good: d[S2[j].tobytes()].append((i, i + int(j)))
    sols = []
    for key, prs in d.items():
        comp = (tgt - np.frombuffer(key, dtype=np.int64)).astype(np.int64).tobytes()
        if comp in d and key <= comp:
            for pr in prs:
                for qr in d[comp]:
                    if key != comp or pr <= qr: sols.append(pr + qr)
    print(len(sols), 'unordered solutions (multisets of two pairs of blocks)')
    out = []
    for s4 in sols:
        ys = [[blocks[i][0][cellidx[q]] for q in range(K)] for i in s4]
        # independent re-check
        tot = [sum(walsh(y)[r] ** 2 for y in ys) for r in range(K)]
        assert tot == [28 if r in Tset else n for r in range(K)], tot
        assert sum(v * v for y in ys for v in y) == ((K - t) * n + 4 * t * t) // K
        out.append(dict(cells=[list(blocks[i][0]) for i in s4], walsh=[walsh(y) for y in ys]))
    json.dump(dict(p=p, T=T, cell_orbits=qorb, char_orbits=corb, solutions=out), open('logs/fable_fresh/q599_solutions.json', 'w'))
    cnt = collections.Counter(tuple(sorted(tuple(sorted(map(abs, c))) for c in s['cells'])) for s in out)
    print('distinct multisets of |cell sums| per block:', len(cnt))
    for k_, v in list(cnt.items())[:12]: print('  ', v, k_)
    if out:
        s0 = out[0]
        print('example cells (per sigma-orbit of Q, sizes', [len(x) for x in qorb], '):', s0['cells'])
        print('Walsh values on T:', [[w[r] for r in T] for w in s0['walsh']])
        print('Walsh values off T:', [[w[r] for r in range(K) if r not in Tset] for w in s0['walsh']])
if __name__ == '__main__':
    main()
