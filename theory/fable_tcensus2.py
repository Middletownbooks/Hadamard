#!/usr/bin/env python3
"""Character borders with a general character set over Z_2^4 x Z_s, multiplier group <sigma> x U:
sigma an odd-order element of GL(4,2) acting on Q = Z_2^4 (cells and characters), U a cyclic group of units of Z_s.
T must be a union of orbits of sigma on the characters of Q; |T| = t = p mod 16.

For every open prime p = 7, 11 (mod 16), every sigma (5 classes of cyclic subgroups), every U with
RMIN <= r = (#sigma-orbits on Q) * (#U-orbits on Z_s) <= RMAX, and every sigma-invariant T of size t that is not
classical (point / flat / flat minus point): decide the Q-level exactly:
   cells y_k(j) in Y (j = sigma-orbit of Q), yhat_k(i) = sum_j P[i][j] y_k(j),
   sum_k yhat_k(i)^2 = 4t (i in T), = 4p (i not in T).
Exact enumeration per block with the bound |yhat_k(i)| <= sqrt(4t) on T, then a 2 + 2 join.  [M]
Usage: fable_tcensus2.py [RMIN RMAX]
"""
import sys, json, itertools, math, time, collections
import numpy as np
sys.path.insert(0, 'theory')
from fable_tcensus import orbits, cyc_subgroups
from fable_tclasses import kind

A = 4; K = 16
def mat_apply(cols, x):
    y = 0
    for i in range(A):
        if (x >> i) & 1: y ^= cols[i]
    return y
def gl_odd_classes():
    """one generator per conjugacy class of cyclic subgroups of odd order > 1 (classes told apart by order and
    number of fixed points), as permutations of F_2^4 (cells) and of the characters (transpose-inverse action)."""
    seen = {}; out = []
    for cols in itertools.product(range(1, K), repeat=A):
        img = [mat_apply(cols, x) for x in range(K)]
        if len(set(img)) != K: continue
        # order
        cur = img[:]; o = 1
        while cur != list(range(K)):
            cur = [img[c] for c in cur]; o += 1
            if o > 15: break
        if o % 2 == 0 or o == 1 or o > 15: continue
        fix = sum(1 for x in range(K) if img[x] == x)
        if (o, fix) in seen: continue
        seen[(o, fix)] = 1
        # dual action on characters m: rho_m(sigma^{-1} q) = rho_{m'}(q), m' = sigma^{-T} m; orbits found directly:
        # m' is characterised by <m', q> = <m, sigma^{-1} q>, i.e. <m', sigma q> = <m, q> for all q
        dual = [None] * K
        dot = lambda m, q: bin(m & q).count('1') & 1
        for m in range(K):
            for m2 in range(K):
                if all(dot(m2, img[q]) == dot(m, q) for q in range(K)): dual[m] = m2; break
        out.append((o, fix, img, dual))
    return out
def perm_orbits(perm):
    seen = set(); orbs = []
    for x in range(len(perm)):
        if x in seen: continue
        o = []; y = x
        while y not in seen:
            seen.add(y); o.append(y); y = perm[y]
        orbs.append(o)
    return orbs

def qlevel(P, osz, Y, p, t, Torb):
    """exact: number of distinct pairs of pair-energy vectors that solve the Q-level (0 = dead), or None if a cap
    is hit.  Blocks are vectors in Y^N; they are enumerated by a meet in the middle on the coordinates, using the
    bound |yhat(i)| <= sqrt(4t) on the character orbits i in T."""
    N = len(osz); n = 4 * p
    tgt = np.array([4 * t if i in Torb else n for i in range(N)], dtype=np.int64)
    Trows = sorted(Torb); lowb = int(math.isqrt(4 * t))
    Yv = np.array(Y, dtype=np.int64)
    h = N // 2
    if len(Y) ** (N - h) > 4_000_000: return 'cap: block space'
    GL = np.array(list(itertools.product(range(len(Y)), repeat=h)), dtype=np.int16)
    GR = np.array(list(itertools.product(range(len(Y)), repeat=N - h)), dtype=np.int16)
    VL = Yv[GL]; VR = Yv[GR]
    WL = VL @ P[:, :h].T; WR = VR @ P[:, h:].T                 # (M, N) contributions to yhat
    # hash the right halves by their T-row vector
    dR = collections.defaultdict(list)
    TR = WR[:, Trows]
    for idx, key in enumerate(map(bytes, TR.astype(np.int64))): dR[key].append(idx)
    rng = range(-lowb, lowb + 1)
    offs = np.array(list(itertools.product(rng, repeat=len(Trows))), dtype=np.int64)
    if len(offs) * len(VL) > 200_000_000: return 'cap: offsets'
    energies = set()
    for il in range(len(VL)):
        need = offs - WL[il, Trows][None, :]                    # required right T-row vectors
        for row in need:
            lst = dR.get(row.tobytes())
            if not lst: continue
            W = WL[il][None, :] + WR[lst]
            E = W * W
            ok = np.all(E <= tgt[None, :], axis=1)
            for e in E[ok]: energies.add(e.tobytes())
            if len(energies) > 60000: return 'cap: more than 60000 admissible block energy vectors'
    if not energies: return 0
    U = np.array([np.frombuffer(e, dtype=np.int64) for e in energies])
    # 2 + 2 join with 64-bit hashes of pair sums (exact re-check of matches)
    rs = np.random.default_rng(7); wts = rs.integers(1, 1 << 40, size=N).astype(np.int64)
    np.seterr(over='ignore')
    hU = U @ wts; ht = int(tgt @ wts)
    sums = []; 
    for i in range(len(U)):
        S2 = U[i][None, :] + U[i:]
        good = np.all(S2 <= tgt[None, :], axis=1)
        if good.any(): sums.append(hU[i] + hU[i:][good])
    if not sums: return 0
    A = np.unique(np.concatenate(sums))
    B = np.int64(ht) - A
    m = np.isin(B, A)
    return int(m.sum())

def main():
    RMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    RMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    op = json.load(open('data/open_orders_lt10000.json'))['open']
    isprime = lambda m: m > 1 and all(m % d for d in range(2, int(m ** .5) + 1))
    primes = [n // 4 for n in op if n % 4 == 0 and isprime(n // 4) and (n // 4) % 16 in (7, 11)]
    sig = gl_odd_classes()
    print('sigma classes (order, fixed points):', [(o, f) for o, f, _, _ in sig], flush=True)
    done = set()
    import os
    if os.path.exists('logs/fable_fresh/tcensus2.jsonl'):
        ps = [json.loads(l)['p'] for l in open('logs/fable_fresh/tcensus2.jsonl')]
        done = set(ps) - ({ps[-1]} if ps else set())            # the last prime may be incomplete
        keep = [l for l in open('logs/fable_fresh/tcensus2.jsonl') if json.loads(l)['p'] in done]
        open('logs/fable_fresh/tcensus2.jsonl', 'w').writelines(keep)
    out = open('logs/fable_fresh/tcensus2.jsonl', 'a')
    for p in primes:
        if p in done: continue
        t = p % K; s = (p - t) // K
        for o, fix, img, dual in sig:
            qorb = perm_orbits(img); corb = perm_orbits(dual)
            N = len(qorb); osz = [len(x) for x in qorb]
            dot = lambda m, q: bin(m & q).count('1') & 1
            P = np.array([[sum((-1) ** dot(co[0], q) for q in qo) for qo in qorb] for co in corb], dtype=np.int64)
            # sigma-invariant T of size t, non-classical
            Ts = []
            for mask in range(1 << len(corb)):
                sel = [i for i in range(len(corb)) if (mask >> i) & 1]
                if sum(len(corb[i]) for i in sel) != t: continue
                T = sorted(x for i in sel for x in corb[i])
                if kind(T, K) == 'NEW': Ts.append((set(sel), T))
            if not Ts: continue
            for g, u in cyc_subgroups(s):
                sizes = orbits(s, g); r = N * len(sizes)
                if not (RMIN <= r <= RMAX): continue
                Y = {0}
                for z in sizes: Y = {y + z for y in Y} | {y - z for y in Y}
                Y = sorted(Y)
                for sel, T in Ts:
                    t0 = time.time()
                    c = qlevel(P, osz, Y, p, t, sel)
                    rec = dict(p=p, t=t, s=s, sigma=[o, fix], U=[g, u], r=r, T=T, nY=len(Y),
                               qlevel=(c if isinstance(c, str) else int(c)), secs=round(time.time() - t0, 1))
                    out.write(json.dumps(rec) + '\n'); out.flush()
                    print(rec, flush=True)
if __name__ == '__main__':
    main()
