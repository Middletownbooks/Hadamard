r"""H(4892) neighbourhood: independent stage-1 enumerator for the b1 slices on G = (Z47 x Z2) x Z13 with U = the
subgroup of order d of Z47^* (d | 46, d >= 23 so that there are at most 6 U-orbits) and M of index s in Z13^* [M].

Same method as theory/h4892_stage1.py (shares no code with the stage-1 tool): U-orbits recomputed, letters
alpha + f sigma (f = 12/s, sigma a sum of s signs), PAF by the orbit quadratic form with brute-force counts, every
letter vector of the 6 (or 4) orbits listed under the energy bound T1(0) - 3*94 and r^2 <= sum T1 = 4, ordered 2+2
join on the PAF vector over the U-orbit lag classes, every solution re-checked on all 94 lags by the definition.
T1 = 4*1222 - 4*12 = 4840 at 0 and -4*13 = -52 elsewhere (it does not depend on s).
Planted: one row-2 and three row-0 enumerated blocks of total energy <= 4840 (so the list is complete for the
planted target); recovery of the exact ordered quadruple.
usage: python3 theory/h4892_stage1_gen.py d s [--planted N] [--out FILE] [--canonical]
(--canonical: only the classes up to block order and sign, as orbit-value 6-vectors; no plants, no ordered list)"""
import sys, json, time, itertools
from collections import defaultdict
import numpy as np

p, n = 47, 94


def setup(d, s):
    f = 12 // s
    gen = next(u for u in range(2, p) if pow(u, d, p) == 1 and all(pow(u, e, p) != 1 for e in range(1, d)))
    H = sorted({pow(gen, e, p) for e in range(d)})
    orbs = []
    seen = set()
    for b in range(2):
        for a in range(p):
            if (a, b) in seen: continue
            O = sorted({(a * h % p, b) for h in H}); seen |= set(O); orbs.append([2 * x + y for x, y in O])
    sig = sorted({sum(t) for t in itertools.product((1, -1), repeat=s)})
    let = sorted({al + f * sg for al in (1, -1) for sg in sig})
    W = np.array([len(o) for o in orbs])
    oid = np.zeros(n, int)
    for j, o in enumerate(orbs): oid[o] = j
    add = lambda g, h: 2 * ((g // 2 + h // 2) % p) + ((g + h) % 2)
    lags = [o[0] for o in orbs]; k = len(orbs)
    N = np.zeros((k, k, k), np.int64)
    for l, h in enumerate(lags):
        for g in range(n): N[l, oid[g], oid[add(g, h)]] += 1
    return dict(orbs=orbs, let=let, W=W, N=N, add=add, k=k)


def main():
    a = sys.argv[1:]; d, s = int(a[0]), int(a[1])
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    outp = a[a.index('--out') + 1] if '--out' in a else f'logs/h4892/stage1_gen_d{d}_s{s}.json'
    t0 = time.time(); S = setup(d, s); W, N, k = S['W'], S['N'], S['k']
    T = np.array([4840] + [-52] * (k - 1), np.int64)
    V = np.array(list(itertools.product(S['let'], repeat=k)), np.int64)
    E = (V * V) @ W; r = V @ W
    V = V[(E <= 4840 - 3 * n) & (r * r <= 4)]
    P = np.einsum('bi,bj,lij->bl', V, V, N)
    def join(Tt):
        dd = defaultdict(list)
        for i, pv in enumerate(map(tuple, P)): dd[pv].append(i)
        keys = list(dd); K = np.array(keys, np.int64).reshape(len(keys), k)
        pair = defaultdict(list)
        for x in range(len(keys)):
            Ssum = K[x] + K
            for y in range(len(keys)): pair[tuple(Ssum[y])].append((x, y))
        out = []; Tt = tuple(int(t) for t in Tt)
        for sm, L in pair.items():
            R = pair.get(tuple(t - z for t, z in zip(Tt, sm)))
            if not R: continue
            for x, y in L:
                for u, w in R:
                    for i in dd[keys[x]]:
                        for j in dd[keys[y]]:
                            for q in dd[keys[u]]:
                                for l in dd[keys[w]]: out.append((i, j, q, l))
        return out, len(keys)
    if '--canonical' in a:                                    # unordered join over sign-normalised blocks
        keep = V[:, 0] > 0; V = V[keep]; P = P[keep]
        dd = defaultdict(list)
        for i, pv in enumerate(map(tuple, P)): dd[pv].append(i)
        keys = list(dd); K = np.array(keys, np.int64).reshape(len(keys), k); nk = len(keys)
        pair = defaultdict(list)
        for x in range(nk):
            Ssum = K[x] + K[x:]
            for y, sm in enumerate(map(tuple, Ssum)): pair[sm].append((x, x + y))
        can = set(); Tt = tuple(int(t) for t in T)
        for sm, L in pair.items():
            R = pair.get(tuple(t - z for t, z in zip(Tt, sm)))
            if not R: continue
            for x, y in L:
                for u, w in R:
                    for t4 in itertools.product(dd[keys[x]], dd[keys[y]], dd[keys[u]], dd[keys[w]]):
                        can.add(tuple(sorted(tuple(int(z) for z in V[i]) for i in t4)))
        rec = dict(script='theory/h4892_stage1_gen.py --canonical', d=d, s=s, U_orbits=k, orbit_sizes=W.tolist(),
                   letters=S['let'], blocks_signnorm=int(len(V)), distinct_block_paf=nk, canonical_classes=len(can),
                   secs=round(time.time() - t0, 1), canonical=sorted(map(list, can)))
        json.dump(rec, open(outp, 'w')); print(json.dumps({x: v for x, v in rec.items() if x != 'canonical'})); return
    sols, nk = join(T)
    add = S['add']
    def exp(v):
        y = [0] * n
        for j, o in enumerate(S['orbs']):
            for g in o: y[g] = int(v[j])
        return y
    paf = lambda y: [sum(y[g] * y[add(g, h)] for g in range(n)) for h in range(n)]
    PF = {}
    for sol in sols:
        for i in sol:
            if i not in PF: PF[i] = paf(exp(V[i]))
        assert [sum(c) for c in zip(*[PF[i] for i in sol])] == [4840] + [-52] * (n - 1)
    can = sorted({tuple(sorted(tuple(y) if y[0] > 0 else tuple(-x for x in y) for y in (exp(V[i]) for i in sol)))
                  for sol in sols})
    rng = np.random.default_rng(1223 + d + s); pl = []
    z = np.nonzero(V @ W == 0)[0]; t2 = np.nonzero(np.abs(V @ W) == 2)[0]; EE = (V * V) @ W
    for t in range(npl):
        if len(z) == 0 or len(t2) == 0: pl.append(dict(skipped='no row-0 or row-2 blocks')); continue
        for _ in range(100000):
            idx = [int(rng.choice(t2))] + [int(rng.choice(z)) for _ in range(3)]
            if EE[idx].sum() <= 4840: break
        else:
            pl.append(dict(skipped='no plant with energy <= 4840 in 1e5 draws')); continue
        idx = [idx[i] for i in rng.permutation(4)]
        sp, _ = join(P[idx].sum(0)); pl.append(dict(recovered=tuple(idx) in set(sp), solutions=len(sp)))
    rec = dict(script='theory/h4892_stage1_gen.py', d=d, s=s, U_orbits=k, orbit_sizes=W.tolist(), letters=S['let'],
               blocks=int(len(V)), distinct_block_paf=nk, ordered_solutions=len(sols), canonical_classes=len(can),
               planted=pl, secs=round(time.time() - t0, 1), canonical=[list(map(list, Y)) for Y in can])
    json.dump(rec, open(outp, 'w'))
    print(json.dumps({x: v for x, v in rec.items() if x != 'canonical'}))


if __name__ == "__main__":
    main()
