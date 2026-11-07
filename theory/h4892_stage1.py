r"""H(4892) slice, step 1(b): an independent complete enumeration of stage 1 [M].

Slice key ["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]]: G = G' x Z13, G' = Z47 x Z2 (index (a, b) -> 2a + b),
S = U x M, U = squares of Z47^* on axis 0, M = {1, 3, 9} on Z13.  Stage 1: four U-invariant y_k on G' with letters
y = alpha + 3 sigma, alpha = +-1, sigma in {-4, -2, 0, 2, 4}, and sum_k PAF_{G'}(y_k) = T1, T1(0) = 4*1222 - 4*12 =
4840, T1(h) = -4*13 = -52 (h != 0).

Shares no code with theory/open_screen_stage1_exact.py, ab_lib, ab_b0_stage1.  Method:
  * the six U-orbits of G' are {0}x{b}, QR x {b}, QNR x {b} (QR computed as {a^2 mod 47});
  * a block is a 6-vector v of letters; PAF_{G'}(y)(h) = sum_{i,j} v_i v_j N_ij(h) with N_ij(h) = #{g in O_i :
    g + h in O_j} counted by brute force over G'; PAF is constant on U-orbits of lags, so it is recorded on one lag
    per orbit (6 numbers; QR and QNR lags agree because PAF(-h) = PAF(h) and -QR = QNR);
  * all 10^6 letter vectors per block are listed (numpy), kept if the lag-0 energy is <= 4840 - 3*94 (each other
    block has energy >= 94) and the row sum r satisfies r^2 <= 4 (sum_k r_k^2 = sum_h T1(h) = 4);
  * every ORDERED quadruple (v_0, ..., v_3) with sum_k PAF = T1 on the 6 lag classes is produced by a 2+2 join on the
    6-dim PAF vectors (dictionary of distinct block PAF vectors -> blocks; ordered pairs of distinct vectors);
  * every solution is re-checked by expanding to sequences on G' (length 94) and computing all 94 PAF values by
    the definition.
Output: data/h4892_stage1.json: the ordered list (each entry four lists of 94 integers), the canonical classes up to
block order and block sign, and the symmetry statement.
usage: python3 theory/h4892_stage1.py [--planted N] [--out data/h4892_stage1.json]"""
import sys, json, time, itertools
from collections import defaultdict
import numpy as np

P47, Q = 47, 13
LET = sorted({a + 3 * s for a in (1, -1) for s in (-4, -2, 0, 2, 4)})
n = 2 * P47
els = [(a, b) for a in range(P47) for b in range(2)]          # index 2a + b
QR = sorted({(x * x) % P47 for x in range(1, P47)})
QNR = sorted(set(range(1, P47)) - set(QR))
ORB = []
for b in range(2):
    ORB.append([2 * 0 + b]); ORB.append([2 * a + b for a in QR]); ORB.append([2 * a + b for a in QNR])
ORB = [ORB[i] for i in (0, 3, 1, 4, 2, 5)]                    # order: (0,0),(0,1),(QR,0),(QR,1),(QNR,0),(QNR,1)
W = np.array([len(o) for o in ORB])
LAGS = [o[0] for o in ORB]                                    # one lag per U-orbit
oid = np.zeros(n, int)
for j, o in enumerate(ORB): oid[o] = j
add = lambda g, h: 2 * ((g // 2 + h // 2) % P47) + ((g + h) % 2)
N = np.zeros((6, 6, 6), np.int64)                             # N[l, i, j]
for l, h in enumerate(LAGS):
    for g in range(n): N[l, oid[g], oid[add(g, h)]] += 1
T1 = np.array([4840] + [-52] * 5, np.int64)                   # on LAGS
T1full = [4840] + [-52] * (n - 1)
EMAX = 4840 - 3 * n


def expand(v):
    y = [0] * n
    for j, o in enumerate(ORB):
        for g in o: y[g] = int(v[j])
    return y


def paf_full(y):
    return [sum(y[g] * y[add(g, h)] for g in range(n)) for h in range(n)]


def blocks():
    V = np.array(list(itertools.product(LET, repeat=6)), np.int64)
    E = (V * V) @ W; r = V @ W
    V = V[(E <= EMAX) & (r * r <= 4)]
    P = np.einsum('bi,bj,lij->bl', V, V, N)
    return V, P


def join(V, P, T):
    d = defaultdict(list)
    for i, p in enumerate(map(tuple, P)): d[p].append(i)
    keys = list(d); K = np.array(keys, np.int64)
    pair = defaultdict(list)
    for a in range(len(keys)):
        S = K[a] + K                                          # ordered pairs (a, b)
        for b in range(len(keys)):
            pair[tuple(S[b])].append((a, b))
    T = tuple(int(t) for t in T); out = []
    for s, L in pair.items():
        need = tuple(t - x for t, x in zip(T, s))
        R = pair.get(need)
        if not R: continue
        for a, b in L:
            for c, e in R:
                for i in d[keys[a]]:
                    for j in d[keys[b]]:
                        for k in d[keys[c]]:
                            for l in d[keys[e]]: out.append((i, j, k, l))
    return out, len(keys)


def canon(Y):
    return tuple(sorted(tuple(y) if y[0] > 0 else tuple(-x for x in y) for y in Y))


def main():
    a = sys.argv[1:]
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    outp = a[a.index('--out') + 1] if '--out' in a else 'data/h4892_stage1.json'
    t0 = time.time()
    assert (N.sum(2) == W[None, :]).all() and len(QR) == 23
    V, P = blocks()
    print('blocks', len(V), 'secs', round(time.time() - t0, 1), flush=True)
    sols, nk = join(V, P, T1)
    print('distinct PAF vectors', nk, 'ordered solutions', len(sols), 'secs', round(time.time() - t0, 1), flush=True)
    Ys = []; PF = {}
    for s in sols:                                            # full 94-lag PAF by the definition, per used block
        for i in s:
            if i not in PF: PF[i] = paf_full(expand(V[i]))
        tot = [sum(c) for c in zip(*[PF[i] for i in s])]
        assert tot == T1full, 'recheck failed'
        Ys.append([expand(V[i]) for i in s])
    can = sorted({canon(Y) for Y in Ys})
    print('canonical classes (block order, block sign)', len(can), 'secs', round(time.time() - t0, 1), flush=True)
    # planted controls: random U-invariant quadruples from the enumerated blocks with sum r^2 = 4, target = their PAF
    rng = np.random.default_rng(4892); pl = []
    rows = V @ W
    z = np.nonzero(rows == 0)[0]; t2 = np.nonzero(np.abs(rows) == 2)[0]
    E = (V * V) @ W
    for t in range(npl):
        # one row-2 block and three row-0 blocks with total energy <= 4840: the target's lag-0 value Tp(0) is then
        # <= 4840, so the enumerated list (energy <= 4840 - 282, r^2 <= 4) contains every block of every solution
        # of the planted problem (each block of it has energy <= Tp(0) - 3*94 and sum r^2 = sum Tp = 4)
        while True:
            idx = [int(rng.choice(t2))] + [int(rng.choice(z)) for _ in range(3)]
            if E[idx].sum() <= 4840: break
        perm = rng.permutation(4); idx = [idx[p] for p in perm]
        Tp = P[idx].sum(0)
        ok = int(Tp[0]) - 3 * n <= EMAX and int((V[idx] @ W) @ (V[idx] @ W)) == 4
        sp, _ = join(V, P, Tp)
        found = tuple(idx) in set(sp)
        pl.append(dict(blocks=idx, T=[int(x) for x in Tp], list_complete_for_target=bool(ok), recovered=bool(found),
                       solutions=len(sp)))
        print('planted', t, pl[-1], flush=True)
    used = sorted({i for s_ in sols for i in s_}); pos = {i: t for t, i in enumerate(used)}
    rec = dict(key=["b1", 1223, [], [47, 2], 13, 4, [[0, 23]]], script='theory/h4892_stage1.py',
               index='G\' = Z47 x Z2, element (a, b) at index 2a + b', T1='4840 at 0, -52 elsewhere',
               blocks=int(len(V)), distinct_block_paf=nk, ordered_solutions=len(Ys), canonical_classes=len(can),
               symmetry=('ordered list: no symmetry modded out (every ordered quadruple of U-invariant y_k that solves '
                         'stage 1); canonical: up to block permutation and block negation (y[0] > 0, sorted)'),
               planted=pl, secs=round(time.time() - t0, 1),
               canonical=[list(map(list, Y)) for Y in can],
               block_table=[[int(x) for x in V[i]] for i in used],
               ordered=[[pos[i] for i in s] for s in sols],
               ordered_format=('each ordered solution is four indices into block_table; a block is its 6 orbit values '
                               'on (0,0),(0,1),(QR,0),(QR,1),(QNR,0),(QNR,1) (expand with ORB of this script)'))
    json.dump(rec, open(outp, 'w'))
    print('written', outp)


if __name__ == "__main__":
    main()
