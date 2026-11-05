r"""[M] Exact decision of stage 1 (the H-level G/Z43) of the S14 sibling at v = 863 (census key fbf4f510b2566552: the
S7 host 2e93be9718324148 plus tau = (q1, q2, u, t) -> (q2, q1, -u, t)).  Set-up of both hosts by OpenAI Codex
(theory/h3452_compress.py, kind s14); construction Kumar / EKLO.  Independent of the Codex scripts and of bt_engine.

Exhaustive: for every Q-level solution z (theory/h863_s14_stage1.q_solutions, up to block order, per-block sign and
the translation by (1,1)), every block's tau-invariant y (11 orbit values in A = {+-1, +-13, +-15, +-27, +-29, +-41,
+-43}) with the block's Q-compressions is listed, subject only to the energy bound that every solution satisfies:
sum over the 4 blocks of sum_cells (y^2 - 1) = 2936 - 80 = 2856, so block i has at most 2856 - sum_{j != i} min_j.
The four lists are joined 2 + 2 on the exact 20-lag H-level PAF vector (target 2936 / -516 / 172).  Every hit is
re-checked directly.  Output: one JSON line per Q-level solution with list sizes and the number of solutions.
usage: python3 theory/h863_s14_stage1_exact.py [first] [last] [out.jsonl]"""
import sys, os, json, time, itertools as it, bisect
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

A = [s * v for v in (1, 13, 15, 27, 29, 41, 43) for s in (1, -1)]
ASET = set(A)
SW = [0, 2, 1, 3]
TAU = [SW[k // 5] * 5 + (-(k % 5)) % 5 for k in range(20)]
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'h863_s14_qsolutions.json')


def q_solutions():
    if os.path.exists(CACHE):
        return [tuple(tuple(b) for b in z) for z in json.load(open(CACHE))['z']]
    from h863_s14_stage1 import q_solutions as qs
    Z = qs()
    json.dump(dict(note='S14 (863) Q-level solutions (a, b, c) = (z(00), z(01) = z(10), z(11)) per block, up to block '
                        'order, per-block sign and the translation by (1,1); theory/h863_s14_stage1.q_solutions',
                   z=[list(map(list, z)) for z in Z]), open(CACHE, 'w'))
    return Z


def g3(a):
    """(f0, p1, p2) with f0 + 2 p1 + 2 p2 = a: the cell values (q,0), (q,+-1), (q,+-2) for q = 00 or 11"""
    out = []
    for p1 in A:
        for p2 in A:
            f0 = a - 2 * p1 - 2 * p2
            if f0 in ASET: out.append(((f0, p1, p2), f0 * f0 + 2 * p1 * p1 + 2 * p2 * p2 - 5))
    return out


def g5(b):
    """(m0..m4) with sum b: m_u = y(01, u) = y(10, -u)"""
    out = []
    for m in it.product(A, repeat=4):
        m4 = b - sum(m)
        if m4 in ASET:
            mm = m + (m4,); out.append((mm, 2 * sum(x * x for x in mm) - 10))
    return out


def cells(x00, x11, m):
    y = [0] * 20
    f0, p1, p2 = x00; y[0], y[1], y[4], y[2], y[3] = f0, p1, p1, p2, p2
    g0, s1, s2 = x11; y[15], y[16], y[19], y[17], y[18] = g0, s1, s1, s2, s2
    for u in range(5):
        y[5 + u] = m[u]; y[10 + (-u) % 5] = m[u]
    return y


LAGS = [((((k // 5) ^ (h // 5)) * 5) + (k % 5 + h % 5) % 5) for h in range(20) for k in range(20)]
LAGS = np.array(LAGS, np.int64).reshape(20, 20)          # LAGS[h0, k] = index of cell k + h0
TGT = np.array([2936 if h == 0 else (-516 if h // 5 == 0 else 172) for h in range(20)], np.int64)


def sigs(Y):
    Y = np.asarray(Y, np.int64)
    return np.stack([(Y * Y[:, LAGS[h]]).sum(1) for h in range(20)], 1)


_G5 = {}
def block_list(a, b, c, cap):
    G00 = [x for x in g3(a) if x[1] <= cap]; G11 = [x for x in g3(c) if x[1] <= cap]
    if b not in _G5: _G5[b] = sorted(g5(b), key=lambda t: t[1])
    G01 = _G5[b]; en01 = [t[1] for t in G01]
    Y = []
    for x, ex in G00:
        for w, ew in G11:
            r = cap - ex - ew
            if r < 0: continue
            for m, em in G01[:bisect.bisect_right(en01, r)]:
                Y.append(cells(x, w, m))
    return np.array(Y, np.int64).reshape(-1, 20)


def min_energy(a, b, c):
    e00 = [x[1] for x in g3(a)]; e11 = [x[1] for x in g3(c)]
    if not e00 or not e11: return None
    if b not in _G5: _G5[b] = sorted(g5(b), key=lambda t: t[1])
    if not _G5[b]: return None
    return min(e00) + min(e11) + _G5[b][0][1]


def join22(S, T, chunk=400):
    """all index tuples (i0, i1, i2, i3) with S[0][i0] + S[1][i1] + S[2][i2] + S[3][i3] == T, exactly.  Linear 64-bit
    keys (key(u + v) = key(u) + key(v) with wraparound), the two smallest lists form the left side (sorted keys), the
    right side is streamed in chunks of its larger list; every key match is re-checked on the integer vectors."""
    R = np.random.default_rng(5).integers(1, 2 ** 62, size=T.shape[0], dtype=np.int64) | 1
    K = [(s * R).sum(1) for s in S]                       # int64 wraparound is fine: equal vectors, equal keys
    tk = int((T * R).sum())
    order = sorted(range(4), key=lambda i: len(S[i]))
    L0, L1, R0, R1 = order[0], order[1], order[2], order[3]
    lk = (K[L0][:, None] + K[L1][None, :]).ravel()
    lo = np.argsort(lk, kind='stable'); lks = lk[lo]; del lk
    out = []
    n1 = len(S[L1])
    for c0 in range(0, len(S[R1]), chunk):
        kk = K[R1][c0:c0 + chunk]
        rk = (np.int64(tk) - K[R0][:, None] - kk[None, :]).ravel()
        pos = np.searchsorted(lks, rk); pos[pos >= len(lks)] = len(lks) - 1
        for r in np.nonzero(lks[pos] == rk)[0]:
            p_ = pos[r]
            while p_ < len(lks) and lks[p_] == rk[r]:
                a, b = divmod(int(lo[p_]), n1); c, d = divmod(int(r), len(kk)); d += c0
                idx = {L0: a, L1: b, R0: c, R1: d}
                if (S[0][idx[0]] + S[1][idx[1]] + S[2][idx[2]] + S[3][idx[3]] == T).all():
                    out.append((idx[0], idx[1], idx[2], idx[3]))
                p_ += 1
    return out


def decide(z, tgt=None):
    tgt = TGT if tgt is None else np.asarray(tgt, np.int64)
    etot = int(tgt[0]) - 80
    mins = [min_energy(*bl) for bl in z]
    if any(m is None for m in mins) or sum(mins) > etot:
        return dict(status='EMPTY (Q-compression or energy)', mins=mins, sizes=None, solutions=[])
    Ls = [block_list(a, b, c, etot - (sum(mins) - mins[i])) for i, (a, b, c) in enumerate(z)]
    sizes = [len(L) for L in Ls]
    if min(sizes) == 0:
        return dict(status='EMPTY (a block list is empty)', mins=mins, sizes=sizes, solutions=[])
    S = [sigs(L) for L in Ls]
    # energy filter inside the join is implied by the lag-0 coordinate
    hits = join22(S, tgt)
    sols = []
    for (a, b, c, d) in hits:
        ys = [Ls[0][a].tolist(), Ls[1][b].tolist(), Ls[2][c].tolist(), Ls[3][d].tolist()]
        # direct re-check: all 20 lags, tau-invariance, alphabet, Q-compressions
        ok = all(sum(sum(y[k] * y[(((k // 5) ^ (h // 5)) * 5) + (k % 5 + h % 5) % 5] for k in range(20)) for y in ys)
                 == int(tgt[h]) for h in range(20))
        ok &= all(y[k] == y[TAU[k]] and y[k] in ASET for y in ys for k in range(20))
        ok &= all(sum(y[q * 5:(q + 1) * 5]) == v for y, (a_, b_, c_) in zip(ys, z) for q, v in ((0, a_), (1, b_), (2, b_), (3, c_)))
        sols.append(dict(y=ys, recheck=ok))
    return dict(status='SOLUTIONS' if sols else 'EMPTY (exhaustive join)', mins=mins, sizes=sizes, solutions=sols)


if __name__ == '__main__':
    Z = q_solutions()
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    last = int(sys.argv[2]) if len(sys.argv) > 2 else len(Z)
    outp = sys.argv[3] if len(sys.argv) > 3 else None
    for n in range(first, min(last, len(Z))):
        t = time.time()
        r = decide(Z[n]); r.update(index=n, z=Z[n], rows=[a + 2 * b + c for a, b, c in Z[n]],
                                   nsol=len(r['solutions']), secs=round(time.time() - t, 1))
        s = json.dumps(r); print(s[:400], flush=True)
        if outp:
            with open(outp, 'a') as fo: fo.write(s + '\n')
