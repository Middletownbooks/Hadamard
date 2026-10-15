r"""Exact stage-1 counts of the c = 6 family (theory/b599_fam.py) for a fixed Q-compression quadruple z [M].

Per block: every block vector with Q-compression z_k (DP over cosets, energy <= T1(0) - the other blocks' minima) and
its class signature (PAF at the 4J lag classes, class 0 = energy).  Ordered 4-tuples with sum of signatures = T1 are
counted by a 2 + 2 join: left pairs (blocks 0, 1) keyed by a 64-bit linear hash of the signature sum, sorted; right
pairs (blocks 2, 3) streamed against T1 - (sum); every hash match is re-checked exactly.  Solutions can be returned.
This is an independent enumerator (no code shared with theory/ab_s1sd.py); the two are cross-checked in main.
usage: python3 theory/b599_exact.py --n 5 --q 11 [--out FILE] [--crosscheck]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam


def block_list(F, z, Emax):
    """all block vectors with Q-compression z and energy <= Emax: (V (N x no), energies)"""
    J = F.J
    TUP = np.array(list(itertools.product(F.LET.tolist(), repeat=J)), np.int64)
    TSUM = TUP @ F.WJ; TEN = (TUP * TUP) @ F.WJ
    per = []
    for q in range(6):
        ix = np.nonzero(TSUM == z[q])[0]
        per.append((TUP[ix], TEN[ix]))
    # combine coset by coset with energy pruning (minimum of remaining cosets)
    mins = [int(e.min()) if len(e) else 10 ** 9 for _, e in per]
    if max(mins) >= 10 ** 9: return np.zeros((0, 6 * J), np.int64), np.zeros(0, np.int64)
    rem = [sum(mins[q + 1:]) for q in range(6)]
    V = np.zeros((1, 0), np.int64); E = np.zeros(1, np.int64)
    for q in range(6):
        T, e = per[q]
        newE = (E[:, None] + e[None, :])
        ok = newE + rem[q] <= Emax
        i, j = np.nonzero(ok)
        V = np.concatenate([V[i], T[j]], 1); E = newE[i, j]
    return V, E


@njit(cache=True)
def _count_join(SA, SB, SC, SD, T, hA, hB, hC, hD, hT, ret, maxret):
    """count (a, b, c, d) with SA[a] + SB[b] + SC[c] + SD[d] == T (rows of int64 signatures); h* = 64-bit linear
    hashes of the rows (wrapping arithmetic), hT the hash of T"""
    nA, nB = SA.shape[0], SB.shape[0]; d = SA.shape[1]
    L = np.empty(nA * nB, np.int64); LI = np.empty(nA * nB, np.int64)
    m = 0
    for a in range(nA):
        for b in range(nB):
            L[m] = hA[a] + hB[b]; LI[m] = a * nB + b; m += 1
    o = np.argsort(L); Ls = L[o]
    cnt = 0; nret = 0
    for c in range(SC.shape[0]):
        for dd in range(SD.shape[0]):
            key = hT - hC[c] - hD[dd]
            lo = np.searchsorted(Ls, key)
            while lo < m and Ls[lo] == key:
                p = LI[o[lo]]; a = p // nB; b = p % nB
                ok = True
                for i in range(d):
                    if SA[a, i] + SB[b, i] + SC[c, i] + SD[dd, i] != T[i]:
                        ok = False; break
                if ok:
                    if nret < maxret:
                        ret[nret, 0] = a; ret[nret, 1] = b; ret[nret, 2] = c; ret[nret, 3] = dd
                    nret += 1; cnt += 1
                lo += 1
    return cnt, nret


def count_z(F, z, maxret=1000, seed=7):
    T0 = int(F.T1SIG[0])
    Ls = []
    mins = []
    for k in range(4):
        V, E = block_list(F, z[k], T0)
        if len(V) == 0: return dict(count=0, sizes=[0] * 4), []
        Ls.append(V); mins.append(int(E.min()))
    # prune each block by the other blocks' minimum energies
    for k in range(4):
        cap = T0 - (sum(mins) - mins[k])
        Ek = F.energy(Ls[k]); Ls[k] = Ls[k][Ek <= cap]
    sizes = [len(x) for x in Ls]
    if min(sizes) == 0: return dict(count=0, sizes=sizes), []
    S = [F.sig(x) for x in Ls]
    rng = np.random.default_rng(seed)
    w = rng.integers(1, 2 ** 62, size=F.nl, dtype=np.int64)
    ret = np.zeros((maxret, 4), np.int64)
    # put the two smallest lists on the left
    order = np.argsort(sizes)
    oL = [int(order[0]), int(order[1])]; oR = [int(order[2]), int(order[3])]
    with np.errstate(over='ignore'):
        H = [x @ w for x in S]; hT = np.int64(F.T1SIG.astype(np.int64) @ w)
    cnt, nret = _count_join(S[oL[0]], S[oL[1]], S[oR[0]], S[oR[1]], F.T1SIG.astype(np.int64),
                            H[oL[0]], H[oL[1]], H[oR[0]], H[oR[1]], hT, ret, maxret)
    sols = []
    for r in ret[:min(nret, maxret)]:
        Y = np.zeros((4, F.no), np.int64)
        Y[oL[0]] = Ls[oL[0]][r[0]]; Y[oL[1]] = Ls[oL[1]][r[1]]; Y[oR[0]] = Ls[oR[0]][r[2]]; Y[oR[1]] = Ls[oR[1]][r[3]]
        sols.append(Y)
    return dict(count=int(cnt), sizes=sizes), sols


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '5')), int(opt('--q', '11'))
    outp = opt('--out', '')
    F = Fam(n, q)
    import b599_flevels as FL
    quads = FL.class_quads(F)
    done = set()
    if outp and os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l); done.add((d['class'], tuple(d['c'])))
            except Exception: pass
    fo = open(outp, 'a') if outp else None
    t0 = time.time(); tot = 0
    for ci, cl, c, z in quads:
        if (ci, tuple(c)) in done: continue
        t1 = time.time()
        r, sols = count_z(F, z)
        ok = all(F.check_stage1(Y) and (F.qcomp(Y) == z).all() for Y in sols)
        rec = dict(n=n, q=q, v=F.v, **{'class': ci}, c=list(c), rows=[int(x) for x in z.sum(1)], count=r['count'],
                   sizes=r['sizes'], verified=bool(ok), secs=round(time.time() - t1, 2))
        if sols: rec['example'] = sols[0].tolist()
        tot += r['count']
        s = json.dumps(rec); print(s, flush=True)
        if fo: fo.write(s + "\n"); fo.flush()
    print(json.dumps(dict(total=tot, secs=round(time.time() - t0, 1))))
