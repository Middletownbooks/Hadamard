r"""Exact decision (and count) of small slices, b = 0 or b >= 3, for validating the screens of theory/ab_sweep.py [M].

All 2^r S-invariant blocks (r = orbits per block <= ~24) are enumerated; per-block PSD filter
|X^(chi)|^2 <= That(chi) (That = 4b on the nontrivial characters trivial on H, else 4v); signatures = PAF at the
<S, -1>-lag-class representatives; a 2+2 join (pair sums keyed exactly) counts ordered quadruples with
sum = T at every lag class.  Normalisation: none (ordered, signed quadruples are counted).
usage: python3 theory/ab_exact_small.py --sweep VMAX [--rmax 22] [--out FILE]
       (every GP config of ab_sweep for b >= 3 with c = 4, 6 and v <= VMAX, and b = 0 composite v <= VMAX)"""
import sys, os, json, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_lib import factorint


def that_mask(sl, cfg):
    G = sl.G
    chis = G.co
    if cfg['b'] >= 3:
        nQ = cfg['nQ']
        # characters trivial on H = those with zero coordinates on the H axes (split G = Q x H, dual coordinates)
        onH = (chis[:, nQ:] == 0).all(1)
        nontriv = (chis != 0).any(1)
        return np.where(onH & nontriv, 4 * cfg['b'], 4 * cfg['v']).astype(float)
    return np.full(G.n, 4.0 * cfg['v'])


def decide(cfg, rmax=22, chunk=1 << 15, pairmax=4 * 10 ** 7):
    sl = W.build(cfg)
    r = sl.r
    if r > rmax: return dict(status="r too large", r=r)
    G = sl.G
    That = that_mask(sl, cfg)
    Up = [G.act(u) for u in sl.S]
    negp = G.idx(-G.co)
    seen = np.zeros(G.n, bool); reps = []
    for g in range(G.n):
        if seen[g]: continue
        reps.append(g)
        for P in Up: seen[P[g]] = True; seen[negp[P[g]]] = True
    T = sl.T[reps]
    keepP = []; keepL = []
    for s0 in range(0, 1 << r, chunk):
        idx = np.arange(s0, min(1 << r, s0 + chunk), dtype=np.int64)
        lab = (1 - 2 * ((idx[:, None] >> np.arange(r)) & 1)).astype(np.int8)
        X = lab[:, sl.orb]
        Pw = G.psd(X)
        ok = (Pw <= That[None] + 1e-6).all(1)
        if not ok.any(): continue
        Pf = G.paf(X[ok])
        keepP.append(Pf[:, reps]); keepL.append(lab[ok])
    if not keepP: return dict(status="no block passes PSD", r=r, solutions=0)
    P = np.concatenate(keepP); L = np.concatenate(keepL); n = len(P)
    # unique signatures with multiplicities, grouped by e = a^2 (determined by the signature)
    up, first, cnt = np.unique(P, axis=0, return_index=True, return_counts=True)
    sizes = np.zeros(len(reps), np.int64)
    for P_ in Up:
        pass
    # class sizes of the lag classes: a^2 = sum_g PAF(g) = sum_classes |class| PAF(rep)
    cls = -np.ones(G.n, np.int64)
    for ci, g in enumerate(reps):
        for P_ in Up:
            cls[P_[g]] = ci; cls[negp[P_[g]]] = ci
    csz = np.bincount(cls, minlength=len(reps))
    e = up @ csz
    groups = defaultdict(list)
    for i, ee in enumerate(e.tolist()): groups[ee].append(i)
    es = sorted(groups)
    tot = 0; witness = None; work = 0
    pairdicts = {}
    def pd(e1, e2):
        key = (e1, e2)
        if key not in pairdicts:
            A = up[groups[e1]]; B = up[groups[e2]]; ca = cnt[groups[e1]]; cb = cnt[groups[e2]]
            d = defaultdict(int); wit = {}
            for i in range(len(A)):
                S = A[i] + B
                for j in range(len(B)):
                    k = S[j].tobytes(); d[k] += int(ca[i] * cb[j])
                    if k not in wit: wit[k] = (groups[e1][i], groups[e2][j])
            pairdicts[key] = (d, wit)
        return pairdicts[key]
    for e1 in es:
        for e2 in es:
            for e3 in es:
                e4 = 4 * cfg['v'] - e1 - e2 - e3
                if e4 not in groups: continue
                work += len(groups[e1]) * len(groups[e2]) + len(groups[e3]) * len(groups[e4])
                if work > pairmax: return dict(status="join too large", r=r, blocks=n, usig=len(up))
                d12, w12 = pd(e1, e2); d34, w34 = pd(e3, e4)
                for k, c12 in d12.items():
                    need = (T - np.frombuffer(k, dtype=np.int64)).tobytes()
                    if need in d34:
                        tot += c12 * d34[need]
                        if witness is None:
                            i1, i2 = w12[k]; i3, i4 = w34[need]
                            witness = [int(first[i1]), int(first[i2]), int(first[i3]), int(first[i4])]
    rec = dict(status="decided", r=r, blocks=n, solutions=int(tot))
    if witness is not None:
        X = L[witness][:, sl.orb].astype(np.int64)
        assert (G.paf(X).sum(0) == sl.T).all()
        rec['witness_labels'] = L[witness].tolist()
    return rec


if __name__ == "__main__":
    a = sys.argv[1:]
    VMAX = int(a[a.index('--sweep') + 1]); rmax = int(a[a.index('--rmax') + 1]) if '--rmax' in a else 22
    out = a[a.index('--out') + 1] if '--out' in a else None
    fo = open(out, 'a') if out else None
    for v in range(7, VMAX + 1, 2):
        cfgs = []
        for c in (4, 6):
            if (v + 1) % c == 0 and (v + 1) // c - 1 >= 3: cfgs += W.configs_bge3(v, (c,), maxU=1)
        fv = factorint(v)
        if len(fv) >= 2: cfgs += W.configs_b0(v, maxU=1)
        for cfg in cfgs:
            sl = W.build(cfg)
            if sl.r > rmax: continue
            t0 = time.time()
            d = decide(cfg, rmax)
            scr = W.run(cfg)
            rec = dict(v=v, c=cfg['c'], b=cfg['b'], Q=cfg['Q'], H1=cfg['H1'], q=cfg['q'], s=cfg['s'], U=scr['U'],
                       r=sl.r, solutions=d.get('solutions'), status=d['status'], A=scr.get('A'), B=scr.get('B'),
                       secs=round(time.time() - t0, 1))
            if 'witness_labels' in d: rec['witness_labels'] = d['witness_labels']
            if d.get('solutions') and (scr.get('A') is False or scr.get('B') is False): rec['FALSE_NEGATIVE'] = True
            s = json.dumps(rec); print(s, flush=True)
            if fo: fo.write(s + "\n"); fo.flush()
