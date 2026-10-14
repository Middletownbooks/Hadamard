r"""Exact stage 1 of a slice of theory/ab_lib.Slice (any G', any target T1), by enumeration + energy-split 2+2 join [M].

Stage 1: four U-invariant integer sequences y_k on G' (y_k(b) in the letter set of b's fibre) with
    sum_k PAF_{G'}(y_k) = T1   (b >= 3: 4m - 4b(q-1) at 0, -4bq on H' \ 0, +4q off H';  b = 0: 4v delta).
Method: per block, all letter assignments to the G'-orbits with energy <= T1(0) - 3 * (minimum block energy) (DFS);
signature = PAF at the <U, -1>-lag-class representatives (lag 0 = energy).  Blocks are normalised: y and -y have the
same PAF, so each block is stored once up to sign (first nonzero orbit value positive).  Join: for every ordered
energy split (E1 <= E2 | E3 <= E4 as unordered pairs) the pair sums of the left lists are hashed and the right pair
sums looked up; every match is re-verified exactly.  No solution with complete = true proves the stage 1 (hence the
slice) empty.
usage: python3 theory/ab_stage1.py CFGJSON [--max-sol N] [--out FILE] [--max-pairs N]
       (CFGJSON: a config dict as produced by theory/ab_sweep.py, e.g. '{"v":659,...}' or @file)"""
import sys, os, json, time, math, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_gauss import lag_reps


def block_lists(sl, Emax, canon=True, cap=3 * 10 ** 7):
    """all U-invariant y on G' with energy <= Emax; returns (Yorb (n, L1) letters per orbit, E, rows)"""
    L = sl.letters; w = sl.bsz.astype(np.int64); nb = len(L)
    minrest = np.zeros(nb + 1, np.int64)
    for j in range(nb - 1, -1, -1): minrest[j] = minrest[j + 1] + w[j] * min(y * y for y in L[j])
    out = []
    cur = [0] * nb
    def rec(j, E):
        if j == nb:
            out.append(list(cur))
            if len(out) > cap: raise RuntimeError("block list too large")
            return
        for y in L[j]:
            E2 = E + int(w[j]) * y * y
            if E2 + minrest[j + 1] > Emax: continue
            cur[j] = y; rec(j + 1, E2)
    rec(0, 0)
    Y = np.array(out, np.int64).reshape(-1, nb)
    if canon:
        # keep one of y, -y: first nonzero orbit value positive (letters are odd so y[0] != 0 when 0 is not a letter)
        nz = np.argmax(Y != 0, axis=1)
        sgn = np.sign(Y[np.arange(len(Y)), nz])
        Y = Y[sgn > 0]
    E = (Y * Y) @ w
    R = Y @ w
    return Y, E, R


def solve(cfg, max_sol=1000, max_pairs=4 * 10 ** 8, log=print, T1=None):
    sl = W.build(cfg)
    Gp = sl.Gp
    Up = [Gp.act(u) for u in sl.Sp]
    reps, cls = lag_reps(Gp, Up)
    T1 = sl.T1 if T1 is None else T1
    T = T1[reps]
    T0 = int(T1[0])
    minE = int(sum(int(sl.bsz[j]) * min(y * y for y in L) for j, L in enumerate(sl.letters)))
    Emax = T0 - 3 * minE
    t0 = time.time()
    Y, E, R = block_lists(sl, Emax)
    full = Y[:, sl.borb]
    P = Gp.paf(full)[:, reps]
    log(json.dumps(dict(stage="lists", blocks=len(Y), Emax=Emax, lag_classes=len(reps), secs=round(time.time() - t0, 1))))
    # group by energy
    order = np.argsort(E, kind='stable'); Y, E, R, P = Y[order], E[order], R[order], P[order]
    Es, starts = np.unique(E, return_index=True)
    ends = np.r_[starts[1:], len(E)]
    rngE = {int(e): (int(a), int(b)) for e, a, b in zip(Es, starts, ends)}
    # hashing weights
    rng = np.random.default_rng(12345)
    wv = rng.integers(1, 2 ** 61, size=P.shape[1], dtype=np.int64)
    H = (P * wv).sum(1)                                   # wraps mod 2^64 (int64 overflow is fine for hashing)
    Th = int((T * wv).sum())
    sols = []; work = 0
    Elist = [int(e) for e in Es]
    # unordered pairs of energies for the left side: (Ea <= Eb), right side (Ec <= Ed) with Ea+Eb+Ec+Ed = T0 and
    # (Ea, Eb) <= (Ec, Ed) lexicographically to avoid double counting of the split
    Eset = set(Elist)
    for Ea in Elist:
        for Eb in Elist:
            if Eb < Ea: continue
            EL = Ea + Eb
            if EL > T0: break
            for Ec in Elist:
                Ed = T0 - EL - Ec
                if Ed < Ec: break
                if Ed not in Eset: continue
                if (Ec, Ed) < (Ea, Eb): continue
                a0, a1 = rngE[Ea]; b0, b1 = rngE[Eb]; c0, c1 = rngE[Ec]; d0, d1 = rngE[Ed]
                nl = (a1 - a0) * (b1 - b0); nr = (c1 - c0) * (d1 - d0)
                work += nl + nr
                if work > max_pairs: return dict(complete=False, reason="pair cap", work=work, solutions=len(sols)), sols
                Lh = (H[a0:a1, None] + H[None, b0:b1]).ravel()
                o = np.argsort(Lh); Ls = Lh[o]
                Rh = (Th - (H[c0:c1, None] + H[None, d0:d1])).ravel()
                lo = np.searchsorted(Ls, Rh, 'left'); hi = np.searchsorted(Ls, Rh, 'right')
                for ri in np.nonzero(hi > lo)[0]:
                    ci, di = divmod(int(ri), d1 - d0); ci += c0; di += d0
                    for p in o[lo[ri]:hi[ri]]:
                        ai, bi = divmod(int(p), b1 - b0); ai += a0; bi += b0
                        if Ea == Eb and bi < ai: continue
                        if Ec == Ed and di < ci: continue
                        if ((P[ai] + P[bi] + P[ci] + P[di]) != T).any(): continue
                        quad = np.array([full[ai], full[bi], full[ci], full[di]])
                        assert (Gp.paf(quad).sum(0) == T1).all()
                        sols.append(dict(y=quad.tolist(), rows=[int(R[i]) for i in (ai, bi, ci, di)],
                                         E=[int(E[i]) for i in (ai, bi, ci, di)]))
                        if len(sols) >= max_sol:
                            return dict(complete=False, reason="max_sol", work=work, solutions=len(sols)), sols
    return dict(complete=True, work=work, solutions=len(sols), secs=round(time.time() - t0, 1)), sols


def planted(cfg, seed=0, log=print):
    """random quadruple from the block lists (energies summing to anything); target = its own PAF; must be found"""
    sl = W.build(cfg); Gp = sl.Gp
    rng = np.random.default_rng(seed)
    T0 = int(sl.T1[0]); minE = int(sum(int(sl.bsz[j]) * min(y * y for y in L) for j, L in enumerate(sl.letters)))
    Y, E, R = block_lists(sl, T0 - 3 * minE)
    for _ in range(100000):
        idx = rng.integers(len(Y), size=4)
        if E[idx].sum() <= T0 - 0: break
    quad = Y[idx][:, sl.borb]
    T1 = Gp.paf(quad).sum(0)
    # the solver uses T1(0) as the energy total, so plant's total energy is the target energy
    res, sols = solve(cfg, max_sol=10 ** 6, log=lambda s: None, T1=T1)
    canon = lambda Qd: sorted(tuple(int(x) * (1 if r[0] > 0 else -1) for x in r) for r in np.asarray(Qd).tolist())
    keyp = canon(quad)
    found = any(canon(s_['y']) == keyp for s_ in sols)
    out = dict(seed=seed, found=bool(found), solutions=len(sols), complete=res.get('complete'))
    log(json.dumps(out)); return out


if __name__ == "__main__":
    a = sys.argv[1:]
    cfgs = a[0]
    cfg = json.load(open(cfgs[1:])) if cfgs.startswith('@') else json.loads(cfgs)
    cfg['mods'] = tuple(cfg['mods']); cfg['Q'] = tuple(cfg['Q']); cfg['H1'] = tuple(cfg['H1'])
    cfg['U'] = [tuple(u) if u[0] != 'neg' else ('neg', tuple(u[1])) for u in cfg['U']]
    max_sol = int(a[a.index('--max-sol') + 1]) if '--max-sol' in a else 1000
    max_pairs = int(float(a[a.index('--max-pairs') + 1])) if '--max-pairs' in a else 4 * 10 ** 8
    if '--plant' in a:
        n = int(a[a.index('--plant') + 1])
        r = [planted(cfg, sd) for sd in range(n)]
        print(json.dumps(dict(planted=n, found=sum(x['found'] for x in r)))); sys.exit()
    res, sols = solve(cfg, max_sol, max_pairs)
    print(json.dumps(res))
    if '--out' in a:
        json.dump(dict(cfg=dict(cfg, U=[list(u) for u in cfg['U']]), result=res, solutions=sols[:2000]),
                  open(a[a.index('--out') + 1], 'w'))
