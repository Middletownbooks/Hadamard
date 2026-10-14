r"""b = 0 at composite v = F q: exact stage 1 on Z_F (direct enumeration) and the stage-2 Gaussian per y [M, E].
Stage 1: four U-invariant integer sequences y_k on Z_F, letters alpha + f sigma, sum_k PAF(y_k) = 4v delta.
Enumeration: per block all letter assignments with energy <= 4v - 3 (min energy), 2+2 join on the exact PAF vector;
solutions reported up to block order and sign.  Stage 2: theory/ab_gauss.stage2 on each solution (exact H2, Gaussian C2).
usage: python3 theory/ab_b0_stage1.py CFGJSON [--max 20] [--N 20000]"""
import sys, os, json, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json

def stage1(sl, cap=5 * 10 ** 6):
    Gp = sl.Gp; L = sl.letters; w = sl.bsz.astype(np.int64); nb = len(L)
    T1 = sl.T1; T0 = int(T1[0])
    minE = int(sum(int(w[j]) * min(y * y for y in L[j]) for j in range(nb)))
    Emax = T0 - 3 * minE
    combos = [list(t) for t in itertools.product(*L) if sum(int(w[j]) * t[j] * t[j] for j in range(nb)) <= Emax]
    if len(combos) > cap: raise RuntimeError("too many blocks")
    Y = np.array(combos, np.int64)
    Y = Y[Y[:, 0] > 0] if (Y[:, 0] != 0).all() else Y
    full = Y[:, sl.borb]
    P = Gp.paf(full)
    key = defaultdict(list)
    for i, row in enumerate(P): key[row.tobytes()].append(i)
    uk = [np.frombuffer(k_, dtype=np.int64) for k_ in key]
    pairs = defaultdict(list)
    for i in range(len(uk)):
        for j in range(i, len(uk)):
            pairs[(uk[i] + uk[j]).tobytes()].append((i, j))
    sols = set()
    for k_, lst in pairs.items():
        need = (T1 - np.frombuffer(k_, dtype=np.int64)).tobytes()
        if need not in pairs: continue
        for (a, b) in lst:
            for (c, d) in pairs[need]:
                for ia in key[uk[a].tobytes()]:
                    for ib in key[uk[b].tobytes()]:
                        for ic in key[uk[c].tobytes()]:
                            for id_ in key[uk[d].tobytes()]:
                                q = tuple(sorted((ia, ib, ic, id_)))
                                sols.add(q)
    out = []
    for q in sorted(sols):
        Yq = full[list(q)]
        assert (Gp.paf(Yq).sum(0) == T1).all()
        out.append(Yq)
    return out, len(Y)

if __name__ == "__main__":
    a = sys.argv[1:]
    cfg = cfg_from_json(a[0]); sl = W.build(cfg)
    mx = int(a[a.index('--max') + 1]) if '--max' in a else 20
    N = int(a[a.index('--N') + 1]) if '--N' in a else 20000
    t0 = time.time()
    sols, nb = stage1(sl)
    print(json.dumps(dict(v=cfg['v'], F=sl.Gp.mods, q=cfg['q'], s=cfg['s'], U=[list(u) for u in cfg['U']], orbits=int(sl.r),
                          blocks=nb, stage1_solutions=len(sols), secs=round(time.time() - t0, 1))), flush=True)
    rng = np.random.default_rng(0)
    for i in rng.permutation(len(sols))[:mx]:
        Y = sols[i]
        r = AGs.stage2(sl, Y, N=N)
        print(json.dumps(dict(y=Y.tolist(), rows=Y.sum(1).tolist(), H2=r.get('H2'), C2=float(r['bits']) if r.get('bits') is not None else None,
                              mahal=r.get('mahal'), stage2_log2E=float(r['stage2_log2E']) if r.get('stage2_log2E') is not None else None)), flush=True)
