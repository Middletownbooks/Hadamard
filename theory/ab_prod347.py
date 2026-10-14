r"""Positive control: the staged pipeline on the published 347 slice (b = 5, c = 6), end to end [M].

Slice (theory/ab_verify_published.py): G = Z_2 x Z_3 x Z_3 x Z_19 (= Z_2 x Z_3 x Z_57), S = <-1 on the two Z_3 (diagonal),
order 3 on Z_19> (index 6, f = 3).  Stage 1 = a y on G' = Z_2 x Z_3 x Z_3 (the published template's own y by default,
so its Q-compression, hence the template's border, is fixed); stage 2 = the hierarchical solver of theory/ab_stage2.py
(S/D on Z_2, middle level along the H-part Z_3, exact D-check).  Every D-solution is assembled on G, verified
(sum PAF = T exactly), mapped to the template coordinates (Z_2 x Z_3 x Z_57) and completed with the template's border
and corner (theory/ab_assemble.py) to a 1388 x 1388 matrix, checked H H^T = 1388 I.
Restart-safe: one JSONL line per middle value (logs/ab/prod347/run_<tag>.jsonl), done middle values are skipped.
usage: python3 theory/ab_prod347.py --tag T [--R 800 800] [--K 4 256] [--seed 1] [--middles 400] [--stop-on-core]"""
import sys, os, json, time, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_stage2 as S2m, ab_sweep as W, ab_assemble as AS
from ab_s1sd import cfg_from_json
from ab_verify_published import recoordinate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')
AMAP = [(0, 2), (1, 3), (2, 3), (2, 19)]


def to_template_coords(X):
    """inverse of recoordinate (2,3,57) -> (2,3,3,19)"""
    from ab_lib import AG
    Go, Gn = AG((2, 3, 57)), AG((2, 3, 3, 19))
    co_n = np.zeros((Go.n, 4), np.int64)
    for j, (oa, mod) in enumerate(AMAP): co_n[:, j] = Go.co[:, oa] % mod
    newidx = Gn.idx(co_n)
    return X[:, newidx]


def main():
    a = sys.argv[1:]
    def opt(name, n=1, default=None):
        if name in a:
            i = a.index(name); return [int(x) for x in a[i + 1:i + 1 + n]] if n > 1 else a[i + 1]
        return default
    tag = opt('--tag', default='t0'); R = opt('--R', 2, [800, 800]); K = opt('--K', 2, [4, 256])
    seed = int(opt('--seed', default='1')); nmid = int(opt('--middles', default='400'))
    outdir = os.path.join(ROOT, 'logs/ab/prod347'); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, f'run_{tag}.jsonl')
    done = set()
    if os.path.exists(logf):
        for l in open(logf):
            try: d = json.loads(l)
            except Exception: continue
            if 'mid' in d: done.add(d['mid'])
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    cfg = cfg_from_json('@' + CFG); sl = W.build(cfg)
    T = AS.load_template347()
    XF = recoordinate(T['F'], (2, 3, 57), (2, 3, 3, 19), AMAP)
    assert sl.is_invariant(XF) and (sl.G.paf(XF).sum(0) == sl.T).all()
    Y = sl.compress(XF)
    t0 = time.time()
    Sv = S2m.FastSolver2(cfg, Y, Laxis=2, log=lambda s: print(s, flush=True))
    Sv.prepare_rows(R=R, K=K, seed=seed)
    cands, ncommon = Sv.choose_middles(top=nmid, seed=seed + 100)
    log(dict(stage="setup", tag=tag, R=R, K=K, seed=seed, middles=len(cands), common=ncommon, secs=round(time.time() - t0, 1)))
    for ci, (c, sc) in enumerate(cands):
        mid = hashlib.sha1(c.tobytes()).hexdigest()[:16]
        if mid in done: continue
        t1 = time.time()
        try: sols, info = Sv.s_join(c)
        except RuntimeError as e:
            log(dict(mid=mid, ci=ci, skipped=str(e))); continue
        t2 = time.time()
        nd = 0; dcand = 0; cores = 0
        for cf in sols:
            try: ds, sizes = Sv.d_check_fast(cf)
            except RuntimeError: continue
            dcand += float(np.prod([float(x) for x in sizes]))
            for d in ds:
                Xa, ok = Sv.assemble(d)
                nd += 1
                if not ok: continue
                F2 = to_template_coords(Xa)
                Z = AS.compressions(F2, (2, 3, 57), 2)
                Cf = AS.forced_corner(Z, (2, 3), T['layout'], T['K'], T['Kc'])
                H = AS.assemble(F2, (2, 3, 57), 2, T['layout'], T['K'], T['Kc'], Cf) if Cf is not None else None
                had = H is not None and AS.is_hadamard(H)
                same = bool((F2 == T['F']).all())
                cores += 1
                tagc = f'{tag}_{mid}_{cores}'
                np.save(os.path.join(outdir, f'core_{tagc}.npy'), F2)
                if had: np.savez_compressed(os.path.join(outdir, f'H1388_{tagc}.npz'), H=np.packbits(H > 0, axis=1))
                log(dict(CORE=tagc, verified_core=True, hadamard_1388=bool(had), identical_to_published=same))
        log(dict(mid=mid, ci=ci, **info, S_solutions=len(sols), D_candidates_log2=round(float(np.log2(dcand)), 2) if dcand else None,
                 D_solutions=nd, cores=cores, join_secs=round(t2 - t1, 1), dcheck_secs=round(time.time() - t2, 1)))
        if cores and '--stop-on-core' in a: break


if __name__ == "__main__":
    main()
