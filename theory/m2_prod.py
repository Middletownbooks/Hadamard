r"""Restart-safe b >= 3 stage-2 production with the compiled D-check, on the 347 control slice [M].

Pipeline: theory/ab_stage2.FastSolver2 (S/D split on Z2, Wagner middle level along the H-part Z3, exact S-join) ->
theory/m2_dcheck.DCheck (compiled exact D-check, validated against ab_stage2._dcheck_fast on the published and the a1
core, 20 planted lifts and 600 real S-solutions) -> assembly on G, exact core check, the 347 template border and its
forced corner (theory/ab_assemble.py), H H^T = 1388 I.
Worker w of a run builds its own row lists (seed + 1000 w) and works through its own middle values.  Per middle value
the S-solutions are saved (logs/m2/prod347/<tag>_w<w>_<mid>.npy, S-value indices) before any D-check; D-check
progress is logged every --ckpt S-solutions, so a restarted worker resumes inside a middle value.  Logs:
logs/m2/prod347/<tag>_w<w>.jsonl; cores: logs/m2/prod347/core_<tag>_w<w>_<mid>_<n>.npy and H1388_*.npz.
Calibration fields per middle value: sides (left, right), S-solutions found, the prediction left*right*2^-F with
F = C_S - C_M at the true target (theory/m2_levelcal.py), D-checks, D-solutions.
usage: python3 theory/m2_prod.py --tag T --worker W [--seed 1] [--R 800 800] [--K 4 128] [--middles 400]
                                 [--nsample 4000000] [--ckpt 2000] [--stop-cores N]"""
import sys, os, json, time, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_stage2 as S2m, ab_sweep as W, ab_assemble as AS
from ab_s1sd import cfg_from_json
from ab_verify_published import recoordinate
from m2_dcheck import DCheck

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')
AMAP = [(0, 2), (1, 3), (2, 3), (2, 19)]
F_FINE = 69.24 - 36.41          # C_S - C_M (middle level along the H-part Z3) at the true target, logs/m2/levelcal_347.json


def to_template_coords(X):
    from ab_lib import AG
    Go, Gn = AG((2, 3, 57)), AG((2, 3, 3, 19))
    co_n = np.zeros((Go.n, 4), np.int64)
    for j, (oa, mod) in enumerate(AMAP): co_n[:, j] = Go.co[:, oa] % mod
    return X[:, Gn.idx(co_n)]


def main():
    a = sys.argv[1:]
    def opt(name, n=1, default=None):
        if name in a:
            i = a.index(name); return [int(x) for x in a[i + 1:i + 1 + n]] if n > 1 else a[i + 1]
        return default
    tag = opt('--tag', default='m2a'); w = int(opt('--worker', default='0'))
    seed = int(opt('--seed', default='1')) + 1000 * w
    R = opt('--R', 2, [800, 800]); K = opt('--K', 2, [4, 128])
    nmid = int(opt('--middles', default='400')); nsample = int(opt('--nsample', default='4000000'))
    ckpt = int(opt('--ckpt', default='2000')); stop_cores = int(opt('--stop-cores', default='0'))
    outdir = os.path.join(ROOT, 'logs/m2/prod347'); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, f'{tag}_w{w}.jsonl')
    done, progress, cores_so_far = set(), {}, 0
    if os.path.exists(logf):
        for l in open(logf):
            try: d = json.loads(l)
            except Exception: continue
            if d.get('stage') == 'middle_done': done.add(d['mid'])
            if d.get('stage') == 'dprogress': progress[d['mid']] = d
            if 'CORE' in d: cores_so_far += 1
    def log(rec):
        rec['t'] = round(time.time(), 1)
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    cfg = cfg_from_json('@' + CFG); sl = W.build(cfg)
    T = AS.load_template347()
    XF = recoordinate(T['F'], (2, 3, 57), (2, 3, 3, 19), AMAP)
    assert sl.is_invariant(XF) and (sl.G.paf(XF).sum(0) == sl.T).all()
    Y = sl.compress(XF)
    t0 = time.time()
    Sv = S2m.FastSolver2(cfg, Y, Laxis=2, log=lambda s: None)
    Sv.prepare_rows(R=R, K=K, seed=seed)
    cands, ncommon = Sv.choose_middles(nsample=nsample, top=nmid, seed=seed + 100)
    DC = DCheck(Sv)
    log(dict(stage="setup", tag=tag, worker=w, seed=seed, R=R, K=K, middles=len(cands), common=ncommon,
             secs=round(time.time() - t0, 1), resumed_done=len(done)))
    for ci, (c, sc) in enumerate(cands):
        mid = hashlib.sha1(c.tobytes()).hexdigest()[:16]
        if mid in done: continue
        fn = os.path.join(outdir, f'{tag}_w{w}_{mid}.npy')
        t1 = time.time()
        if os.path.exists(fn):
            S = np.load(fn); info = progress.get(mid, {}).get('info', {}); jsecs = 0.0
        else:
            try: sols, info = Sv.s_join(c)
            except RuntimeError as e:
                log(dict(stage="middle_done", mid=mid, ci=ci, skipped=str(e))); continue
            S = np.array([np.stack(cf) for cf in sols], np.int16).reshape(-1, 4, len(Sv.S2.pairs))
            np.save(fn, S); jsecs = time.time() - t1
            pred = info['left'] * info['right'] * 2.0 ** -F_FINE if info.get('left') else 0.0
            log(dict(stage="sjoin", mid=mid, ci=ci, info=info, S_solutions=len(S), join_secs=round(jsecs, 2),
                     predicted=round(pred, 1),
                     excess_bits=round(float(np.log2(len(S) / pred)), 3) if len(S) and pred else None))
        start = progress.get(mid, {}).get('next', 0)
        t2 = time.time(); nd = 0; cores = 0; skipped = 0; ncand = 0; logcfg = 0.0
        for i in range(start, len(S)):
            cf = [S[i, k].astype(np.int64) for k in range(4)]
            try: ds, sizes, nc = DC.check(cf)
            except RuntimeError: skipped += 1; continue
            ncand += nc; logcfg += float(np.log2(np.prod(np.array(sizes, np.float64))))
            for d in ds:
                Xa, ok = Sv.assemble(d); nd += 1
                if not ok: continue
                F2 = to_template_coords(Xa)
                Z = AS.compressions(F2, (2, 3, 57), 2)
                Cf = AS.forced_corner(Z, (2, 3), T['layout'], T['K'], T['Kc'])
                H = AS.assemble(F2, (2, 3, 57), 2, T['layout'], T['K'], T['Kc'], Cf) if Cf is not None else None
                had = H is not None and AS.is_hadamard(H)
                cores += 1; tagc = f'{tag}_w{w}_{mid}_{cores}'
                np.save(os.path.join(outdir, f'core_{tagc}.npy'), F2)
                if had: np.savez_compressed(os.path.join(outdir, f'H1388_{tagc}.npz'), H=np.packbits(H > 0, axis=1))
                log(dict(CORE=tagc, mid=mid, s_index=i, verified_core=True, hadamard_1388=bool(had),
                         identical_to_published=bool((F2 == T['F']).all()),
                         sha=hashlib.sha256(F2.astype(np.int8).tobytes()).hexdigest()[:16]))
            if (i + 1) % ckpt == 0:
                log(dict(stage="dprogress", mid=mid, next=i + 1, of=len(S), D_solutions=nd, cores=cores,
                         secs=round(time.time() - t2, 1), info=info))
        n_checked = len(S) - start - skipped
        log(dict(stage="middle_done", mid=mid, ci=ci, S_solutions=len(S), checked_from=start, D_checked=n_checked,
                 skipped=skipped, D_solutions=nd, cores=cores, key_candidates=ncand,
                 mean_log2_Dlists=round(logcfg / max(1, n_checked), 2), dcheck_secs=round(time.time() - t2, 1)))
        cores_so_far += cores
        if stop_cores and cores_so_far >= stop_cores: break
    log(dict(stage="worker_end", worker=w, cores=cores_so_far))


if __name__ == "__main__":
    main()
