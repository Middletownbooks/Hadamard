r"""Stage 2 production for a GP slice (b >= 3, c = 4 or 6) on given stage-1 solutions: the hierarchical solver of
theory/ab_stage2.py (S/D on axis 0 of G', middle level along LAXIS, exact D-check), every D-solution assembled on G
and verified (sum_k PAF_G = T), then bordered (c = 4: theory/ab_border_c4.border_c4; c = 6: the monomial Singer border
theory/ab_mono_border.border_for, else the 347 template) into H(4v) and checked H H^T = 4v I [M].
Restart-safe: one JSONL line per (y, middle value); finished ones are skipped.
usage: python3 theory/ab_prod_generic.py '<cfg json>' YFILE.jsonl --tag T [--ys 0,1,2] [--laxis 2] [--R 800 800 ...]
                                         [--K 4 256 ...] [--middles 400] [--seed 1] [--stop-on-core] [--memgb 3.2]
(YFILE: JSONL with key "y" per line, e.g. the output of a stage-1 run)"""
import sys, os, json, time, hashlib, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_stage2 as S2m, ab_sweep as W, ab_assemble as AS
from ab_s1sd import cfg_from_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def border_and_matrix(cfg, X):
    c = cfg['c']; mods = tuple(cfg['mods']); nQ = cfg['nQ']
    if c == 4:
        from ab_border_c4 import border_c4
        r = border_c4(X, mods, nQ)
        return (r[0], 'c4-frames') if r is not None else (None, None)
    import ab_mono_border as MB
    Z = AS.compressions(X, mods, nQ)
    br = MB.border_for(Z, tuple(cfg['Q']))
    if br is not None:
        Cf = AS.forced_corner(Z, tuple(cfg['Q']), MB.LAYOUT, br['K'], br['Kc'])
        if Cf is not None:
            H = AS.assemble(X, mods, nQ, MB.LAYOUT, br['K'], br['Kc'], Cf)
            if AS.is_hadamard(H): return H, br['template']
    return None, None


def main():
    a = sys.argv[1:]
    def opt(name, n=1, default=None):
        if name in a:
            i = a.index(name); return [int(x) for x in a[i + 1:i + 1 + n]] if n > 1 else a[i + 1]
        return default
    memgb = float(opt('--memgb', default='3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    cfgj = a[0]; yfile = a[1]
    cfg = cfg_from_json(cfgj); sl = W.build(cfg)
    tag = opt('--tag', default='t0'); Lax = int(opt('--laxis', default='2'))
    seed = int(opt('--seed', default='1')); nmid = int(opt('--middles', default='400'))
    ys = [json.loads(l)['y'] for l in open(yfile) if l.strip()]
    pick = [int(x) for x in opt('--ys', default=','.join(map(str, range(len(ys))))).split(',')]
    outdir = os.path.join(ROOT, 'logs/ab/prod', tag); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, 'run.jsonl')
    done = set()
    if os.path.exists(logf):
        for l in open(logf):
            try: d = json.loads(l)
            except Exception: continue
            if 'mid' in d: done.add((d['yi'], d['mid']))
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    cores = 0
    for yi in pick:
        Y = np.array(ys[yi], np.int64)
        t0 = time.time()
        Sv = S2m.FastSolver2(cfg, Y, Laxis=Lax, log=lambda s: None)
        nrows = len(Sv.rows)
        R = opt('--R', nrows, None) if '--R' in a else [800] * nrows
        K = opt('--K', nrows, None) if '--K' in a else [256] * nrows
        R = R if isinstance(R, list) else [int(R)]; K = K if isinstance(K, list) else [int(K)]
        Sv.prepare_rows(R=R, K=K, seed=seed)
        cands, ncommon = Sv.choose_middles(top=nmid, seed=seed + 100)
        log(dict(stage="setup", tag=tag, yi=yi, rows=nrows, R=R, K=K, middles=len(cands), common=ncommon,
                 parents=[len(p) for p in Sv.pidx], secs=round(time.time() - t0, 1)))
        for ci, (c, sc) in enumerate(cands):
            mid = hashlib.sha1(c.tobytes()).hexdigest()[:16]
            if (yi, mid) in done: continue
            t1 = time.time()
            try: sols, info = Sv.s_join(c)
            except RuntimeError as e:
                log(dict(yi=yi, mid=mid, ci=ci, skipped=str(e)[:200])); continue
            t2 = time.time(); nd = 0; dcand = 0.0; found = 0
            for cf in sols:
                try: ds, sizes = Sv.d_check_fast(cf)
                except RuntimeError: continue
                dcand += float(np.prod([float(x) for x in sizes]))
                for d in ds:
                    Xa, ok = Sv.assemble(d); nd += 1
                    if not ok: continue
                    found += 1; cores += 1
                    ct = f'{tag}_{yi}_{mid}_{found}'
                    np.save(os.path.join(outdir, f'core_{ct}.npy'), Xa.astype(np.int8))
                    H, how = border_and_matrix(cfg, Xa)
                    if H is not None:
                        np.savez_compressed(os.path.join(outdir, f'H{4 * cfg["v"]}_{ct}.npz'), H=np.packbits(H > 0, axis=1))
                    log(dict(CORE=ct, verified_core=True, hadamard=H is not None, border=how, order=4 * cfg['v']))
            log(dict(yi=yi, mid=mid, ci=ci, **info, S_solutions=len(sols),
                     D_candidates_log2=round(float(np.log2(dcand)), 2) if dcand else None, D_solutions=nd, cores=found,
                     join_secs=round(t2 - t1, 1), dcheck_secs=round(time.time() - t2, 1)))
            if cores and '--stop-on-core' in a: return


if __name__ == "__main__":
    main()
