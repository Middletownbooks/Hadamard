r"""Monomial-class c = 6 slices end to end: exact stage 1, monomial border, exact stage 2, assembly [M].

For a c = 6 GP slice of an open order: exact stage 1 restricted to monomial Q-compressions (theory/ab_mono6.py
mono_allowed; S/D solver theory/ab_s1sd.py), every solution's border from a Singer template
(theory/ab_mono_border.py border_for: forced corner exactly +-1), then exact stage 2 per bordered y by the one-shot
2+2 join of the full PAF over G (theory/ab_small_exact.solve_y, small fibres only).  Every core is re-checked
(sum_k PAF_G = T), assembled with its border into the 4v x 4v matrix (theory/ab_assemble.assemble, 347 layout) and
checked H H^T = 4v I.  Cores and matrices are saved under logs/ab/mono6/cores/ (the first NSAVE; all are verified and counted).
Restart-safe: one JSONL line per stage-1 solution (key: its index and a hash of y); finished ones are skipped.
usage: python3 theory/ab_mono6_s2.py '<cfg json>' OUT.jsonl [--maxsol 100000] [--smax 21] [--max-fibre 65536]
                                    [--memgb 3.2] [--stop-on-core] [--nsave 20]"""
import sys, os, json, time, hashlib, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_mono_border as MB, ab_gauss as AGs, ab_assemble as AS
from ab_mono6 import mono_allowed
from ab_small_exact import solve_y

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(cfgj, outp, maxsol=100000, smax=21, cap=65536, nsave=20, stop_on_core=False):
    cfg = SD.cfg_from_json(cfgj); sl = W.build(cfg)
    tag = f"v{cfg['v']}_q{cfg['q']}_s{cfg['s']}_" + hashlib.sha1(json.dumps(json.loads(cfgj) if not cfgj.startswith('@') else cfgj).encode()).hexdigest()[:8]
    coredir = os.path.join(ROOT, 'logs/ab/mono6/cores'); os.makedirs(coredir, exist_ok=True)
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l)
            except Exception: continue
            if 'yhash' in d: done.add(d['yhash'])
    fo = open(outp, 'a')
    def log(rec):
        fo.write(json.dumps(rec) + "\n"); fo.flush(); print(json.dumps(rec), flush=True)
    t0 = time.time()
    Z = mono_allowed(json.loads(cfgj) if not cfgj.startswith('@') else json.load(open(cfgj[1:])), smax)
    P = SD.S1SD(cfg, zallow=Z, sl=sl)
    rq, A = SD.row_quads_of(cfg, sl=sl)
    if not A.get('A'):
        log(dict(stage='stage1', tag=tag, A=False)); fo.close(); return 0
    res, sols = P.solve(rq, max_sol=maxsol, log=lambda s: None)
    log(dict(stage='stage1', tag=tag, **{k: v for k, v in res.items() if k not in ('stage', 'secs')}, secs=round(time.time() - t0, 1)))
    reps, _ = AGs.lag_reps(sl.G, sl.perms)
    gw = np.random.default_rng(5).integers(1, 2 ** 62, size=len(reps), dtype=np.int64)
    cq = P.zallow.shape[1]; Qmods = tuple(cfg['Q']); nQ = cfg['nQ']
    ncore = 0
    for i, s_ in enumerate(sols):
        Y = np.array(s_['y'], np.int64)
        yh = hashlib.sha1(Y.tobytes()).hexdigest()[:16]
        if yh in done: continue
        z = np.zeros((4, cq), np.int64)
        for k in range(4):
            for qq in range(cq): z[k, qq] = Y[k][P.qfull == qq].sum()
        br = MB.border_for(z, Qmods)
        rec = dict(i=i, yhash=yh, bordered=br is not None, z=z.tolist())
        if br is not None:
            t1 = time.time()
            hits, m = solve_y(sl, Y, reps, gw, cap)
            rec.update(fibres=m, cores=None if hits is None else len(hits), s2_secs=round(time.time() - t1, 2))
            for Xc in (hits or []):
                ok = bool((sl.G.paf(Xc).sum(0) == sl.T).all())
                Cf = AS.forced_corner(AS.compressions(Xc, cfg['mods'], nQ), Qmods, MB.LAYOUT, br['K'], br['Kc'])
                H = AS.assemble(Xc, cfg['mods'], nQ, MB.LAYOUT, br['K'], br['Kc'], Cf) if Cf is not None else None
                had = H is not None and AS.is_hadamard(H)
                ncore += 1; ct = f'{tag}_{yh}_{ncore}'
                if ncore > nsave:
                    if not had: log(dict(CORE=ct, core_ok=ok, hadamard=False))
                    continue
                np.save(os.path.join(coredir, f'core_{ct}.npy'), Xc.astype(np.int8))
                if had:
                    np.savez_compressed(os.path.join(coredir, f'H{4 * cfg["v"]}_{ct}.npz'), H=np.packbits(H > 0, axis=1),
                                        K=br['K'], Kc=br['Kc'], C=Cf, template=br['template'])
                log(dict(CORE=ct, core_ok=ok, hadamard=bool(had), order=4 * cfg['v'], template=br['template']))
        log(rec)
        if ncore and stop_on_core: break
    log(dict(stage='end', tag=tag, solutions=len(sols), cores=ncore, secs=round(time.time() - t0, 1)))
    fo.close()
    return ncore


def main():
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    memgb = float(opt('--memgb', '3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    run(a[0], a[1], int(opt('--maxsol', '100000')), int(opt('--smax', '21')), int(opt('--max-fibre', '65536')),
        int(opt('--nsave', '20')), '--stop-on-core' in a)


if __name__ == "__main__":
    main()
