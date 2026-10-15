r"""Exact stage 1 for c = 6 GP slices of the open orders in the monomial (Singer) border class [M].

Allowed per-block Q-compressions: z = s*1 + a*delta_q with s = h (mod 2), a in {0, +-2, +-4}, q in Q (|z| <= h).
The stage-1 solver (theory/ab_s1sd.py, zallow) enumerates every stage-1 solution whose four compressions lie in that
set; the Q-level equation then forces the profile {4, 2, 0, 0}.  Every solution is bordered exactly
(theory/ab_mono_border.py: a Singer template carried to the solution's arrangement, signs and spike positions; the
forced corner is checked to be +-1).  Restart-safe: one JSONL line per configuration.
Configurations are taken smallest first (by G'-orbits); each runs under a time cap and the process under a memory cap
(RLIMIT_AS), so a blow-up is recorded as an error instead of killing the run.
usage: python3 theory/ab_mono6.py OUT.jsonl [--vs 863,719,...] [--secs 1200] [--L1max 40] [--maxsol 20000] [--smax 21]
                                           [--memgb 3.5]"""
import sys, os, json, time, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_mono_border as MB
from ab_cbatch import configs

OPEN6 = [599, 659, 719, 839, 863, 917, 947, 971]


def mono_allowed(cfg, smax=21):
    c = 6; h = cfg['h']; par = h % 2; out = []
    for s in range(-smax, smax + 1):
        if (s - par) % 2: continue
        base = np.full(c, s, np.int64)
        out.append(base.tolist())
        for q in range(c):
            for a in (2, -2, 4, -4):
                z = base.copy(); z[q] += a
                if np.abs(z).max() <= h: out.append(z.tolist())
    return out


class TO(Exception): pass


def run_cfg(cfg, secs, maxsol, smax, keep=50):
    c2 = SD.cfg_from_json(json.dumps(cfg))
    sl = W.build(c2)
    rec = dict(orbits=int(sl.r), L1=len(sl.breps))
    Z = mono_allowed(cfg, smax)
    P = SD.S1SD(c2, zallow=Z, sl=sl)
    rq, A = SD.row_quads_of(c2, sl=sl)
    rec['A'] = A.get('A')
    if not A.get('A'): return rec, []
    res, sols = P.solve(rq, max_sol=maxsol, log=lambda s: None)
    rec.update(alpha_quads=res['alpha_quads'], S_solutions=res['S_solutions'], stage1_solutions=res['stage1_solutions'],
               complete=res['complete'])
    bordered = []
    cq = P.zallow.shape[1]
    for s_ in sols:
        Yb = np.array(s_['y']); z = np.zeros((4, cq), np.int64)
        for k in range(4):
            for qq in range(cq): z[k, qq] = Yb[k][P.qfull == qq].sum()
        br = MB.border_for(z, tuple(cfg['Q']))
        if br is not None:
            bordered.append(dict(y=Yb.tolist(), z=z.tolist(), a=br['a'], q=br['q'], template=br['template']))
    rec['bordered'] = len(bordered)
    return rec, bordered[:keep]


def main():
    a = sys.argv[1:]; outp = a[0]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    vs = [int(x) for x in opt('--vs', ','.join(map(str, OPEN6))).split(',')]
    secs = int(opt('--secs', '1200')); L1max = int(opt('--L1max', '40')); maxsol = int(opt('--maxsol', '20000'))
    smax = int(opt('--smax', '21'))
    import resource
    memgb = float(opt('--memgb', '3.5'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30), int(memgb * 2 ** 30)))
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    soldir = os.path.join(os.path.dirname(os.path.abspath(outp)), 'mono6_sols'); os.makedirs(soldir, exist_ok=True)
    fo = open(outp, 'a')
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for v in vs:
        cl = []
        for cfg in configs(v, 6):
            try: cl.append((len(W.build(SD.cfg_from_json(json.dumps(cfg))).breps), cfg))
            except Exception: cl.append((10 ** 6, cfg))
        cl.sort(key=lambda t: t[0])
        for _, cfg in cl:
            key = [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]
            if json.dumps(key) in done: continue
            rec = dict(cfgkey=key, v=v, q=cfg['q'], s=cfg['s'], f=cfg['f'], H1=cfg['H1'], U=cfg['U'])
            t0 = time.time(); signal.alarm(secs)
            try:
                sl = W.build(SD.cfg_from_json(json.dumps(cfg)))
                if len(sl.breps) > L1max:
                    rec.update(orbits=int(sl.r), L1=len(sl.breps), skipped=f'L1 > {L1max}')
                else:
                    r, bordered = run_cfg(cfg, secs, maxsol, smax)
                    rec.update(r)
                    if bordered:
                        fn = os.path.join(soldir, f"v{v}_q{cfg['q']}_s{cfg['s']}_{abs(hash(json.dumps(key))) % 10 ** 8}.json")
                        json.dump(dict(cfg=cfg, solutions=bordered), open(fn, 'w'))
                        rec['saved'] = os.path.relpath(fn, os.path.dirname(os.path.dirname(os.path.abspath(outp))))
            except TO: rec['timeout'] = secs
            except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
            signal.alarm(0)
            rec['secs'] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec) + "\n"); fo.flush()
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
