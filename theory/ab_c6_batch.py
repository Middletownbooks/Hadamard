r"""Batch: exact stage 1 restricted to the 347 border class for every c = 6 GP slice of the open orders [M].
Configurations: Q = Z_2 x Z_3 (axis order (2, 3) so that the S/D axis is Z_2), H = H' x Z_q split, q a prime with
q || h, index s of M on Z_q (f >= 2), U = trivial / one cyclic subgroup on one axis of H' / -1 on Q's Z_3 / -1 on
Q's Z_3 diagonal with -1 on H' / combinations; G'-orbits <= 26.  Restart-safe (one JSONL line per config).
usage: python3 theory/ab_c6_batch.py OUT.jsonl [--vs 863,917,...] [--secs 1800]"""
import sys, os, json, time, itertools, math, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_borderfit as BF
from ab_s1_bordered import allowed
from ab_lib import factorint, primitive_root

def configs(v):
    c = 6
    if (v + 1) % c: return []
    h = (v + 1) // c - 1; out = []
    for q, e in factorint(h).items():
        if e != 1 or q < 5: continue
        hp = h // q
        for Ht in W.abelian_types(hp):
            mods = (2, 3) + tuple(Ht) + (q,)
            for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                f = (q - 1) // s
                if f < 2: continue
                Hax = list(range(2, 2 + len(Ht)))
                gens_opts = [[]]
                for ai in Hax:
                    for (o, u) in W.unit_subgroup_gens(mods[ai]): gens_opts.append([[ai, o, u]])
                base = list(gens_opts)
                for g in base:
                    gens_opts.append([["neg", [1]]] + g)
                    if Hax: gens_opts.append([["neg", [1] + Hax]] + g)
                for U in gens_opts:
                    cfg = dict(v=v, c=c, b=5, h=h, m=c * h, Q=[2, 3], H1=list(Ht), q=q, s=s, f=f, mods=list(mods), nQ=2, U=U)
                    out.append(cfg)
    return out

class TO(Exception): pass

if __name__ == "__main__":
    a = sys.argv[1:]; outp = a[0]
    vs = [int(x) for x in a[a.index('--vs') + 1].split(',')] if '--vs' in a else [599, 659, 719, 839, 863, 917, 947, 971]
    secs = int(a[a.index('--secs') + 1]) if '--secs' in a else 1800
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    fo = open(outp, 'a')
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for v in vs:
        for cfg in configs(v):
            key = [cfg['v'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]
            if json.dumps(key) in done: continue
            c2 = SD.cfg_from_json(json.dumps(cfg))
            try:
                sl = W.build(c2)
            except AssertionError as e:
                continue
            L1 = len(sl.breps)
            rec = dict(cfgkey=key, v=v, q=cfg['q'], s=cfg['s'], f=cfg['f'], H1=cfg['H1'], U=cfg['U'], orbits=int(sl.r), L1=L1)
            if L1 > 26: rec['skipped'] = 'L1 > 26'
            else:
                t0 = time.time(); signal.alarm(secs)
                try:
                    Z = allowed(c2, 1388)
                    P = SD.S1SD(c2, zallow=Z, sl=sl)
                    rq, A = SD.row_quads_of(c2, sl=sl)
                    rec['A'] = A.get('A')
                    if A.get('A'):
                        res, sols = P.solve(rq, max_sol=20000, log=lambda s: None)
                        nfit = 0; ex = None
                        for s_ in sols:
                            Yb = np.array(s_['y']); z = np.zeros((4, 6), np.int64)
                            for k in range(4):
                                for qq in range(6): z[k, qq] = Yb[k][P.qfull == qq].sum()
                            if BF.fits(z, 1388):
                                nfit += 1
                                if ex is None: ex = Yb.tolist()
                        rec.update(alpha_quads=res['alpha_quads'], S_solutions=res['S_solutions'],
                                   stage1_solutions=res['stage1_solutions'], complete=res['complete'], bordered=nfit, example_y=ex)
                except TO:
                    rec['timeout'] = secs
                except (Exception, MemoryError) as e:
                    rec['error'] = repr(e)[:200]
                signal.alarm(0)
                rec['secs'] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec) + "\n"); fo.flush()
            print(json.dumps({k: v_ for k, v_ in rec.items() if k != 'example_y'}), flush=True)
