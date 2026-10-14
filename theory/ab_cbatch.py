r"""Batch exact stage 1 for GP slices of the open orders with c = 4, 6 or 8 [M].
c = 4: Q = Z_2^2, unrestricted (every c = 4 compression has a border, report/bge3_answer.md 4.1 and the 479 transfer);
c = 6: Q = Z_2 x Z_3, restricted to the 347 border class;  c = 8: Q = Z_2^3, restricted to the 359 border class.
Stage-1 solutions are checked exactly against the frame (theory/ab_borderfit.py) for c = 6, 8.
usage: python3 theory/ab_cbatch.py C OUT.jsonl [--vs ...] [--secs 1200] [--L1max 26] [--maxsol 20000]"""
import sys, os, json, time, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_borderfit as BF
from ab_s1_bordered import allowed
from ab_lib import factorint

QOF = {4: [2, 2], 6: [2, 3], 8: [2, 2, 2]}
TPL = {4: None, 6: 1388, 8: 1436}
OPEN = [571, 599, 643, 659, 719, 739, 751, 823, 839, 859, 863, 883, 907, 917, 919, 947, 955, 971, 991]

def configs(v, c):
    if (v + 1) % c: return []
    h = (v + 1) // c - 1; out = []; Q = QOF[c]; nQ = len(Q)
    for q, e in factorint(h).items():
        if e != 1 or q < 5: continue
        hp = h // q
        for Ht in W.abelian_types(hp):
            mods = tuple(Q) + tuple(Ht) + (q,)
            for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                f = (q - 1) // s
                if f < 2: continue
                Hax = list(range(nQ, nQ + len(Ht)))
                gens_opts = [[]]
                for ai in Hax:
                    for (o, u) in W.unit_subgroup_gens(mods[ai]): gens_opts.append([[ai, o, u]])
                if c == 6:
                    base = list(gens_opts)
                    for g in base:
                        gens_opts.append([["neg", [1]]] + g)
                        if Hax: gens_opts.append([["neg", [1] + Hax]] + g)
                for U in gens_opts:
                    out.append(dict(v=v, c=c, b=c - 1, h=h, m=c * h, Q=Q, H1=list(Ht), q=q, s=s, f=f, mods=list(mods), nQ=nQ, U=U))
    return out

class TO(Exception): pass

if __name__ == "__main__":
    a = sys.argv[1:]; c = int(a[0]); outp = a[1]
    vs = [int(x) for x in a[a.index('--vs') + 1].split(',')] if '--vs' in a else OPEN
    secs = int(a[a.index('--secs') + 1]) if '--secs' in a else 1200
    L1max = int(a[a.index('--L1max') + 1]) if '--L1max' in a else 26
    maxsol = int(a[a.index('--maxsol') + 1]) if '--maxsol' in a else 20000
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    fo = open(outp, 'a')
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    N = TPL[c]
    for v in vs:
        for cfg in configs(v, c):
            key = [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]
            if json.dumps(key) in done: continue
            c2 = SD.cfg_from_json(json.dumps(cfg))
            try: sl = W.build(c2)
            except AssertionError: continue
            L1 = len(sl.breps)
            rec = dict(cfgkey=key, v=v, c=c, q=cfg['q'], s=cfg['s'], f=cfg['f'], H1=cfg['H1'], U=cfg['U'], orbits=int(sl.r), L1=L1)
            if L1 > L1max: rec['skipped'] = f'L1 > {L1max}'
            else:
                t0 = time.time(); signal.alarm(secs)
                try:
                    Z = allowed(c2, N) if N else None
                    P = SD.S1SD(c2, zallow=Z, sl=sl)
                    rq, A = SD.row_quads_of(c2, sl=sl)
                    rec['A'] = A.get('A')
                    if A.get('A'):
                        res, sols = P.solve(rq, max_sol=maxsol, log=lambda s: None)
                        nfit = 0; ex = None
                        for s_ in sols:
                            Yb = np.array(s_['y'])
                            if N:
                                cq = len(P.zallow[0]); z = np.zeros((4, cq), np.int64)
                                for k in range(4):
                                    for qq in range(cq): z[k, qq] = Yb[k][P.qfull == qq].sum()
                                ok = bool(BF.fits(z, N))
                            else: ok = True
                            if ok:
                                nfit += 1
                                if ex is None: ex = Yb.tolist()
                        rec.update(alpha_quads=res['alpha_quads'], S_solutions=res['S_solutions'],
                                   stage1_solutions=res['stage1_solutions'], complete=res['complete'], bordered=nfit, example_y=ex)
                except TO: rec['timeout'] = secs
                except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
                signal.alarm(0)
                rec['secs'] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec) + "\n"); fo.flush()
            print(json.dumps({k: v_ for k, v_ in rec.items() if k != 'example_y'}), flush=True)
