r"""Batch of typical-y stage-2 level budgets for the c = 4, 6, 8 GP configurations of the open orders [E].

Per configuration: orbits, G'-orbits; the energy screen [P] (4 * min block energy > T1(0) => stage 1 empty);
the whole-slice and stage-1 Gaussians (ab_gauss.whole, ab_gauss.stage1); the typical-y level budget
(theory/ab_typical_levels.budget: H2, C2, S-level bits, h_D, D_def).  Configurations with more than ORBMAX stage-2
orbits per block are recorded and skipped.  Restart-safe JSONL (one line per configuration), time cap per
configuration (SIGALRM), memory cap per process (RLIMIT_AS).
usage: python3 theory/ab_levels_batch.py OUT.jsonl [--part i/W] [--vs ...] [--cs 4,6,8] [--orbmax 130] [--secs 900]
                                         [--memgb 3.2] [--wmin -30]"""
import sys, os, json, time, signal, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json
from ab_cbatch import configs, OPEN
from ab_typical_levels import budget


class TO(Exception): pass


def main():
    a = sys.argv[1:]; outp = a[0]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    vs = [int(x) for x in opt('--vs', ','.join(map(str, OPEN))).split(',')]
    cs = [int(x) for x in opt('--cs', '4,6,8').split(',')]
    orbmax = int(opt('--orbmax', '130')); secs = int(opt('--secs', '900'))
    memgb = float(opt('--memgb', '3.2')); wmin = float(opt('--wmin', '-30'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    allc = []
    for v in vs:
        for c in cs:
            if (v + 1) % c: continue
            allc += configs(v, c)
    allc = [cfg for i, cfg in enumerate(allc) if i % nW == part]
    fo = open(outp, 'a')
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for cfg in allc:
        key = [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]
        if json.dumps(key) in done: continue
        rec = dict(cfgkey=key, v=cfg['v'], c=cfg['c'], q=cfg['q'], s=cfg['s'], f=cfg['f'], H1=cfg['H1'], U=cfg['U'])
        t0 = time.time(); signal.alarm(secs)
        try:
            c2 = cfg_from_json(json.dumps(cfg)); sl = W.build(c2)
            rec.update(orbits=int(sl.r), L1=len(sl.breps))
            st2 = (int(sl.r) - len(sl.breps))
            w = sl.bsz.astype(np.int64)
            emin = int(sum(int(w[j]) * min(y * y for y in L) for j, L in enumerate(sl.letters)))
            rec['energy_ok'] = bool(4 * emin <= int(sl.T1[0]))
            if not rec['energy_ok']: rec['status'] = 'stage 1 empty (energy) [P]'
            elif st2 > orbmax: rec['skipped'] = f'stage-2 orbits {st2} > {orbmax}'
            else:
                wh = AGs.whole(sl, N=4000); rec['whole_log2E'] = None if wh.get('log2E') is None else float(wh['log2E'])
                if rec['whole_log2E'] is None or rec['whole_log2E'] < wmin:
                    rec['skipped'] = f'whole < {wmin}'
                    raise StopIteration
                try:
                    s1 = AGs.stage1(sl, N=3000)
                    rec['stage1_log2E'] = None if s1.get('stage1_log2E') is None else float(s1['stage1_log2E'])
                    rec['stage1_mahal'] = s1.get('mahal'); rec['stage1_log2Ny'] = s1.get('log2Ny')
                except (Exception, MemoryError) as e: rec['stage1_error'] = repr(e)[:120]
                r = budget(c2, seed=3, N=4000)
                for k in ('H2', 'C2', 'rank', 'h_D', 'C_D_given_S', 'D_def', 'stage2_log2E'): rec[k] = r.get(k)
                rec['C_S'] = r['levels']['S(axis0)']['bits']
        except StopIteration: pass
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
        signal.alarm(0)
        rec['secs'] = round(time.time() - t0, 1)
        fo.write(json.dumps(rec) + "\n"); fo.flush()
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
