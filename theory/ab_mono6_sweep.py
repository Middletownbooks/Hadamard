r"""Sweep of the c = 6 GP slices of the open orders in the monomial border class, end to end [M].

Every configuration of theory/ab_cbatch.configs(v, 6) (v in 599, 659, 719, 839, 863, 917, 947, 971), smallest first
(G'-orbits, then index s), runs theory/ab_mono6_s2.run: screen A, exact stage 1 in the monomial class, the Singer
border of every solution, exact stage 2 of every bordered y whose block fibres are <= MAXFIB, assembly and
H H^T = 4v I for every core.  Configurations decided by the earlier probe (logs/ab/mono6/mono6_*.jsonl: screen A
fails, or stage 1 complete with 0 solutions) are skipped.  Per configuration: a time cap (SIGALRM) and the process
memory cap (RLIMIT_AS); a summary line in OUT.jsonl; per-y lines in OUTDIR/<tag>.jsonl (restart-safe).
usage: python3 theory/ab_mono6_sweep.py OUT.jsonl OUTDIR [--part i/W] [--secs 1800] [--L1max 30] [--memgb 3.2]"""
import sys, os, json, time, signal, resource, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD
from ab_cbatch import configs
from ab_mono6 import OPEN6
from ab_mono6_s2 import run

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TO(Exception): pass


def probe_decided():
    out = {}
    for f in glob.glob(os.path.join(ROOT, 'logs/ab/mono6/mono6_*.jsonl')):
        for l in open(f):
            try: d = json.loads(l)
            except Exception: continue
            if d.get('A') is False: out[json.dumps(d['cfgkey'])] = 'A fails'
            elif d.get('complete') and d.get('stage1_solutions') == 0: out[json.dumps(d['cfgkey'])] = 'stage 1 empty'
    return out


def main():
    a = sys.argv[1:]; outp, outdir = a[0], a[1]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    secs = int(opt('--secs', '1800')); L1max = int(opt('--L1max', '30')); memgb = float(opt('--memgb', '3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    os.makedirs(outdir, exist_ok=True)
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    prev = probe_decided()
    cl = []
    for v in OPEN6:
        for cfg in configs(v, 6):
            try: L1 = len(W.build(SD.cfg_from_json(json.dumps(cfg))).breps)
            except Exception: L1 = 10 ** 6
            cl.append((L1, cfg['s'], cfg['v'], cfg))
    cl.sort(key=lambda t: t[:3])
    cl = [t for i, t in enumerate(cl) if i % nW == part]
    fo = open(outp, 'a')
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for L1, s, v, cfg in cl:
        key = [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]
        k = json.dumps(key)
        if k in done: continue
        rec = dict(cfgkey=key, v=v, q=cfg['q'], s=s, f=cfg['f'], H1=cfg['H1'], U=cfg['U'], L1=L1)
        if k in prev: rec['status'] = prev[k] + ' (probe)'
        elif L1 > L1max: rec['status'] = f'skipped: L1 > {L1max}'
        else:
            fn = os.path.join(outdir, f"v{v}_q{cfg['q']}_s{s}_{abs(hash(k)) % 10 ** 8}.jsonl")
            rec['ylog'] = fn
            t0 = time.time(); signal.alarm(secs)
            try:
                rec['cores'] = run(json.dumps(cfg), fn); rec['status'] = 'done'
            except TO: rec['status'] = f'timeout {secs}'
            except (Exception, MemoryError) as e: rec['status'] = 'error'; rec['error'] = repr(e)[:200]
            signal.alarm(0)
            rec['secs'] = round(time.time() - t0, 1)
            try:
                lines = [json.loads(l) for l in open(fn)]
                st = [d for d in lines if d.get('stage') == 'stage1']
                if st: rec['stage1'] = {kk: st[-1].get(kk) for kk in ('A', 'alpha_quads', 'S_solutions', 'stage1_solutions', 'complete')}
                ys = [d for d in lines if 'yhash' in d]
                rec['ys_done'] = len(ys); rec['ys_bordered'] = sum(1 for d in ys if d.get('bordered'))
                rec['ys_s2_undecided'] = sum(1 for d in ys if d.get('bordered') and d.get('cores') is None)
                rec['cores_total'] = sum(d.get('cores') or 0 for d in ys)
            except Exception: pass
        fo.write(json.dumps(rec) + "\n"); fo.flush()
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
