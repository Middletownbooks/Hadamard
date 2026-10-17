r"""Typical-y stage-2 level budget for the b >= 3 configurations the exact batches skipped (more than 26 G'-orbits) [E].

The 117 entries of logs/ab/c4_batch.jsonl, c6_bordered.jsonl, c8_bordered.jsonl with `skipped: "L1 > 26"` are 112
distinct configurations (five 599 c = 6 entries appear twice, once with the older key without c).  Per configuration:
the energy screen [P], then theory/ab_typical_levels.budget (H2, C2, S-level bits C_S, exact D-freedom h_D, D_def)
with no stage-2 orbit cap, and the cost model of report/paley_523_answer.md 3 (report/b599_answer.md 5):
cost per core = 2^(D_def + C_S - H_S/2), H_S = H2 - h_D.  Restart-safe JSONL; time and memory cap per configuration.
usage: python3 theory/s1z_budget_batch.py OUT.jsonl [--part i/W] [--secs 1800] [--memgb 3.5] [--N 4000]"""
import sys, os, json, time, signal, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_s1sd import cfg_from_json
from ab_typical_levels import budget

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = ['logs/ab/c4_batch.jsonl', 'logs/ab/c6_bordered.jsonl', 'logs/ab/c8_bordered.jsonl']
QOF = {4: [2, 2], 6: [2, 3], 8: [2, 2, 2]}


class TO(Exception): pass


def skipped_configs():
    """the distinct configurations skipped for L1 > 26, in the order of the batch logs"""
    out = []; seen = set()
    for f, cf in zip(SRC, (4, 6, 8)):
        for l in open(os.path.join(ROOT, f)):
            z = json.loads(l)
            if 'skipped' not in z or 'L1 >' not in z['skipped']: continue
            ck = z['cfgkey']
            if len(ck) == 6: v, c, q, s, H1, U = ck
            else: (v, q, s, H1, U), c = ck, cf
            key = json.dumps([v, c, q, s, H1, U])
            if key in seen: continue
            seen.add(key)
            Q = QOF[c]; h = (v + 1) // c - 1
            out.append(dict(v=v, c=c, b=c - 1, h=h, m=c * h, Q=Q, H1=H1, q=q, s=s, f=(q - 1) // s,
                            mods=Q + H1 + [q], nQ=len(Q), U=U, L1=z['L1'], orbits=z['orbits']))
    return out


def cfgkey(cfg):
    return [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]


def main():
    a = sys.argv[1:]; outp = a[0]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    secs = int(opt('--secs', '1800')); memgb = float(opt('--memgb', '3.5')); N = int(opt('--N', '4000'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    cfgs = sorted(skipped_configs(), key=lambda c: (c['orbits'], c['L1']))
    cfgs = [c for i, c in enumerate(cfgs) if i % nW == part]
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for cfg in cfgs:
        key = cfgkey(cfg)
        if json.dumps(key) in done: continue
        rec = dict(cfgkey=key, L1=cfg['L1'], orbits=cfg['orbits'], f=cfg['f'])
        t0 = time.time(); signal.alarm(secs)
        try:
            c2 = cfg_from_json(json.dumps(cfg)); sl = W.build(c2)
            w = sl.bsz.astype(np.int64)
            emin = int(sum(int(w[j]) * min(y * y for y in L) for j, L in enumerate(sl.letters)))
            rec['energy_ok'] = bool(4 * emin <= int(sl.T1[0]))
            if rec['energy_ok']:
                r = budget(c2, seed=3, N=N)
                for k in ('H2', 'C2', 'rank', 'h_D', 'C_D_given_S', 'D_def', 'stage2_log2E'): rec[k] = r.get(k)
                rec['C_S'] = r['levels']['S(axis0)']['bits']
                if None not in (rec['D_def'], rec['C_S'], rec['H2'], rec['h_D']):
                    HS = rec['H2'] - rec['h_D']
                    rec['H_S'] = round(HS, 2)
                    rec['log2_cost_per_core'] = round(max(rec['D_def'], 0) + rec['C_S'] - HS / 2, 2)
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
        signal.alarm(0)
        rec['secs'] = round(time.time() - t0, 1)
        with open(outp, 'a') as fo: fo.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
