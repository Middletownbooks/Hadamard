r"""Batch z-conditioned stage-1 screen over the b >= 3 configurations the exact batches skipped (more than 26
G'-orbits; logs/ab/c4_batch.jsonl, c6_bordered.jsonl, c8_bordered.jsonl) [counts exact per z; totals E, calibrated on
347 at about 2.3 bits low, report/prompt_flatpool_stage1.md].

Per configuration (theory/ab_s1_zcount.py): compressions within the Parseval bound, all ordered z-quadruples meeting
the Q-level equation, and the stage-1 total estimated from K random z-quadruples (exact Q-level/lag-0 counts, lattice
Gaussian for the remaining lags), scaled to all z-quadruples.  A configuration with no compression or no z-quadruple
has an empty stage 1 [M].  Restart-safe (configurations already in OUT are skipped); a per-configuration time limit.
usage: python3 theory/ab_s1_zscreen.py OUT.jsonl [--K 12] [--N 400] [--secs 900] [--maxL1 200]"""
import sys, os, json, time, math, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from ab_s1_zcount import ZCount

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = ['logs/ab/c4_batch.jsonl', 'logs/ab/c6_bordered.jsonl', 'logs/ab/c8_bordered.jsonl']
QOF = {4: [2, 2], 6: [2, 3], 8: [2, 2, 2]}


class TO(Exception): pass


def skipped_configs():
    out = []
    for f, cf in zip(SRC, (4, 6, 8)):
        for l in open(os.path.join(ROOT, f)):
            z = json.loads(l)
            if 'skipped' not in z or 'L1 >' not in z['skipped']: continue
            ck = z['cfgkey']
            if len(ck) == 6: v, c, q, s, H1, U = ck
            else: (v, q, s, H1, U), c = ck, cf                      # c6_bordered keys omit c
            Q = QOF[c]; h = (v + 1) // c - 1
            out.append(dict(v=v, c=c, b=c - 1, h=h, m=c * h, Q=Q, H1=H1, q=q, s=s, f=(q - 1) // s,
                            mods=Q + H1 + [q], nQ=len(Q), U=U, L1=z['L1']))
    return out


def main():
    a = sys.argv[1:]; outp = a[0]
    K = int(a[a.index('--K') + 1]) if '--K' in a else 12
    N = int(a[a.index('--N') + 1]) if '--N' in a else 400
    secs = int(a[a.index('--secs') + 1]) if '--secs' in a else 900
    maxL1 = int(a[a.index('--maxL1') + 1]) if '--maxL1' in a else 200
    done = set()
    if os.path.exists(outp):
        for l in open(outp): done.add(json.dumps(json.loads(l)['cfgkey']))
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    cfgs = sorted(skipped_configs(), key=lambda c: c['L1'])
    for cfg in cfgs:
        key = [cfg['v'], cfg['c'], cfg['q'], cfg['s'], cfg['H1'], cfg['U']]     # c always included here
        if json.dumps(key) in done or cfg['L1'] > maxL1: continue
        rec = dict(cfgkey=key, L1=cfg['L1']); t0 = time.time(); signal.alarm(secs)
        try:
            c2 = SD.cfg_from_json(json.dumps(cfg)); Zc = ZCount(c2)
            rq, A = SD.row_quads_of(c2, sl=Zc.sl)
            rows = sorted({int(s_) * int(x) for qd in rq for x in qd for s_ in (1, -1)})
            dev = float(sum(abs(x) for x in np.real(np.fft.fftn(Zc.TQ.reshape(Zc.Q.mods))).ravel()[1:])) / Zc.c
            Z = Zc.block_z(dev, rows) if rows else {}
            rec.update(screenA=A.get('A'), row_values=len(rows), compressions=len(Z))
            Q4 = Zc.zquads(Z) if Z else []
            rec['zquads'] = len(Q4)
            if Q4:
                rng = np.random.default_rng(1); tot = 0.0; ests = []
                for i in rng.choice(len(Q4), size=min(K, len(Q4)), replace=False):
                    r = Zc.estimate(Q4[int(i)], Z, N=N, seed=int(i) + 1)
                    ests.append(r.get('log2E'))
                    if r.get('log2E') is not None: tot += 2.0 ** r['log2E']
                rec['est_per_zquad'] = ests
                rec['log2_stage1_total'] = round(math.log2(tot * len(Q4) / min(K, len(Q4))), 2) if tot > 0 else None
            else:
                rec['stage1_empty'] = True                      # no admissible compression quadruple [M]
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
        signal.alarm(0); rec['secs'] = round(time.time() - t0, 1)
        with open(outp, 'a') as f: f.write(json.dumps(rec) + '\n')
        print(json.dumps(rec), flush=True)


if __name__ == '__main__':
    main()
