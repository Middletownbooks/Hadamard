r"""Cross-checks of the <61> = Paley(19) stage-1 sweep (report section 7), on random samples of classes that model A
(theory/p571_p19cpsat.py: lags via the A/X/s algebra, plus Lemma 7.4 at psi) proved INFEASIBLE:
  y     - model y (pure S1 written directly in y = A + 18 W; no X/s algebra, no Lemma 7.4), 1 worker;
  dp    - the solver-free single-lag DP of theory/p571_p19lagdp.py: fraction of classes with a one-lag certificate;
  pass2 - a 4-worker CP-SAT pass on a random sample of the unresolved (UNKNOWN) classes.
usage: python3 theory/p571_p19crosscheck.py {y|dp|pass2} nsample [secs]"""
import sys, os, json, pickle, glob, random, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p571_p19cpsat as C
import p571_p19lagdp as D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(ROOT, 'theory/logs/p571/p19')

def load():
    cl = pickle.load(open(os.path.join(DIR, 'gadget_classes.pkl'), 'rb'))
    keys = sorted(cl, key=lambda c: (cl[c]['K'], str(cl[c]['dist']), str(c)))
    verdict = {}
    import gzip
    files = sorted(glob.glob(os.path.join(DIR, 'cpsat_s1_chunk*.jsonl'))) or [os.path.join(DIR, 'verdicts_pass1.jsonl.gz')]
    for f in files:
        for l in (gzip.open(f, 'rt') if f.endswith('.gz') else open(f)):
            r = json.loads(l)
            if r['st'] != 'UNKNOWN' or r['i'] not in verdict: verdict[r['i']] = r['st']
    return cl, keys, verdict

def gadget_s(W):
    return np.array([0] + [sum(int((W[k] * np.roll(W[k], -h)).sum()) for k in range(4)) for h in range(1, 16)])

if __name__ == "__main__":
    mode, ns = sys.argv[1], int(sys.argv[2]); secs = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
    cl, keys, verdict = load(); random.seed(7)
    pool = [i for i, v in verdict.items() if v == ('UNKNOWN' if mode == 'pass2' else 'INFEASIBLE')]
    samp = random.sample(pool, min(ns, len(pool)))
    out = open(os.path.join(DIR, f'crosscheck_{mode}.jsonl'), 'w'); tally = {}
    t0 = time.time()
    for idx in samp:
        info = cl[keys[idx]]; W = np.array(info['W']); E = np.array(info['eps'])
        if mode == 'y':
            st, wt = C.solve_stage1_y(W, E, tlim=secs, workers=1); rec = dict(i=idx, K=info['K'], st=st, t=round(wt, 2))
        elif mode == 'dp':
            h = D.class_certificate(W, E); st = 'cert' if h is not None else 'none'; rec = dict(i=idx, K=info['K'], cert_lag=h)
        else:
            st, sols, wt = C.solve_stage1(W, E, gadget_s(W), tlim=secs, workers=4); rec = dict(i=idx, K=info['K'], dist=info['dist'], st=st, t=round(wt, 2))
            if sols: rec['A'] = sols[0].tolist()
        tally[st] = tally.get(st, 0) + 1
        out.write(json.dumps(rec) + '\n'); out.flush()
    print(json.dumps(dict(summary=True, mode=mode, sample=len(samp), pool=len(pool), tally=tally, secs=round(time.time() - t0, 1))), flush=True)
