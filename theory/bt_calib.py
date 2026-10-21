r"""Calibration of the stratified counts on structures with exact core counts [M counts, E estimates].

At small orders n = 3 mod 4 (k = 4 family, G = P x O, K = O, Q = P of order 4), every structure of the census
enumeration (theory/bt_census_run.structures) with rmin <= r <= rmax is
1. evaluated as in the census (coset-Parseval, Q-level gate, Q-level stratified count log2E_Q, read at the true
   target; theory/bt_census_run.evaluate);
2. if it passes, given the level-A stratified count log2E_A at the coarsest invariant N_A < K whose lists fit
   (theory/bt_level1.run_level);
3. decided by the exhaustive chain (theory/bt_search.Run with the levels of theory/bt_decide.choose_levels, no
   border step), which lists every core.  The exact count N counts ordered quadruples of S-invariant +-1 sequences:
   sum over cores (sorted nonnegative row sums) of the pattern weight (bt_gates.pattern_weight).
One JSON line per structure (restart-safe).  The estimators are calibrated by comparing log2E_Q, log2E_A with log2 N.
usage: python3 theory/bt_calib.py OUT.jsonl --ns 63,87,135 [--rmin 10] [--rmax 34] [--hours 0.25] [--part i/W]
"""
import sys, os, json, math, time, argparse, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_census_run as R, bt_gates as BG, bt_level1 as L1, bt_decide as D, bt_search as BS


def exact_count(rec, workdir, hours):
    st = R.build(rec)
    lv = D.choose_levels(st)
    if lv is None:
        return dict(exact='no levels')
    spec = D.spec_from_record(rec, lv)
    d = os.path.join(workdir, rec['key'])
    if os.path.exists(d):
        shutil.rmtree(d)                      # a partial run is redone from scratch (counts must be complete)
    run = BS.Run(spec, d, seed=1, border=False)
    run.run(hours=hours)
    ev = [json.loads(x) for x in open(run.logf)]
    last = [e for e in ev if e['event'] in ('exhausted', 'stop')][-1]
    cores = []
    if os.path.exists(run.coref):
        cores = [json.loads(x) for x in open(run.coref)]
    N = sum(BG.pattern_weight(tuple(sorted(abs(v) for v in c['rows']))) for c in cores)
    out = dict(levels=lv, levelA_solutions=sum(e.get('solutions', 0) for e in ev if e['event'] == 'levelA'),
               cores_sorted=len(cores), N_exact=int(N), complete=last['event'] == 'exhausted',
               patterns_with_cores=sorted({tuple(sorted(abs(v) for v in c['rows'])) for c in cores}))
    shutil.rmtree(d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--ns', default='63,87,135')
    ap.add_argument('--rmin', type=int, default=10); ap.add_argument('--rmax', type=int, default=34)
    ap.add_argument('--hours', type=float, default=0.25); ap.add_argument('--part', default='0/1')
    ap.add_argument('--N', type=int, default=1500)
    a = ap.parse_args()
    pi, pw = map(int, a.part.split('/'))
    done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            done.add(json.loads(x)['key'])
    workdir = os.path.join(os.path.dirname(a.out) or '.', 'calib_work')
    os.makedirs(workdir, exist_ok=True)
    k = 0
    for n in map(int, a.ns.split(',')):
        recs = sorted(R.structures(n, a.rmin, a.rmax), key=lambda d: d['r'])
        for rec in recs:
            k += 1
            if (k - 1) % pw != pi:
                continue
            key = R.struct_key(rec)
            rec['key'] = key
            if key in done:
                continue
            t0 = time.time()
            try:
                res = R.evaluate(rec, N=a.N)
            except Exception as e:
                res = dict(rec, key=key, verdict=f'error: {type(e).__name__}: {e}')
            res['rule'] = round(4 * res.get('r', 0) - 6.5 * res.get('e', 0), 1)
            if res.get('verdict', '').startswith('passes'):
                st = R.build(rec)
                cands = L1.level_candidates(st)
                for NA in cands[::-1]:                       # coarsest first
                    try:
                        la = L1.run_level(st, NA, N=a.N, maxA=40)
                    except Exception as e:
                        la = dict(cap=f'{type(e).__name__}: {e}')
                    if 'cap' not in la:
                        res['levelA'] = dict(NA=la['NA'], solutions=la['levelA_solutions_weighted'],
                                             log2E_A=la['log2E_A'])
                        break
                try:
                    res.update(exact_count(rec, workdir, a.hours))
                except MemoryError as e:
                    res['exact'] = f'cap: {e}'
            else:
                res['N_exact'] = 0
                res['complete'] = True
            res['secs_total'] = round(time.time() - t0, 1)
            res.pop('q_terms', None)
            with open(a.out, 'a') as f:
                f.write(json.dumps(res) + '\n')
            print(json.dumps({x: res.get(x) for x in ('n', 'Pmods', 'Omods', 'order', 'r', 'e', 'rule', 'verdict',
                                                      'log2E_Q', 'levelA', 'N_exact', 'complete', 'secs_total')}),
                  flush=True)


if __name__ == '__main__':
    main()
