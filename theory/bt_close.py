r"""Exact gates for the structures the census prefilter skipped [P screens, M Q-level gate].

theory/bt_census2.py skips structures with 4r - 6.5e < -45 (rule of thumb, an estimate).  To close a family exactly,
this re-evaluates every skipped record with the prefilter off (theory/bt_census2.evaluate(..., prefilter=None):
Theorem alpha, coset-Parseval, Q-level gate, Q-level stratified count).  One JSON line per record, restart-safe;
theory/bt_summary.py reads the verdicts.
usage: python3 theory/bt_close.py OUT.jsonl CENSUS.jsonl[,...] [--ns 719,863]
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census2 as C2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('census'); ap.add_argument('--ns', default=None)
    a = ap.parse_args()
    ns = set(map(int, a.ns.split(','))) if a.ns else None
    done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            done.add(json.loads(x)['key'])
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x)
            if not d.get('verdict', '').startswith('skipped') or d['key'] in done or (ns and d['n'] not in ns):
                continue
            base = {k: v for k, v in d.items() if k not in ('verdict', 'q_solutions', 'log2E_Q', 'log2E_whole', 'secs',
                                                            'rule', 'orbit_sizes', 'key', 'r', 'e')}
            base['r'] = d['r']
            t = time.time()
            try:
                res = C2.evaluate(base, prefilter=None)
            except Exception as ex:
                res = dict(base, verdict=f'error: {type(ex).__name__}: {ex}')
            res['key'] = d['key']
            res.pop('q_terms', None)
            res['secs'] = round(time.time() - t, 1)
            done.add(d['key'])
            with open(a.out, 'a') as fo:
                fo.write(json.dumps(res) + '\n')
            print(json.dumps({k: res.get(k) for k in ('n', 'fam', 'order', 'r', 'e', 'verdict', 'q_solutions',
                                                      'log2E_Q', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
