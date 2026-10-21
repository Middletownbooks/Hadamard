r"""Calibration of the census estimators and gates against exact core counts at small orders [M counts, E estimates].

For every k = 4 structure of the census enumeration (theory/bt_census_run.structures) at the small orders given, with
r <= rmax: the census evaluation (coset-Parseval, Theorem alpha, Q-level gate, Q-level stratified count log2E_Q read
at the true target, whole-structure Gaussian log2E_whole; theory/bt_census_run.evaluate) and the exact number of
ordered cores by the independent brute-force counter (theory/bt_brute_count.count_cores, shares no code with the
census).  Every structure the gates kill must have 0 cores; the estimates are compared with log2 of the counts.
usage: python3 theory/bt_calib_small.py OUT.jsonl --ns 7,11,15,... [--rmax 18] [--N 2000]
"""
import sys, os, json, math, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census_run as R, bt_brute_count as BF


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--ns', default='7,11,15,19,23,27,31,35,39,43,47')
    ap.add_argument('--rmin', type=int, default=2); ap.add_argument('--rmax', type=int, default=18)
    ap.add_argument('--N', type=int, default=2000)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            done.add(json.loads(x)['key'])
    for n in map(int, a.ns.split(',')):
        for rec in sorted(R.structures(n, a.rmin, a.rmax), key=lambda d: d['r']):
            key = R.struct_key(rec)
            if key in done:
                continue
            t = time.time()
            res = R.evaluate(rec, N=a.N)
            res.pop('q_terms', None)
            bf = BF.count_cores(rec, rmax=a.rmax + 2)
            res['N_exact'] = bf.get('cores_ordered')
            res['hashes_agree'] = bf.get('hashes_agree')
            res['log2N'] = round(math.log2(res['N_exact']), 2) if res['N_exact'] else None
            res['secs_total'] = round(time.time() - t, 1)
            with open(a.out, 'a') as f:
                f.write(json.dumps(res) + '\n')
            print(json.dumps({x: res.get(x) for x in ('n', 'Pmods', 'Omods', 'order', 'r', 'e', 'verdict', 'log2E_Q',
                                                      'log2E_whole', 'N_exact', 'log2N', 'secs_total')}), flush=True)


if __name__ == '__main__':
    main()
