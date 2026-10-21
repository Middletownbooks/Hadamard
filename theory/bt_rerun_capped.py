r"""Re-run the census evaluation of structures whose Q-level gate hit a cap ("Q-level cap (...)"), with larger caps [M].

The product family at 669 (theory/bt_census2.py --fam prod12) left 17 structures undecided: 9 ran out of memory under
the earlier 6 GB address-space limit (the 2 + 2 join output array) and 8 exceeded the two-step Q-level list cap of
2*10^6 (theory/bt_chain.q_level_list_two_step).  This driver raises both caps (--maxlen, --maxout) and evaluates the
structures again with the unchanged census code (bt_census2.evaluate); each result is appended to the census file it
came from, where it supersedes the capped record (theory/bt_summary.py keeps the last record per key).
usage: python3 theory/bt_rerun_capped.py CENSUS.jsonl [--maxlen 1e7] [--maxout 3e7] [--N 1200]
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_chain as C, bt_census2 as C2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('census'); ap.add_argument('--maxlen', type=float, default=1e7)
    ap.add_argument('--maxout', type=float, default=3e7); ap.add_argument('--N', type=int, default=1200)
    a = ap.parse_args()
    last = {}
    for x in open(a.census):
        d = json.loads(x)
        last[d['key']] = d
    todo = [d for d in last.values() if d.get('verdict', '').startswith('Q-level cap')]
    print(f'{len(todo)} capped structures', flush=True)
    f2s = C.q_level_list_two_step
    C.q_level_list_two_step = lambda st, qQ, i, row, chars, cvals=None, maxlen=int(a.maxlen): f2s(st, qQ, i, row, chars, cvals, maxlen)
    j4 = C.join4
    C.join4 = lambda P, T, maxout=int(a.maxout), chunk=4 * 10 ** 7: j4(P, T, maxout, chunk)
    for d in sorted(todo, key=lambda d: d['r']):
        rec = {k: d[k] for k in d if k not in ('verdict', 'secs', 'key', 'r', 'e', 'orbit_sizes', 'rule', 'cp',
                                                'inv_subgroups', 'q_patterns', 'q_solutions')}
        t0 = time.time()
        try:
            res = C2.evaluate(rec, N=a.N)
        except MemoryError as e:
            res = dict(d, verdict=f'Q-level cap ({e})')
        res['rerun'] = dict(maxlen=a.maxlen, maxout=a.maxout)
        res['secs'] = round(time.time() - t0, 1)
        with open(a.census, 'a') as fo:
            fo.write(json.dumps(res) + '\n')
        print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'e', 'verdict', 'q_solutions', 'log2E_Q',
                                                  'secs')}), flush=True)


if __name__ == '__main__':
    main()
