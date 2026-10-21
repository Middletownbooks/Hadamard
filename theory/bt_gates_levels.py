r"""Exact gates at chosen levels, with a time limit instead of a cost cap [M].

The pair-count bound of theory/bt_cost.py can overestimate a gate by many orders of magnitude (at 669 a level bounded
by 2^37 pairs ran in 8 s).  This runner takes the levels by the order of N (every S-invariant N <= K of that order)
and the structures by key or by family, and runs theory/bt_gates_run.run_level under an alarm.  A level with no
solution closes the structure: no core [M].  One JSON line per (structure, level), restart-safe; output in the format
of bt_gates_run (read by theory/bt_summary.py).
usage: python3 theory/bt_gates_levels.py OUT.jsonl CENSUS.jsonl[,...] --levels 111,74 (--keys k1,k2|@FILE | --fam n:fam)
       [--tmax 900] [--maxlen 2e7]
"""
import sys, os, json, time, math, argparse, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census_run as R, bt_census2 as C2, bt_cost as CO, bt_gates_run as GR
import bt_summary as SU


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('census')
    ap.add_argument('--levels', required=True); ap.add_argument('--keys', default=None)
    ap.add_argument('--fam', default=None, help='n:fam, every alive structure of that order and family')
    ap.add_argument('--tmax', type=float, default=900); ap.add_argument('--maxlen', type=float, default=2e7)
    a = ap.parse_args()
    levels = [int(x) for x in a.levels.split(',')]
    recs = {}
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x)
            if 'key' in d:
                recs[d['key']] = d
    if a.keys:
        ks = open(a.keys[1:]).read().replace('\n', ',') if a.keys.startswith('@') else a.keys
        keys = [k for k in ks.split(',') if k]
    else:
        n, fam = a.fam.split(':')
        D = SU.Data()
        keys = [d['key'] for d in sorted(D.groups().get((int(n), fam), []), key=lambda d: -SU.e_of(d))
                if D.status(d) == 'alive']
    done = set(); dead = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            d = json.loads(x)
            done.add((d['key'], d.get('Nkey')))
            if d.get('verdict', '').startswith('no core'):
                dead.add(d['key'])
    signal.signal(signal.SIGALRM, _alarm)
    for key in keys:
        if key in dead or key not in recs:
            continue
        rec = recs[key]
        st = C2.build(rec) if 'fam' in rec else R.build(rec)
        for Nsz in levels:
            if key in dead:
                break
            for NA in [N for N in st.invariant_subgroups() if len(N) == Nsz and np.isin(N, st.Kel).all()]:
                nk = CO.nkey(NA)
                if (key, nk) in done:
                    continue
                t = time.time()
                signal.alarm(int(a.tmax))
                try:
                    res = GR.run_level(st, NA, cap_pairs=2.0 ** 40, maxlen=int(a.maxlen))
                    res['verdict'] = 'no core [M] (level empty)' if res['solutions_weighted'] == 0 else 'passes'
                except Timeout:
                    res = dict(verdict=f'capped (time limit {a.tmax:.0f} s)')
                except MemoryError as e:
                    res = dict(verdict=f'capped ({e})')
                finally:
                    signal.alarm(0)
                res.update(key=key, n=rec['n'], fam=rec.get('fam', 'k4'), order=rec['order'], r=rec['r'], e=rec['e'],
                           log2E_Q=rec.get('log2E_Q'), N=int(Nsz), Nkey=nk, log2_pairs=None,
                           secs=round(time.time() - t, 1))
                with open(a.out, 'a') as fo:
                    fo.write(json.dumps(res) + '\n')
                done.add((key, nk))
                print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'log2E_Q', 'N', 'verdict',
                                                          'solutions_weighted', 'secs')}), flush=True)
                if res['verdict'].startswith('no core'):
                    dead.add(key)
                    break


if __name__ == '__main__':
    main()
