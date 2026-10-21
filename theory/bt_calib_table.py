r"""Calibration tables for report/b3_push_answer.md section 6 from logs/bt/calib/calib_small.jsonl [M counts, E estimates].

Per bin of the Q-level stratified count log2 E_Q (theory/bt_census_run.evaluate, read at the true target): the number of
structures that pass every census gate, how many have cores (exact count, theory/bt_brute_count.py), and for those the
spread of log2 N - log2 E_Q.  Also: structures the gates kill that have cores (must be 0).
usage: python3 theory/bt_calib_table.py [FILE]
"""
import sys, os, json, math, collections

f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs', 'bt',
                                                       'calib', 'calib_small.jsonl')
rows = [json.loads(x) for x in open(f)]
killed = [d for d in rows if not d['verdict'].startswith('passes')]
bad = [d for d in killed if d['N_exact']]
print(f"structures: {len(rows)} at orders {sorted(set(d['n'] for d in rows))}; killed by the gates: {len(killed)} "
      f"(coset-Parseval {sum('coset' in d['verdict'] for d in killed)}, Theorem alpha "
      f"{sum('alpha' in d['verdict'] for d in killed)}, Q-level {sum('Q-level' in d['verdict'] for d in killed)}); "
      f"killed but with cores: {len(bad)}; hash disagreements: {sum(d.get('hashes_agree') is False for d in rows)}")
print()
P = [d for d in rows if d['verdict'].startswith('passes')]
print('| log2 E_Q | structures passing the gates | with cores | log2 N − log2 E_Q (median; min, max) over those with cores |')
print('|---|---|---|---|')
tab = collections.defaultdict(list)
for d in P:
    e = d.get('log2E_Q')
    k = 'none (−∞)' if e is None else (f'< 0' if e < 0 else f'[{int(e // 10) * 10}, {int(e // 10) * 10 + 10})')
    tab[k].append(d)
order = lambda k: -1e9 if k.startswith('none') else (-1e8 if k == '< 0' else int(k[1:].split(',')[0]))
for k in sorted(tab, key=order):
    L = tab[k]
    wc = [d for d in L if d['N_exact']]
    diffs = sorted(math.log2(d['N_exact']) - d['log2E_Q'] for d in wc if d.get('log2E_Q') is not None)
    s = '—' if not diffs else f"{diffs[len(diffs) // 2]:+.2f}; {diffs[0]:+.2f}, {diffs[-1]:+.2f}"
    print(f"| {k} | {len(L)} | {len(wc)} | {s} |")
print()
S1 = [d for d in P if d['order'] == 1]
print('S = 1 (no symmetry): ' + ', '.join(f"n = {d['n']} {d['Pmods']}x{d['Omods']}: log2 E_Q {d['log2E_Q']}, log2 N "
                                          f"{d['log2N']}" for d in S1))
