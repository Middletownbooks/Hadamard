"""Recompute the open skew census list (prompt §1) against Kumar's repository.

Census: Kumar's search/data/skew_census_undecided.txt (Dokovic-Seberry, May 2025, 737 odd v < 2500).
Settled: every directory skew/<n> (and hadamard/3732, which README says is skew) with n = 4v, v odd and in the census;
plus Schumer's 919 (Zenodo 2026). Orders 8v, 16v do not settle v.
  python3 theory/skew_census_recompute.py /home/user/vatskishan03/hadamard
"""
import json, os, sys, re, subprocess
K = sys.argv[1] if len(sys.argv) > 1 else '/home/user/vatskishan03/hadamard'
cen = [int(l) for l in open(f'{K}/search/data/skew_census_undecided.txt') if l.strip() and not l.startswith('#')]
assert len(cen) == 737, len(cen)
orders = [int(d) for d in os.listdir(f'{K}/skew') if d.isdigit()] + [3732]
settled = {}
for n in orders:
    if n % 4 == 0 and (n // 4) % 2 == 1 and (n // 4) in cen:
        settled[n // 4] = 'Kumar skew/%d' % n
settled[919] = settled.get(919, '') + '; Schumer 2026 (Zenodo)'
not4v = sorted(n for n in orders if not (n % 4 == 0 and (n // 4) % 2 == 1))
openv = sorted(set(cen) - set(settled))
isprime = lambda n: n > 1 and all(n % p for p in range(2, int(n ** .5) + 1))
commit = subprocess.run(['git', '-C', K, 'log', '-1', '--format=%h %ci'], capture_output=True, text=True).stdout.strip()
rec = dict(kumar_commit=commit, census=len(cen), settled=len(settled), open=len(openv),
           open_prime=sum(map(isprime, openv)), open_le600=sum(v <= 600 for v in openv),
           orders_not_4v_odd=not4v, smallest_open=openv[:20])
print(json.dumps(rec, indent=1))
old = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'skew_open_census.json')))
print('agrees with data/skew_open_census.json:', old['open_v'] == openv, sorted(old['settled_census_v']) == sorted(settled))
