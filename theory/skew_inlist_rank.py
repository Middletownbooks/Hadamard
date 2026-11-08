"""Rank alive screen slices by the k-tree cost model (skew_lib.inlist_estimate) [E].
For each alive (slice, pattern) of logs/skew/screen.jsonl with C <= CMAX and log10E >= EMIN (best pattern per slice):
lists L = (min(N_A, 2^24), 5e7, 3e7, 3e7); cost = log10(pairs per solution) = log10 W_lists - log10 E_in.
Output logs/skew/inlist_rank.jsonl.   python3 theory/skew_inlist_rank.py CMAX EMIN"""
import json, os, sys, math, ast
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skew_lib import *
from skew_enum import all_slices
CMAX, EMIN = int(sys.argv[1]), float(sys.argv[2])
HERE = os.path.dirname(os.path.abspath(__file__)); LOG = os.path.join(HERE, '..', 'logs', 'skew')
best = {}
for l in open(os.path.join(LOG, 'screen.jsonl')):
    r = json.loads(l)
    if not r.get('status', '').startswith('alive') or r.get('log10E') is None: continue
    if r['C'] > CMAX or r['log10E'] < EMIN: continue
    k = (r['v'], r['gens'], str(r['group']))
    if k not in best or r['log10E'] > best[k]['log10E']: best[k] = r
done = set()
OUT = os.path.join(LOG, 'inlist_rank.jsonl')
if os.path.exists(OUT):
    for l in open(OUT): x = json.loads(l); done.add((x['v'], x['gens'], str(x['group'])))
out = open(OUT, 'a')
cache = {}
for k, r in sorted(best.items(), key=lambda kv: kv[1]['v']):
    if k in done: continue
    v = r['v']
    if v not in cache: cache = {v: all_slices(v)}
    st = next((s for s in cache[v] if str(s.gdesc['gens']) == r['gens'] and s.gdesc['group'] == r['group']), None)
    if st is None: continue
    pat = tuple(r['pattern'][1:])
    LA = min(2.0 ** len(st.pairs), 2 ** 24)
    try:
        e, W, d = inlist_estimate(st, pat, (LA, 5e7, 3e7, 3e7), rng=np.random.default_rng(v))
    except Exception as ex:
        print('ERR', k, ex); continue
    rec = dict(v=v, group=r['group'], gens=r['gens'], H=r['H'], r=r['r'], C=r['C'], pattern=r['pattern'], gates=r.get('gates'),
               log10E=r['log10E'], log10E_in=round(float(e), 2), log10W_lists=round(W, 2), cost=round(W - float(e), 2), logN=[round(x, 1) for x in d['N']])
    out.write(json.dumps(rec) + '\n'); out.flush()
    print(v, r['gens'][:50], 'C', r['C'], 'E', r['log10E'], 'E_in', rec['log10E_in'], 'cost', rec['cost'], flush=True)
