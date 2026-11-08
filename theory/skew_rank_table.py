"""Step 2 ranking: merge our screen (logs/skew/screen.jsonl) with the estimate-only run of Kumar's cp_hunt.py
(logs/skew/cp_hunt_est.jsonl, copied from $HADAMARD_WORK) and rank the candidates [E].

Ranking key (estimated cost on 4 cores / 15 GB):
- cp_hunt rows (cyclotomic products Z_m x GF(q) with K_m, searched by cp_lowe.py / kt4.c): logwork = log10 of the pair
  work AB + CD per solution inside low-energy lists of the listed sizes (Kumar's own estimate), at 4e7 pairs/s;
- screen rows (all odd-order multiplier slices, 24..130 orbits): log10 W - log10 E, with W the pair work of one full
  exhaustive 2+2 join (a full join is the only search the screen costs; it is pessimistic for large slices).
Writes logs/skew/rank_table.jsonl and prints the top N.
  python3 theory/skew_rank_table.py [N]
"""
import json, math, os, sys
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); LOG = os.path.join(HERE, '..', 'logs', 'skew')
N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
RATE = math.log10(4e7)
rows = []
for l in open(os.path.join(LOG, 'screen.jsonl')):
    r = json.loads(l)
    if not r.get('status', '').startswith('alive') or r.get('log10E') is None: continue
    cost = r['log10W'] - min(r['log10E'], 0) if r['log10E'] < 0 else r['log10W'] - r['log10E']
    rows.append(dict(src='screen', v=r['v'], shape=f"G={r['group']} H={r['gens']} |H|={r['H']} r={r['r']} C={r['C']}",
                     pattern=r['pattern'], gates=r.get('gates'), log10E=r['log10E'], log10W=r['log10W'],
                     cost_log10_pairs=round(cost, 2), cost_hours_4c=round(10 ** (cost - RATE) / 3600, 3)))
p = os.path.join(LOG, 'cp_hunt_est.jsonl')
if os.path.exists(p):
    for l in open(p):
        r = json.loads(l)
        if 'logwork' not in r: continue
        rows.append(dict(src='cp_hunt', v=r['v'], shape=f"Z{r['m']} x GF({r['q']}) e={r['e']} F={r['F']} K=<{r['k']}> |K|={r['K']} twist={r['twist']}",
                         pattern=r.get('rows'), gates=dict(stage1_classes=r.get('classes')), log10E=r.get('logE'),
                         log10E_lists=r.get('logE_lists'), L=r.get('L'), sizes=r.get('sizes'),
                         cost_log10_pairs=r['logwork'], cost_hours_4c=round(10 ** (r['logwork'] - RATE) / 3600, 3)))
rows.sort(key=lambda r: r['cost_log10_pairs'])
with open(os.path.join(LOG, 'rank_table.jsonl'), 'w') as f:
    for r in rows: f.write(json.dumps(r) + '\n')
best = {}
for r in rows:
    best.setdefault((r['v'], r['src']), r)
for r in sorted(best.values(), key=lambda r: r['cost_log10_pairs'])[:N]:
    print(r['v'], r['src'], r['shape'], r['pattern'], 'E', r['log10E'], 'E_lists', r.get('log10E_lists'), 'cost', r['cost_log10_pairs'], r['cost_hours_4c'], 'h')
