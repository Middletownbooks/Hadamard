"""Tables for report/skew_census_answer.md from the logs (screen, in-list ranking, exact gates, prime pass).
  python3 theory/skew_report_tables.py > logs/skew/report_tables.md"""
import json, os
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); LOG = os.path.join(HERE, '..', 'logs', 'skew')
rd = lambda f: [json.loads(l) for l in open(os.path.join(LOG, f))]
S = rd('screen.jsonl'); IR = rd('inlist_rank.jsonl'); GC = rd('gates_candidates.jsonl')
PI = rd('prime_inlist.jsonl') if os.path.exists(os.path.join(LOG, 'prime_inlist.jsonl')) else []

# ---- screen totals
vals = defaultdict(list)
for r in S: vals[r['v']].append(r)
st = Counter(); vst = Counter()
for v, rs in vals.items():
    ss = [r.get('status', '') for r in rs]
    if any(s.startswith('skipped') for s in ss): vst['skipped (Kumar lead / Schumer)'] += 1
    elif any(s.startswith('alive') for s in ss): vst['some slice alive'] += 1
    elif ss == ['no slice in range']: vst['no slice with 24-130 orbits'] += 1
    else: vst['every slice dead (exact gates)'] += 1
    for s in ss: st[s.split(' (')[0] if not s.startswith('dead') else ('dead ' + s.split()[1][0] if s.split()[1][0] in 'RL' else 'dead Z_d')] += 1
print('### Screen totals\n')
print('| values | count |\n|---|---:|')
for k, n in vst.most_common(): print(f'| {k} | {n} |')
print('\n| (slice, pattern) records | count |\n|---|---:|')
for k, n in st.most_common(): print(f'| {k} | {n} |')

# ---- gates of the candidates
gate = defaultdict(list)
for g in GC: gate[g['name']].append(g)
print('\n### Exact gates of the leading slices [M]\n')
for name, gs in gate.items():
    pats = len(gs); dead = sum(1 for g in gs if any(x[0] == 'dead' for x in g['gates'].values()))
    print(f'- **{name}**: {pats} row patterns, {dead} dead; ' + '; '.join(f"{g['pattern']}: " + ', '.join(f'{k} {x[0]}' + (f' ({x[1]})' if x[1] is not None else '') for k, x in g['gates'].items()) for g in gs[:2]) + (' ...' if pats > 2 else ''))

# ---- ranked table
rows = []
for r in IR:
    rows.append(dict(v=r['v'], shape=f"G={'x'.join('Z%d' % m for m in r['group'])}, H={r['gens']}", H=r['H'], r=r['r'], C=r['C'],
                     pattern=r['pattern'], E=r['log10E'], Ein=r['log10E_in'], cost=r['cost'], gates=r.get('gates')))
for r in PI:
    if r['v'] in (761, 1123, 613): continue
    rows.append(dict(v=r['v'], shape=f"G=Z{r['v']}, H=<{r['gen']}>", H=r['H'], r=2 * r['C'] + 1, C=r['C'], pattern=r['pattern'],
                     E=None, Ein=r['log10E_in'], cost=r['cost'], gates={}))
dead_names = {}
for name, gs in gate.items():
    dead_names[name] = all(any(x[0] == 'dead' for x in g['gates'].values()) for g in gs)
def verdict(r):
    g = r['gates'] or {}
    if r['v'] == 549 and (r['shape'].endswith('[[22]]') or r['shape'].endswith('[[16]]') or "('unip', [[1, 1], [0, 1]], 12)]" in r['shape']): return 'dead at Z3 [M] (all patterns)'
    if r['v'] == 361 and '[4, 0], [0, 7]' in r['shape']: return 'alive: both Z19 compressions pass [M]'
    if r['v'] == 325 and 'singer' in r['shape']: return 'alive: Z13 passes; Z5^2 undecided (cap)'
    if any(x[0] == 'dead' for x in g.values()): return 'dead'
    if g and all(x[0] == 'pass' for x in g.values()): return 'gates pass'
    return 'gates undecided/none' if g else 'prime: no Z_d gate'
seen = set(); out = []
for r in sorted(rows, key=lambda r: r['cost']):
    k = (r['v'], r['C'], tuple(r['pattern']))
    if k in seen: continue
    seen.add(k); out.append(r)
print('\n### Ranked screen table (top 30 by k-tree cost) [E]\n')
print('Cost = log10 of pair work per solution in a four-list k-tree with lists (min(N_A, 2^24), 5e7, 3e7, 3e7), from the analytic in-list estimate (skew_lib.inlist_estimate). On this scale the solved controls are 613: 10.2, 1123: 11.4, 761: 12.9.\n')
print('| # | v | slice | \\|H\\| | orbits | classes | rows | log10 E (slice) | log10 E in lists | cost | gates | next step |')
print('|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|---|')
for i, r in enumerate(out[:30], 1):
    vd = verdict(r)
    nxt = 'searched (§5)' if vd.startswith('alive') else ('none (dead)' if 'dead' in vd else ('k-tree beyond budget (cost >= 14)' if r['cost'] >= 14 else 'gate first'))
    print(f"| {i} | {r['v']} | {r['shape'][:60]} | {r['H']} | {r['r']} | {r['C']} | {r['pattern']} | {r['E'] if r['E'] is not None else '-'} | {r['Ein']} | {r['cost']} | {vd} | {nxt} |")
