r"""Gate 1 summary for report/prompt_blockswap_863.md: compares, tuple by tuple, the primary exhaustive method
(theory/blockswap_stage1_exact.py, logs/blockswap/stage1_exact_w*.jsonl) with the second exhaustive method
(theory/blockswap_stage1_check2.py, logs/blockswap/check2_w*.jsonl), collects the planted tests, and runs the primary
method also on the 8 Gate 0 tuples that fail the joint energy gate (consistency check).
Credit: ansatz by the OpenAI model (codex/order-3452-compression, 5b26dae); 863 set-up by OpenAI Codex; bordered
construction Kumar / EKLO Theorem 7.   usage: python3 theory/blockswap_stage1_summary.py"""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
L = lambda pat: [json.loads(l) for f in sorted(glob.glob(os.path.join(ROOT, 'logs', 'blockswap', pat))) for l in open(f)]
P = {r['index']: r for r in L('stage1_exact_w*.jsonl')}
S = {r['index']: r for r in L('check2_w*.jsonl')}
rows, agree = [], True
for i in range(40):
    p, s = P.get(i), S.get(i)
    ok = (p is not None and s is not None and p['sizes'] == s['sizes'] and p['mins'] == s['mins']
          and p['hits'] == 0 and s['nsol'] == 0)
    agree &= ok
    rows.append(dict(index=i, z=p['z'][:3] if p else None, sizes=p['sizes'] if p else None, mins=p['mins'] if p else None,
                     primary=p['status'] if p else None, primary_key_matches=p.get('key_matches') if p else None,
                     second=s['status'] if s else None, sizes_equal=bool(p and s and p['sizes'] == s['sizes']),
                     secs_primary=p['secs'] if p else None, secs_second=s['secs'] if s else None))
pl = L('stage1_planted.jsonl') + L('stage1_planted_real.jsonl')
c2 = L('check2_planted.jsonl')
import blockswap_stage1_exact as X
d = json.load(open(os.path.join(ROOT, 'logs', 'blockswap', 'gate0_m5.json')))
non = []
for t in d['tuples']:
    if not t['survives']:
        r = X.decide(*[tuple(z) for z in t['z'][:3]])
        non.append(dict(z=t['z'][:3], status=r['status'], mins=r['mins'], sizes=r['sizes'], hits=r['hits']))
out = dict(n_survivors=40, all_empty_both_methods=agree, decided_primary=len(P), decided_second=len(S),
           total_secs_primary=round(sum(r['secs'] for r in P.values()), 1),
           total_secs_second=round(sum(r['secs'] for r in S.values()), 1),
           max_sizes=[max(r['sizes'][k] for r in P.values()) for k in range(3)],
           planted_primary=dict(n=len(pl), recovered=sum(r['planted_recovered'] for r in pl)),
           planted_second=dict(n=len(c2), recovered=sum(r['planted_recovered'] for r in c2),
                               hit_counts_equal=all(r['nsol'] == r['primary_hits'] for r in c2)),
           non_survivors_primary=non, tuples=rows)
json.dump(out, open(os.path.join(ROOT, 'logs', 'blockswap', 'stage1_summary.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('tuples', 'non_survivors_primary')}))
print([(n['status'], n['mins']) for n in non])
