r"""Cross-check for report/blockswap_answer.md section 2.3: every CP-SAT solution of a planted target
(logs/blockswap/cpsat_planted.jsonl, theory/blockswap_stage1_cpsat.py) must lie in the complete hit set of the
primary exhaustive method (theory/blockswap_stage1_exact.py) for the same target, because that enumeration is
complete [P].  Writes logs/blockswap/cpsat_membership.json.
Credit: ansatz by the OpenAI model (codex/order-3452-compression, 5b26dae); 863 set-up by OpenAI Codex;
bordered construction Kumar / EKLO Theorem 7.   usage: python3 theory/blockswap_cpsat_membership.py"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blockswap_stage1_exact as X
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
L = lambda f: [json.loads(l) for l in open(os.path.join(ROOT, 'logs', 'blockswap', f))]
plants = {}
for f in ('stage1_planted.jsonl', 'stage1_planted_real.jsonl'):
    for p in L(f):
        plants[(p.get('index'), p['seed'], p['trial'])] = p
out = []
for r in L('cpsat_planted.jsonl'):
    if 'y' not in r:
        out.append(dict(trial=r['trial'], index=r.get('index'), status=r['status'], member=None))
        continue
    # the plant this record came from: same (trial, index) and the same Q row sums
    match = [p for p in plants.values() if p['trial'] == r['trial'] and p.get('index') == r.get('index')
             and [sum(y[q * 5:(q + 1) * 5]) for y in r['y'][:3] for q in range(4)] ==
             [v for z in p['z'][:3] for v in z]]
    if len(match) != 1:
        out.append(dict(trial=r['trial'], index=r.get('index'), status=r['status'], member=None, note='plant not identified'))
        continue
    p = match[0]
    zA, zB, zC = [tuple(z) for z in p['z'][:3]]
    d = X.decide(zA, zB, zC, tgt=p['target'], etot=p['target'][0] - 80)
    hits = {(tuple(s['A']), tuple(s['B']), tuple(s['C'])) for s in d['solutions']}
    mem = (tuple(r['y'][0]), tuple(r['y'][1]), tuple(r['y'][2])) in hits
    out.append(dict(trial=r['trial'], index=r.get('index'), status=r['status'], cpsat_secs=r['secs'],
                    recheck=r.get('recheck'), member=mem, n_hits=len(hits)))
    print(out[-1], flush=True)
json.dump(out, open(os.path.join(ROOT, 'logs', 'blockswap', 'cpsat_membership.json'), 'w'), indent=1)
print('members:', sum(1 for o in out if o['member']), 'of', len(out))
