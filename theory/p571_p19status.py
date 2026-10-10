r"""Status of the <61> = Paley(19) stage-1 sweep: per gadget distribution, the CP-SAT verdicts from all pass logs
(theory/logs/p571/p19/cpsat_s1_*.jsonl; a later INFEASIBLE/OPTIMAL verdict overrides an earlier UNKNOWN), and the list
of unresolved class indices (written to theory/logs/p571/p19/unresolved.txt).
usage: python3 theory/p571_p19status.py"""
import json, collections, pickle, glob, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'theory/logs/p571/p19')
cl = pickle.load(open(os.path.join(D, 'gadget_classes.pkl'), 'rb'))
keys = sorted(cl, key=lambda c: (cl[c]['K'], str(cl[c]['dist']), str(c)))
verdict = {}
order = sorted(glob.glob(os.path.join(D, 'cpsat_s1_chunk*.jsonl'))) + sorted(glob.glob(os.path.join(D, 'cpsat_s1_pass*.jsonl')))
import gzip
if not glob.glob(os.path.join(D, 'cpsat_s1_chunk*.jsonl')): order = [os.path.join(D, 'verdicts_pass1.jsonl.gz')] + order
for f in order:
    for l in (gzip.open(f, 'rt') if f.endswith('.gz') else open(f)):
        r = json.loads(l)
        if r['st'] != 'UNKNOWN' or r['i'] not in verdict: verdict[r['i']] = r['st']
byd = collections.defaultdict(collections.Counter); tot = collections.Counter()
for idx, k in enumerate(keys):
    v = verdict.get(idx, 'TODO'); byd[(cl[k]['K'], tuple(cl[k]['dist']))][v] += 1; tot[(cl[k]['K'], v)] += 1
for d, c in sorted(byd.items()): print(d, dict(c))
print(json.dumps({f'K{K}:{v}': n for (K, v), n in sorted(tot.items())}))
unres = [idx for idx in range(len(keys)) if verdict.get(idx) == 'UNKNOWN']
open(os.path.join(D, 'unresolved.txt'), 'w').write('\n'.join(map(str, unres)) + '\n')
print('unresolved', len(unres), 'todo', sum(1 for i in range(len(keys)) if i not in verdict))
