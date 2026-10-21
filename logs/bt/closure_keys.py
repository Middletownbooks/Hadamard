"""Alive structures of the families judged 'effectively closed' [E] (and the low-estimate rest), for the exact
closure pass (logs/bt/closure.sh): prints one comma-separated key list per priority tier.  Within a tier: structures
without a closure record first, then the cheap families, then by orbit count."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'theory'))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import bt_summary as S
D = S.Data()
eff = [(739, 'k4'), (947, 'k4'), (863, 'k8'), (719, 'k8'), (971, 'k6'), (599, 'k6'), (599, 'k8'), (669, 'prod12')]
done = set()
if os.path.exists('logs/bt/gates_run_closure.jsonl'):
    done = {json.loads(x)['key'] for x in open('logs/bt/gates_run_closure.jsonl')}
rr = lambda d: d['r'] if isinstance(d['r'], int) else sum(d['r'])
tier1, tier2 = [], []
for (n, fam), L in sorted(D.groups().items()):
    if fam == 'k4split':
        continue
    for d in L:
        if D.status(d) != 'alive':
            continue
        e = d.get('log2E_Q')
        if (n, fam) in eff:
            tier1.append(d)
        elif e is None or e < 0:
            tier2.append(d)
key = lambda d: (d['key'] in done, eff.index((d['n'], d['fam'])) if (d['n'], d['fam']) in eff else 99, rr(d))
for t in (tier1, tier2):
    print(','.join(d['key'] for d in sorted(t, key=key)))
print(len(tier1), len(tier2), sum(d['key'] not in done for d in tier1), file=sys.stderr)
