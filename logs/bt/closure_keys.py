"""Alive structures of the families judged 'effectively closed' [E] (and the low-estimate rest), for the exact
closure pass (logs/bt/closure.sh): prints one comma-separated key list per priority tier.
  tier 1: the families judged effectively closed (list eff);
  tier 2: the other families with no alive structure at log2 E_Q >= 0, which the pass can still close exactly;
  tier 3: the remaining alive structures with log2 E_Q < 0 or no finite estimate.
Tiers 1 and 2 go smallest family first, one family at a time; tier 3 puts structures without a closure record first,
then the cheap families, then sorts by orbit count."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'theory'))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import bt_summary as S
D = S.Data()
eff = [(739, 'k4'), (947, 'k4'), (863, 'k8'), (719, 'k8'), (971, 'k6'), (599, 'k6'), (599, 'k8'), (669, 'prod12'),
       (751, 'k4one')]
done = set()
if os.path.exists('logs/bt/gates_run_closure.jsonl'):
    done = {json.loads(x)['key'] for x in open('logs/bt/gates_run_closure.jsonl')}
rr = lambda d: d['r'] if isinstance(d['r'], int) else sum(d['r'])
tier1, tier2, tier3 = [], [], []
nalive = {}
for (n, fam), L in sorted(D.groups().items()):
    if fam == 'k4split':
        continue
    alive = [d for d in L if D.status(d) == 'alive']
    nalive[(n, fam)] = len(alive)
    rich = any(d.get('log2E_Q') is not None and d['log2E_Q'] >= 0 for d in alive)
    for d in alive:
        e = d.get('log2E_Q')
        if (n, fam) in eff:
            tier1.append(d)
        elif not rich:
            tier2.append(d)
        elif e is None or e < 0:
            tier3.append(d)
key = lambda d: (d['key'] in done, eff.index((d['n'], d['fam'])) if (d['n'], d['fam']) in eff else 99, rr(d))
key2 = lambda d: (nalive[(d['n'], d['fam'])], d['n'], d['fam'], d['key'] in done, rr(d))
for t, k in ((tier1, key2), (tier2, key2), (tier3, key)):
    print(','.join(d['key'] for d in sorted(t, key=k)))
print(len(tier1), len(tier2), len(tier3), sum(d['key'] not in done for d in tier1), file=sys.stderr)
