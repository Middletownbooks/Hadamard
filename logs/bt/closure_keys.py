"""Alive structures of the families judged 'effectively closed' [E] (and the low-estimate rest), for the exact
closure pass (logs/bt/closure.sh): prints one comma-separated key list per priority tier."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'theory'))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import bt_summary as S
D = S.Data()
eff = {(599, 'k6'), (719, 'k8'), (739, 'k4'), (863, 'k8'), (947, 'k4'), (971, 'k6'), (599, 'k8'), (669, 'prod12')}
tier1, tier2 = [], []
for (n, fam), L in sorted(D.groups().items()):
    if fam == 'k4split':
        continue
    for d in L:
        if D.status(d) != 'alive':
            continue
        e = d.get('log2E_Q')
        if (n, fam) in eff:
            tier1.append(d['key'])
        elif e is None or e < 0:
            tier2.append(d['key'])
print(','.join(tier1))
print(','.join(tier2))
print(len(tier1), len(tier2), file=sys.stderr)
