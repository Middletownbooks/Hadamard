"""H(4892) neighbourhood: the screen's deficit re-read at REAL stage-1 solutions [E]
(open_screen_deficit.analyse(sl, Y, real=True)), for a slice key and a stage-1 list file.
usage: python3 theory/h4892_reread.py 'KEY' LISTFILE NSAMPLE [--out FILE]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import open_screen_deficit as D, ab_sweep as SW
a = sys.argv[1:]; key = json.loads(a[0]); L = json.load(open(a[1])); ns = int(a[2])
outp = a[a.index('--out') + 1] if '--out' in a else None
sols = L['solutions'] if 'solutions' in L else L['canonical']
def full(Y):
    if 'orbits' in L:
        out = []
        for y in Y:
            z = [0] * sum(len(o) for o in L['orbits'])
            for v, o in zip(y, L['orbits']):
                for g in o: z[g] = v
            out.append(z)
        return out
    return Y
sl = SW.build(D.cfg_from_key(key))
rng = np.random.default_rng(5); pick = rng.choice(len(sols), size=min(ns, len(sols)), replace=False)
for i in pick:
    r = D.analyse(sl, np.array(full(sols[int(i)])), real=True)
    best = min((s for s in r.get('splits', []) if s.get('cost') is not None), key=lambda s: s['cost'], default=None)
    rec = dict(key=key, y=int(i), H2=r.get('H2'), C2=r.get('C2'), E_per_y=r.get('E_per_y'), mahal2=r.get('mahal2'),
               best_split=best)
    print(json.dumps(rec), flush=True)
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")
