"""Regression test for the b_setup fix in theory/open_screen_stage1_exact.py (3-entry U keys [axis, order, gen]):
every key of the screen that the old b_setup accepted gives an identical set-up (orbits, letters, T1, lag map) with
the new code; the 5892 key now sets up, and its orbit partition agrees with ab_sweep.build (method A's set-up).
The old version is read from git (commit d9595f9)."""
import sys, os, json, glob, subprocess, types
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import open_screen_stage1_exact as NEW
src = subprocess.run(['git', 'show', 'd9595f9:theory/open_screen_stage1_exact.py'], capture_output=True, text=True).stdout
OLD = types.ModuleType('old'); OLD.__file__ = NEW.__file__; exec(compile(src, 'old', 'exec'), OLD.__dict__)
keys = set()
for f in glob.glob('logs/open_screen/screens_s*.jsonl'):
    for l in open(f):
        r = json.loads(l)
        if r.get('r', 999) <= 60: keys.add(r['key'])
keys = sorted(keys)
same = diff = olderr = 0
for k in keys:
    K = json.loads(k)
    try: a = OLD.b_setup(K)
    except Exception: olderr += 1; continue
    b = NEW.b_setup(K)
    if all(a[x] == b[x] for x in ('orbits', 'letters', 'T1', 'lagmap', 'axes')): same += 1
    else: diff += 1; print('DIFF', k)
K = ["b1", 1473, [], [32, 2], 23, 2, [[0, 8, 3]]]
st = NEW.b_setup(K)
import open_screen_deficit as D, ab_sweep as SW, numpy as np
sl = SW.build(D.cfg_from_key(K))
agree = sorted(sorted(o) for o in st['orbits']) == sorted(sorted(np.nonzero(sl.borb == j)[0].tolist()) for j in range(len(sl.breps)))
out = dict(keys_tested=len(keys), identical=same, different=diff, old_errors=olderr, key5892_orbits=len(st['orbits']),
           key5892_agrees_with_ab_sweep=bool(agree))
print(json.dumps(out)); assert diff == 0 and agree
