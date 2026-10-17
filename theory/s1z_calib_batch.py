r"""Run theory/s1z_calib.py --exact over the calibration rungs of logs/s1z/calib/rungs.json (restart-safe: a rung whose
output already has a summary line is skipped).
usage: python3 theory/s1z_calib_batch.py [RUNGS.json] [--N 3000] [--part i/W]"""
import sys, os, json, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
a = sys.argv[1:]
src = a[0] if a and not a[0].startswith('--') else os.path.join(ROOT, 'logs/s1z/calib/rungs.json')
N = a[a.index('--N') + 1] if '--N' in a else '3000'
part, nW = (int(x) for x in (a[a.index('--part') + 1] if '--part' in a else '0/1').split('/'))
for i, r in enumerate(json.load(open(src))):
    if i % nW != part: continue
    out = os.path.join(ROOT, f"logs/s1z/calib/calib_{r['tag']}.jsonl")
    if os.path.exists(out) and any('"summary"' in l for l in open(out)): continue
    subprocess.run([sys.executable, os.path.join(ROOT, 'theory/s1z_calib.py'), json.dumps(r['cfg']), out, '--N', N, '--exact',
                    '--max-pairs', '4e8'], stdout=open(out.replace('.jsonl', '.out'), 'w'), stderr=subprocess.STDOUT)
