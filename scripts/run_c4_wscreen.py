"""Exact w-screen (theory/paley_wscreen.py one, with Lemma 6) on every c = 4 (b = 3) Paley-coset option of the
3000-4000 list; appends to logs/c4_wscreen_3000_4000.jsonl, skips options already there; per-option time cap."""
import sys, json, subprocess, os
sys.path.insert(0, 'scripts')
LOG = 'logs/c4_wscreen_3000_4000.jsonl'
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        r = json.loads(l); done.add((r['v'], r['p']))
L = [751, 789, 823, 839, 859, 863, 883, 907, 917, 919, 933, 947, 955, 971, 991]
opts = json.loads(subprocess.run(['python3', '-c', 'import sys,json; sys.argv=["x"]; sys.path.insert(0,"scripts"); import paley_bge3 as P; print(json.dumps([o for v in %r for o in P.options(v) if o["c"]==4]))' % L],
                                 capture_output=True, text=True).stdout)
opts.sort(key=lambda o: o['n'])
cap = int(sys.argv[1]) if len(sys.argv) > 1 else 3600
for o in opts:
    if (o['v'], o['p']) in done: continue
    n = o['n']; H = n // 4
    try:
        out = subprocess.run(['python3', 'theory/paley_wscreen.py', 'one', str(o['v']), '3', str(o['p']), str(n), str(H)],
                             capture_output=True, text=True, timeout=cap).stdout.strip()
        w = json.loads(out); alive = w.get('alive'); det = w.get('witness') if alive else w.get('KN')
    except subprocess.TimeoutExpired:
        alive = 'timeout'; det = None
    except Exception as e:
        alive = 'ERR'; det = repr(e)[:200]
    rec = dict(v=o['v'], p=o['p'], n=n, H=H, p_mod4=o['p'] % 4, wscreen_alive=alive, detail=det)
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    print(json.dumps(rec)[:200], flush=True)
