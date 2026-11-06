r"""Open-order screen, step 2 (supplement): the exact Z_2 screen of theory/gp_z2screen.py (lag 0 + rows + the order-2
character, sum_k A_k^2 = 4v) on every b1 survivor with U trivial [P, M].
gp_z2screen.screen(F, q, s) needs only the number F = |G'| of cells and one character of order 2 with halves of size
F/2, so it applies to any G' of even order with U = 1.  The slice is exactly empty when no composition survives.
Also the w-magnitude screen of theory/paley_wscreen.py (Theorem 5 / Lemma 6 of report/paley_design_answer.md,
the Paley-coset ansatz, which is the s = 2 slice with U = 1) on every b1 and k4 survivor with s = 2 and U = 1 [P, M].
usage: python3 theory/open_screen_z2gate.py OUT.jsonl [--secs 600]"""
import sys, os, json, glob, signal, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_z2screen import screen
import math

class TO(Exception): pass
def _al(*a): raise TO()

if __name__ == "__main__":
    out = sys.argv[1]; secs = int(sys.argv[sys.argv.index('--secs') + 1]) if '--secs' in sys.argv else 600
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    done = set()
    if os.path.exists(out): done = {json.loads(l)['key'] for l in open(out)}
    R = [json.loads(l) for f in sorted(glob.glob(os.path.join(root, 'logs/open_screen/screens_s*.jsonl'))) for l in open(f)]
    R = [r for r in R if r['fam'] in ('b1', 'k4') and not r['U'] and r.get('A') and r.get('B') is not False
         and (r['fam'] == 'b1' or r['s'] == 2)]
    R.sort(key=lambda r: r['sGeff'])
    signal.signal(signal.SIGALRM, _al)
    fo = open(out, 'a')
    for r in R:
        if r['key'] in done: continue
        F = math.prod(r['H1']); t0 = time.time(); rec = dict(key=r['key'], v=r['v'], F=F, q=r['q'], s=r['s'])
        signal.alarm(secs)
        try:
            if r['fam'] == 'b1':
                d = screen(F, r['q'], r['s'])
                surv = [x['comp'] for x in d['detail'] if x['rows_and_Z2_ok']]
                rec.update(comps=d['compositions'], survivors=len(surv), surv=surv[:6], Z2=bool(surv))
            if r['s'] == 2:
                import subprocess
                F4 = F * math.prod(json.loads(r['key'])[2] or [1])
                args = ['python3', os.path.join(root, 'theory/paley_wscreen.py'), 'one', str(r['v']),
                        '1' if r['fam'] == 'b1' else '3', str(r['q']), str(F4)] + ([] if r['fam'] == 'b1' else [str(F)])
                o = subprocess.run(args, capture_output=True, text=True, timeout=secs).stdout.strip().splitlines()
                w = json.loads(o[-1]); rec['wscreen'] = w.get('alive'); rec['wscreen_why'] = w.get('why')
            signal.alarm(0)
        except TO: rec['timeout'] = secs
        except Exception as e: signal.alarm(0); rec['error'] = repr(e)[:150]
        rec['secs'] = round(time.time() - t0, 1)
        fo.write(json.dumps(rec) + "\n"); fo.flush()
