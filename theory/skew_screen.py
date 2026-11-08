"""Step 2: screen of skew multiplier slices at the open census values (prompt §4 step 2).

For each open v (minus Kumar's leads and the two values Schumer may have settled), every slice of skew_enum with
RMIN..RMAX orbits per block:
  gate R  [P/M]: skew row sum 1; 4v - 1 = b^2 + c^2 + d^2 with b, c, d attainable by H-invariant blocks (exact DP);
  gate Zd [M]:   exact compression to Z_d for every divisor 1 < d < v of the cyclic factor (d <= DMAX), with the lag-0
                 energy 4v and every lag class of <H_d, -1> on Z_d (2+2 join; 'undecided' at a cap, never a zero);
                 for d = v/q (q a prime factor) this is the stage-1 level of a staged slice;
  estimate [E]:  lattice-Gaussian log10 E of ordered quadruples with positive row sums (skew_lib.gauss_estimate);
  cost [E]:      log10 W, W = min over pairings of N_X N_Y + N_Z N_T (pair work of one exhaustive 2+2 join).
Output: logs/skew/screen.jsonl (appended; resumes by skipping values already done).
  python3 theory/skew_screen.py WORKERS [v ...]
"""
import json, math, os, signal, sys, time, traceback
from multiprocessing import Pool
import numpy as np
from sympy import divisors, isprime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skew_lib import *
from skew_enum import all_slices

RMIN, RMAX, DMAX, TBUDGET, RGATE, GATESECS = 24, 130, 60, 150, 65, 30
class GateTimeout(Exception): pass
def _alarm(sig, frm): raise GateTimeout()
signal.signal(signal.SIGALRM, _alarm)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'logs', 'skew', 'screen.jsonl')
SKIP = {1951: 'Kumar lead (spin, §14)', 2269: 'Kumar lead (spin, §14)', 2437: 'Kumar lead (spin, §14)',
        661: 'Kumar running (§115)', 883: 'Kumar running (§115)', 1661: 'Schumer 6644?', 1963: 'Schumer 7852?'}

def screen_v(v):
    recs = []
    t0 = time.time()
    try:
        sl = all_slices(v, RMIN, RMAX)
    except Exception as e:
        return [dict(v=v, error='enum ' + repr(e))]
    for st in sl:
        base = dict(v=v, group=st.gdesc['group'], gens=str(st.gdesc['gens']), H=int(st.Hord), r=int(st.r), C=int(st.C),
                    skew_ok=bool(st.skew_ok))
        if not st.skew_ok:
            recs.append(dict(base, status='no skew block')); continue
        if not st.gateL:
            recs.append(dict(base, status='dead L', maxforced=int(st.forced.max()))); continue
        pats, rs = admissible_patterns(st, v)
        if not pats:
            recs.append(dict(base, status='dead R')); continue
        mods = st.G.mods; axis = len(mods) - 1
        ds = [d for d in divisors(mods[axis]) if 1 < d < v and (isprime(d) or d <= DMAX)] if mods[axis] > 1 else []
        for pat in pats:
            pat = tuple(int(x) for x in pat)
            rec = dict(base, pattern=[1] + list(pat))
            gates = {}; dead = False
            for d in ds:
                if st.r > RGATE: gates[str(d)] = ['not run (r > %d)' % RGATE, None]; continue
                if time.time() - t0 > TBUDGET: gates[str(d)] = ['skipped (time)', None]; continue
                try:
                    signal.alarm(GATESECS)
                    stt, n, sz = zp_gate(st, d, axis, pat, cap=200_000, joincap=3 * 10 ** 6)
                    signal.alarm(0)
                except GateTimeout:
                    stt, n, sz = 'undecided (gate time %ds)' % GATESECS, None, []
                except Exception as e:
                    signal.alarm(0); stt, n, sz = 'error', None, [repr(e)]
                gates[str(d)] = [stt, n]
                if stt == 'dead': dead = True; break
            rec['gates'] = gates
            nA = len(st.pairs)
            logN = [nA * math.log10(2)] + [math.log10(rs[r]) for r in pat]
            rec['logN'] = [round(x, 2) for x in logN]
            a, b, c, d_ = logN
            W = min(math.log10(10 ** (a + b - 30) + 10 ** (c + d_ - 30)) + 30,
                    math.log10(10 ** (a + c - 30) + 10 ** (b + d_ - 30)) + 30,
                    math.log10(10 ** (a + d_ - 30) + 10 ** (b + c - 30)) + 30)
            rec['log10W'] = round(W, 2)
            if dead:
                rec['status'] = 'dead Z' + [k for k, g in gates.items() if g[0] == 'dead'][0]
            else:
                rec['status'] = 'alive' if all(g[0] == 'pass' for g in gates.values()) else 'alive (some gate undecided)'
                try:
                    g = gauss_estimate(st, pat, nsamp=2000, rng=np.random.default_rng(v))
                    rec['log10E'] = None if g is None else round(float(g[0]), 2)
                except Exception as e:
                    rec['log10E'] = None; rec['est_error'] = repr(e)
            recs.append(rec)
    for r in recs: r['secs_v'] = round(time.time() - t0, 1)
    return recs or [dict(v=v, status='no slice in range')]

def main():
    W = int(sys.argv[1])
    vs = [int(x) for x in sys.argv[2:]] or json.load(open(os.path.join(HERE, '..', 'data', 'skew_open_census.json')))['open_v']
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT): done.add(json.loads(l)['v'])
    todo = [v for v in vs if v not in done and v not in SKIP]
    with open(OUT, 'a') as f:
        for v in vs:
            if v in SKIP and v not in done: f.write(json.dumps(dict(v=v, status='skipped: ' + SKIP[v])) + '\n')
    with Pool(W) as P, open(OUT, 'a') as f:
        for recs in P.imap_unordered(screen_v, todo):
            for r in recs: f.write(json.dumps(r) + '\n')
            f.flush()
            print(recs[0]['v'], len(recs), sum(r.get('status', '').startswith('alive') for r in recs), flush=True)

if __name__ == '__main__':
    main()
