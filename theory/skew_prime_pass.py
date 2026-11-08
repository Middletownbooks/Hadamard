"""Open census primes v: every odd d | v - 1 with e = (v - 1)/d in [EMIN, EMAX]; H = the d-th roots (e-th powers).
Exact row-sum gate and the lattice-Gaussian class estimate of skew_lib (as theory/skew_screen.py, without Z_d gates,
which do not exist for prime v). Output logs/skew/prime_pass.jsonl.
  python3 theory/skew_prime_pass.py EMIN EMAX"""
import json, math, os, sys
import numpy as np
from sympy import isprime, primitive_root, divisors
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skew_lib import *
EMIN, EMAX = int(sys.argv[1]), int(sys.argv[2])
HERE = os.path.dirname(os.path.abspath(__file__))
openv = json.load(open(os.path.join(HERE, '..', 'data', 'skew_open_census.json')))['open_v']
SKIP = {1951, 2269, 2437, 661, 883}
out = open(os.path.join(HERE, '..', 'logs', 'skew', 'prime_pass.jsonl'), 'w')
for v in openv:
    if not isprime(v) or v in SKIP: continue
    g = primitive_root(v)
    for d in divisors(v - 1):
        if d % 2 == 0: continue
        e = (v - 1) // d
        if not (EMIN <= e <= EMAX): continue
        h = pow(g, e, v)
        G = Group([v]); st = Struct(G, [G.perm(aut_diag([v], [h]))])
        pats, rs = admissible_patterns(st, v)
        rec = dict(v=v, H=d, gen=h, e=e, C=st.C, patterns=[[1] + [int(x) for x in p] for p in pats])
        ests = []
        for p in pats:
            p = tuple(int(x) for x in p)
            ge = gauss_estimate(st, p, nsamp=3000, rng=np.random.default_rng(v))
            logN = [len(st.pairs) * math.log10(2)] + [math.log10(rs[r]) for r in p]
            ests.append(dict(pattern=[1] + list(p), log10E=None if ge is None else round(float(ge[0]), 2),
                             logN=[round(x, 2) for x in logN]))
        rec['est'] = ests
        out.write(json.dumps(rec) + '\n'); out.flush()
        if ests: print(v, d, e, [(x['pattern'], x['log10E']) for x in ests], flush=True)
