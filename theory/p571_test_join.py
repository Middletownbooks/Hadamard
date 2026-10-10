"""Unit test of the hierarchical solver's 2+2 join (theory/p571_hsolve.join_all): planted quadruples must be found.
usage: python3 theory/p571_test_join.py"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import join_all

if __name__ == "__main__":
    rng = np.random.default_rng(0); ok = True
    for L, D, tb in ((1000, 8, 16), (20000, 8, 20), (40000, 8, 22)):
        sigs = [rng.integers(-50, 50, (L, D)) for _ in range(4)]
        pick = rng.integers(0, L, 4)
        T = sum(s[p] for s, p in zip(sigs, pick))
        t = time.time(); hits, st = join_all(sigs, T, seed=3, tbits=tb)
        found = any((h == pick).all() for h in hits)
        ok &= found
        print(json.dumps(dict(L=L, D=D, planted=pick.tolist(), found=bool(found), n_solutions=len(hits),
                              secs=round(time.time() - t, 1), **st)), flush=True)
    print(json.dumps(dict(all_planted_found=bool(ok))))
