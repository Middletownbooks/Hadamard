"""Validation of tmat_sat: (1) known T(103) satisfies the CNF; (2) all-solution counts equal brute force
(numpy over 8^m labelings, same row-sum normalisation) on small (t, k)."""
import sys, json, itertools
import numpy as np
sys.path.insert(0, "scripts")
from tmat_sat import TEnc, TProp, orbit_data, rowsum_cases, solve, verify_gs_product
from pysat.solvers import Cadical195

H = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])
PAT = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])  # label & 3 -> GS pattern
PAT = np.vstack([PAT, -PAT])

def known():
    from cw_build import t103
    Q = H.T @ t103(); case = tuple(int(x) for x in Q.sum(1))
    e = TEnc(103, 3, case)
    assum = [(e.v[X][o] if Q[X, orb[0]] == -1 else -e.v[X][o]) for X in range(4) for o, orb in enumerate(e.orbits)]
    with Cadical195(bootstrap_with=e.cnf.clauses) as S:
        return {"known_T103_case": case, "cnf_sat_under_known": S.solve(assumptions=assum)}

def brute(t, k):
    K, orbits, oid, reps = orbit_data(t, k); m = len(orbits)
    cases = set(rowsum_cases(t)); cnt = 0
    sizes = np.array([len(o) for o in orbits])
    for L in itertools.product(range(8), repeat=m):
        Ql = PAT[list(L)].T                       # 4 x m
        rs = tuple(int(x) for x in Ql @ sizes)
        if rs not in cases: continue
        Q = Ql[:, oid]
        if verify_gs_product(Q): cnt += 1
    return cnt

if __name__ == "__main__":
    print(json.dumps(known()), flush=True)
    for t, k in [(7, 3), (13, 3), (19, 3), (31, 5), (37, 9), (43, 7)]:
        b = brute(t, k)
        s = solve(t, k, all_sols=True, log=False)
        print(json.dumps({"t": t, "k": k, "m": len(orbit_data(t, k)[1]), "brute": b, "sat": len(s["solutions"]),
                          "match": b == len(s["solutions"])}), flush=True)
