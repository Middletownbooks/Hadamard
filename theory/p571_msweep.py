r"""Sweep ALL multiplier subgroups M <= (Z/m)^* for b = 1 slices on Z_m (571: m = 570, (Z/570)^* = Z2 x Z4 x Z18).

The repo sweep (scripts/slice_sweep.py) only covered "Paley(p) + <mu>": M = squares mod p on the Z_p factor times a
cyclic <mu> on G' (for p = 3 the squares are trivial, so only mu = 1 mod 3).  Here every subgroup is listed (cyclic or
not, acting on any factor, including t = -1 mod 3), with its orbit count (+-1 variables per block), the number of
equations (orbits of <M, -1> on lags != 0), the residues of M mod each prime factor, and a Gaussian log2E
(theory/p571_slicelib.estimate) for subgroups whose naive bit balance is not hopeless.
usage: python3 theory/p571_msweep.py m [--omin 8] [--omax 90] [--N 20000] [--list]"""
import sys, os, json, math, itertools as it
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_slicelib import mgroup, orbits, lag_reps, estimate

def units(m): return [t for t in range(1, m) if gcd(t, m) == 1]

def all_subgroups(m, maxgens=3):
    U = units(m); seen = {}
    for r in range(0, maxgens + 1):
        for gens in it.combinations(U[1:], r):
            G = frozenset(mgroup(m, gens) if gens else [1])
            if G not in seen: seen[G] = gens
    return seen

def factor(m):
    f, p = [], 2
    while p * p <= m:
        while m % p == 0: f.append(p); m //= p
        p += 1
    if m > 1: f.append(m)
    return sorted(set(f))

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--omin": "8", "--omax": "90", "--N": "20000"}; pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(int(a[i])); i += 1
    m = pos[0]; P = factor(m)
    subs = all_subgroups(m)
    rows = []
    for G, gens in subs.items():
        orbs = orbits(m, sorted(G)); reps, _ = lag_reps(m, sorted(G))
        res = {str(p): sorted({t % p for t in G}) for p in P if p > 2}
        rows.append(dict(order=len(G), gens=list(gens), orbits=len(orbs), equations=len(reps),
                         minus1=((m - 1) in G), res=res, naive=4 * len(orbs) - 5.4 * len(reps)))
    rows.sort(key=lambda r: (-r['naive'], r['orbits']))
    print(json.dumps(dict(m=m, subgroups=len(rows))), flush=True)
    for r in rows:
        if not (int(opt["--omin"]) <= r['orbits'] <= int(opt["--omax"])): continue
        if '--list' in flags or r['naive'] < -40:
            print(json.dumps(r), flush=True); continue
        try: e = estimate(m, r['gens'], N=int(opt["--N"]))
        except Exception as ex: e = dict(error=repr(ex))
        print(json.dumps({**r, **{k: e[k] for k in e if k in ('log2E', 'constraint_bits', 'rowsum_bits', 'eqs', 'feasible')}}), flush=True)
