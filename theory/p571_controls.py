r"""Scan for positive controls of the hierarchical solver: v = 6q + 1, gcd(q, 6) = 1, slices M <= (Z/q)^* lifted to
Z_m (t = 1 mod 6, so Z6 is fixed pointwise).  Per slice: orbits R on Z_q, Gaussian log2E (p571_slicelib), and the
exact L1 compression test at every prime p | q onto which M projects surjectively ([P]: theta_p is rational; with
A = sum of X1 over t = 0 mod p, B = (r - A)/(p - 1) the value on each nonzero class, the L1 equation at theta_p is
sum_k (A_k - B_k)^2 = 4v, with (A_k, B_k) reachable for the block's row sum r_k).
usage: python3 theory/p571_controls.py vmin vmax [--rmax 12]"""
import sys, os, json, math, itertools as it
from math import gcd
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_slicelib import mgroup, orbits, estimate

def subgroups_mod(q):
    U = [t for t in range(1, q) if gcd(t, q) == 1]; seen = {}
    for r in range(0, 3):
        for gens in it.combinations(U[1:], r):
            G = frozenset(mgroup(q, gens) if gens else [1])
            if G not in seen: seen[G] = gens
    return seen

def lift(m, q, t):
    """t mod q -> T mod m with T = 1 mod 6"""
    for T in range(t, m, q):
        if T % 6 == 1: return T
    raise ValueError

def primes(n):
    out, p = [], 2
    while p * p <= n:
        while n % p == 0: out.append(p); n //= p
        p += 1
    if n > 1: out.append(n)
    return sorted(set(out))

def theta_p_test(q, G, v, p, rows=(2, 0, 0, 0)):
    """exact: does sum_k (A_k - B_k)^2 = 4v have a solution with (A_k) reachable?  (L1 level, X1 in {0,+-2,+-4,+-6})"""
    orbs = orbits(q, sorted(G)); sz = [len(o) for o in orbs]; c0 = [sum(1 for t in o if t % p == 0) for o in orbs]
    reach = {(0, 0)}
    for s, c in zip(sz, c0):
        reach = {(a + s * x, b + c * x) for (a, b) in reach for x in (-6, -4, -2, 0, 2, 4, 6)}
    vals = []
    for r in rows:
        vv = set()
        for (a, A) in reach:
            if a == r and (r - A) % (p - 1) == 0: vv.add(A - (r - A) // (p - 1))
        vals.append(sorted(vv))
    T = 4 * v; sq = [{w * w for w in V} for V in vals]
    for s0 in sq[0]:
        for s1 in sq[1]:
            for s2 in sq[2]:
                if T - s0 - s1 - s2 in sq[3]: return True
    return False

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--rmax": "12", "--N": "10000"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    for v in range(pos[0], pos[1] + 1):
        m = v - 1
        if m % 6 or gcd(m // 6, 6) != 1: continue
        q = m // 6
        for G, gens in subgroups_mod(q).items():
            R = len(orbits(q, sorted(G)))
            if R > int(opt["--rmax"]) or R < 3: continue
            tests = {}
            for p in primes(q):
                if len({t % p for t in G}) == p - 1: tests[p] = theta_p_test(q, G, v, p)
            if not all(tests.values()):
                print(json.dumps(dict(v=v, q=q, gens=list(gens), R=R, theta_tests=tests, dead=True)), flush=True); continue
            gl = [lift(m, q, g) for g in gens]
            e = estimate(m, gl, N=int(opt["--N"]))
            print(json.dumps(dict(v=v, q=q, gens=list(gens), lifted=gl, order=len(G), R=R, theta_tests=tests,
                                  log2E=e.get('log2E'), constraint_bits=e.get('constraint_bits'))), flush=True)
