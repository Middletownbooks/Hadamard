"""Generalised theta_q criterion (report/paley_571_answer.md Theorem 4, generalised): an exact, slice-free screen.

Setting: b = 1 core on Z_m, m = v - 1 (rows (+-2, 0, 0, 0)), invariant under a multiplier group M <= (Z/m)^*.
For an odd prime q | m with M mapping ONTO (Z/q)^*, the Z_q-compression W_k(c) = sum_{i = c mod q} x_k(i) is constant on
c != 0: W_k = (W0, W1).  For a character chi of Z_m of order q,  x^_k(chi) = W0 - W1  (an integer), and the core condition
at chi reads  sum_k (W0_k - W1_k)^2 = 4v,  with row sums r_k = W0_k + (q - 1) W1_k.
W0, W1 are sums of m/q entries +-1: W = m/q (mod 2), |W| <= m/q.
If no (W0, W1)_k satisfy this, EVERY b = 1 slice whose multiplier group maps onto (Z/q)^* is empty (for that q).
Necessary conditions only (orbit granularity is ignored, so passing proves nothing).
usage: theta_screen.py v [v ...]       -> one JSON line per (v, q)"""
import sys, json, itertools as it

def odd_prime_factors(m):
    f, p = set(), 2
    while p * p <= m:
        while m % p == 0:
            if p > 2: f.add(p)
            m //= p
        p += 1
    if m > 2: f.add(m)
    return sorted(f)

def screen(v, q):
    m = v - 1; f = m // q
    vals = range(-f, f + 1, 2) if f % 2 == 0 else range(-f, f + 1, 2)
    vals = [w for w in range(-f, f + 1) if (w - f) % 2 == 0]
    by_r = {}
    for W0 in vals:
        for W1 in vals:
            r = W0 + (q - 1) * W1
            by_r.setdefault(r, set()).add((W0 - W1) ** 2)
    T = 4 * v
    s0 = by_r.get(2, set()) | by_r.get(-2, set()); s = by_r.get(0, set())
    if not s0 or not s: return False, None
    sl = sorted(s)
    for a in s0:
        for b in sl:
            if a + b > T: break
            for c in sl:
                if a + b + c > T: break
                if T - a - b - c in s: return True, (a, b, c, T - a - b - c)
    return False, None

if __name__ == "__main__":
    for v in map(int, sys.argv[1:]):
        for q in odd_prime_factors(v - 1):
            ok, wit = screen(v, q)
            print(json.dumps(dict(v=v, q=q, passes=ok, witness_squares=wit)), flush=True)
