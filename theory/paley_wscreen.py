r"""Exact "w-magnitude" screen for the Paley-coset ansatz, any odd prime p >= 5 (report/paley_design_answer.md section 4).
Stage 1 only involves the Z_p-coset sums y, so it holds verbatim for p = 1 mod 4 as well.

Notation: x_k(g', .) = alpha delta + u chi + w eps on every Z_p-coset; y = alpha + (p-1) w; K = #constant cosets.
Stage 1 (exact, verified in theory/paley_algebra.py):  sum_k PAF(y_k)(h') = T_J(h') with
    T_J = 4p off H', -4bp on H' \ 0 (b >= 3);   T_J = -4p off 0 (b = 1);   T_J = 0 off 0 (b = 0),
and PAF(y) = PAF(alpha) + (p-1) [CC(alpha,w) + CC(w,alpha)] + (p-1)^2 PAF(w).  Hence for every h' != 0
    (p-1)^2 s(h') + a(h') + (p-1) c(h') = T_J(h'),       s = sum_k PAF(w_k)(h'),
    a = sum_k PAF(alpha_k)(h') in [-4n', 4n'], a = 0 mod 4   (PAF of a full +-1 sequence = n' mod 4),
    c = sum_k [CC(alpha_k,w_k) + CC(w_k,alpha_k)](h') in [-2K, 2K], c even.
This confines s(h') to a small set S(h').  Summing PAF(w) over G' gives the exact identity
    K + sum_{h' != 0} s(h') = sum_k (sum_g w_k(g))^2 =: sum_k t_k^2,
while the row sums r_k = (sum alpha_k) + (p-1) t_k with |sum alpha_k| <= n', sum alpha_k = n' (mod 2), and the row-sum
profile (+-2,0,0,0) for b = 1 or sum r_k^2 = 4v for b != 1, confine the t_k.  The (K, N_-, N_+) come from the h' = 0
equation sum y^2 = T_J(0).  An option is killed if no (K, t_1..t_4, s) satisfies all of this.
usage: python3 theory/paley_wscreen.py b1                 (all b = 1 options of the open orders, p = 3 mod 4)
       python3 theory/paley_wscreen.py bge3               (b >= 3 options, p = 3 mod 4)
       python3 theory/paley_wscreen.py one v b p n' [|H'|]"""
import sys, os, json, itertools as it, math
from functools import lru_cache

OPEN = [523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739]

def factor(n):
    f, q = {}, 2
    while q * q <= n:
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def KN_options(v, b, p, n):
    """(K, N-, N+) with N_tw + (p-2)^2 N- + p^2 N+ = T_J(0) = 4m + (p-1) Ts(0), N_tw = 4n - K"""
    m = n * p; TJ0 = 4 * m + (p - 1) * (-4 * b)
    out = []
    for Nm in range(0, 4 * n + 1):
        for Np in range(0, 4 * n + 1 - Nm):
            K = Nm + Np
            if (4 * n - K) + (p - 2) ** 2 * Nm + p ** 2 * Np == TJ0: out.append((K, Nm, Np))
    return out

def S_allowed(TJ, p, n, K):
    out = []
    smax = (abs(TJ) + 4 * n + 2 * (p - 1) * K) // (p - 1) ** 2 + 1
    for s in range(-smax, smax + 1):
        rest = TJ - (p - 1) ** 2 * s                       # = a + (p-1) c
        ok = False
        for c in range(-2 * K, 2 * K + 1, 2):
            a = rest - (p - 1) * c
            if abs(a) <= 4 * n and a % 4 == 0: ok = True; break
        if ok: out.append(s)
    return out

def screen(v, b, p, n, hsz=1, verbose=False):
    """hsz = |H'| for b >= 3 (H' of index c = b + 1 in G'); lag classes: n - 1 nonzero lags, of which hsz - 1 in H'"""
    res = {"v": v, "b": b, "p": p, "n'": n, "H'": hsz}
    KNs = KN_options(v, b, p, n)
    if not KNs: res.update(alive=False, why="no (K, N-, N+) from h'=0"); return res
    if b == 1: profiles = [(2, 0, 0, 0)]
    else:
        R = [r for r in range(0, 2 * int(math.isqrt(v)) + 2) if (r - n * p) % 2 == 0]   # r = m mod 2
        profiles = [q for q in it.combinations_with_replacement(R, 4) if sum(x * x for x in q) == 4 * v]
    alive = []
    for (K, Nm, Np) in KNs:
        # allowed s on lags off H' (n - hsz of them) and on H' \ 0 (hsz - 1)
        if b == 1: Toff, TH = -4 * p, None
        elif b == 0: Toff, TH = 0, None
        else: Toff, TH = 4 * p, -4 * b * p
        Soff = S_allowed(Toff, p, n, K)
        SH = S_allowed(TH, p, n, K) if TH is not None else []
        noff = n - hsz; nH = hsz - 1
        if not Soff or (nH and not SH): continue
        lo = noff * min(Soff) + (nH * min(SH) if nH else 0); hi = noff * max(Soff) + (nH * max(SH) if nH else 0)
        # achievable sums: every value between lo and hi in steps of the gcd of the allowed differences is NOT
        # guaranteed; compute exactly by DP over lag classes (pairs h, -h have equal s; small sets -> interval DP)
        @lru_cache(None)
        def sums(cnt, which):
            S = Soff if which == 0 else SH
            if S == list(range(min(S), max(S) + 1)):                # contiguous: all integers in [cnt*min, cnt*max]
                return frozenset(range(cnt * min(S), cnt * max(S) + 1))
            acc = {0}
            for _ in range(cnt): acc = {x + y for x in acc for y in S}
            return frozenset(acc)
        # lag multiplicity: s(h') = s(-h'); a self-inverse h' appears once, other h' in pairs. Use the loose version
        # (every lag independent) which is a superset -> still an exact necessary condition.
        tot = {x + y for x in sums(noff, 0) for y in (sums(nH, 1) if nH else {0})}
        for prof in profiles:
            for signs in it.product(*[(1, -1) if r else (1,) for r in prof]):
                rs = [s_ * r for s_, r in zip(signs, prof)]
                # t_k options per block
                topts = []
                for r in rs:
                    ts = [t for t in range(-K, K + 1) if abs(r - (p - 1) * t) <= n and (r - (p - 1) * t - n) % 2 == 0]
                    topts.append(ts)
                for ts in it.product(*topts):
                    if sum(abs(t) for t in ts) > K: continue
                    if (sum(abs(t) for t in ts) - K) % 2: continue          # sum |supp w_k| = K, t_k = |supp| mod 2
                    need = sum(t * t for t in ts) - K
                    if need in tot:
                        alive.append(dict(K=K, Nminus=Nm, Nplus=Np, profile=rs, t=list(ts), Soff=Soff, SH=SH, sum_s=need))
                        if not verbose: break
                if alive and not verbose: break
            if alive and not verbose: break
    res.update(alive=bool(alive), witness=alive[0] if alive else None, KN=KNs)
    return res

def b1_options():
    out = []
    for v in OPEN:
        m = v - 1
        for p, e in factor(m).items():
            if e != 1 or p < 7: continue
            out.append((v, 1, p, m // p, 1))
    return out

def bge3_options():
    out = []
    for v in OPEN:
        for c in range(4, v, 2):
            if (v + 1) % c: continue
            h = (v + 1) // c - 1; b = c - 1; m = c * h
            if h < 3: continue
            for p, e in factor(h).items():
                if p < 7 or factor(m)[p] != 1 or c % p == 0: continue
                n = m // p
                if n % c: continue
                out.append((v, b, p, n, n // c))
    return out

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd in ("b1", "bge3"):
        for o in (b1_options() if cmd == "b1" else bge3_options()):
            r = screen(*o)
            print(json.dumps({k: r[k] for k in ("v", "b", "p", "n'", "H'", "alive")} | ({"witness": r["witness"]} if r["alive"] else {"KN": r.get("KN")})), flush=True)
    elif cmd == "one":
        a = list(map(int, sys.argv[2:])); r = screen(*a, verbose=False); print(json.dumps(r))
