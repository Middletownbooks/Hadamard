r"""Stage-2 rigidity of 599, c = 6, q = 11, s = 2, H' = Z3^2, U = -1 on one Z3 axis of H' ("Z3^2-a") [P, M, E].

Slice (logs/s1z/budget_seeds_599rich.log, key [599, 6, 11, 2, [3, 3], [[3, 2, 2]]]): G = Q x Z3' x Z3'' x Z11,
Q = Z2 x Z3, H' = Z3' x Z3'', U = -1 on Z3''.  Write L = Q x Z3' x Z11 (198 points), r = x(., 0, .), s = x(., 1, .) =
x(., 2, .), d = r - s in {0, +-2}.

(1) [P] For chi = (chi_L, phi'') with phi'' != 1, x^(chi) = d^(chi_L) (phi''(1) + phi''(2) = -1).  These chi are
    nontrivial on H, so sum_k PAF_L(d_k) = 4v delta: at lag 0, W := sum_{k, l} (d_k(l)/2)^2 = v = 599.  Its
    Z11-compression (stage 1) gives sum_{k, l'} (delta_k(l')/2)^2 = v, delta = y(., 0) - y(., 1) on L' = Q x Z3'.
(2) [P] Per fibre (k, l') with d/2 = (e0 at z = 0, e1 on QR, e2 on NQR):
        (delta/2)^2 - sum_z (d/2)^2 = 20 e1^2 + 20 e2^2 + 10 e0 (e1 + e2) + 50 e1 e2 =: cross,
    and, by the letters (y = alpha + 5 (p + n)), cross = -10 [free signs differ] for a pair of small letters (+-1, +-1),
    and is fixed by y for every other pair.  Hence a core needs
        D11 = D*(y),   D11 = #{(k, l'): both letters +-1 and their free signs differ},   D* = sum of the fixed parts / 10.
    Under uniform lifts D11 ~ Bin(N11, 1/2), N11 = #{(k, l'): both letters +-1}.
(3) [M] The per-pair values are tabulated here by brute force over the fibre patterns (no algebra used), the identity is
    checked on random lifts, and an exact DP over all letter-pair compositions with the stage-1 lag-0 energy T1(0) = 2176
    and the d-energy 599 (both necessary for every stage-1 solution) gives the smallest -log2 P(Bin(N11, 1/2) = D*) that
    any stage-1 solution can have, and its excess over the modal value (what an at-mean Gaussian charges) [E for the
    reading as a cost: one equation of the stage-2 system, correlations with the rest ignored].
(4) --case negQ3: the two 36-orbit slices with U = -1 on Q's Z3 (H' = Z9 or Z3^2).  The same split along Q's Z3 gives
    L = Z2 x H' x Z11 and characters (chi_L, phi) with phi != 1 on Q's Z3; those with chi_L trivial on H need 4b, so
    W = (196 v + 2b)/198 = 593 and sum (delta/2)^2 = (16 v + 2b)/18 = 533: a core needs D11 = D* + 6.
usage: python3 theory/b599x4_z3a_rigidity.py [--case z3a|negQ3] [--out FILE]"""
import sys, os, json, math, itertools
from collections import defaultdict
import numpy as np

LET = (-11, -9, -1, 1, 9, 11)
V = 599


def patterns(y):
    """fibre patterns (alpha, p, n) in {+-1}^3 with alpha + 5 p + 5 n = y"""
    return [(a, p, n) for a in (1, -1) for p in (1, -1) for n in (1, -1) if a + 5 * p + 5 * n == y]


def pair_table():
    """per letter pair (y_r, y_s): energy y_r^2 + 2 y_s^2, (delta/2)^2, the multiset of cross values over the lifts"""
    T = {}
    for yr in LET:
        for ys in LET:
            cr = []
            for (ar, pr, nr) in patterns(yr):
                for (as_, ps, ns) in patterns(ys):
                    e0, e1, e2 = (ar - as_) // 2, (pr - ps) // 2, (nr - ns) // 2
                    dd = (e0 + 5 * e1 + 5 * e2) ** 2
                    w = e0 * e0 + 5 * e1 * e1 + 5 * e2 * e2
                    assert dd == ((yr - ys) // 2) ** 2
                    cr.append(dd - w)
            T[(yr, ys)] = dict(E=yr * yr + 2 * ys * ys, dd=((yr - ys) // 2) ** 2, cross=sorted(cr),
                               small=(abs(yr) == 1 and abs(ys) == 1))
    return T


def check_identity(T, trials=200, seed=1):
    """random y (any letters) and random lifts: sum (delta/2)^2 - W = 10 (D* - D11)"""
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        tot_dd = 0; tot_w = 0; D11 = 0; Dstar = 0
        for _ in range(72):
            yr, ys = LET[rng.integers(6)], LET[rng.integers(6)]
            Pr, Ps = patterns(yr), patterns(ys)
            (ar, pr, nr) = Pr[rng.integers(len(Pr))]; (as_, ps, ns) = Ps[rng.integers(len(Ps))]
            d = np.array([ar - as_] + [pr - ps] * 5 + [nr - ns] * 5) // 2
            tot_dd += int(d.sum()) ** 2; tot_w += int((d * d).sum())
            if T[(yr, ys)]['small']:
                D11 += int(pr != ps)
            else:
                Dstar += T[(yr, ys)]['cross'][0] // 10
        assert tot_dd - tot_w == 10 * (Dstar - D11)
    return True


def log2binom(n, k):
    if k < 0 or k > n: return -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2) - n


def best_penalty(T, npairs=72, T0=2176, DD=V, shift=0):
    """exact DP over the multiset of non-small pairs: excess energy (E - 3), dd, D*; then A small pairs with
    alpha_r != alpha_s carry the rest of the d-energy (dd 1 each, energy 3 each).  A core needs D11 = D* + shift."""
    big = [(k, t) for k, t in T.items() if not t['small']]
    for k, t in big: assert len(set(t['cross'])) == 1 and t['cross'][0] % 10 == 0
    opts = sorted({(t['E'] - 3, t['dd'], t['cross'][0] // 10) for k, t in big})
    exc_tot = T0 - 3 * npairs
    # states: (excess, dd, D*, n) reachable; multiset DP in a fixed option order
    states = {(0, 0, 0, 0)}
    for (ex, dd, ds) in opts:
        new = set(states)
        frontier = set(states)
        while frontier:
            nxt = set()
            for (e, d_, s_, n) in frontier:
                e2, d2, s2, n2 = e + ex, d_ + dd, s_ + ds, n + 1
                if e2 > exc_tot or d2 > DD or n2 > npairs: continue
                st = (e2, d2, s2, n2)
                if st not in new: new.add(st); nxt.add(st)
            frontier = nxt
        states = new
    res = []
    for (e, d_, s_, n) in states:
        if e != exc_tot: continue
        n11 = npairs - n; A = DD - d_
        if A < 0 or A > n11: continue
        lp = log2binom(n11, s_ + shift)
        if lp == -math.inf: continue
        mode = log2binom(n11, n11 // 2)
        res.append(dict(n_nonsmall=n, N11=n11, A=A, Dstar=s_, D11_needed=s_ + shift, log2P=round(lp, 2),
                        excess_over_mode=round(mode - lp, 2)))
    res.sort(key=lambda r: -r['log2P'])
    return res


def main(out=None, case='z3a'):
    T = pair_table()
    small_ok = all(sorted(set(t['cross'])) == [-10, 0] for t in T.values() if t['small'])
    fixed_ok = all(len(set(t['cross'])) == 1 for t in T.values() if not t['small'])
    assert small_ok and fixed_ok
    check_identity(T)
    if case == 'z3a':
        name = "599 c=6 q=11 s=2 H'=Z3^2, -1 on one axis (Z3^2-a)"; DD, W = V, V
    else:
        b = 5; name = "599 c=6 q=11 s=2 H'=Z9 or Z3^2, -1 on Q's Z3"
        assert (196 * V + 2 * b) % 198 == 0 and (16 * V + 2 * b) % 18 == 0
        W = (196 * V + 2 * b) // 198; DD = (16 * V + 2 * b) // 18
    assert (DD - W) % 10 == 0
    rows = best_penalty(T, DD=DD, shift=-(DD - W) // 10)
    summary = dict(slice=name, W_target=W, stage1_dd_target=DD, D11_minus_Dstar=-(DD - W) // 10, identity_checked=True,
                   small_pair_cross=[-10, 0], fixed_cross_values={f"{k}": T[k]['cross'][0] for k in T if not T[k]['small']},
                   feasible_compositions=len(rows),
                   best=rows[:8],
                   worst_of_best10_excess=max(r['excess_over_mode'] for r in rows[:10]) if rows else None,
                   min_excess_over_mode=min(r['excess_over_mode'] for r in rows) if rows else None)
    print(json.dumps(summary, indent=1))
    if out: json.dump(dict(summary, all=rows), open(out, 'w'), indent=1)
    return summary


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[a.index('--out') + 1] if '--out' in a else None, a[a.index('--case') + 1] if '--case' in a else 'z3a')
