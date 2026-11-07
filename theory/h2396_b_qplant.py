r"""Part B controls [M]: the two Q-level methods on planted (non-empty) targets.  A uniform H-invariant quadruple on G
is drawn (only its Q-compressions matter: y(0) = e0 + 9 * (4 signs), y(x != 0) = e0 + 3 * (12 signs)); the target is
its own Q-level PAF sum.  Method F (Fourier classes, as theory/h2396_b_qlevel.py) and method S (spatial dictionary join
on all 16 lags, as theory/h2396_b_qcheck.py) are re-run as functions of the target, counting ordered quadruples without
listing them; they must agree and be >= 1.  Per-block values are capped (|y| <= VMAX) in both, so the counts refer to
the same capped space.
usage: python3 theory/h2396_b_qplant.py [trials] [VMAX]"""
import sys, itertools
from collections import Counter
import numpy as np
E = [(a, b) for a in range(4) for b in range(4)]; ix = {e: i for i, e in enumerate(E)}
def mul(a, x): return ((a[0]*x[0] - a[1]*x[1]) % 4, (a[0]*x[1] + a[1]*x[0] - a[1]*x[1]) % 4)
def ad(a, x): return ((a[0]+x[0]) % 4, (a[1]+x[1]) % 4)
def tr(x): return (2*x[0] - x[1]) % 4
xi = (0, 1)
orb = []; seen = set()
for x in E:
    if x in seen: continue
    o = sorted({x, mul(xi, x), mul(xi, mul(xi, x))}); seen |= set(o); orb.append(o)
CH = {a: np.array([1j ** tr(mul(a, x)) for x in E]) for a in E}
A_T, A_2, A_U = (1, 0), (2, 0), (1, 2)
Y0 = sorted({e + 9 * s for e in (1, -1) for s in (-4, -2, 0, 2, 4)})
Y1 = sorted({e + 3 * s for e in (1, -1) for s in range(-12, 13, 2)})
def expand(vals):
    y = [0] * 16
    for o, v in zip(orb, vals):
        for x in o: y[ix[x]] = v
    return y
def paf(y): return tuple(sum(y[ix[x]] * y[ix[ad(x, h)]] for x in E) for h in E)
def fourier_key(y):
    y = np.array(y); return (int(round(y.sum() ** 2)), int(round(abs((y * CH[A_T]).sum()) ** 2)),
                             int(round(((y * CH[A_2]).sum().real) ** 2)), int(round(abs((y * CH[A_U]).sum()) ** 2)))
def blocks(vmax, E0, rb):
    for vals in itertools.product([v for v in Y0 if abs(v) <= vmax], *([[v for v in Y1 if abs(v) <= vmax]] * 5)):
        if vals[0] ** 2 + 3 * sum(v * v for v in vals[1:]) > E0: continue
        y = expand(vals)
        if sum(y) ** 2 > rb: continue
        yield y
def count_F(tgtF, vmax, E0, rb):
    C = Counter(fourier_key(y) for y in blocks(vmax, E0, rb))
    AB = Counter()
    for k1, c1 in C.items():
        for k2, c2 in C.items(): AB[tuple(a + b for a, b in zip(k1, k2))] += c1 * c2
    return sum(c * AB.get(tuple(t - x for t, x in zip(tgtF, s)), 0) for s, c in AB.items())
def count_S(tgtS, vmax, E0, rb):
    C = Counter(paf(y) for y in blocks(vmax, E0, rb))
    AB = Counter()
    for k1, c1 in C.items():
        for k2, c2 in C.items():
            s = tuple(a + b for a, b in zip(k1, k2))
            if s[0] <= tgtS[0]: AB[s] += c1 * c2
    return sum(c * AB.get(tuple(t - x for t, x in zip(tgtS, s)), 0) for s, c in AB.items())
if __name__ == "__main__":
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 3; vmax = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    rng = np.random.default_rng(5); agree = 0
    for t in range(trials):
        while True:
            Yb = []
            for _ in range(4):
                vals = [int(rng.choice((1, -1))) + 9 * int(rng.choice((1, -1), 4).sum())]
                vals += [int(rng.choice((1, -1))) + 3 * int(rng.choice((1, -1), 12).sum()) for _ in range(5)]
                Yb.append(expand(vals))
            if max(abs(v) for y in Yb for v in y) <= vmax: break
        tgtS = tuple(sum(c) for c in zip(*[paf(y) for y in Yb]))
        tgtF = tuple(sum(c) for c in zip(*[fourier_key(y) for y in Yb]))
        rb = max(sum(y) ** 2 for y in Yb)
        nF = count_F(tgtF, vmax, tgtS[0], rb); nS = count_S(tgtS, vmax, tgtS[0], rb)
        ok = (nF == nS and nS >= 1); agree += ok
        print(f"planted {t}: Fourier-class count {nF}, spatial count {nS}, agree and >= 1: {ok}", flush=True)
    print(f"{agree}/{trials} planted targets: both methods agree on a non-zero count")
