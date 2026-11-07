r"""Part B control (shares no code with theory/h2396_b_qlevel.py): the Q-level system of the 599 Z4^2 x Z37 structure
solved in the spatial domain.  Q = GR(4,2) as Z4^2 with the multiplication of Z4[xi]/(xi^2 + xi + 1); target on Q:
Y(0) = 4|G| - 36 * 4t, Y(x != 0) = -37 Gamma(x), Gamma = 4 sum_{chi in T} chi; blocks: xi-invariant y with
y(0) in e0 + 9 * (4 signs), y(x != 0) in e0 + 3 * (12 signs); every block with sum_x y^2 <= Y(0) and |PAF_y(x)| bounds is
listed and a dictionary 2+2 join on the full PAF vector (16 lags) counts ordered quadruples.  Also the variant with
H acting trivially on Q (y any function with these letters) for comparison is NOT run (2^16 x ... too large).
usage: python3 theory/h2396_b_qcheck.py"""
import itertools, cmath
from collections import Counter
E = [(a, b) for a in range(4) for b in range(4)]
ix = {e: i for i, e in enumerate(E)}
def mul(a, x): return ((a[0]*x[0] - a[1]*x[1]) % 4, (a[0]*x[1] + a[1]*x[0] - a[1]*x[1]) % 4)
def ad(a, x): return ((a[0]+x[0]) % 4, (a[1]+x[1]) % 4)
def ng(a): return ((-a[0]) % 4, (-a[1]) % 4)
def tr(x): return (2*x[0] - x[1]) % 4
xi = (0, 1)
teich = [(1, 0), xi, mul(xi, xi)]
T = [(0, 0)] + teich + [ng(a) for a in teich]
Gam = [round(4 * sum((1j ** tr(mul(a, x))).real for a in T)) for x in E]
Y = [4 * 592 - 36 * 28 if x == (0, 0) else -37 * Gam[ix[x]] for x in E]
print("Gamma:", dict(zip(E, Gam))); print("Y:", Y, "sum", sum(Y))
orb = []; seen = set()
for x in E:
    if x in seen: continue
    o = sorted({x, mul(xi, x), mul(xi, mul(xi, x))}); seen |= set(o); orb.append(o)
Y0 = sorted({e + 9 * s for e in (1, -1) for s in (-4, -2, 0, 2, 4)})
Y1 = sorted({e + 3 * s for e in (1, -1) for s in range(-12, 13, 2)})
def paf(y): return tuple(sum(y[ix[x]] * y[ix[ad(x, h)]] for x in E) for h in E)
blocks = Counter()
for vals in itertools.product(Y0, *([Y1] * 5)):
    if vals[0] ** 2 + 3 * sum(v * v for v in vals[1:]) > Y[0] - 48: continue
    y = [0] * 16
    for o, v in zip(orb, vals):
        for x in o: y[ix[x]] = v
    if sum(y) ** 2 > 28: continue
    blocks[paf(y)] += 1
print("blocks (distinct PAF vectors):", len(blocks), "total", sum(blocks.values()))
keys = list(blocks.items()); AB = Counter()
for ka, ca in keys:
    for kb, cb in keys:
        s = tuple(a + b for a, b in zip(ka, kb))
        if s[0] <= Y[0]: AB[s] += ca * cb
tot = 0
for s, c in AB.items():
    need = tuple(t - x for t, x in zip(Y, s))
    tot += c * AB.get(need, 0)
print("ordered Q-level quadruples:", tot)
