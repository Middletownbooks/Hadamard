r"""Part B, gate B1: Gaussian estimate of the number of solutions of the 599 Z4^2 x Z37 structure (Kumar's 10^24.6,
search/NOTES.md section 103), by this repository's lattice Gaussian (theory/q523_gauss.gauss_bits) [E].
G = GR(4,2) x Z37, H = <(xi, x7)> (70 orbits), T = {1} u +-Teich^*, target sum_i PAF_i(g) = -Gamma(g + K) for g != 0.
Uniform H-invariant blocks (2^70 each); lag classes = H-orbits on G \ 0 modulo g ~ -g; E = 2^(4 * 70) * 2^(-bits).
Also the same estimate with the Q-level equations removed (the 15 lag classes of Q, i.e. the coarse part), to show how
much of the estimate is the Q-level, which is exactly empty (theory/h2396_b_qlevel.py, h2396_b_qcheck.py).
usage: python3 theory/h2396_b_gauss.py [N]"""
import sys, os, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
E = [(a, b) for a in range(4) for b in range(4)]
def mul(a, x): return ((a[0]*x[0] - a[1]*x[1]) % 4, (a[0]*x[1] + a[1]*x[0] - a[1]*x[1]) % 4)
def ng(a): return ((-a[0]) % 4, (-a[1]) % 4)
def tr(x): return (2*x[0] - x[1]) % 4
xi = (0, 1); p = 37; m = 7
teich = [(1, 0), xi, mul(xi, xi)]; T = [(0, 0)] + teich + [ng(a) for a in teich]
Gam = {x: round(4 * sum((1j ** tr(mul(a, x))).real for a in T)) for x in E}
pts = [(x, k) for x in E for k in range(p)]; ix = {P: i for i, P in enumerate(pts)}; n = len(pts)
orb = [-1] * n; reps = []
for P in pts:
    if orb[ix[P]] >= 0: continue
    cur = P; o = len(reps); reps.append(P)
    while orb[ix[cur]] < 0: orb[ix[cur]] = o; cur = (mul(xi, cur[0]), cur[1] * m % p)
no = len(reps); assert no == 70
# lag classes modulo negation
lagcls = {}
for o, P in enumerate(reps):
    if o == orb[ix[((0, 0), 0)]]: continue
    Q = (ng(P[0]), (-P[1]) % p); key = min(o, orb[ix[Q]]); lagcls.setdefault(key, P)
lags = list(lagcls.values())
tgt = np.array([-Gam[h[0]] for h in lags])
addi = np.array([[ix[(((x[0] + h[0][0]) % 4, (x[1] + h[0][1]) % 4), (k + h[1]) % p)] for (x, k) in pts] for h in lags])
orbA = np.array(orb)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
rng = np.random.default_rng(1)
V = np.zeros((N, len(lags)), np.int64)
for b in range(4):
    S = rng.choice([-1, 1], size=(N, no))
    X = S[:, orbA]                                     # N x 592
    for j in range(len(lags)):
        V[:, j] += (X * X[:, addi[j]]).sum(1)
g = gauss_bits(V, tgt)
qlev = [j for j, h in enumerate(lags) if h[1] == 0]   # lags in Q x {0}: not the Q-level, kept for information
print(dict(orbits=no, lag_classes=len(lags), bits=g.get('bits'), rank=g.get('rank'), in_lattice=g.get('in_lattice'),
           log2E=round(280 - g['bits'], 2) if g.get('bits') is not None else None,
           log10E=round((280 - g['bits']) * math.log10(2), 2) if g.get('bits') is not None else None))
