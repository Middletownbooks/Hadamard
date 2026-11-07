r"""Part B of report/prompt_599_plugin.md: Kumar's bordered array for n = 599 over Z4^2 x Z37 with a partial
Teichmueller character set (search/NOTES.md section 103, "599 over Z4^2 x Z37 (t = 7, 70 orbits, 10^24.6): no border in
any symmetric class tested"; the family of section 18).  Credit: the family and the border reduction are K. Kumar's
(vatskishan03/hadamard, Apache-2.0); this code is new.

Reconstruction (gate B1) [P, M]:
  Q = GR(4, 2) = Z4[xi]/(xi^2 + xi + 1) (additively Z4^2), G = Q x Z37, K = 0 x Z37, |G| = 592 = 599 - 7.
  T = {1} u {chi_a : a in +-Teich^*}, Teich^* = {1, xi, xi^2}, chi_a(x) = i^{Tr(ax)}, Tr(x0 + x1 xi) = 2 x0 - x1 (mod 4):
      the trivial character and one orbit of 6 order-4 characters under <xi, -1>; the other such orbit,
      (1 + 2 xi) Teich^* u its negative, is its image under the automorphism x -> (1 + 2 xi) x, which commutes with xi.
  H = <(xi, x m)>, m of order 9 in Z37^*: orbits 1 + 5 + 4 + 60 = 70.
  Target: sum_i PAF_i(g) = -Gamma(g + K) (g != 0), Gamma = 4 sum_{chi in T} chi, border of size 4t = 28.
Q-level (coset sums y_i(x) = sum_k x_i(x, k)): y_i is xi-invariant; y(0) = e0 + 9 (4 signs), y(x != 0) = e0 + 3 (12 signs).
  Fourier form [P]: sum_i |y_i^(chi)|^2 = 28 for chi in T and 2396 = 4n for chi not in T.  The 16 characters fall into
  four xi/conjugation classes (trivial; T \ 1; order 2; the other units), so the Q-level system is four equations.
usage: python3 theory/h2396_b_qlevel.py classes [OUT.json]     (all Q-level solutions, and their (r, z) classes)
       python3 theory/h2396_b_qlevel.py orbits                (check the 70 orbits)"""
import sys, json, itertools, math
from collections import defaultdict, Counter
import numpy as np

ELS = [(a, b) for a in range(4) for b in range(4)]
IDX = {e: i for i, e in enumerate(ELS)}
def mul(a, x): return ((a[0] * x[0] - a[1] * x[1]) % 4, (a[0] * x[1] + a[1] * x[0] - a[1] * x[1]) % 4)
def add(a, x): return ((a[0] + x[0]) % 4, (a[1] + x[1]) % 4)
def neg(a): return ((-a[0]) % 4, (-a[1]) % 4)
def tr(x): return (2 * x[0] - x[1]) % 4
XI = (0, 1); ONE = (1, 0)
TEICH = [ONE, XI, mul(XI, XI)]
TSET = [(0, 0)] + TEICH + [neg(a) for a in TEICH]
CHI = {a: np.array([1j ** tr(mul(a, x)) for x in ELS]) for a in ELS}

def xi_orbits():
    seen = set(); orbs = []
    for x in ELS:
        if x in seen: continue
        o = [x, mul(XI, x), mul(XI, mul(XI, x))]
        o = sorted(set(o)); seen |= set(o); orbs.append(o)
    return orbs
ORBS = xi_orbits()                    # [ [0], 2Teich*, Teich*, ... ] in some order

def orbit_check(m=None):
    p = 37
    ords = {}
    for g in range(2, p):
        o = 1; x = g
        while x != 1: x = x * g % p; o += 1
        ords[g] = o
    m = m or next(g for g in range(2, p) if ords[g] == 9)
    pts = [(x, k) for x in ELS for k in range(p)]
    seen = set(); norb = 0
    for P in pts:
        if P in seen: continue
        norb += 1; cur = P
        while cur not in seen:
            seen.add(cur); cur = (mul(XI, cur[0]), cur[1] * m % p)
    return dict(m=m, orbits=norb)

def letters():
    Y0 = sorted({e + 9 * s for e in (1, -1) for s in (-4, -2, 0, 2, 4)})
    Y1 = sorted({e + 3 * s for e in (1, -1) for s in range(-12, 13, 2)})
    return Y0, Y1

# representative characters of the four classes
A_T = ONE                 # in T
A_2 = (2, 0)              # order 2
A_U = add(ONE, mul((2, 0), XI))   # 1 + 2 xi, a unit outside +-Teich^*
assert A_U not in TSET and A_T in TSET and A_2 not in TSET

def block_fourier(yvals):
    """yvals: dict orbit-index -> value; returns (r, z_T, y_2, z_U) and the full vector"""
    y = np.zeros(16)
    for oi, o in enumerate(ORBS):
        for x in o: y[IDX[x]] = yvals[oi]
    r = y.sum(); zT = (y * CHI[A_T]).sum(); y2 = (y * CHI[A_2]).sum().real; zU = (y * CHI[A_U]).sum()
    return r, zT, y2, zU, y

def classes():
    Y0, Y1 = letters()
    zero_oi = next(i for i, o in enumerate(ORBS) if o == [(0, 0)])
    others = [i for i in range(len(ORBS)) if i != zero_oi]
    Ebound = 1360 - 3 * 16
    blocks = []
    for y0 in Y0:
        if y0 * y0 > Ebound: continue
        for vals in itertools.product([v for v in Y1 if 3 * v * v <= Ebound], repeat=5):
            E = y0 * y0 + 3 * sum(v * v for v in vals)
            if E > Ebound: continue
            yv = {zero_oi: y0}
            for oi, v in zip(others, vals): yv[oi] = v
            r, zT, y2, zU, y = block_fourier(yv)
            if abs(r) not in (2, 4): continue
            if abs(zT) ** 2 > 28 + 1e-9 or y2 ** 2 > 2396 + 1e-9 or abs(zU) ** 2 > 2396 + 1e-9: continue
            blocks.append((int(round(r)), complex(round(zT.real), round(zT.imag)), int(round(y2)),
                           int(round(abs(zU) ** 2)), tuple(int(v) for v in y), E))
    # 2+2 join on (r^2, |zT|^2, y2^2, |zU|^2)
    key = lambda b: (b[0] ** 2, int(round(abs(b[1]) ** 2)), b[2] ** 2, b[3])
    by = defaultdict(list)
    for b in blocks: by[key(b)].append(b)
    ks = list(by)
    pair = defaultdict(list)
    for k1 in ks:
        for k2 in ks:
            s = tuple(a + c for a, c in zip(k1, k2))
            if s[0] <= 28 and s[1] <= 28 and s[2] <= 2396 and s[3] <= 2396: pair[s].append((k1, k2))
    sols = []
    tgt = (28, 28, 2396, 2396)
    for s, lst in pair.items():
        need = tuple(t - x for t, x in zip(tgt, s))
        if need in pair:
            for (k1, k2) in lst:
                for (k3, k4) in pair[need]:
                    sols.append((k1, k2, k3, k4))
    nq = 0; rz = Counter(); examples = {}
    for (k1, k2, k3, k4) in sols:
        cnt = len(by[k1]) * len(by[k2]) * len(by[k3]) * len(by[k4])
        nq += cnt
    # (r, z) classes: multiset of (|r|, |z|^2) per block, and the realised (r, z) combinations
    for (k1, k2, k3, k4) in sols:
        for combo in itertools.product(by[k1], by[k2], by[k3], by[k4]):
            sig = tuple(sorted((abs(b[0]), int(round(abs(b[1]) ** 2))) for b in combo))
            rz[sig] += 1
            if sig not in examples: examples[sig] = [list(b[4]) for b in combo]
    return dict(blocks=len(blocks), key_classes=len(ks), ordered_quadruples=nq,
                rz_classes={str(k): v for k, v in rz.items()}, examples={str(k): v for k, v in examples.items()})

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'classes'
    if cmd == 'orbits':
        print(json.dumps(dict(T=TSET, xi_orbits=ORBS, **orbit_check())))
    elif cmd == 'classes':
        out = classes()
        print(json.dumps({k: v for k, v in out.items() if k != 'examples'}, indent=1))
        if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1)
