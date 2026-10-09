"""Finite-field seed tests on the published cores (v prime = 3 mod 4).
b=0 (core on Z_v = GF(v)+):     agreement of each block with +-Legendre chi(u*x + t), u in {1, nonresidue}, t in Z_v.
b=1 (core on Z_{v-1} = GF(v)*): agreement of each block with +-Sidelnikov s(i) = chi(g^(u*i + t) + 1), all units u, all t
                                 (covers chi(g^i + c) for every c != 0 and every primitive root).
Null: the same maximum over random +-1 sequences of the same length (200 draws)."""
import sys, json
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/scripts')
from hadlib import all_cores
from sympy import primitive_root
rng = np.random.default_rng(0)

def legendre(v):
    sq = np.zeros(v, dtype=int)
    sq[(np.arange(1, v) ** 2) % v] = 1
    chi = np.where(sq == 1, 1, -1); chi[0] = 0
    return chi

def best_agree_family(x, fam):
    # fam: array (F, n) of +-1/0 sequences; agreement counted on nonzero entries, both signs
    num = (fam != 0).sum(1)
    dot = (fam * x[None, :]).sum(1)
    agree = (num + np.abs(dot)) / 2 / num
    return agree.max()

def fam_b0(v):
    chi = legendre(v); nr = next(a for a in range(2, v) if chi[a] == -1)
    fam = []
    for u in (1, nr):
        for t in range(v):
            fam.append(chi[(u * (np.arange(v) + t)) % v])
    return np.array(fam)

def fam_b1(v):
    n = v - 1; g = primitive_root(v); chi = legendre(v)
    pw = np.array([pow(g, i, v) for i in range(n)])
    sid = chi[(pw + 1) % v]  # chi(g^i + 1); zero where g^i = -1
    units = [u for u in range(1, n) if np.gcd(u, n) == 1]
    fam = []
    for u in units:
        base = sid[(u * np.arange(n)) % n]
        for t in range(n):
            fam.append(np.roll(base, -t))
    return np.array(fam)

for c in all_cores():
    if c.b not in (0, 1) or len(c.mods) != 1 and c.b == 0:
        pass
    if c.b == 0 and c.m == c.v:
        fam = fam_b0(c.v); kind = 'b0 Legendre'
    elif c.b == 1 and c.m == c.v - 1:
        fam = fam_b1(c.v); kind = 'b1 Sidelnikov'
    else:
        continue
    # cores are stored in mixed-radix coordinates of their group; convert to the cyclic index if needed
    if len(c.mods) > 1:
        # CRT: element with coords (x1..xk) <-> integer i with i = xj mod mods[j]
        n = c.m; coords = c.coords(); idx_cyc = np.zeros(n, dtype=int)
        for i in range(n):
            idx_cyc[i] = c.idx(np.array([i % mm for mm in c.mods]))
        X = c.X[:, idx_cyc]
    else:
        X = c.X
    real = [round(float(best_agree_family(X[k], fam)), 3) for k in range(4)]
    null = [best_agree_family(rng.choice([-1, 1], size=c.m), fam) for _ in range(60)]
    print(json.dumps(dict(core=c.name, v=c.v, test=kind, family_size=len(fam), real_max_agreement=real,
                          null_mean=round(float(np.mean(null)), 3), null_max=round(float(np.max(null)), 3))), flush=True)
