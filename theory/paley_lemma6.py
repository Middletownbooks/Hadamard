r"""Lemma 6 (report section 4): S2 at the real characters of G' fixes the u row sums by stage-1 data.

For p = 3 (mod 4) and every core of the Paley-coset ansatz (any b), at every character psi of G' with psi^2 = 1:
    sum_k |beta^_k(psi)|^2 + p sum_k |u^_k(psi)|^2 = 4v,      beta = alpha - w = y - p w,
so with S1 at psi (sum_k y^_k(psi)^2 = T^_J(psi)):
    sum_k u^_k(psi)^2 = (4v - T^_J(psi)) / p + 2 sum_k y^_k(psi) w^_k(psi) - p sum_k w^_k(psi)^2,
and u^_k(psi) = |T_k| (mod 2), |u^_k(psi)| <= |T_k| (T_k = twisted cosets of block k).  At psi = 1: T^_J = sum r_k^2,
y^_k = r_k, w^_k = t_k.  Corollary (b = 1, psi = 1, mod 8): the four blocks cannot all contain an odd number of constant
cosets.
This script verifies the identity exactly at every real character on the published 359 core and on every
data/core_v*_b1_paley.txt control, and prints the parity data.
usage: python3 theory/paley_lemma6.py"""
import sys, os, json, glob, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_algebra import load_core, components, paf

def real_chars(shape):
    """all characters with values +-1 of Z_{n1} x ... : psi(g) = prod (-1)^{e_i g_i} with e_i in {0,1}, e_i = 0 if n_i odd"""
    out = []
    for e in it.product(*[(0, 1) if n % 2 == 0 else (0,) for n in shape]):
        grids = np.meshgrid(*[np.arange(n) for n in shape], indexing='ij')
        out.append((e, (-1) ** sum(ei * g for ei, g in zip(e, grids))))
    return out

def check(v, b, X, p):
    comps = [components(x, p) for x in X]
    assert all(c[3] for c in comps)
    A = np.array([c[0] for c in comps]); U = np.array([c[1] for c in comps]); W = np.array([c[2] for c in comps])
    Y = A + (p - 1) * W; B = A - W
    shape = A.shape[1:]
    TJ = paf(Y[0]) * 0
    for k in range(4): TJ = TJ + paf(Y[k])                     # = T_J exactly (S1 holds on a core)
    res = []
    for e, psi in real_chars(shape):
        yh = (Y * psi).sum(axis=tuple(range(1, Y.ndim))); wh = (W * psi).sum(axis=tuple(range(1, W.ndim)))
        uh = (U * psi).sum(axis=tuple(range(1, U.ndim))); bh = (B * psi).sum(axis=tuple(range(1, B.ndim)))
        TJh = int((TJ * psi).sum())
        lhs = int((bh ** 2).sum() + p * (uh ** 2).sum())
        rhs_u = (4 * v - TJh) // p + 2 * int((yh * wh).sum()) - p * int((wh ** 2).sum())
        Tk = [int((W[k] == 0).sum()) for k in range(4)]
        res.append(dict(char=list(e), S2_at_psi=lhs == 4 * v, u_hat=uh.tolist(), sum_u2=int((uh ** 2).sum()), predicted=rhs_u,
                        identity=int((uh ** 2).sum()) == rhs_u and (4 * v - TJh) % p == 0,
                        parity_ok=all((int(uh[k]) - Tk[k]) % 2 == 0 for k in range(4))))
    Kk = [int((W[k] != 0).sum()) for k in range(4)]
    return dict(v=v, b=b, p=p, K_k=Kk, all_K_k_odd=all(x % 2 for x in Kk), t=W.reshape(4, -1).sum(1).tolist(),
                rho=U.reshape(4, -1).sum(1).tolist(), chars=len(res), all_identities=all(r['identity'] and r['S2_at_psi'] and r['parity_ok'] for r in res),
                detail=res[:2])

def main():
    out = []
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt'))
    out.append(check(v, b, X, mods[-1]))
    for f in sorted(glob.glob(os.path.join(ROOT, 'data/core_v*_b1_paley.txt'))):
        L = open(f).read().split('\n'); hd = dict(t.split('=') for t in L[0].split(';')[0].split() if '=' in t)
        n, p, v = int(hd['n']), int(hd['p']), int(hd['v'])
        X = [np.array([1 if c == '+' else -1 for c in l.split('=')[1]], np.int64).reshape(n, p) for l in L[1:5]]
        out.append(check(v, 1, X, p))
    for r in out:
        print(json.dumps({k: r[k] for k in r if k != 'detail'}))
    return out

if __name__ == "__main__":
    main()
