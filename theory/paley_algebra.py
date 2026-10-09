r"""Exact algebraic reformulation of the Paley-coset ansatz (p = 3 mod 4), checked on the published v = 359 core.

Setting: G = G' x Z_p, each block x_k(g', .) = alpha*delta + u*chi + w*eps on every Z_p-coset (exactly one of u, w is
nonzero).  As G-developed matrices the blocks lie in Z[G'] (x) span{I, Q, J}, Q = Paley core of Z_p (Q^T = -Q,
Q^2 = J - pI, QJ = 0).  The algebra span{I, Q, J} (x) Q splits by the idempotents J/p and I - J/p into
    Q  (J -> p, Q -> 0)              : x  ->  y = alpha + (p-1) w        (Z_p-coset sums, real sequence on G')
    Q(sqrt(-p)) (J -> 0, Q -> i*sqrt(p)): x  ->  Z = beta + i*sqrt(p) u,  beta = alpha - w
so sum_k PAF(x_k) = T (target) holds iff BOTH
    (S1)  sum_k PAF(y_k) = T_J                      (T_J = T0 + (p-1) Ts on G'),
    (S2)  sum_k PAF^C(Z_k) = 4v * delta_0           (complex periodic autocorrelation, conj on the shifted copy),
where for the b >= 3 target T = 4m on 0, -4b on H \ 0, +4 off H (H = H' x Z_p) one has T_J(0) = 4v + 4p(n' - n') = 4v,
T_J = 4b - 4v + 4v... (computed below), and the S2 target is 4v * delta for every b (b = 0, 1, >= 3).
Writing pi = (1 + i sqrt p)/2 (norm q = (p+1)/4):  Z/2 = pi*u + e,  e = (alpha - u - w)/2 in {0, +-1}.

usage: python3 theory/paley_algebra.py check359      (all identities on the published core; exact integer arithmetic)
       python3 theory/paley_algebra.py anatomy359    (components in the (y, beta, u, Z) coordinates)"""
import sys, os, json, itertools as it
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def legendre(p):
    chi = np.zeros(p, dtype=np.int64)
    for z in range(1, p): chi[z] = 1 if pow(z, (p - 1) // 2, p) == 1 else -1
    return chi

def load_core(path):
    L = open(path).read().split('\n')
    hdr = L[0]
    mods = tuple(int(t) for t in hdr.split('group=')[1].split()[0].split('x'))
    X = [np.array([1 if ch == '+' else -1 for ch in l.split('=')[1].strip()], np.int64) for l in L if '=' in l and l[0] in 'ABCD']
    b = int(hdr.split('b=')[1].split()[0]); v = int(hdr.split('v=')[1].split()[0])
    return v, b, mods, [x.reshape(mods) for x in X]

def paf(x):
    """periodic autocorrelation over the abelian group given by x.shape: P(g) = sum_y x(y) x(y + g) (exact for ints)"""
    F = np.fft.fftn(x)
    return np.rint(np.fft.ifftn(np.abs(F) ** 2).real).astype(np.int64)

def cc(f, g):
    """CC(f, g)(h) = sum_x f(x) g(x + h)"""
    F, Gf = np.fft.fftn(f), np.fft.fftn(g)
    return np.rint(np.fft.ifftn(np.conj(F) * Gf).real).astype(np.int64)

def paf_exact(x):
    """direct integer PAF over Z_{n1} x ... (no FFT), for verification"""
    shp = x.shape; out = np.zeros(shp, np.int64)
    for h in it.product(*[range(n) for n in shp]):
        out[h] = int((x * np.roll(x, [-t for t in h], axis=tuple(range(len(shp))))).sum())
    return out

def components(x, p):
    """x: array over G' x Z_p (last axis Z_p).  Returns alpha, u, w on G' (int arrays) and ok flag (in ansatz)."""
    chi = legendre(p)
    a = x[..., 0]
    rest = x[..., 1:]
    f1 = x[..., 1]                                    # value at z = 1 (a QR)
    nq = next(z for z in range(1, p) if chi[z] == -1)
    fN = x[..., nq]
    u = (f1 - fN) // 2; w = (f1 + fN) // 2
    recon = a[..., None] * (np.arange(p) == 0) + u[..., None] * chi[None, :] + w[..., None] * (np.arange(p) != 0)
    ok = bool((recon == x).all()) and bool(((u == 0) ^ (w == 0)).all())
    return a.astype(np.int64), u.astype(np.int64), w.astype(np.int64), ok

def subgroup_mask_from_target(P, b):
    """H' x Z_p from the sum PAF: lags (g', 0) with value -4b (plus 0)"""
    Hm = (P[..., 0] == -4 * b); Hm[(0,) * Hm.ndim] = True
    return Hm

def check359():
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt'))
    p = mods[-1]; Gp = mods[:-1]; n = int(np.prod(Gp)); m = n * p
    out = {"v": v, "b": b, "G'": list(Gp), "p": p}
    P = sum(paf(x) for x in X)
    assert (P == sum(paf_exact(x) for x in X)).all()
    Hm = subgroup_mask_from_target(P, b)
    out["core_ok"] = bool(P[(0,) * len(mods)] == 4 * m and all(
        P[g] == (-4 * b if Hm[g[:-1]] else 4) for g in it.product(*[range(t) for t in mods]) if any(g)))
    comps = [components(x, p) for x in X]
    out["in_ansatz"] = all(c[3] for c in comps)
    A = [c[0] for c in comps]; U = [c[1] for c in comps]; W = [c[2] for c in comps]
    beta = [a - w for a, w in zip(A, W)]
    y = [a + (p - 1) * w for a, w in zip(A, W)]
    # (S2 real part)  sum PAF(beta) + p sum PAF(u) = 4v delta
    I1 = sum(paf_exact(bb) for bb in beta) + p * sum(paf_exact(uu) for uu in U)
    d0 = np.zeros(Gp, np.int64); d0[(0,) * len(Gp)] = 1
    out["S2_real: sum PAF(beta) + p sum PAF(u) == 4v delta"] = bool((I1 == 4 * v * d0).all())
    # (S2 imag part)  sum [CC(beta, u) - CC(u, beta)] = 0
    I2 = sum(cc(bb, uu) - cc(uu, bb) for bb, uu in zip(beta, U))
    out["S2_imag: sum CC(beta,u) - CC(u,beta) == 0"] = bool((I2 == 0).all())
    # complex form: Z = beta + i sqrt(p) u ; sum |FFT Z|^2 = 4v at every character of G'
    spec = sum(np.abs(np.fft.fftn(bb + 1j * np.sqrt(p) * uu)) ** 2 for bb, uu in zip(beta, U))
    out["S2_complex: max |sum |Zhat|^2 - 4v|"] = float(np.abs(spec - 4 * v).max())
    # (S1) sum PAF(y) = T_J ; T_J(0) = 4v - ... computed from the target: T0 + (p-1) Ts
    TJ = np.where(Hm, -4 * b * p, 4 * p).astype(np.int64); TJ[(0,) * len(Gp)] = 4 * m - 4 * b * (p - 1)
    I3 = sum(paf_exact(yy) for yy in y)
    out["S1: sum PAF(y) == T0 + (p-1) Ts"] = bool((I3 == TJ).all())
    # spectral form of S1: sum |yhat|^2 = 4v off H'^perp, 4b on H'^perp \ 1, (sum r)^2 = 4v at 1
    specy = sum(np.abs(np.fft.fftn(yy)) ** 2 for yy in y)
    # H'^perp: characters trivial on H'
    Hel = [g for g in it.product(*[range(t) for t in Gp]) if Hm[g]]
    perp = np.ones(Gp, bool)
    for chi_idx in it.product(*[range(t) for t in Gp]):
        val = sum(np.exp(2j * np.pi * sum(c * h / t for c, h, t in zip(chi_idx, hh, Gp))) for hh in Hel)
        perp[chi_idx] = abs(val - len(Hel)) < 1e-6
    tgt = np.where(perp, 4 * b, 4 * v).astype(float); tgt[(0,) * len(Gp)] = 4 * v
    out["S1_spectral: max |sum |yhat|^2 - (4b on H'perp\\1, 4v else)|"] = float(np.abs(specy - tgt).max())
    # pi-decomposition: Z/2 = pi*u + e, e = (alpha - u - w)/2 in {0, +-1}
    E = [(a - uu - w) // 2 for a, uu, w in zip(A, U, W)]
    out["e in {0,+-1} and alpha-u-w even"] = all(((a - uu - w) % 2 == 0).all() and np.isin(e, [-1, 0, 1]).all()
                                                for a, uu, w, e in zip(A, U, W, E))
    q = (p + 1) // 4
    # sum |pi*u + e|^2-spectrum = v
    pi = (1 + 1j * np.sqrt(p)) / 2
    spec2 = sum(np.abs(np.fft.fftn(pi * uu + e)) ** 2 for uu, e in zip(U, E))
    out["sum |(pi u + e)^|^2 == v everywhere: max dev"] = float(np.abs(spec2 - v).max())
    # norm count: q * #twisted + #(constant with w = -alpha) = v
    ntw = sum(int((uu != 0).sum()) for uu in U)
    nminus = sum(int(((w != 0) & (w == -a)).sum()) for a, w in zip(A, W))
    nplus = sum(int(((w != 0) & (w == a)).sum()) for a, w in zip(A, W))
    out.update(twisted=ntw, const_w_eq_minus_alpha=nminus, const_w_eq_alpha=nplus,
               norm_identity=f"q*N_tw + N_- = {q}*{ntw} + {nminus} = {q * ntw + nminus} (v = {v})")
    print(json.dumps(out, indent=1))
    return out

def anatomy359():
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt'))
    p = mods[-1]; Gp = mods[:-1]
    for k, x in enumerate(X):
        a, u, w, ok = components(x, p)
        y = a + (p - 1) * w; beta = a - w
        print(f"block {'ABCD'[k]}: rowsum {int(x.sum())}")
        for name, arr in (("alpha", a), ("u", u), ("w", w), ("y", y), ("beta", beta)):
            print(f"  {name:5s}", ' | '.join(' '.join(f"{int(t):+d}" if t else ' 0' for t in arr[i, j]) for i in range(Gp[0]) for j in range(Gp[1])))

if __name__ == "__main__":
    {"check359": check359, "anatomy359": anatomy359}[sys.argv[1]]()
