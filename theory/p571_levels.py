r"""Exact Z6-level decomposition of b = 1 GS cores on Z_m = Z2 x Z3 x Z_q (q = m/6 coprime to 6), for multiplier
slices M that fix Z6 pointwise (t = 1 mod 6; the 571 x27 slice: q = 95, M = <27> on Z95; its control analogue
v = 211: q = 35, M = <27> on Z35, same orbit skeleton).

[P] Write x(a, b, t) (a in Z2, b in Z3, t in Z_q; cyclic index i <-> (i mod 2, i mod 3, i mod q)).  For each character
lambda = (eps, eta) of Z6 put c_lambda(t) = sum_{a,b} eps^a eta^b x(a, b, t).  Then x^(lambda, theta) = c^_lambda(theta),
so sum_k PAF_k = 4m delta - 4 (b = 1) is equivalent to four independent "level" conditions on Z_q
    L1  (eps, eta) = (+1, 1):  X1 = y0 + y1            sum_k PAF(X1_k) = 4v delta - (4v - 4)/q
    L2  (-1, 1):               Xm = y0 - y1            sum_k PAF(Xm_k) = 4v delta
    L3  (+1, w):               S  = Z0 + Z1 (complex)  sum_k PAF_C(S_k) = 4v delta
    L4  (-1, w):               D  = Z0 - Z1 (complex)  sum_k PAF_C(D_k) = 4v delta
(y_a = sum_b x(a, b, .), Z_a = sum_b w^b x(a, b, .); (eta = w-bar gives the conjugate condition)), plus row sums.
L1 + L2 = S1 on y (Z3-compression, Z_2q = Z2 x Z_q), L3 + L4 = S2 on Z' = Z / -2 (p571_ternary.py).
All four levels are M-invariant sequences on Z_q (R orbits); they are coupled only through the +-1 symbol
x(., ., O) in {+-1}^6 per (block, orbit).

Stage-2 hierarchy used by the solver (p571_hsolve.py): given y, S-level classes first (L3), then D | S (L4).
This module: skeleton, components, exact check on real cores, per-level Gaussian constraint bits [E], and the
per-block S/D class entropies for a given y.
usage: python3 theory/p571_levels.py check | python3 theory/p571_levels.py est m gen [--N 20000] [--ny 8]"""
import sys, os, json, math, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_slicelib import mgroup, orbits, lag_reps, gauss_logp

W = np.exp(2j * np.pi / 3)

class Skel:
    def __init__(self, m, gens):
        assert m % 6 == 0 and math.gcd(m // 6, 6) == 1
        self.m = m; self.q = q = m // 6; self.v = m + 1
        G = mgroup(m, gens); assert all(t % 6 == 1 for t in G), "slice must fix Z6 pointwise"
        self.G = G; self.Gq = sorted({t % q for t in G})
        self.orbs = orbits(q, self.Gq); self.R = len(self.orbs); self.sizes = np.array([len(o) for o in self.orbs])
        self.oid = np.zeros(q, np.int64)
        for j, o in enumerate(self.orbs): self.oid[o] = j
        self.reps, _ = lag_reps(q, self.Gq)                     # nonzero lag reps of <M, -1> on Z_q
        i = np.arange(m); self.ia, self.ib, self.it = i % 2, i % 3, i % q
    def to_cube(self, x):
        """x (..., m) -> (..., 2, 3, q)"""
        X = np.zeros(x.shape[:-1] + (2, 3, self.q), x.dtype); X[..., self.ia, self.ib, self.it] = x; return X
    def from_orbits(self, sym):
        """sym (..., 2, 3, R) orbit values -> x (..., m)"""
        return sym[..., self.ia, self.ib, self.oid[self.it]]

def comps(Xc):
    """Xc (..., 2, 3, L) -> X1, Xm (real), S, D (complex), each (..., L)"""
    y = Xc.sum(-2); Z = (Xc * W ** np.arange(3)[:, None]).sum(-2)
    return y[..., 0, :] + y[..., 1, :], y[..., 0, :] - y[..., 1, :], Z[..., 0, :] + Z[..., 1, :], Z[..., 0, :] - Z[..., 1, :]

def paf_c(c):
    """sum over blocks of PAF_C along the last axis; c (..., 4, q) -> (..., q), P(h) = sum_g c(g + h) conj c(g)"""
    F = np.fft.fft(c, axis=-1)
    return np.fft.ifft(np.abs(F) ** 2, axis=-1).sum(-2)

def targets(S):
    q, v = S.q, S.v
    t1 = np.full(q, -(4 * v - 4) // q); t1[0] = 4 * v - (4 * v - 4) // q
    t0 = np.zeros(q); t0[0] = 4 * v
    return t1, t0

def check_core(x, S):
    """exact level check of a b = 1 core x (4, m) (any x, the slice is not needed)"""
    Xc = S.to_cube(x); X1, Xm, Ss, Dd = comps(Xc)
    t1, t0 = targets(S)
    P = [paf_c(c.astype(complex)) for c in (X1, Xm, Ss, Dd)]
    full = np.rint(np.fft.ifft(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
    return dict(core=bool(full[0] == 4 * S.m and (full[1:] == -4).all()),
                L1=bool(np.allclose(P[0], t1)), L2=bool(np.allclose(P[1], t0)),
                L3=bool(np.allclose(P[2], t0)), L4=bool(np.allclose(P[3], t0)))

def check():
    from p571_ternary import load_core
    out = []
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    x, m = load_core(os.path.join(root, 'data', 'core_v43_b1_paley.txt'))
    S = Skel(m, []); out.append(dict(file='v43 (m = 42 = 2*3*7)', **check_core(x, S)))
    rng = np.random.default_rng(1)
    for m, g in ((210, [97]), (570, [103])):
        S = Skel(m, g); sym = rng.choice([1, -1], size=(4, 2, 3, S.R)); x = S.from_orbits(sym)
        inv = all(np.array_equal(x[:, (np.arange(m) * t) % m], x) for t in S.G)
        Xc = S.to_cube(x); X1, Xm, Ss, Dd = comps(Xc)
        full = np.fft.ifft(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=1).real.sum(0)
        F = np.zeros((2, 3, S.q)); F[S.ia, S.ib, S.it] = full                   # full PAF at lag (a, b, t)
        ok = True
        for (eps, eta), c in zip(((1, 1), (-1, 1), (1, W), (-1, W)), (X1, Xm, Ss, Dd)):
            lam = np.array([[eps ** a * eta ** b for b in range(3)] for a in range(2)])
            # sum over (a, b) of conj(lambda)(a, b) * PAF(a, b, t) equals sum_k PAF_C(c_k)(t) up to conjugation
            red = (np.conj(lam)[:, :, None] * F).sum((0, 1)); red2 = (lam[:, :, None] * F).sum((0, 1))
            P = paf_c(c.astype(complex))
            ok &= bool(np.allclose(P, red) or np.allclose(P, red2))
        out.append(dict(slice=(m, g), invariant=inv, level_split_exact=ok))
    return out

# ---------------------------------------------------------------- symbol tables
SYMS = np.array(list(it.product((1, -1), repeat=6))).reshape(64, 2, 3)
def sym_comps():
    X1, Xm, Ss, Dd = comps(SYMS[..., None].astype(complex))
    return X1[:, 0].real.round().astype(int), Xm[:, 0].real.round().astype(int), Ss[:, 0], Dd[:, 0]

def cond_entropy():
    """per-orbit entropies (uniform symbol): H(X1), H(Xm|X1), H(S|y), H(D|S,y) in bits [exact]"""
    X1, Xm, Ss, Dd = sym_comps()
    y = SYMS.sum(2)                                       # (64, 2)
    def H(keys):
        from collections import Counter
        c = Counter(keys); p = np.array(list(c.values())) / 64.0
        return float(-(p * np.log2(p)).sum())
    k1 = [(a,) for a in X1]; k12 = [(a, b) for a, b in zip(X1, Xm)]
    rS = lambda z: (round(z.real, 6), round(z.imag, 6))
    k123 = [(a, b, rS(s)) for a, b, s in zip(X1, Xm, Ss)]; k1234 = [(a, b, rS(s), rS(d)) for a, b, s, d in zip(X1, Xm, Ss, Dd)]
    return dict(H_X1=H(k1), H_Xm_given_X1=H(k12) - H(k1), H_S_given_y=H(k123) - H(k12), H_D_given_Sy=H(k1234) - H(k123),
                total=H(k1234))

# ---------------------------------------------------------------- per-level Gaussian constraint bits [E]
def level_bits(S, N=20000, seed=1, rows=(2, 0, 0, 0)):
    """sample uniform slice x with row sums; per level: Gaussian -log2 density of the level equations at the
    target (L1, L2 real: lags = reps of <M, -1>; L3, L4 complex: Re/Im (a, b) coordinates, lattice index 3 per lag)"""
    from p571_slicelib import RowSumSampler
    rng = np.random.default_rng(seed)
    # orbit sizes of the full slice: (a, b, O) has size |O|
    sizes = [int(S.sizes[j]) for a in range(2) for b in range(3) for j in range(S.R)]
    rs = RowSumSampler(sizes)
    lab = np.stack([rs.sample(rng, r, N) for r in rows], 1).reshape(N, 4, 2, 3, S.R).astype(np.int64)
    xo = lab                                                   # orbit values
    Xq = xo[..., S.oid]                                        # (N, 4, 2, 3, q)
    X1, Xm, Ss, Dd = comps(Xq.astype(complex))
    t1, t0 = targets(S); reps = S.reps
    out = {}
    for name, c, t in (('L1', X1, t1), ('L2', Xm, t0)):
        P = np.rint(paf_c(c).real).astype(np.int64)
        cb, eqs = gauss_logp(P[:, reps], t[reps].astype(np.int64)); out[name] = (round(cb, 1), eqs)
    for name, c in (('L3', Ss), ('L4', Dd)):
        P = paf_c(c) / 4.0                                     # = PAF of Z' combos (Z = -2 Z'), Eisenstein integers
        b = np.rint(2 * P.imag / math.sqrt(3)); a = np.rint(P.real + b / 2)
        cols = [a[:, reps]] + [b[:, [h]] for h in reps if not (b[:, h] == 0).all()]
        M_ = np.concatenate(cols, 1).astype(np.int64)
        from p571_ternary import gauss_cell
        cb, eqs = gauss_cell(M_, np.zeros(M_.shape[1]), len(reps) * math.log2(3))
        out[name] = (round(cb, 1), eqs)
    # lag-0 of each level is an exact (integer) constraint too: report its entropy separately
    return out

if __name__ == "__main__":
    if sys.argv[1] == 'check':
        for r in check(): print(json.dumps(r))
        print(json.dumps(dict(per_orbit_entropy=cond_entropy())))
    elif sys.argv[1] == 'est':
        a = sys.argv[2:]; opt = {"--N": "20000"}; pos = []; i = 0
        while i < len(a):
            if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
            else: pos.append(int(a[i])); i += 1
        S = Skel(pos[0], pos[1:])
        print(json.dumps(dict(m=S.m, q=S.q, R=S.R, sizes=S.sizes.tolist(), lag_reps=len(S.reps),
                              level_bits=level_bits(S, N=int(opt["--N"])))), flush=True)

# ---------------------------------------------------------------- stage-2 class structure for a given y
def sd_options(y0, y1):
    """per orbit pair: list of (S', D') options (in units of Z' = Z / -2: S' = Z'0 + Z'1, D' = Z'0 - Z'1) for
    y0, y1 in {+-1, +-3}; Z'_a = y_a w^j (j in Z3) if |y_a| = 1, 0 if |y_a| = 3"""
    za = [y0 * W ** j for j in range(3)] if abs(y0) == 1 else [0j]
    zb = [y1 * W ** j for j in range(3)] if abs(y1) == 1 else [0j]
    return [(a + b, a - b) for a in za for b in zb]

def key(z): return (round(z.real, 6), round(z.imag, 6))

def class_entropy(y):
    """y (4, 2, R) -> per block: log2 #S-classes (h3), log2 #configs (h3 + h4), and the mean log2 #D per S-class"""
    out = []
    for k in range(y.shape[0]):
        nS = 1.0; nAll = 1.0
        for j in range(y.shape[2]):
            opts = sd_options(int(y[k, 0, j]), int(y[k, 1, j]))
            nS *= len({key(s) for s, d in opts}); nAll *= len(opts)
        out.append((math.log2(nS), math.log2(nAll)))
    return out
