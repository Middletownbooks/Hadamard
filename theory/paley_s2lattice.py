r"""Why do b = 1, p = 11 stage-1 solutions never complete?  Lattice (congruence) test of stage 2 [M].

For every exact stage-1 solution (alpha, w) of theory/paley_b1pipe.py, stage 2 asks for u in {+-1} on the twisted cosets
with  sum_k PAF(u_k)(h) = tau(h) := -sum_k PAF(beta_k)(h) / p   and   X(h) = sum_k [CC(beta_k,u_k) - CC(u_k,beta_k)](h) = 0
for h != 0.  The vector F(u) = (sum PAF(u)(h), X(h))_{h != 0} lies in  F(u0) + L, where L is the lattice spanned by the
differences F(u) - F(u0).  If (tau, 0) - F(u0) is not in L, stage 2 is empty for that y by a congruence alone.
L is computed exactly from random samples (integer row reduction; more samples never enlarge it once saturated), and
membership is decided exactly.  Also reported: the Gaussian log2 E of stage 2 for this y (same estimator as
theory/paley_split_b1.py) and the coordinates whose residue class excludes the target.
usage: python3 theory/paley_s2lattice.py v p n [samples]"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import paley_b1pipe as BP
from paley_compfix import gauss_logp

def paf(x):
    F = np.fft.fft(x, axis=-1)
    return np.rint(np.fft.ifft(np.abs(F) ** 2, axis=-1).real).astype(np.int64)

def cc(f, g):
    return np.rint(np.fft.ifft(np.conj(np.fft.fft(f, axis=-1)) * np.fft.fft(g, axis=-1), axis=-1).real).astype(np.int64)

def hnf_rows(V):
    """row-style Hermite-type reduction of the integer row lattice spanned by V (exact, Python ints)"""
    rows = [list(map(int, r)) for r in V if any(r)]
    d = len(rows[0]) if rows else 0
    basis = []; col = 0
    while rows and col < d:
        piv = [r for r in rows if r[col] != 0]
        if not piv: col += 1; continue
        rest = [r for r in rows if r[col] == 0]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[col]))
            p0 = piv[0]; new = [p0]
            for r in piv[1:]:
                q = r[col] // p0[col]; r2 = [a - q * b for a, b in zip(r, p0)]
                (new if r2[col] != 0 else rest).append(r2) if any(r2) else None
            piv = new
        p0 = piv[0]
        if p0[col] < 0: p0 = [-a for a in p0]
        basis.append((col, p0)); rows = rest; col += 1
    return basis

def in_lattice(basis, t):
    t = list(map(int, t))
    for col, b in basis:
        if t[col] % b[col]: return False, col
        q = t[col] // b[col]; t = [a - q * c for a, c in zip(t, b)]
    return (not any(t)), (None if not any(t) else next(i for i, a in enumerate(t) if a))

def analyse(A, W, p, n, S=4000, seed=0):
    rng = np.random.default_rng(seed)
    beta = A - W; tw = (W == 0)
    PB = paf(beta).sum(0)
    if (PB[1:] % p).any(): return dict(error='sum PAF(beta) not divisible by p')
    tau = -PB[1:] // p
    def F(U):
        s = paf(U).sum(1)[:, 1:]
        X = (cc(np.broadcast_to(beta, U.shape), U) - cc(U, np.broadcast_to(beta, U.shape))).sum(1)[:, 1:]
        return np.concatenate([s, X], 1)
    U = np.where(tw[None], rng.choice((1, -1), (S, 4, n)), 0)
    Fs = F(U); base = Fs[0]; D = Fs[1:] - base
    basis = hnf_rows(D[:600])
    # saturation check: adding more differences must not change membership verdicts
    basis2 = hnf_rows(np.concatenate([D[:600], D[600:2400]]))
    target = np.concatenate([tau, np.zeros(n - 1, np.int64)])
    ok1, c1 = in_lattice(basis, target - base); ok2, c2 = in_lattice(basis2, target - base)
    # per-coordinate residues
    steps = [int(np.gcd.reduce(np.abs(D[:, j]))) for j in range(D.shape[1])]
    bad = [j for j in range(D.shape[1]) if steps[j] and (target[j] - base[j]) % steps[j]]
    names = [f'PAFu({h})' for h in range(1, n)] + [f'X({h})' for h in range(1, n)]
    lp, r, d = gauss_logp(Fs[:, list(range(n // 2)) + list(range(n - 1, n - 1 + n // 2 - 1))],
                          np.concatenate([tau[:n // 2], np.zeros(n // 2 - 1, np.int64)]), 'S2')
    return dict(u_bits=int(tw.sum()), lattice_rank=len(basis), target_in_lattice=bool(ok1), stable=bool(ok1 == ok2),
                first_bad_coord=(names[c1] if c1 is not None else None), per_coord_bad=[names[j] for j in bad][:8],
                steps_PAFu=steps[:n - 1], steps_X=steps[n - 1:], gauss_stage2_log2E=round(int(tw.sum()) - lp, 1))

def main(v, p, n, S=4000):
    rec = []
    BP.run(v, p, n, max_cores=10 ** 9, record=rec, lemma6=False, tag='_nolemma6')
    out = []
    for i, (A, W, nU) in enumerate(rec):
        r = analyse(A, W, p, n, S)
        r.update(v=v, p=p, n=n, y_index=i, stage2_exact_solutions=nU, K=int((W != 0).sum()),
                 t=W.sum(1).tolist(), K_k=(W != 0).sum(1).tolist(), lemma6_trivial_char=bool(
                     (4 * v - int((((A + (p - 1) * W).sum(1) - p * W.sum(1)) ** 2).sum())) % p == 0 and
                     BP.sq4((4 * v - int((((A + (p - 1) * W).sum(1) - p * W.sum(1)) ** 2).sum())) // p,
                            [int((W[k] == 0).sum()) for k in range(4)])),
                 lemma6_real_chars=bool(BP.lemma6_real_chars(A, W, p, n)))
        print(json.dumps(r), flush=True); out.append(r)
    return out

if __name__ == "__main__":
    v, p, n = map(int, sys.argv[1:4]); S = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
    main(v, p, n, S)
