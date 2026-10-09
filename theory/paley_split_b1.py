r"""Stage split of the b = 1 Paley-coset ansatz over a cyclic G' = Z_n (p = 3 mod 4), as a Gaussian estimate [E].

For each (K, N-, N+) allowed by the norm identity (report section 2), sample random cores of the ansatz with the exact
row sums (2, 0, 0, 0): K constant cosets placed uniformly over the 4n cosets, N- of them J' (w = -alpha), w signs
uniform, alpha uniform on the twisted cosets with the required block sums, u uniform.  Report
  * stage-1 entropy  H1 = log2 #(gadget placements, J/J' choice, w signs, alpha with the right sums),
  * stage-1 bits     C1 = -log2 Pr[sum PAF(y)(h) = T_J(h), h = 1..n/2]          (Gaussian, lattice-corrected),
  * stage-2 entropy  H2 = #twisted cosets (u bits),
  * stage-2 bits     C2 = -log2 Pr[S2-real(h) = 0, h = 1..n/2;  X(h) = 0, h = 1..n/2-1 | y]   (X(n/2) = 0 always),
so that log2 E(stage 1) = H1 - C1, log2 E(stage 2 | y) = H2 - C2 and log2 E(core) = H1 + H2 - C1 - C2.
Same density estimator as theory/paley_compfix.py (calibrated on 359).
usage: python3 theory/paley_split_b1.py v p n [N] [Kmax_rows]"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_compfix import gauss_logp
from paley_wscreen import KN_options

def paf(x):            # rows, cyclic PAF via FFT (exact for |entries| small)
    F = np.fft.fft(x, axis=-1)
    return np.rint(np.fft.ifft(np.abs(F) ** 2, axis=-1).real).astype(np.int64)

def cc(f, g):          # CC(f, g)(h) = sum_x f(x) g(x + h), rows
    return np.rint(np.fft.ifft(np.conj(np.fft.fft(f, axis=-1)) * np.fft.fft(g, axis=-1), axis=-1).real).astype(np.int64)

def sample_stage1(rng, N, n, p, K, Nm):
    """returns alpha, w (N, 4, n) and the log2 count of alpha choices per sample"""
    rs = np.array([2, 0, 0, 0])
    A = np.zeros((N, 4, n), np.int64); W = np.zeros((N, 4, n), np.int64); la = np.zeros(N); ok = np.zeros(N, bool)
    for i in range(N):
        pos = rng.choice(4 * n, K, replace=False)
        jm = np.zeros(4 * n, bool); jm[pos[:Nm]] = True                     # first Nm are J'
        w = np.zeros(4 * n, np.int64); w[pos] = rng.choice((1, -1), K)
        a = np.zeros(4 * n, np.int64); a[pos] = np.where(jm[pos], -w[pos], w[pos])
        w = w.reshape(4, n); a = a.reshape(4, n); good = True; lg = 0.0
        for k in range(4):
            free = np.nonzero(w[k] == 0)[0]; mfree = len(free)
            need = rs[k] - (p - 1) * int(w[k].sum()) - int(a[k].sum())
            if abs(need) > mfree or (mfree - need) % 2: good = False; break
            nminus = (mfree - need) // 2
            lg += math.log2(math.comb(mfree, nminus))
            vals = np.ones(mfree, np.int64); vals[rng.choice(mfree, nminus, replace=False)] = -1
            a[k, free] = vals
        if good: A[i], W[i], la[i], ok[i] = a, w, lg, True
    return A[ok], W[ok], la[ok]

def main(v, p, n, N=20000, seed=1):
    b = 1; m = n * p; assert m == v - 1
    rng = np.random.default_rng(seed)
    TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m - 4 * (p - 1)
    half = n // 2
    out = []
    for (K, Nm, Np) in KN_options(v, b, p, n):
        A, W, la = sample_stage1(rng, N, n, p, K, Nm)
        if len(A) < 1000: out.append(dict(K=K, Nm=Nm, note='row sums rarely feasible')); continue
        Y = A + (p - 1) * W
        S1 = paf(Y).sum(1)[:, 1:half + 1]
        c1, r1, d1 = gauss_logp(S1, TJ[1:half + 1], 'S1')
        # row-sum feasibility fraction (the sampler discards infeasible gadgets)
        feas = len(A) / N
        H1 = (math.log2(math.comb(4 * n, K)) + math.log2(math.comb(K, Nm)) + K + math.log2(feas) + float(np.mean(la)))
        # stage 2 given one sampled y (statistics depend weakly on y)
        a0, w0 = A[0], W[0]; beta = a0 - w0; tw = (w0 == 0)
        U = np.where(tw[None], rng.choice((1, -1), (N, 4, n)), 0)
        S2 = (paf(beta[None]).sum(1) + p * paf(U).sum(1))[:, 1:half + 1]
        X = (cc(U, np.broadcast_to(beta, U.shape)) - cc(np.broadcast_to(beta, U.shape), U)).sum(1)[:, 1:half]
        # S2 target: 0 (the y-dependent shift sum PAF(beta) is included in S2 itself)
        c2, r2, d2 = gauss_logp(np.concatenate([S2, X], 1), np.zeros(S2.shape[1] + X.shape[1], np.int64), 'S2X')
        H2 = int(tw.sum())
        rec = dict(v=v, p=p, n=n, K=K, Nminus=Nm, Nplus=Np, H1=round(H1, 1), C1=round(c1, 1), stage1_log2E=round(H1 - c1, 1),
                   H2=H2, C2=round(c2, 1), stage2_log2E_per_y=round(H2 - c2, 1), core_log2E=round(H1 + H2 - c1 - c2, 1),
                   eqs1=r1, eqs2=r2, samples=len(A))
        print(json.dumps(rec), flush=True); out.append(rec)
    return out

if __name__ == "__main__":
    v, p, n = map(int, sys.argv[1:4]); N = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    main(v, p, n, N)
