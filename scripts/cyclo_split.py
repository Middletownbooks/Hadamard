"""Stage split [E] for "generalised Paley" slices: G = Z_F x Z_q (gcd(F, q) = 1, q prime), multiplier group M fixing
Z_F pointwise and acting on Z_q by the index-s subgroup H of (Z/q)^*.  Per Z_q-coset g in Z_F and block k the slice has
1 + s free signs: alpha at z = 0 and c_j on the cyclotomic class C_j (|C_j| = f = (q - 1)/s); 32 configs for s = 4.

Stage 1 = the Z_q-coset sums y = alpha + f * sum_j c_j on Z_F.  They satisfy, exactly (b = 1):
   sum_k PAF_{Z_F}(y_k)(h) = T_J(h),   T_J(0) = 4m - 4(q - 1),  T_J(h != 0) = -4q,   row sums (+-2, 0, 0, 0).
log2 #(stage-1 solutions) = log2 #(y vectors) + log2 P_uniform-y(S1)   (Gaussian, row sums exact via 4 coordinates)
Given the whole-slice estimate log2E_total (theory/p571_msweep.py), the mean stage-2 yield per y is
   log2E_total - log2 #(stage-1 solutions).
The split is what decided 573 (stage 2 about 2^-12 per y: solved) versus 571 (2^-39 or worse: out of reach).
usage: cyclo_split.py F q s total_log2E [N]"""
import sys, math, itertools as it, json
import numpy as np

def main(F, q, s, total, N=200000, seed=1):
    f = (q - 1) // s; v = F * q + 1; m = F * q
    yvals = sorted({a + f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=s)})
    rng = np.random.default_rng(seed)
    Y = rng.choice(yvals, size=(N, 4, F))
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(Y, axis=-1)) ** 2, axis=-1).real).astype(np.int64).sum(1)   # (N, F)
    TJ = np.full(F, -4 * q, np.int64); TJ[0] = 4 * m - 4 * (q - 1)
    lags = list(range(0, F // 2 + 1))
    R = Y.sum(2)                                                   # row sums per block
    best = None; tot = 0.0
    # row-sum profile (2, 0, 0, 0) up to the 4 positions and sign; use block 0 = 2 and multiply by 8
    S = np.concatenate([P[:, lags], R], 1).astype(float)
    T = np.concatenate([TJ[lags], [2, 0, 0, 0]]).astype(float)
    mu = S.mean(0); C = np.cov(S.T); w, V = np.linalg.eigh(C); kw = w > 1e-9 * w.max()
    d = S - S[0]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:2000, j]).astype(np.int64)))) for j in range(S.shape[1])])
    pr = V.T @ (T - mu)
    lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum() - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
    ny = 4 * F * math.log2(len(yvals))
    s1 = ny + lp + math.log2(8)
    return dict(v=v, F=F, q=q, s=s, f=f, y_values=yvals, y_vector_bits=round(ny, 1), S1_constraint_bits=round(-lp, 1),
                stage1_log2=round(s1, 1), total_log2E=total, stage2_per_y=round(total - s1, 1))

if __name__ == "__main__":
    F, q, s = map(int, sys.argv[1:4]); total = float(sys.argv[4]); N = int(sys.argv[5]) if len(sys.argv) > 5 else 200000
    print(json.dumps(main(F, q, s, total, N)))
