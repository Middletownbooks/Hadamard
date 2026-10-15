r"""Conditional level budgets in the index-s, U = <2> slice of v = 7q [M fibres, E Gaussians].

Given a configuration X (4 blocks, drawn uniformly given y and the dual energy by theory/b917_levels.LevelModel):
  cross | w   the dual data w (types, x0 row, sym signs) fixed; the fibre is the anti-column signs e with
              sum_anti e = (sigma1 - sigma2)/2 per block: F = sum_k log2 C(n_a, (n_a + d)/2) exactly.
              Lattice Gaussian (theory/q523_gauss.gauss_bits) of the cross lag classes over uniform fibre samples,
              evaluated at the typical target (log-det part, ab_gauss.gauss_bits_at_mean): C(cross | w).
              Split into the symmetric (C+) and antisymmetric (C-) parts of the cross PAF.
  dual | u    u fixed (all but the signs of |w| = 7 columns); fibre = those signs with fixed sums.
The deficit of a dual-first order is C(cross | w) - F per dual solution (2^deficit dual solutions per core).
usage: python3 theory/b917_cond.py [--q 131 --s 26] [--ys 0] [--n 6] [--N 4000]"""
import sys, os, math, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b917_slice import Slice917
from b917_levels import LevelModel, stage1_ys, lag_class_reps, alpha_sigma
from ab_gauss import gauss_bits_at_mean

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cross_coords(S, Xs_batch, reps_cross):
    """sum_k PAF_G at the cross lag classes (1, gamma^j), j < s, for a batch (n, 4, 3, s+1)"""
    x = S.expand(Xs_batch).astype(np.float64)
    F = np.fft.fft2(x)
    P = np.rint(np.fft.ifft2((np.abs(F) ** 2).sum(1)).real).astype(np.int64)
    hh = np.array([r[0] for r in reps_cross]); zz = np.array([r[1] for r in reps_cross])
    return P[:, hh, zz]


def dual_coords(S, Xs_batch, reps_dual):
    x = S.expand(Xs_batch).astype(np.float64)
    F = np.fft.fft2(x)
    P = np.rint(np.fft.ifft2((np.abs(F) ** 2).sum(1)).real).astype(np.int64)
    hh = np.array([r[0] for r in reps_dual]); zz = np.array([r[1] for r in reps_dual])
    return P[:, hh, zz]


def fibre_e(X, rng, n):
    """n resamples of the anti-column signs e (x1 = -x2 = e) of block X (3, s+1), keeping sum_anti e"""
    cols = X[:, 1:]
    anti = np.nonzero(cols[1] == -cols[2])[0]
    e = cols[1, anti].copy()
    out = np.repeat(X[None], n, 0)
    for t in range(n):
        p = rng.permutation(len(e))
        out[t, 1, 1 + anti] = e[p]; out[t, 2, 1 + anti] = -e[p]
    na = len(anti); d = int(e.sum())
    return out, (math.log2(math.comb(na, (na + d) // 2)) if na else 0.0), na


def fibre_tau(X, rng, n):
    """n resamples of the signs of the |w| = 7 columns (+-(1,1,1)) keeping their sum"""
    cols = X[:, 1:]
    s7 = np.nonzero((cols[0] == cols[1]) & (cols[1] == cols[2]))[0]
    t = cols[0, s7].copy()
    out = np.repeat(X[None], n, 0)
    for i in range(n):
        p = rng.permutation(len(t))
        out[i, :, 1 + s7] = t[p][:, None]
    n7 = len(s7); d = int(t.sum())
    return out, (math.log2(math.comb(n7, (n7 + d) // 2)) if n7 else 0.0), n7


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    q = int(opt('--q', '131')); s = int(opt('--s', '26')); N = int(opt('--N', '4000')); nconf = int(opt('--n', '6'))
    S = Slice917(q, s)
    Ys = stage1_ys(q, s)
    reps = lag_class_reps(S); h = s // 2
    reps_dual = reps[1:1 + h]; reps_cross = reps[2 + h:]
    out = os.path.join(ROOT, 'logs/b917/cond.jsonl')
    for yi in [int(x) for x in opt('--ys', '0').split(',')]:
        lm = LevelModel(S, Ys[yi])
        rng = np.random.default_rng(100 + yi)
        draw = lm.sampler(rng)
        for c in range(nconf):
            X = draw()
            # cross | w: resample e in every block independently
            samples = np.zeros((N, 4, 3, s + 1), np.int64); F = 0.0; nas = []
            for k in range(4):
                fb, fk, na = fibre_e(X[k], rng, N); samples[:, k] = fb; F += fk; nas.append(na)
            V = cross_coords(S, samples, reps_cross)
            Vp = V[:, :h] + V[:, h:]; Vm = V[:, :h] - V[:, h:]
            g = gauss_bits_at_mean(V); gp = gauss_bits_at_mean(Vp); gm = gauss_bits_at_mean(Vm)
            # dual | u
            samples2 = np.zeros((N, 4, 3, s + 1), np.int64); F2 = 0.0; n7s = []
            for k in range(4):
                fb, fk, n7 = fibre_tau(X[k], rng, N); samples2[:, k] = fb; F2 += fk; n7s.append(n7)
            Vd = dual_coords(S, samples2, reps_dual)
            gd = gauss_bits_at_mean(Vd) if Vd.std(0).max() > 0 else dict(bits_at_mean=0.0, rank=0)
            rec = dict(q=q, s=s, y_index=yi, conf=c, n_anti=nas, F_e=round(F, 2), C_cross_given_w=g.get('bits_at_mean'),
                       rank=g.get('rank'), C_plus=gp.get('bits_at_mean'), C_minus=gm.get('bits_at_mean'),
                       deficit=round(g['bits_at_mean'] - F, 2) if g.get('bits_at_mean') is not None else None,
                       n7=n7s, F_tau=round(F2, 2), C_dual_given_u=gd.get('bits_at_mean'), rank_dual_u=gd.get('rank'),
                       std_cross=round(float(V.std(0).mean()), 1))
            with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
