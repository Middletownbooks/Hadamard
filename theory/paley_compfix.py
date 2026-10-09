r"""Gaussian estimate of the Paley-coset ansatz at v = 359 with the PUBLISHED Q-compression and w-gadget held fixed
(the border forces the compression; report section 5).  alpha is sampled uniformly among +-1 functions with the given
coset sums (number of -1's per H'-coset fixed), the gadget (positions, w, type) is the published one, u uniform +-1 on
twisted positions.  Reports log2 E = log2(#alpha configs) + #u bits - constraint bits, and the same split into
stage 1 (y: S1 only) and stage 2 (u given y).
usage: python3 theory/paley_compfix.py [N]"""
import sys, os, json, math, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_algebra import load_core, components
import paley_ls as LSm
from paley_ansatz import LABELS, LC

def gauss_logp(S, T, label):
    """S: (N, d) samples of an integer vector; T: target.  Returns -log2 density at T times lattice cell volume."""
    cols, seen = [], {}
    for j in range(S.shape[1]):
        kk = S[:3000, j].tobytes()
        if kk in seen:
            assert T[j] == T[seen[kk]], (label, j); continue
        if S[:, j].std() == 0:
            assert S[0, j] == T[j], (label, 'const coord mismatch', j, S[0, j], T[j]); continue
        seen[kk] = j; cols.append(j)
    S = S[:, cols].astype(float); T = T[cols].astype(float)
    mu = S.mean(0); C = np.cov(S.T); w, V = np.linalg.eigh(np.atleast_2d(C)); kw = w > 1e-6 * w.max()
    pr = V.T @ (T - mu)
    d = S - S[0]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j]).astype(np.int64)))) for j in range(S.shape[1])])
    lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum() - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
    return -lp, int(kw.sum()), len(cols)

def main(N=200000, seed=1):
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt'))
    p = 11; n = 32
    P = LSm.Prob(359, 7, 11, (2, 2, 8), [(0, 0, 2)])
    comps = [components(x, p) for x in X]
    A0 = [c[0].reshape(-1) for c in comps]; W0 = [c[2].reshape(-1) for c in comps]
    co = P.G.co; cos = co[:, 0] * 4 + co[:, 1] * 2 + (co[:, 2] % 2)          # H'-coset index
    rng = np.random.default_rng(seed)
    labs = np.zeros((N, 4, n), np.int64); nalpha_bits = 0.0
    for k in range(4):
        w = W0[k]; a0 = A0[k]
        for q in range(8):
            pos = np.nonzero(cos == q)[0]
            nminus = int((a0[pos] == -1).sum())
            nalpha_bits += math.log2(math.comb(len(pos), nminus))
            # random alpha on this coset with nminus -1's
            r = np.argsort(rng.random((N, len(pos))), axis=1)[:, :nminus]
            alpha = np.ones((N, len(pos)), np.int64)
            np.put_along_axis(alpha, r, -1, axis=1)
            for j, g in enumerate(pos):
                if w[g] != 0:
                    labs[:, k, g] = [LABELS.index((int(a), 0, int(w[g]))) for a in alpha[:, j]]
                else:
                    u = rng.choice(np.array([1, -1]), N)
                    labs[:, k, g] = [LABELS.index((int(a), 1, int(s))) for a, s in zip(alpha[:, j], u)]
    nu_bits = sum(int((W0[k] == 0).sum()) for k in range(4))
    # reduced PAF of each sample (sum over blocks)
    R = np.zeros((N, n, 3), np.int64)
    for k in range(4):
        Ak = LC[labs[:, k]].astype(np.int64)
        for h in range(n):
            B = Ak[:, P.add[:, h]]
            R[:, h, :] += np.einsum('ngi,ngj,ijc->nc', Ak, B, P.C)
    T = P.T
    # full system: coordinates (h, c) for all h, c except (0,0)
    S_full = R.reshape(N, -1)[:, 1:]; T_full = T.reshape(-1)[1:]
    cb_full, eq_full, _ = gauss_logp(S_full, T_full, 'full')
    # stage 1 only: R0 + (p-1) * Y  where Y = (R_QR + R_NQR)/2
    Y = (R[:, :, 1] + R[:, :, 2]) // 2
    S1 = (R[:, :, 0] + (p - 1) * Y)[:, 1:]; T1 = (T[:, 0] + (p - 1) * T[:, 1])[1:]
    cb1, eq1, _ = gauss_logp(S1, T1, 'stage1')
    out = dict(alpha_bits=round(nalpha_bits, 1), u_bits=nu_bits, free_bits=round(nalpha_bits + nu_bits, 1),
               constraint_bits_full=round(cb_full, 1), eqs_full=eq_full, log2E_full=round(nalpha_bits + nu_bits - cb_full, 1),
               constraint_bits_stage1=round(cb1, 1), eqs_stage1=eq1, log2E_stage1=round(nalpha_bits - cb1, 1),
               implied_log2E_stage2_given_y=round(nu_bits - (cb_full - cb1), 1))
    print(json.dumps(out))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200000)
