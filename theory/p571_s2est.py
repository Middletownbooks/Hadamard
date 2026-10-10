r"""Per-y stage-2 estimate for the hierarchical solver (theory/p571_hsolve.py) [E], on REAL stage-1 solutions.

For each y: per block the exact S-class count 2^{h3_k} (product of per-orbit S' options) and configuration count
3^{n_free,k}; Gaussian constraint bits of L3 (S' = Z'0 + Z'1) and L4 (D' = Z'0 - Z'1) from uniform z* given y.
Lattice [P]: S' = X1, D' = Xm (mod 1 - w) on every orbit (Z' = y mod (1 - w), and 3 = 0 mod (1 - w)), so each lag's
value P(h) in Z[w] lies in a residue class of index 3 that S1 already fixes correctly; cell volume 3 per lag
(lag 0 included).  Predictions:
  S-solutions per y   = 2^{h3 - C3}          (h3 = sum_k h3_k)
  D-survival per S-sol = 2^{h4 - C4}          (h4 = H2 - h3, average)
  yield per y          = 2^{H2 - C3 - C4}
  S-join pairs         = 2^{(h3_a + h3_b)} per side (balanced pairing)
usage: python3 theory/p571_s2est.py y.npz m gen [--ny 20] [--N 20000]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst, Stage2Block, _sig_eis1, EIS
from p571_ternary import gauss_cell

def level_bits_y(inst, y, N=20000, seed=0):
    rng = np.random.default_rng(seed); R = inst.R
    SA = np.zeros((4, N, R), np.int64); SB = np.zeros((4, N, R), np.int64)
    DA = np.zeros((4, N, R), np.int64); DB = np.zeros((4, N, R), np.int64)
    for k in range(4):
        for i in range(R):
            z = []
            for a in range(2):
                ya = int(y[k, a, i])
                if abs(ya) == 3: z.append((np.zeros(N, np.int64), np.zeros(N, np.int64)))
                else:
                    j = rng.integers(0, 3, N); z.append((ya * EIS[j, 0], ya * EIS[j, 1]))
            SA[k, :, i] = z[0][0] + z[1][0]; SB[k, :, i] = z[0][1] + z[1][1]
            DA[k, :, i] = z[0][0] - z[1][0]; DB[k, :, i] = z[0][1] - z[1][1]
    out = {}
    for name, A, B in (('C3', SA, SB), ('C4', DA, DB)):
        sig = sum(_sig_eis1(A[k], B[k], inst.N) for k in range(4))
        H = len(inst.lags)
        cb, eqs = gauss_cell(sig, np.array([inst.v] + [0] * (sig.shape[1] - 1)), H * math.log2(3))
        out[name] = (round(cb, 1), eqs)
    return out

def estimate_y(inst, y, N=20000, seed=0):
    blocks = [Stage2Block(inst, y[k]) for k in range(4)]
    h3k = [b.log2n for b in blocks]
    nfree = [int((np.abs(y[k]) == 1).sum()) for k in range(4)]
    H2k = [n * math.log2(3) for n in nfree]
    lb = level_bits_y(inst, y, N=N, seed=seed)
    h3 = sum(h3k); H2 = sum(H2k); h4 = H2 - h3; C3 = lb['C3'][0]; C4 = lb['C4'][0]
    srt = sorted(h3k)
    pairs = max(srt[3] + srt[0], srt[2] + srt[1])
    both_free = int(sum(((np.abs(y[k, 0]) == 1) & (np.abs(y[k, 1]) == 1)).sum() for k in range(4)))
    return dict(nfree=sum(nfree), both_free_pairs=both_free, H2=round(H2, 1), h3=round(h3, 1), h4=round(h4, 1),
                h3_blocks=[round(x, 1) for x in h3k], C3=C3, C4=C4,
                s_solutions_log2=round(h3 - C3, 1), d_survival_log2=round(h4 - C4, 1),
                yield_log2=round(H2 - C3 - C4, 1), s_join_pairs_log2=round(pairs, 1))

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--ny": "20", "--N": "20000"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    Y = np.load(pos[0])['Y']; inst = Inst(int(pos[1]), [int(g) for g in pos[2:]])
    for q in range(min(len(Y), int(opt['--ny']))):
        print(json.dumps(dict(y=q, **estimate_y(inst, Y[q], N=int(opt['--N']), seed=q))), flush=True)
