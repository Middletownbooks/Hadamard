r"""Is the Gaussian tail of the D-level right?  (report/paley_523_answer.md, sections 3.4 and 6) [M against E].

The D-success per S-solution (theory/q523_dsuccess_real.py) is (D-candidates per S-solution) x P(full D-level match |
D3 match), the second factor from the lattice-exact Gaussian (theory/q523_gauss.py) fitted to the residuals of sampled
candidates: about 2^-37.2 at 523, far below anything that can be counted.  This script checks the same Gaussian where
counting is possible: for coordinate subsets of the D-level signature (lags of the Z9 sequence, s coordinates each) it
compares the exact fraction of candidates that hit the target on the subset with the Gaussian prediction for that
subset, on every D-candidate of the given production S-solutions.
usage: python3 theory/q523_dtail_check.py F q s y.json P4.npy [--n 1000]"""
import sys, os, json, math, time, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_dcheck import DChecker, hjoin4
from q523_gauss import gauss_bits
from q523_sdsolve import SDInst

def all_candidate_sums(dc, pidx4, cap=200000):
    """full D-level signature sums of every D-candidate (configurations inside matched D3 groups) of one S-solution"""
    B = [dc.block(k, pidx4[k]) for k in range(4)]
    h3, n3 = hjoin4(B[0]['k3'], B[1]['k3'], B[2]['k3'], B[3]['k3'], dc.k3T, 16, 1 << 20)
    out = []
    for q in h3:
        if not (B[0]['sig3'][q[0]] + B[1]['sig3'][q[1]] + B[2]['sig3'][q[2]] + B[3]['sig3'][q[3]] == dc.T3).all(): continue
        idx = [B[k]['order'][B[k]['start'][q[k]]:B[k]['start'][q[k] + 1]] for k in range(4)]
        S = [B[k]['sigf'][idx[k]] for k in range(4)]
        tot = (S[0][:, None, None, None, :] + S[1][None, :, None, None, :] + S[2][None, None, :, None, :]
               + S[3][None, None, None, :, :]).reshape(-1, S[0].shape[1])
        out.append(tot)
        if sum(len(x) for x in out) > cap: break
    return np.concatenate(out) if out else np.zeros((0, dc.TF.shape[0]), np.int64)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--n": "1000"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y); P4 = np.load(pos[4], mmap_mode='r'); n = min(int(opt["--n"]), len(P4))
    dc = DChecker(inst, 1); t0 = time.time(); V = []
    for j in range(n):
        V.append(all_candidate_sums(dc, [np.asarray(P4[j, k], np.int64) for k in range(4)]))
    V = np.concatenate(V); TF = dc.TF; L = TF.shape[0]; nl = L // s
    R = V - TF[None, :]
    print(json.dumps(dict(stage="samples", s_solutions=n, candidates=len(V), log2=round(math.log2(len(V)), 2), lags=nl,
                          secs=round(time.time() - t0, 1))), flush=True)
    rng = np.random.default_rng(5); fit = V[rng.choice(len(V), size=min(len(V), 200000), replace=False)]
    subsets = [[0], [1], [2], [0, 1], [0, 3], [1, 2], [0, 1, 2], [0, 1, 3], [0, 1, 2, 3], list(range(nl))]
    for lg in subsets:
        cols = [l * s + t for l in lg for t in range(s)]
        hit = int((R[:, cols] == 0).all(1).sum())
        g = gauss_bits(fit[:, cols], TF[cols])
        emp = -math.log2(hit / len(V)) if hit else None
        print(json.dumps(dict(lags=lg, hits=hit, empirical_bits=round(emp, 2) if emp is not None else None,
                              gaussian_bits=round(g['bits'], 2), mahal=round(g.get('mahal', float('nan')), 2),
                              in_lattice=g.get('in_lattice'))), flush=True)
