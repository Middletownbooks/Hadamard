r"""Stage-2 Gaussian for ANY multiplier slice given a stage-1 solution y on the compressed group Z_F [E].

Slice: x_k constant on the orbits of M = <mu> on Z_m; compression along the Hall factor Z_K onto Z_F (F = m/K);
y_k(b) = sum_{g = b mod F} x_k(g) is fixed (a verified stage-1 solution).  Stage 2 = the signs of the M-orbits given y.
 * Freedom H2 (exact): per block and per M-orbit of Z_F (representative b), the number of sign patterns of the
   M_b-orbits partitioning the fibre over b with signed size-sum y(b).
 * Constraint: the full core equations sum_k PAF(x_k)(h) = -4 (h != 0) at one lag per <M, -1>-orbit of Z_m, and 4m at
   h = 0.  Configurations are drawn uniformly given y (independently per fibre-orbit and block), their summed PAF
   vectors V are computed by FFT, and the lattice-coordinate Gaussian of theory/q523_gauss.py (exact difference
   lattice, so the stage-1 implied combinations drop out as rank deficiency) gives -log2 P(V = T).
 * Output: log2 (expected stage-2 completions of this y) = H2 - bits.  As in the 523 work this ignores nothing exact
   except congruences outside the sampled lattice (in_lattice false => an exact obstruction to be proved separately).
usage: python3 theory/nd_stage2_est.py m mu K Yfile [--idx i] [--N 40000] [--seed 1]
       (Yfile: json with "solutions": [{"y": [[...] x4]}] as written by theory/nd_s1_sd.py --out, or {"y": ...})"""
import sys, os, json, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
from nd_s1_sd import mgroup, lag_reps

class Slice:
    def __init__(self, m, mu, K):
        self.m, self.mu, self.K = m, mu, K; F = m // K; self.F = F
        Mg = mgroup(m, mu); self.Mg = Mg
        orb = -np.ones(m, np.int64); reps = []
        for g in range(m):
            if orb[g] >= 0: continue
            i = len(reps)
            for a in Mg: orb[g * a % m] = i
            reps.append(g)
        self.orb = orb; self.nor = len(reps)
        Fo = -np.ones(F, np.int64); freps = []
        for b in range(F):
            if Fo[b] >= 0: continue
            j = len(freps)
            for a in Mg: Fo[b * a % F] = j
            freps.append(b)
        self.Fo, self.freps = Fo, freps
        # parts of each fibre representative: the M-orbits meeting it, with their sizes inside the fibre
        self.parts = []
        for b in freps:
            P = {}
            for g in range(b, m, F): P[int(orb[g])] = P.get(int(orb[g]), 0) + 1
            self.parts.append(sorted(P.items()))
        self.lags = lag_reps(m, Mg)                     # starts with 0
        self.T = np.array([4 * m] + [-4] * (len(self.lags) - 1), np.int64)

    def patterns(self, j, target):
        """all sign vectors (one per part) with sum sign*size = target"""
        sizes = [s for _, s in self.parts[j]]
        out = []
        if len(sizes) > 22: raise ValueError("too many parts for enumeration")
        for msk in range(1 << len(sizes)):
            t = sum(s if (msk >> i) & 1 else -s for i, s in enumerate(sizes))
            if t == target: out.append([1 if (msk >> i) & 1 else -1 for i in range(len(sizes))])
        return np.array(out, np.int8).reshape(-1, len(sizes))

    def estimate(self, Y, N=40000, seed=1, chunk=2000):
        rng = np.random.default_rng(seed); Y = np.array(Y)
        pats = [[self.patterns(j, int(Y[k][b])) for j, b in enumerate(self.freps)] for k in range(4)]
        H2 = sum(math.log2(len(p)) for pk in pats for p in pk if len(p))
        if any(len(p) == 0 for pk in pats for p in pk): return dict(H2=None, note="a coset sum is not realisable")
        V = np.zeros((N, len(self.lags)), np.int64); done = 0
        lags = np.array(self.lags)
        while done < N:
            c = min(chunk, N - done)
            acc = np.zeros((c, self.m))
            for k in range(4):
                lab = np.zeros((c, self.nor), np.int8)
                for j in range(len(self.freps)):
                    P = pats[k][j]; pick = P[rng.integers(len(P), size=c)]
                    for i, (o, s) in enumerate(self.parts[j]): lab[:, o] = pick[:, i]
                x = lab[:, self.orb].astype(np.float64)
                fx = np.fft.rfft(x, axis=1); acc += np.fft.irfft(fx * fx.conj(), n=self.m, axis=1)
            V[done:done + c] = np.rint(acc[:, lags]).astype(np.int64); done += c
        gb = gauss_bits(V, self.T)
        return dict(H2=round(H2, 2), lag_classes=len(self.lags), **gb,
                    stage2_log2E=round(H2 - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--idx": "-1", "--N": "40000", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    m, mu, K = map(int, pos[:3]); d = json.load(open(pos[3]))
    sols = d['solutions'] if 'solutions' in d else [d]
    S = Slice(m, mu, K)
    idxs = range(len(sols)) if int(opt["--idx"]) < 0 else [int(opt["--idx"])]
    for ix in idxs:
        r = S.estimate(sols[ix]['y'], int(opt["--N"]), int(opt["--seed"]))
        print(json.dumps(dict(m=m, mu=mu, K=K, idx=ix, **r)), flush=True)
