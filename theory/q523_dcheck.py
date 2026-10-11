r"""Fast exact D-check for one S-solution (numba), used by theory/q523_wagner.py (report/paley_523_answer.md, section 3).

Given S (one S-value per pair, every block), each pair has its D-options; D runs over the product of the options.
The D-level condition is sum_k P(D_k) = T_D on Z_F2 (lags 0..F2/2, eta-coordinates).  Exact algorithm:
  1. per block enumerate all D configurations, their D3 = row sums (rows g = r mod 3) and full D-signatures;
  2. group by D3; D3 class split: 2+2 join of the distinct D3 values on the Z3-compressed D-level signature;
  3. per D3-solution, 2+2 join of the configurations in the four groups on the full D-level signature.
All hashes are linear 64-bit hashes of the integer signatures; every reported solution is re-checked exactly.
Validated against theory/q523_sdsolve.d_check (same solutions) by `python3 theory/q523_dcheck.py test`."""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_sdsolve import seq_sig

@numba.njit(cache=True)
def enum_opts(opts, m):
    """opts (G, M, s), m (G,) -> all configurations (n, G, s) (mixed radix, last pair fastest)"""
    G = m.shape[0]; s = opts.shape[2]; n = 1
    for g in range(G): n *= m[g]
    out = np.zeros((n, G, s), np.int64)
    for l in range(n):
        rem = l
        for g in range(G - 1, -1, -1):
            t = rem % m[g]; rem //= m[g]
            for j in range(s): out[l, g, j] = opts[g, t, j]
    return out

@numba.njit(cache=True)
def hjoin4(kA, kB, kC, kD, kT, tbits, maxhits):
    """all (a, b, c, d) with kA[a] + kB[b] + kC[c] + kD[d] == kT (uint64), chained table of AB keys"""
    nab = kA.shape[0] * kB.shape[0]
    size = 1 << tbits
    while size < 2 * nab + 2: size <<= 1
    mask = size - 1
    tk = np.zeros(size, np.uint64); tv = np.full(size, -1, np.int64)
    for a in range(kA.shape[0]):
        for b in range(kB.shape[0]):
            key = kA[a] + kB[b]; slot = int(key >> np.uint64(20)) & mask
            while tv[slot] >= 0: slot = (slot + 1) & mask
            tk[slot] = key; tv[slot] = a * kB.shape[0] + b
    hits = np.zeros((maxhits, 4), np.int64); nh = 0
    for c in range(kC.shape[0]):
        for d in range(kD.shape[0]):
            need = kT - kC[c] - kD[d]; slot = int(need >> np.uint64(20)) & mask
            while tv[slot] >= 0:
                if tk[slot] == need:
                    if nh < maxhits:
                        hits[nh, 0] = tv[slot] // kB.shape[0]; hits[nh, 1] = tv[slot] % kB.shape[0]
                        hits[nh, 2] = c; hits[nh, 3] = d
                    nh += 1
                slot = (slot + 1) & mask
    return hits[:min(nh, maxhits)], nh

class DChecker:
    def __init__(self, inst, seed=7):
        self.inst = inst; rng = np.random.default_rng(seed); s = inst.s
        self.w3 = rng.integers(1, 2 ** 63, 2 * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.wf = rng.integers(1, 2 ** 63, inst.nl * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.wd = rng.integers(1, 2 ** 63, 3 * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.T3 = inst.tsigD3(); self.TF = inst.tsigD()
        self.k3T = (self.T3.astype(np.uint64) * self.w3).sum(); self.kfT = (self.TF.astype(np.uint64) * self.wf).sum()

    def block(self, k, pidx):
        inst = self.inst; F2 = inst.F2; s = inst.s
        m = np.array([len(inst.pairs[k][g]['D'][pidx[g]]) for g in range(F2)], np.int64)
        opts = np.zeros((F2, int(m.max()), s), np.int64)
        for g in range(F2):
            o = inst.pairs[k][g]['D'][pidx[g]]; opts[g, :len(o)] = o
        Dc = enum_opts(opts, m)
        D3 = np.stack([Dc[:, inst.rows[r]].sum(1) for r in range(3)], 1)
        kd = (D3.reshape(len(D3), -1).astype(np.uint64) * self.wd[None, :]).sum(1)
        u, first, inv = np.unique(kd, return_index=True, return_inverse=True)
        order = np.argsort(inv, kind='stable'); counts = np.bincount(inv.ravel(), minlength=len(u))
        start = np.zeros(len(u) + 1, np.int64); start[1:] = np.cumsum(counts)
        sig3 = seq_sig(D3[first], inst.cy.T, 2)
        k3 = (sig3.astype(np.uint64) * self.w3[None, :]).sum(1)
        sigf = seq_sig(Dc, inst.cy.T, inst.nl)
        kf = (sigf.astype(np.uint64) * self.wf[None, :]).sum(1)
        return dict(D=Dc, k3=k3, sig3=sig3, start=start, order=order.ravel(), sigf=sigf, kf=kf, n=len(Dc))

    def check(self, pidx4, maxD3=1 << 22):
        B = [self.block(k, pidx4[k]) for k in range(4)]
        h3, n3 = hjoin4(B[0]['k3'], B[1]['k3'], B[2]['k3'], B[3]['k3'], self.k3T, 16, maxD3)
        sols = []; nd3 = 0
        for q in h3:
            if not (B[0]['sig3'][q[0]] + B[1]['sig3'][q[1]] + B[2]['sig3'][q[2]] + B[3]['sig3'][q[3]] == self.T3).all(): continue
            nd3 += 1
            idx = [B[k]['order'][B[k]['start'][q[k]]:B[k]['start'][q[k] + 1]] for k in range(4)]
            hf, _ = hjoin4(B[0]['kf'][idx[0]], B[1]['kf'][idx[1]], B[2]['kf'][idx[2]], B[3]['kf'][idx[3]], self.kfT, 8, 1 << 16)
            for h in hf:
                cfg = [idx[k][h[k]] for k in range(4)]
                tot = sum(B[k]['sigf'][cfg[k]] for k in range(4))
                if (tot == self.TF).all(): sols.append([B[k]['D'][cfg[k]] for k in range(4)])
        return sols, dict(d3_solutions=nd3, d_options_log2=[round(math.log2(b['n']), 2) for b in B], d3_hits=int(n3))

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    # agreement with the reference d_check on planted S-solutions (v = 79 instance, planted target)
    from q523_sdsolve import SDInst, planted_target, d_check
    Y = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'logs', 'q523', 'stage1_v79_s1.json')))['y']
    inst0 = SDInst(6, 13, 4, Y); rng = np.random.default_rng(1); agree = 0; tot = 0
    for t in range(20):
        Pt, AC = planted_target(inst0, rng); inst = SDInst(6, 13, 4, Y, T=Pt)
        from q523_lib import zcoords
        Z = zcoords(AC)
        pidx4 = []
        for k in range(4):
            row = []
            for g in range(inst.F2):
                S = tuple((Z[k, g] + Z[k, g + inst.F2]).tolist())
                row.append(next(i for i, sv in enumerate(inst.pairs[k][g]['S'].tolist()) if tuple(sv) == S))
            pidx4.append(np.array(row))
        a, _ = d_check(inst, pidx4); b, _ = DChecker(inst).check(pidx4)
        ka = sorted(tuple(np.concatenate([x.ravel() for x in s_]).tolist()) for s_ in a)
        kb = sorted(tuple(np.concatenate([x.ravel() for x in s_]).tolist()) for s_ in b)
        tot += 1; agree += int(ka == kb and len(kb) >= 1)
    print(json.dumps(dict(test="DChecker vs d_check on planted S-solutions (v=79)", agree=agree, total=tot)))
