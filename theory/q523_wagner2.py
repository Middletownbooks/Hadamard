r"""Large middle level for the S-step (report/paley_523_answer.md, section 3.3) [M: same S-solutions as
theory/q523_wagner.py on shared parents; `python3 theory/q523_wagner2.py test`].

theory/q523_wagner.py keeps every parent (row-sum triple) of a block in a table (64 row sums per row, 2^18 parents),
which caps a Z3 middle value at a few hundred parent pairs and the second-level sides at about 2^35.5.  Here a
block's parents are never tabulated: with nrow row sums per row (up to 256, so 2^24 parents) the Z3 signature of a
parent (i0, i1, i2) is the sum of per-row self terms (lag 0) and per-row-pair cross terms (lag 1), and so is its
64-bit linear key.  Keys of blocks 1 and 3 are computed in one compiled pass and sorted; parent pairs of a middle
value c are found by scanning blocks 0 / 2 and binary-searching the needed key (every hit re-checked exactly).
The second level, the D-check and the checkpoints are those of theory/q523_wagner.py.
usage: python3 theory/q523_wagner2.py test"""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_wagner as W
from q523_sdsolve import SDInst, seq_sig, next_prime, planted_target

def row_tables(inst, rows_k, w):
    """per row: row sums (n_r, s), self terms (n_r, s), log2 fibre sizes, log2 masses; per row pair: cross terms;
    keys: per-row self keys and per-pair cross keys (uint64)"""
    T = inst.cy.T; s = inst.s
    RS = [np.array([e[0] for e in rows_k[r]], np.int64).reshape(len(rows_k[r]), s) for r in range(3)]
    def tc(x, y):                                                    # T-contraction of x (x) y, vectorised
        return np.einsum('ijl,...i,...j->...l', T, x, y)
    selfv = [tc(RS[r], RS[r]) for r in range(3)]
    c01 = tc(RS[1][None, :, :], RS[0][:, None, :])                   # [i0, i1] : X1 (x) X0
    c12 = tc(RS[2][None, :, :], RS[1][:, None, :])                   # [i1, i2] : X2 (x) X1
    c20 = tc(RS[0][None, :, :], RS[2][:, None, :])                   # [i2, i0] : X0 (x) X2
    w0 = w[:s]; w1 = w[s:]
    ks = [(selfv[r].astype(np.uint64) * w0[None, :]).sum(1) for r in range(3)]
    k01 = (c01.astype(np.uint64) * w1[None, None, :]).sum(2); k12 = (c12.astype(np.uint64) * w1[None, None, :]).sum(2)
    k20 = (c20.astype(np.uint64) * w1[None, None, :]).sum(2)
    lsz = [np.array([math.log2(len(e[1])) for e in rows_k[r]]) for r in range(3)]
    lms = [np.array([e[2] for e in rows_k[r]]) for r in range(3)]
    return dict(RS=RS, selfv=selfv, c01=c01, c12=c12, c20=c20, ks=ks, k01=k01, k12=k12, k20=k20, lsz=lsz, lms=lms)

def sig_of(tb, i0, i1, i2):
    s = tb['selfv'][0].shape[1]
    lag0 = tb['selfv'][0][i0] + tb['selfv'][1][i1] + tb['selfv'][2][i2]
    lag1 = tb['c01'][i0, i1] + tb['c12'][i1, i2] + tb['c20'][i2, i0]
    return np.concatenate([lag0, lag1])

@numba.njit(cache=True)
def _all_keys(ks0, ks1, ks2, k01, k12, k20):
    n0 = ks0.shape[0]; n1 = ks1.shape[0]; n2 = ks2.shape[0]
    K = np.empty(n0 * n1 * n2, np.uint64); t = 0
    for a in range(n0):
        for b in range(n1):
            base = ks0[a] + ks1[b] + k01[a, b]
            for c in range(n2):
                K[t] = base + ks2[c] + k12[b, c] + k20[c, a]; t += 1
    return K

@numba.njit(cache=True)
def _pairs_scan(ks0, ks1, ks2, k01, k12, k20, ckey, sk, so, out, cap):
    """all (parent of block X, parent of block Y) with key_X + key_Y == ckey: scan X's parents, binary-search the
    needed key in Y's sorted keys sk (order so); returns (count, pairs (cap, 2) of flat parent indices)"""
    n0 = ks0.shape[0]; n1 = ks1.shape[0]; n2 = ks2.shape[0]; n = 0; N = sk.shape[0]
    for a in range(n0):
        for b in range(n1):
            base = ks0[a] + ks1[b] + k01[a, b]
            for c in range(n2):
                need = ckey - (base + ks2[c] + k12[b, c] + k20[c, a])
                lo = np.searchsorted(sk, need)
                while lo < N and sk[lo] == need:
                    if n < cap:
                        out[n, 0] = (a * n1 + b) * n2 + c; out[n, 1] = so[lo]
                    n += 1; lo += 1
    return n

class BigMiddle:
    """lazy parents of the four blocks for one run (rows from W.block_candidates with large nrow)"""
    def __init__(self, inst, rows, seed=12345):
        self.inst = inst; self.rows = rows; s = inst.s
        self.w = np.random.default_rng(seed).integers(1, 2 ** 62, 2 * s, dtype=np.int64).astype(np.uint64)
        self.tb = [row_tables(inst, rows[k], self.w) for k in range(4)]
        self.n = [tuple(len(rows[k][r]) for r in range(3)) for k in range(4)]
        self.sorted = {}
        for k in (1, 3):
            t = self.tb[k]; K = _all_keys(t['ks'][0], t['ks'][1], t['ks'][2], t['k01'], t['k12'], t['k20'])
            o = np.argsort(K); self.sorted[k] = (K[o], o.astype(np.int64)); del K

    def unflat(self, k, f):
        n0, n1, n2 = self.n[k]; return f // (n1 * n2), (f // n2) % n1, f % n2

    def key_of_sig(self, sig):
        return (np.asarray(sig, np.int64).astype(np.uint64) * self.w).sum()

    def pairs(self, kx, ky, c, cap=1 << 22):
        """exact parent pairs (flat indices) of blocks kx, ky with Z3sig sum == c"""
        t = self.tb[kx]; sk, so = self.sorted[ky]
        out = np.zeros((cap, 2), np.int64)
        n = _pairs_scan(t['ks'][0], t['ks'][1], t['ks'][2], t['k01'], t['k12'], t['k20'], self.key_of_sig(c), sk, so, out, cap)
        out = out[:min(n, cap)]; good = []
        for fx, fy in out:
            if (sig_of(self.tb[kx], *self.unflat(kx, fx)) + sig_of(self.tb[ky], *self.unflat(ky, fy)) == c).all():
                good.append((int(fx), int(fy)))
        return good

    def lsz(self, k, f):
        i = self.unflat(k, f); return sum(self.tb[k]['lsz'][r][i[r]] for r in range(3))

    def sample_parent(self, k, rng, beta=1.0):
        idx = []
        for r in range(3):
            m = self.tb[k]['lms'][r]; p = np.exp2(beta * (m - m.max())); p /= p.sum(); idx.append(int(rng.choice(len(m), p=p)))
        return tuple(idx)

def choose_middle_big(bm, rng, T3, ncand=12, log=None):
    """middle value maximising min(log2 left, log2 right configuration pairs) among ncand sampled left parent pairs"""
    best = None
    for _ in range(ncand):
        a = bm.sample_parent(0, rng); b = bm.sample_parent(1, rng)
        c = sig_of(bm.tb[0], *a) + sig_of(bm.tb[1], *b)
        R = bm.pairs(2, 3, T3 - c)
        if not R: continue
        L = bm.pairs(0, 1, c)
        if not L: continue
        wl = math.log2(sum(2.0 ** (bm.lsz(0, i) + bm.lsz(1, j)) for i, j in L))
        wr = math.log2(sum(2.0 ** (bm.lsz(2, i) + bm.lsz(3, j)) for i, j in R))
        if log: log(dict(stage="middle-candidate", left_pairs=len(L), right_pairs=len(R), wl=round(wl, 2), wr=round(wr, 2)))
        if best is None or min(wl, wr) > best[0]: best = (min(wl, wr), c, L, R, wl, wr)
    return best

def grid_row(bm, k, f):
    return np.array(bm.unflat(k, f), np.int64)

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    # 1. decomposed signatures agree with z3sig_of; 2. pair search agrees with W.parent_pairs on the same rows
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    Y = json.load(open(os.path.join(ROOT, 'logs/q523/y523_best.json')))['y']; inst = SDInst(18, 29, 4, Y)
    rng = np.random.default_rng(3)
    rows = [W.block_candidates(inst, k, 1.0, 24, 24, rng) for k in range(4)]
    bm = BigMiddle(inst, rows)
    C = [W.cand_table(inst, rows[k]) for k in range(4)]
    ok = 0
    for k in range(4):
        for t in range(50):
            f = int(rng.integers(len(C[k][0]))); g = C[k][0][f]
            ok += int((sig_of(bm.tb[k], *g) == C[k][1][f]).all())
    T3 = inst.tsig3(); agree = 0; tot = 0
    for t in range(5):
        a = bm.sample_parent(0, rng); b = bm.sample_parent(1, rng); c = sig_of(bm.tb[0], *a) + sig_of(bm.tb[1], *b)
        L = sorted(bm.pairs(0, 1, c))
        w = np.random.default_rng(12345).integers(1, 2 ** 62, C[0][1].shape[1], dtype=np.int64).astype(np.uint64)
        h1 = W.SigHash(C[1][1], w); k0 = (C[0][1].astype(np.uint64) * w[None, :]).sum(1)
        L2 = W.parent_pairs(C[0][1], k0, h1, c)
        # map table indices (grid rows) to flat indices
        n0, n1, n2 = bm.n[0]; m0, m1, m2 = bm.n[1]
        L2f = sorted(((int(C[0][0][i][0]) * n1 + int(C[0][0][i][1])) * n2 + int(C[0][0][i][2]),
                      (int(C[1][0][j][0]) * m1 + int(C[1][0][j][1])) * m2 + int(C[1][0][j][2])) for i, j in L2)
        agree += int(L == L2f); tot += 1
    print(json.dumps(dict(sig_agree=f"{ok}/200", pairs_agree=f"{agree}/{tot}")))
