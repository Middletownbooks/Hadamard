r"""S-join kernels on one full-size centre (report/paley_523_answer.md, section 3.3) [M: same hit set, timing].

Builds the join sides of one centre exactly as theory/q523_wagner.run does (lists from seed, middle value from
default_rng([seed, ci]), all parent pairs up to logN), then runs mjoin_rad and mjoin_rad2 (G residue classes per pass,
2^rbits partitions, bitset of 2^bb bits) over the same residue classes p0 .. p0 + n - 1 and compares the hit sets.
usage: python3 theory/q523_joinbench.py [--seed 4242] [--centre 0] [--logN 38] [--K 32] [--classes 64] [--p0 1000]
          [--configs 2,5,23 ...]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_wagner as W
from q523_sdsolve import SDInst, seq_sig, next_prime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def build(inst, seed, ci, logN, K, pass_log2=23.0):
    rng = np.random.default_rng(seed)
    rows = [W.block_candidates(inst, k, 1.0, 64, K, rng) for k in range(4)]
    C = [W.cand_table(inst, rows[k]) for k in range(4)]
    T = inst.cy.T; tS = inst.tsigS(); D = tS.shape[0]
    wts = rng.integers(1, 2 ** 63, D, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    rho = rng.integers(1, 1 << 30, D, dtype=np.int64)
    sc, c, L, R, wl, wr = W.choose_middle(inst, C, np.random.default_rng([seed, ci]))
    def trim(Pp, ka, kb):
        Pp = sorted(Pp, key=lambda ij: -(C[ka][2][ij[0]] + C[kb][2][ij[1]])); keep = []; tot = 0.0
        for i, j in Pp:
            keep.append((i, j)); tot += 2 ** (C[ka][2][i] + C[kb][2][j])
            if math.log2(tot) >= logN: break
        return keep, math.log2(tot)
    Lk, wL = trim(L, 0, 1); Rk, wR = trim(R, 2, 3)
    need = [sorted({p[0] for p in Lk}), sorted({p[1] for p in Lk}), sorted({p[0] for p in Rk}), sorted({p[1] for p in Rk})]
    P = next_prime(max(2, 2 ** wL / 2.0 ** pass_log2))
    def fibres_of(k):
        for x in need[k]:
            pidx, S = W.materialise(inst, k, rows[k], C[k][0][x], K)
            sig = seq_sig(S, T, inst.nl); del S
            yield pidx, sig
    sides = [W.Side(fibres_of(k), wts, rho, P) for k in range(4)]
    fid = [{x: i for i, x in enumerate(need[k])} for k in range(4)]
    lp = [np.array([fid[0][i] for i, j in Lk], np.int64), np.array([fid[1][j] for i, j in Lk], np.int64),
          np.array([fid[2][i] for i, j in Rk], np.int64), np.array([fid[3][j] for i, j in Rk], np.int64)]
    kT = (tS.astype(np.int64).astype(np.uint64) * wts).sum(); rT = int(((tS % P) * (rho % P) % P).sum() % P)
    A, B, Cc, Dd = sides
    args = (A.key, A.res, A.off, B.key, B.rst, B.off, Cc.key, Cc.res, Cc.off, Dd.key, Dd.rst, Dd.off, *lp, P, rT, kT)
    return args, dict(parents=(len(Lk), len(Rk)), side_log2=(round(wL, 2), round(wR, 2)), P=P), wL, wR, P

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--seed": "4242", "--centre": "0", "--logN": "38", "--K": "32", "--classes": "64", "--p0": "1000"}
    cfg = []; i = 0
    while i < len(a):
        if a[i] == "--configs":
            i += 1
            while i < len(a) and not a[i].startswith('--'): cfg.append(tuple(map(int, a[i].split(',')))); i += 1
        else: opt[a[i]] = a[i + 1]; i += 2
    cfg = cfg or [(2, 5, 23)]
    Y = json.load(open(os.path.join(ROOT, 'logs/q523/y523_best.json')))['y']; inst = SDInst(18, 29, 4, Y)
    t0 = time.time(); args, info, wL, wR, P = build(inst, int(opt["--seed"]), int(opt["--centre"]), float(opt["--logN"]), int(opt["--K"]))
    info['build_secs'] = round(time.time() - t0, 1); print(json.dumps(info), flush=True)
    n = int(opt["--classes"]); p0 = int(opt["--p0"])
    bs, hk, hv, bs2 = np.zeros(1 << 17, np.uint64), np.zeros(1 << 19, np.uint64), np.full(1 << 19, -1, np.int64), np.zeros(1 << 15, np.uint64)
    def bufs(G):
        cL = int(G * 2 ** wL / P * 1.06) + (1 << 18); cR = int(G * 2 ** wR / P * 1.06) + (1 << 18)
        return np.zeros(cL, np.uint64), np.zeros(cL, np.int64), np.zeros(cR, np.uint64), np.zeros(cR, np.int64)
    Lk, Li, Rk, Ri = bufs(1)
    W.mjoin_rad(*args, 4, Lk, Li, Rk, Ri, bs, 23, hk, hv, 1 << 20, bs2, 21, 0, 1)
    t = time.time(); h0, nh0, nl, nr, _ = W.mjoin_rad(*args, 4, Lk, Li, Rk, Ri, bs, 23, hk, hv, 1 << 20, bs2, 21, p0, p0 + n)
    t_old = time.time() - t; H0 = set(map(tuple, h0.tolist())); del Lk, Li, Rk, Ri
    print(json.dumps(dict(kernel="mjoin_rad", ms_per_class=round(1e3 * t_old / n, 1), ns_per_pair=round(1e9 * t_old / (nl + nr), 1),
                          hits=int(nh0))), flush=True)
    for G, rb, bb in cfg:
        Lk, Li, Rk, Ri = bufs(G)
        W.mjoin_rad2(*args, rb, Lk, Li, Rk, Ri, bs, bb, hk, hv, 1 << 20, bs2, 21, G, 0, G)
        t = time.time(); h1, nh1, nl1, nr1, _ = W.mjoin_rad2(*args, rb, Lk, Li, Rk, Ri, bs, bb, hk, hv, 1 << 20, bs2, 21, G, p0, p0 + n)
        t_new = time.time() - t
        print(json.dumps(dict(kernel="mjoin_rad2", G=G, rbits=rb, bb=bb, ms_per_class=round(1e3 * t_new / n, 1),
                              ns_per_pair=round(1e9 * t_new / (nl1 + nr1), 1), speedup=round(t_old / t_new, 2), hits=int(nh1),
                              same_hits=set(map(tuple, h1.tolist())) == H0, same_pairs=bool(nl1 == nl and nr1 == nr),
                              buffers_clean=bool((bs == 0).all() and (bs2 == 0).all() and (hv == -1).all()))), flush=True)
        del Lk, Li, Rk, Ri
