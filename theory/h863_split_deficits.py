r"""[M sizes, E bits] Deficits of candidate splits of the v = 863 lift (S7 host, Codex certificate y; see h863_lib.py).

523-style S/D split along an involution z of Q = Z2^2 (the U/W split of Kumar's 2148 route is the same split):
  S(g) = x(g) + x(g + z) (the compression to G/<z>, characters trivial on z), D = x(g) - x(g + z).
  h_D = H(D | S) exactly: per pair of cells (h, h + z) the two class-sign rows are uniform with fixed sums; H(D|S) =
        log2(n0 n1) - H(r0 + r1), summed over the 40 cell pairs (the singletons are fixed by y).
  C_S = Gaussian bits of the S-level equations (class PAF summed over {g, g + z}); C(D|S) = C_total - C_S.
  D_def = C(D|S) - h_D  (523 form): log2 of the S-solutions needed per core.
Inside the fibre of a G/Z5 solution (y, w) (--inner k, k = index in data/h863_g5_solutions.json): the same split
of the 16 doubly fixed 5 x 6 tables; H(S | y, w) and h_D = H(D | S, y, w) exactly, by enumerating every table of
both cells of each pair (q, q ^ z) and the distribution of their sum table; C over the mixed equations, uniform
sampling in the (y, w) fibre.
usage: python3 theory/h863_split_deficits.py [--N samples] [--out F]            (given y)
       python3 theory/h863_split_deficits.py --inner 0,1,2 [--N samples] [--out F] (inside (y, w))
"""
import sys, os, json, math, time, argparse, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, q523_gauss as QG


def pair_classes(zmask):
    """for each lag class c, the class of (q ^ zmask, u, t); returns the list of unordered pairs {c, c'}"""
    partner = np.empty(72, np.int64)
    for c, (q, u, t) in enumerate(L.LAGREPS):
        partner[c] = L.LAGCLS[(q ^ zmask, u, t)]
    pairs = sorted({tuple(sorted((c, int(partner[c])))) for c in range(72)})
    return partner, pairs


def row_pair_entropy(s0, s1):
    """exact H(D | S) for two uniform 6-sign rows with sums s0, s1: log2(n0 n1) - H(r0 + r1)"""
    R0 = [r for r in it.product((-1, 1), repeat=6) if sum(r) == s0]
    R1 = [r for r in it.product((-1, 1), repeat=6) if sum(r) == s1]
    cnt = {}
    for a in R0:
        for b in R1:
            k = tuple(x + y for x, y in zip(a, b)); cnt[k] = cnt.get(k, 0) + 1
    tot = len(R0) * len(R1)
    HS = -sum(v / tot * math.log2(v / tot) for v in cnt.values())
    return math.log2(tot) - HS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cert', default='data/h3452_stage1.json')
    ap.add_argument('--N', type=int, default=16000)
    ap.add_argument('--seed', type=int, default=9)
    ap.add_argument('--out', default=None)
    ap.add_argument('--inner', default=None, help='comma-separated w indices')
    a = ap.parse_args()
    Y, E, S = L.load_cert(a.cert)
    rng = np.random.default_rng(a.seed)
    if a.inner is not None:
        W = json.load(open('data/h863_g5_solutions.json'))['w']
        fo = open(a.out, 'a') if a.out else None
        for k in [int(v) for v in a.inner.split(',')]:
            for zmask, name in ((1, 'z=(0,1)'), (2, 'z=(1,0)'), (3, 'z=(1,1)')):
                r = inner_split(S, W[k], zmask, a.N, rng, E); r.update(w_index=k, z=name)
                print(json.dumps(r), flush=True)
                if fo: fo.write(json.dumps(r) + '\n'); fo.flush()
        return
    c = L.sample_y_fibre(S, a.N, rng)
    x = L.assemble(np.broadcast_to(E, (a.N, 4, 20)), c)
    V = L.class_paf(x).sum(1)
    tot = QG.gauss_bits(V, L.FULL_TARGET)
    fo = open(a.out, 'a') if a.out else None
    def log(d):
        s = json.dumps(d); print(s, flush=True)
        if fo: fo.write(s + '\n'); fo.flush()
    log(dict(what='given y', total=tot, log2_fibre=round(sum(math.log2(math.comb(6, (int(s) + 6) // 2))
                                                         for s in S.ravel()), 2)))
    for zmask, name in ((1, 'z=(0,1)'), (2, 'z=(1,0)'), (3, 'z=(1,1)')):
        partner, pairs = pair_classes(zmask)
        Vs = np.stack([V[:, p] + V[:, q] if p != q else 2 * V[:, p] for p, q in pairs], 1)
        Ts = np.array([L.FULL_TARGET[p] + L.FULL_TARGET[q] if p != q else 2 * L.FULL_TARGET[p] for p, q in pairs])
        gs = QG.gauss_bits(Vs, Ts)
        hD = 0.0
        for i in range(4):
            for q in range(4):
                if q < q ^ zmask:
                    for u in range(5):
                        hD += row_pair_entropy(int(S[i, q * 5 + u]), int(S[i, (q ^ zmask) * 5 + u]))
        CDS = tot['bits'] - gs['bits']
        log(dict(what='S/D split given y', z=name, C_S=gs['bits'], rank_S=gs['rank'], in_lattice_S=gs['in_lattice'],
                 C_D_given_S=round(CDS, 2), h_D=round(hD, 2), D_def=round(CDS - hD, 2),
                 log2_S_solutions_per_y=round(sum(math.log2(math.comb(6, (int(s) + 6) // 2)) for s in S.ravel())
                                              - hD - gs['bits'], 2)))



# ---------------------------------------------------------------- inside a (y, w) fibre
def all_tables(rowsums, colsums):
    """explicit list of all +-1 tables (R x C) with the given margins, as an (n, R, C) int8 array"""
    R, C = len(rowsums), len(colsums)
    rp = [(s + C) // 2 for s in rowsums]; cp = [(s + R) // 2 for s in colsums]
    out = []
    rows = {n: [r for r in it.product((0, 1), repeat=C) if sum(r) == n] for n in set(rp)}
    def rec(r, rem, acc):
        if r == R:
            if not any(rem): out.append(list(acc))
            return
        left = R - r - 1
        for row in rows[rp[r]]:
            nr = [a - b for a, b in zip(rem, row)]
            if min(nr) < 0 or max(nr) > left: continue
            acc.append(row); rec(r + 1, nr, acc); acc.pop()
    rec(0, cp, [])
    return 2 * np.array(out, np.int8) - 1


def inner_split(S, w4, zmask, N, rng, E):
    import h863_yw_estimate as YW
    ts, lf = YW.yw_fibre(S, w4)
    # exact H(S | y, w) per table pair
    HSsum = 0.0; hD = 0.0
    for i in range(4):
        for q in range(4):
            qq = q ^ zmask
            if q > qq: continue
            A = all_tables([int(S[i, q * 5 + u]) for u in range(5)], w4[i][q])
            B = all_tables([int(S[i, qq * 5 + u]) for u in range(5)], w4[i][qq])
            assert len(A) == ts[i, q].total and len(B) == ts[i, qq].total
            ka = (A.reshape(len(A), -1).astype(np.int64) + 1) // 2
            kb = (B.reshape(len(B), -1).astype(np.int64) + 1) // 2
            # sum table entries in {0, 1, 2} (number of +1): key = sum over cells of (a + b) * 3^cell  (as python ints)
            # keys fit in 3^30 < 2^48
            pw64 = np.array([3 ** k for k in range(30)], dtype=np.int64)
            sa = ka @ pw64; sb = kb @ pw64
            keys = (sa[:, None] + sb[None, :]).ravel()
            _, cnt = np.unique(keys, return_counts=True)
            p = cnt / cnt.sum()
            Hs = float(-(p * np.log2(p)).sum())
            HSsum += Hs; hD += math.log2(len(A) * len(B)) - Hs
    # Gaussian bits of S-level mixed equations and of all mixed equations, uniform in the (y, w) fibre
    x = YW.sample_yw(E, ts, N, rng)
    V = L.class_paf(x).sum(1)
    tot = QG.gauss_bits(V, L.FULL_TARGET)
    partner, pairs = pair_classes(zmask)
    Vs = np.stack([V[:, p] + V[:, q] if p != q else 2 * V[:, p] for p, q in pairs], 1)
    Ts = np.array([L.FULL_TARGET[p] + L.FULL_TARGET[q] if p != q else 2 * L.FULL_TARGET[p] for p, q in pairs])
    gs = QG.gauss_bits(Vs, Ts)
    return dict(log2_fibre_yw=round(lf, 2), H_S=round(HSsum, 2), h_D=round(hD, 2), C_mixed=tot['bits'],
                rank_mixed=tot['rank'], C_S=gs['bits'], rank_S=gs['rank'], C_D_given_S=round(tot['bits'] - gs['bits'], 2),
                D_def=round(tot['bits'] - gs['bits'] - hD, 2), log2_S_solutions_per_yw=round(HSsum - gs['bits'], 2))


if __name__ == '__main__':
    main()
