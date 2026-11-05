r"""[M] Exact solutions of the G/Z5 level (Z2^2 x Z43) inside the fibre of a Z43-compressed solution y at v = 863
(S7 host 2e93be9718324148; the certificate data/h3452_stage1.json and its set-up are by OpenAI Codex, directed by the
repository owner; the bordered construction is Kumar's / EKLO's).

Unknowns: per block i and Z2^2-cell q, the column-sum vector w_i(q) in Z^6 of a 5 x 6 sign table whose row sums
S_i(q, u) are fixed by y (h863_lib.cell_lists).  Equations: the 12 independent G/Z5 lag equations, split by the
characters psi = (a1, a2) of Z2^2 into
  stage alpha (a2 = 0: psi00, psi10): only A0_i = w_i(00) + w_i(01) and B0_i = w_i(10) + w_i(11) enter;
  stage beta  (a2 = 1: psi01, psi11): A1_i = w_i(00) - w_i(01), B1_i = w_i(10) - w_i(11), given A0, B0.
Each stage is a 2 + 2 join over the four blocks on a 6-dimensional integer signature (Z43 PAF of the psi-twisted
sequences at the lags 1, 3, 7; the lag-0 value is then implied by y).  Every hit is re-checked on all 172 lags by
h863_lib.g5_check (explicit sequences).  Independent verification: theory/h863_verify_g5.py.

usage: python3 theory/h863_g5_solve.py [--cert C] [--na N] [--nb N] [--seed s] [--want k] [--out F] [--split q2|q1|q12]
"""
import sys, os, json, time, math, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L

KEYR = None


def keys_of(V):
    """64-bit hash keys of integer signature rows (exactness re-checked after matching)"""
    global KEYR
    if KEYR is None or len(KEYR) != V.shape[-1]:
        KEYR = np.random.default_rng(12345).integers(1, 2 ** 62, size=V.shape[-1], dtype=np.int64) | 1
    return (V.astype(np.int64) * KEYR).sum(-1)


def join22(s0, s1, s2, s3, T, cap=10 ** 6):
    """all (a, b, c, d) with s0[a] + s1[b] + s2[c] + s3[d] == T (exact), s_k: (n_k, D) int64"""
    lk = keys_of(s0[:, None, :] + s1[None, :, :]).ravel()
    rk = keys_of(T[None, None, :] - s2[:, None, :] - s3[None, :, :]).ravel()
    lo = np.argsort(lk, kind='stable'); lks = lk[lo]
    pos = np.searchsorted(lks, rk)
    pos[pos >= len(lks)] = len(lks) - 1
    hit = np.nonzero(lks[pos] == rk)[0]
    out = []
    n1, n3 = len(s1), len(s3)
    for r in hit:
        # all left entries with that key
        k = rk[r]; p = pos[r]
        while p < len(lks) and lks[p] == k:
            li = lo[p]; a, b = divmod(int(li), n1); c, d = divmod(int(r), n3)
            if (s0[a] + s1[b] + s2[c] + s3[d] == T).all():
                out.append((a, b, c, d))
                if len(out) >= cap: return out
            p += 1
    return out


class Problem:
    def __init__(self, cert, split='q2'):
        self.Y, self.E, self.S = L.load_cert(cert)
        L.check_y(self.Y)
        self.CL = L.cell_lists(self.S)
        self.E4 = np.array([[self.E[i, q * 5:(q + 1) * 5].sum() for q in range(4)] for i in range(4)], np.int64)
        # cell permutation so that the stage-alpha pairing is {perm[0], perm[1]} | {perm[2], perm[3]}
        # q2: halves q1 = 0 / q1 = 1 (compress along q2); q1: halves q2 = 0 / 1; q12: halves {00, 11} / {01, 10}
        self.perm = {'q2': [0, 1, 2, 3], 'q1': [0, 2, 1, 3], 'q12': [0, 3, 1, 2]}[split]
        self.split = split
        self.lists = {}
        for i in range(4):
            for q in range(4):
                W, m = self.CL[i, q]
                mf = np.array([float(x) for x in m])
                self.lists[i, q] = (W, mf / mf.sum(), {tuple(w): k for k, w in enumerate(W.tolist())})
        # characters, in the permuted cell order c = (c0, c1, c2, c3): alpha uses (c0 + c1) +- (c2 + c3)
        p = self.perm
        self.Ep = self.E4[:, p]

    def alpha_sig(self, i, A, B):
        e = self.Ep[i]
        z0p = e[0] + e[1] + e[2] + e[3]; z0m = e[0] + e[1] - e[2] - e[3]
        return np.concatenate([L.z43_paf3(np.full(len(A), z0p), A + B), L.z43_paf3(np.full(len(A), z0m), A - B)], -1)

    def beta_sig(self, i, A1, B1):
        e = self.Ep[i]
        z0p = e[0] - e[1] + e[2] - e[3]; z0m = e[0] - e[1] - e[2] + e[3]
        return np.concatenate([L.z43_paf3(np.full(len(A1), z0p), A1 + B1), L.z43_paf3(np.full(len(A1), z0m), A1 - B1)], -1)

    def targets(self):
        # alpha: characters trivial on the pairing-internal difference: psi with psi(c0) = psi(c1), psi(c2) = psi(c3)
        # in terms of the original psi-index the targets are psi0 -> (0,0,0) and the others -> (-80,)*3
        return np.array([0, 0, 0, -80, -80, -80], np.int64), np.array([-80] * 6, np.int64)

    def sample_halves(self, i, n, rng):
        p = self.perm
        out = []
        for h in (0, 1):
            qa, qb = p[2 * h], p[2 * h + 1]
            Wa, pa, _ = self.lists[i, qa]; Wb, pb, _ = self.lists[i, qb]
            ia = rng.choice(len(Wa), size=n, p=pa); ib = rng.choice(len(Wb), size=n, p=pb)
            out.append(Wa[ia] + Wb[ib])
        return out

    def splits(self, i, h, A):
        """all (w_a, w_b) with w_a + w_b = A for half h of block i, and their weights"""
        p = self.perm
        qa, qb = p[2 * h], p[2 * h + 1]
        Wa, pa, _ = self.lists[i, qa]; Wb, pb, db = self.lists[i, qb]
        res = []
        for ka, wa in enumerate(Wa.tolist()):
            kb = db.get(tuple(int(x) - y for x, y in zip(A, wa)))
            if kb is not None: res.append((ka, kb))
        res = np.array(res, np.int64).reshape(-1, 2)
        return Wa[res[:, 0]], Wb[res[:, 1]], pa[res[:, 0]] * pb[res[:, 1]]

    def assemble(self, wq):
        """wq: dict (i, permuted position) -> w; returns w4 (4, 4, 6) in the original cell order"""
        w4 = np.zeros((4, 4, 6), np.int64)
        for (i, pos), w in wq.items(): w4[i, self.perm[pos]] = w
        return w4


def solve(pb, na=1024, nb=512, seed=1, want=4, log=print, max_alpha=50, tlimit=3600):
    rng = np.random.default_rng(seed)
    Ta, Tb = pb.targets()
    sols = []
    t0 = time.time(); nalpha = 0; rounds = 0
    while len(sols) < want and time.time() - t0 < tlimit:
        rounds += 1
        H = [pb.sample_halves(i, na, rng) for i in range(4)]
        sig = [pb.alpha_sig(i, H[i][0], H[i][1]) for i in range(4)]
        hits = join22(sig[0], sig[1], sig[2], sig[3], Ta)
        log(json.dumps(dict(stage='alpha', round=rounds, na=na, hits=len(hits), secs=round(time.time() - t0, 1))))
        for (a, b, c, d) in hits[:max_alpha]:
            nalpha += 1
            pick = [a, b, c, d]
            AB = [(H[i][0][pick[i]], H[i][1][pick[i]]) for i in range(4)]
            # stage beta: per block, all splits of A0 and B0, subsampled to nb
            opts = []
            for i in range(4):
                Wa, Wb, pa = pb.splits(i, 0, AB[i][0]); Wc, Wd, pc = pb.splits(i, 1, AB[i][1])
                ka = rng.choice(len(Wa), size=nb, p=pa / pa.sum()); kc = rng.choice(len(Wc), size=nb, p=pc / pc.sum())
                opts.append((Wa[ka], Wb[ka], Wc[kc], Wd[kc], len(Wa), len(Wc)))
            bs = [pb.beta_sig(i, o[0] - o[1], o[2] - o[3]) for i, o in enumerate(opts)]
            bh = join22(bs[0], bs[1], bs[2], bs[3], Tb, cap=8)
            log(json.dumps(dict(stage='beta', alpha=nalpha, nsplits=[[o[4], o[5]] for o in opts], hits=len(bh),
                                secs=round(time.time() - t0, 1))))
            for (a2, b2, c2, d2) in bh:
                pk = [a2, b2, c2, d2]
                wq = {}
                for i in range(4):
                    o = opts[i]; k = pk[i]
                    wq[i, 0], wq[i, 1], wq[i, 2], wq[i, 3] = o[0][k], o[1][k], o[2][k], o[3][k]
                w4 = pb.assemble(wq)
                bad = L.g5_check(pb.E4, w4)
                assert not bad, bad[:5]
                key = w4.tobytes()
                if all(key != s.tobytes() for s in sols):
                    sols.append(w4)
                    log(json.dumps(dict(stage='solution', n=len(sols), w=w4.tolist())))
            if len(sols) >= want: break
    return sols


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--cert', default='data/h3452_stage1.json')
    ap.add_argument('--na', type=int, default=1024)
    ap.add_argument('--nb', type=int, default=512)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--want', type=int, default=4)
    ap.add_argument('--split', default='q2')
    ap.add_argument('--tlimit', type=float, default=3600)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    pb = Problem(a.cert, a.split)
    sols = solve(pb, a.na, a.nb, a.seed, a.want, tlimit=a.tlimit)
    res = dict(cert=a.cert, split=a.split, seed=a.seed, na=a.na, nb=a.nb, nsol=len(sols),
               E4=pb.E4.tolist(), w=[s.tolist() for s in sols])
    print(json.dumps(dict(nsol=len(sols))))
    if a.out:
        json.dump(res, open(a.out, 'w'))
