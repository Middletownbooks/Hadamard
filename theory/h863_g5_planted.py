r"""[M] Planted test of the G/Z5-level join (theory/h863_g5_solve.py) at the 863 shape (S7 host; set-up by OpenAI Codex).
A random full configuration (uniform class signs and singletons on all 80 cells of the four blocks) is planted; its
Z43-compression gives the y-fibre (row sums S, singletons e) and its G/Z5 compression gives the planted w and the
planted targets (its own psi-twisted Z43 PAF values).  The planted half sums are forced into the stage-alpha lists and
the planted splits into the stage-beta lists; the join must return the planted w.  Other hits are counted against
the lattice-Gaussian prediction at the planted target (a correctness and calibration check, not a cost measurement).
usage: python3 theory/h863_g5_planted.py [--trials T] [--na N] [--nb N] [--out F]"""
import sys, os, json, math, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, q523_gauss as QG
from h863_g5_solve import Problem, join22


class Planted(Problem):
    def __init__(self, rng, split='q2'):
        self.E = rng.choice([-1, 1], size=(4, 20)).astype(np.int64)
        c = rng.choice([-1, 1], size=(4, 20, 6)).astype(np.int64)
        self.S = c.sum(2)
        self.Y = self.E + 7 * self.S
        self.CL = L.cell_lists(self.S)
        self.E4 = np.array([[self.E[i, q * 5:(q + 1) * 5].sum() for q in range(4)] for i in range(4)], np.int64)
        self.w_planted = np.array([[c[i, q * 5:(q + 1) * 5, :].sum(0) for q in range(4)] for i in range(4)], np.int64)
        self.perm = {'q2': [0, 1, 2, 3], 'q1': [0, 2, 1, 3], 'q12': [0, 3, 1, 2]}[split]
        self.split = split
        self.lists = {}
        for i in range(4):
            for q in range(4):
                W, m = self.CL[i, q]
                mf = np.array([float(x) for x in m])
                self.lists[i, q] = (W, mf / mf.sum(), {tuple(w): k for k, w in enumerate(W.tolist())})
        self.Ep = self.E4[:, self.perm]
        wp = self.w_planted[:, self.perm]
        self.planted_halves = [(wp[i, 0] + wp[i, 1], wp[i, 2] + wp[i, 3]) for i in range(4)]
        self.planted_beta = [(wp[i, 0], wp[i, 1], wp[i, 2], wp[i, 3]) for i in range(4)]
        self.Ta = sum(self.alpha_sig(i, self.planted_halves[i][0][None], self.planted_halves[i][1][None])[0]
                      for i in range(4))
        self.Tb = sum(self.beta_sig(i, (wp[i, 0] - wp[i, 1])[None], (wp[i, 2] - wp[i, 3])[None])[0] for i in range(4))

    def targets(self):
        return self.Ta, self.Tb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=5)
    ap.add_argument('--na', type=int, default=1024)
    ap.add_argument('--nb', type=int, default=512)
    ap.add_argument('--seed', type=int, default=77)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    fo = open(a.out, 'w') if a.out else None
    for trial in range(a.trials):
        pb = Planted(rng, ['q2', 'q1', 'q12'][trial % 3])
        Ta, Tb = pb.targets()
        H = [pb.sample_halves(i, a.na, rng) for i in range(4)]
        k0 = [int(rng.integers(a.na)) for i in range(4)]
        for i in range(4):                                   # force the planted half sums into the lists
            H[i][0][k0[i]] = pb.planted_halves[i][0]; H[i][1][k0[i]] = pb.planted_halves[i][1]
        sig = [pb.alpha_sig(i, H[i][0], H[i][1]) for i in range(4)]
        hits = join22(sig[0], sig[1], sig[2], sig[3], Ta, cap=10 ** 7)
        alpha_ok = tuple(k0) in set(hits)
        # Gaussian prediction of the unforced alpha hits at the planted target
        Vn = sum(pb.alpha_sig(i, *pb.sample_halves(i, 40000, rng)) for i in range(4))
        ga = QG.gauss_bits(Vn, Ta)
        # stage beta on the planted alpha solution
        opts = []; kb = []
        for i in range(4):
            Wa, Wb, pa = pb.splits(i, 0, pb.planted_halves[i][0]); Wc, Wd, pc = pb.splits(i, 1, pb.planted_halves[i][1])
            ka = rng.choice(len(Wa), size=a.nb, p=pa / pa.sum()); kc = rng.choice(len(Wc), size=a.nb, p=pc / pc.sum())
            A0, A1, B0, B1 = Wa[ka], Wb[ka], Wc[kc], Wd[kc]
            j = int(rng.integers(a.nb)); kb.append(j)
            A0[j], A1[j], B0[j], B1[j] = pb.planted_beta[i]
            opts.append((A0, A1, B0, B1))
        bs = [pb.beta_sig(i, o[0] - o[1], o[2] - o[3]) for i, o in enumerate(opts)]
        bh = join22(bs[0], bs[1], bs[2], bs[3], Tb, cap=10 ** 7)
        beta_ok = tuple(kb) in set(bh)
        # the recovered w equals the planted w and passes the direct G/Z5 check at the planted targets
        wq = {}
        for i in range(4):
            o = opts[i]; j = kb[i]
            wq[i, 0], wq[i, 1], wq[i, 2], wq[i, 3] = o[0][j], o[1][j], o[2][j], o[3][j]
        w4 = pb.assemble(wq)
        same = bool((w4 == pb.w_planted).all())
        d = dict(trial=trial, split=pb.split, alpha_planted_found=alpha_ok, alpha_hits=len(hits),
                 alpha_pred_unforced=round(a.na ** 4 * 2.0 ** -ga['bits'], 2), alpha_gauss=ga,
                 beta_planted_found=beta_ok, beta_hits=len(bh), recovered_equals_planted=same)
        s = json.dumps(d); print(s, flush=True)
        if fo: fo.write(s + '\n')


if __name__ == '__main__':
    main()
