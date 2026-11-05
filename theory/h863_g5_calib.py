r"""[E vs M] Calibration of the lattice Gaussian (theory/q523_gauss.gauss_bits) on the G/Z5 level at v = 863, by exact
hit counts of the two-stage join of theory/h863_g5_solve.py (certificate by OpenAI Codex, see h863_lib.py).

For each Z2-split of Z2^2 (q2, q1, q12):
  stage alpha: predicted hits per round = na^4 * 2^-C_alpha, with C_alpha from the Gaussian of the same sampling
               distribution (table-weighted half sums); counted exactly by join22 over R rounds;
  stage beta:  for each of the first B alpha solutions, C_beta|alpha from the Gaussian of the weighted splits actually
               drawn, predicted = nb^4 * 2^-C_beta, counted exactly (uncapped join).
Writes one JSON line per measurement.
usage: python3 theory/h863_g5_calib.py [--rounds R] [--na N] [--nb N] [--nbeta B] [--out F]
"""
import sys, os, json, time, math, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, q523_gauss as QG
from h863_g5_solve import Problem, join22


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cert', default='data/h3452_stage1.json')
    ap.add_argument('--rounds', type=int, default=10)
    ap.add_argument('--na', type=int, default=1024)
    ap.add_argument('--nb', type=int, default=512)
    ap.add_argument('--nbeta', type=int, default=4)
    ap.add_argument('--nest', type=int, default=40000)
    ap.add_argument('--seed', type=int, default=5)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    fo = open(a.out, 'a') if a.out else None
    def log(d):
        s = json.dumps(d); print(s, flush=True)
        if fo: fo.write(s + '\n'); fo.flush()
    for split in ('q2', 'q1', 'q12'):
        pb = Problem(a.cert, split)
        Ta, Tb = pb.targets()
        rng = np.random.default_rng(a.seed)
        # alpha Gaussian
        V = 0
        for i in range(4):
            A, B = pb.sample_halves(i, a.nest, rng)
            V = V + pb.alpha_sig(i, A, B)
        ga = QG.gauss_bits(V, Ta)
        pred_a = a.na ** 4 * 2.0 ** -ga['bits']
        tot = 0; t0 = time.time(); alpha_sols = []
        for r in range(a.rounds):
            H = [pb.sample_halves(i, a.na, rng) for i in range(4)]
            sig = [pb.alpha_sig(i, H[i][0], H[i][1]) for i in range(4)]
            hits = join22(sig[0], sig[1], sig[2], sig[3], Ta)
            tot += len(hits)
            for (p, q, s, t) in hits[:2]:
                if len(alpha_sols) < a.nbeta:
                    pk = [p, q, s, t]
                    alpha_sols.append([(H[i][0][pk[i]], H[i][1][pk[i]]) for i in range(4)])
        log(dict(split=split, stage='alpha', na=a.na, rounds=a.rounds, gauss=ga, pred_per_round=round(pred_a, 2),
                 found=tot, found_per_round=round(tot / a.rounds, 2),
                 log2_found_minus_pred=round(math.log2(max(tot, 0.5) / (pred_a * a.rounds)), 2),
                 secs=round(time.time() - t0, 1)))
        for n, AB in enumerate(alpha_sols):
            opts = []
            for i in range(4):
                Wa, Wb, pa = pb.splits(i, 0, AB[i][0]); Wc, Wd, pc = pb.splits(i, 1, AB[i][1])
                opts.append((Wa, Wb, pa / pa.sum(), Wc, Wd, pc / pc.sum()))
            def draw(n_, rg):
                out = []
                for i, o in enumerate(opts):
                    ka = rg.choice(len(o[0]), size=n_, p=o[2]); kc = rg.choice(len(o[3]), size=n_, p=o[5])
                    out.append(pb.beta_sig(i, o[0][ka] - o[1][ka], o[3][kc] - o[4][kc]))
                return out
            Vb = sum(draw(a.nest, rng))
            gb = QG.gauss_bits(Vb, Tb)
            pred_b = a.nb ** 4 * 2.0 ** -gb['bits']
            bs = draw(a.nb, rng)
            hits = join22(bs[0], bs[1], bs[2], bs[3], Tb, cap=10 ** 7)
            log(dict(split=split, stage='beta', alpha_sol=n, nb=a.nb, gauss=gb, pred=round(pred_b, 2), found=len(hits),
                     log2_found_minus_pred=round(math.log2(max(len(hits), 0.5) / pred_b), 2),
                     nsplits=[[len(o[0]), len(o[3])] for o in opts]))


if __name__ == '__main__':
    main()
