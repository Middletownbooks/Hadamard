r"""Calibration of the cross-level Gaussian of theory/ab_b0z3.py against exact core counts [E vs M].

For q, y and the dual solutions of ab_b0z3 (same selection K, same order), per dual solution:
  predicted  log2 E(cores) = log2(m0 m1 m2 m3) - bits,  bits = lattice-exact Gaussian (theory/q523_gauss.gauss_bits)
             fitted to N random 4-tuples of fibre elements (sum_k PAF_G at every lag class of G);
  exact      number of cores = the exact cross join (ab_b0z3.cross_join), when --exact.
usage: python3 theory/ab_b0z3_calib.py --q 71 --y 0 [--K 50] [--n 10] [--N 20000] [--exact]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_b0z3 as B
from ab_b0_stage1 import stage1
from q523_gauss import gauss_bits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    q = int(opt('--q', '71')); yi = int(opt('--y', '0')); K = int(opt('--K', '50')); nd = int(opt('--n', '10'))
    N = int(opt('--N', '20000')); exact = '--exact' in a; skip = int(opt('--skip', '0'))
    out = os.path.join(ROOT, 'logs/ab/b0z3/calib.jsonl')
    S = B.SliceZ3(q); sols, _ = stage1(S.sl); Y = sols[yi]
    D = S.dual_lists(Y); idx = B.select(D, K)
    duals, _ = B.dual_join(D, idx, S.Tq[S.qreps])
    score = np.array([sum(np.log2(D[k]['mult'][q4[k]]) for k in range(4)) for q4 in duals])
    order = np.argsort(-score, kind='stable')
    rng = np.random.default_rng(3)
    for r, oi in enumerate(order[skip:skip + nd]):
        q4 = duals[oi]; t0 = time.time()
        Xs = [S.fibre(D[k]['alpha'], D[k]['sigma'], D[k]['vecs'][q4[k]]) for k in range(4)]
        V = np.zeros((N, len(S.greps)), np.int64)
        for c0 in range(0, N, 2000):
            c = min(2000, N - c0); acc = np.zeros((c, S.G.n), np.int64)
            for k in range(4):
                acc += S.G.paf(Xs[k][rng.integers(len(Xs[k]), size=c)].astype(np.int64))
            V[c0:c0 + c] = acc[:, S.greps]
        gb = gauss_bits(V, S.T)
        rec = dict(q=q, y=yi, K=K, rank=skip + r, fibres=[len(x) for x in Xs], log2_fibre_product=round(float(score[oi]), 2),
                   gauss_bits=gb.get('bits'), lattice_rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'),
                   pred_log2E=round(float(score[oi]) - gb['bits'], 2) if gb.get('bits') is not None and np.isfinite(gb['bits']) else None)
        if exact:
            hits, _ = B.cross_join(S, Xs); rec['exact_cores'] = len(hits)
        rec['secs'] = round(time.time() - t0, 1)
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
