r"""[M sizes, E counts] The lift (y, w) -> G at v = 863: exact fibre of a G/Z5-level solution w inside the Codex fibre,
and the lattice-Gaussian bits of the remaining (mixed) class equations, by uniform sampling of the doubly fixed
5 x 6 tables (row sums from y, column sums from w; h863_lib.TableSampler).
Certificate y by OpenAI Codex (data/h3452_stage1.json); w from data/h863_g5_solutions.json (theory/h863_g5_solve.py).

Prints per w: log2 fibre (exact), bits, rank, lattice membership (a False would be a congruence obstruction),
Mahalanobis term, log2 expected cores per (y, w).
usage: python3 theory/h863_yw_estimate.py [--sols F] [--n K] [--N samples] [--out F]
"""
import sys, os, json, time, math, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, q523_gauss as QG


def yw_fibre(S, w4):
    """TableSampler per (i, q) and the exact log2 fibre size"""
    ts = {}
    lf = 0.0
    for i in range(4):
        for q in range(4):
            t = L.TableSampler([int(S[i, q * 5 + u]) for u in range(5)], [int(v) for v in w4[i][q]])
            assert t.total > 0
            ts[i, q] = t; lf += math.log2(t.total)
    return ts, lf


def sample_yw(E, ts, N, rng):
    """uniform full configurations in the (y, w) fibre: (N, 4, 860) int8"""
    c = np.empty((N, 4, 20, 6), np.int8)
    for (i, q), t in ts.items():
        c[:, i, q * 5:(q + 1) * 5, :] = t.sample(N, rng)
    return L.assemble(np.broadcast_to(E, (N, 4, 20)), c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cert', default='data/h3452_stage1.json')
    ap.add_argument('--sols', default='data/h863_g5_solutions.json')
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--N', type=int, default=16000)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    Y, E, S = L.load_cert(a.cert)
    W = json.load(open(a.sols))['w'][:a.n]
    rng = np.random.default_rng(a.seed)
    fo = open(a.out, 'a') if a.out else None
    for n, w4 in enumerate(W):
        t0 = time.time()
        ts, lf = yw_fibre(S, w4)
        x = sample_yw(E, ts, a.N, rng)
        V = L.class_paf(x).sum(1)
        # sanity: every sample has the right y and w levels (lags fixed by (y, w) are constant)
        r = QG.gauss_bits(V, L.FULL_TARGET)
        d = dict(w_index=n, log2_fibre=round(lf, 2), N=a.N, bits=r['bits'], rank=r['rank'],
                 in_lattice=r['in_lattice'], mahal=r.get('mahal'),
                 log2E_cores=round(lf - r['bits'], 2) if r['in_lattice'] else None, secs=round(time.time() - t0, 1))
        s = json.dumps(d); print(s, flush=True)
        if fo: fo.write(s + '\n'); fo.flush()


if __name__ == '__main__':
    main()
