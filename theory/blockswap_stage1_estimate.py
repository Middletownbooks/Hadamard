r"""[E] Post-hoc estimate for report/blockswap_answer.md: how many stage-1 (H-level) solutions the (A, B, C, tau C)
slice at v = 863 was expected to hold, per Q-tuple class, against the exact count 0 (Gate 1).

Credit.  Ansatz: the OpenAI model (codex/order-3452-compression, 5b26dae).  863 set-up: OpenAI Codex.  Bordered
construction: Kumar's form of EKLO Theorem 7.  Uses the list builders of theory/blockswap_stage1_exact.py and the
lattice Gaussian theory/q523_gauss.gauss_bits (read at the true target).

For a class representative (zA, zB, zC): N triples drawn uniformly from the energy-consistent set
{(A, B, C) in LA x LB x LC : exc(A) + exc(B) + 2 exc(C) = 2856}; the 20-lag H-level vector of (A, B, C, tau C) is fitted
by the lattice Gaussian; log2 E[solutions] = log2 |energy-consistent set| - bits.  Tail check: the exact fraction of
samples that hit the target on small lag subsets against the Gaussian of the same subset.
usage: python3 theory/blockswap_stage1_estimate.py N SEED   -> logs/blockswap/stage1_estimate.jsonl
"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blockswap_stage1_exact as X
from q523_gauss import gauss_bits

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
REPS = [0, 8, 16, 24, 32]          # one survivor per class (Gate 1 summary)


def run(idx, N, seed):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    zA, zB, zC, _ = [tuple(z) for z in X.survivors()[idx]]
    mA, mB, mC = X.min_fixed(zA), X.min_fixed(zB), X.min_free(zC)
    LA = X.fixed_list(zA, X.EXC - mB - mC); LB = X.fixed_list(zB, X.EXC - mA - mC)
    eA = (LA * LA - 1).sum(1); eB = (LB * LB - 1).sum(1)
    capC = X.EXC - mA - mB
    Rs, Es = [], []
    for s in zC:
        r = X.rows_free(s)
        e = np.array([2 * X.excess(x) for x in r], np.int64)
        o = np.argsort(e, kind='stable')
        Rs.append(np.array(r, np.int64)[o]); Es.append(e[o])
    RK = np.ones(X.NC, np.int64)
    nC = X._enum_C(*Rs, *Es, capC, X.SH, X.TAULAG, RK, np.zeros(0, np.int64), np.zeros((0, 4), np.int32), np.zeros(0, np.int64))
    kC = np.zeros(nC, np.int64); iC = np.zeros((nC, 4), np.int32); eC = np.zeros(nC, np.int64)
    X._enum_C(*Rs, *Es, capC, X.SH, X.TAULAG, RK, kC, iC, eC)
    del kC
    ca = np.bincount(eA); cb = np.bincount(eB); cc = np.bincount(eC)
    trip = []
    for a in np.nonzero(ca)[0]:
        for b in np.nonzero(cb)[0]:
            c = X.EXC - a - b
            if 0 <= c < len(cc) and cc[c]:
                trip.append((a, b, c, float(ca[a]) * float(cb[b]) * float(cc[c])))
    W = np.array([t[3] for t in trip])
    total = W.sum()
    pick = rng.choice(len(trip), size=N, p=W / total)
    byA = {a: np.nonzero(eA == a)[0] for a in np.unique(eA)}
    byB = {b: np.nonzero(eB == b)[0] for b in np.unique(eB)}
    byC = {c: np.nonzero(eC == c)[0] for c in np.unique(eC)}
    V = np.zeros((N, X.NC), np.int64)
    for n, k in enumerate(pick):
        a, b, c, _ = trip[k]
        ia = rng.choice(byA[a]); ib = rng.choice(byB[b]); ic = rng.choice(byC[c])
        Cy = np.concatenate([Rs[q][iC[ic, q]] for q in range(4)])
        Y = np.stack([LA[ia], LB[ib], Cy, Cy[X.TAU]])
        V[n] = X.sig_block(Y).sum(0)
    assert (V[:, 0] == 2936).all()
    g = gauss_bits(V, X.TGT, nlat=600)
    tails = []
    for lags in ([1], [5], [6], [1, 5], [1, 6], [1, 2, 5]):
        hits = int((V[:, lags] == X.TGT[lags]).all(1).sum())
        gb = gauss_bits(V[:, lags], X.TGT[lags], nlat=400)
        tails.append(dict(lags=lags, hits=hits, exact_bits=round(-math.log2(hits / N), 2) if hits else None,
                          gauss_bits=gb['bits']))
    res = dict(index=idx, z=[list(zA), list(zB), list(zC)], N=N, seed=seed,
               log2_energy_consistent=round(math.log2(total), 2), gauss=g,
               log2_expected_solutions=round(math.log2(total) - g['bits'], 2) if g.get('in_lattice') else None,
               tails=tails, secs=round(time.time() - t0, 1), read_at='true target')
    return res


if __name__ == '__main__':
    N, seed = int(sys.argv[1]), int(sys.argv[2])
    for i in REPS:
        r = run(i, N, seed + i)
        s = json.dumps(r)
        print(s, flush=True)
        with open(os.path.join(ROOT, 'logs', 'blockswap', 'stage1_estimate.jsonl'), 'a') as f:
            f.write(s + '\n')
