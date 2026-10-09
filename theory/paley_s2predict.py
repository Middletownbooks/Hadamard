r"""Predicted stage-2 yield [E] for recorded stage-1 solutions y of the slice pipeline (theory/paley_slicepipe.py guided1):
for each y, u is sampled uniformly (constant on the twisted multiplier orbits), and the Gaussian density of
(S2 real, X) at (0, 0) over the lag representatives gives log2 E(u | y) = (#twisted orbits) - constraint bits.
Also reports the exact u-list sizes per block (2^#twisted orbits) and the Lemma 6 sums Q1 (psi = 1), Q2 (psi = (-1)^g).
usage: python3 theory/paley_s2predict.py file.npz v p n mu [N]"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_compfix import gauss_logp
from paley_split_b1 import paf, cc
from paley_slicepipe import orbits

def main(path, v, p, n, mu, N=20000, seed=0):
    d = np.load(path); As, Ws = d['A'], d['W']; orbs = orbits(n, mu); rng = np.random.default_rng(seed)
    half = n // 2; out = []
    for i, (A, W) in enumerate(zip(As, Ws)):
        beta = A - W; tw = [[o for o in orbs if W[k, o[0]] == 0] for k in range(4)]
        U = np.zeros((N, 4, n), np.int64)
        for k in range(4):
            for o in tw[k]: U[:, k, o] = rng.choice((1, -1), N)[:, None]
        S2 = (paf(beta[None]).sum(1) + p * paf(U).sum(1))[:, 1:half + 1]
        X = (cc(np.broadcast_to(beta, U.shape), U) - cc(U, np.broadcast_to(beta, U.shape))).sum(1)[:, 1:half]
        c2, r2, _ = gauss_logp(np.concatenate([S2, X], 1), np.zeros(S2.shape[1] + X.shape[1], np.int64), 'S2X')
        bits = sum(len(t) for t in tw)
        rec = dict(y=i, u_orbit_bits=bits, u_lists_log2=[len(t) for t in tw], C2=round(c2, 1), stage2_log2E=round(bits - c2, 1), eqs=r2)
        print(json.dumps(rec), flush=True); out.append(rec)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; main(a[0], *map(int, a[1:5]), *(map(int, a[5:6])))
