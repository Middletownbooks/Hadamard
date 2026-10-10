r"""Exact lattice of achievable stage-2 values for a given stage-1 y (report/paley_523_answer.md, section 5) [P, M].

V(c) = sum_k P_k(c_k) (eta-coordinates of sum_k P(Z_k) on Z_F) is a quadratic function of the per-coset choices, and
blocks do not interact.  Hence, with a base configuration c0,
    Lambda := span_Z { V(c) - V(c0) : c }  =  span_Z { first differences D_i(z') } + span_Z { cross terms X_ij(z', z'') },
D_i(z') = V(c0 with coset i -> z') - V(c0),  X_ij = V(c0; i -> z', j -> z'') - V(c0; i -> z') - V(c0; j -> z'') + V(c0)
(i, j in the same block).  [P]: V(c) - V(c0) = sum_i D_i + sum_{i<j} X_ij for the changed cosets, and every
generator is itself in Lambda.  Exact consequence: if T - V(c0) is not in Lambda, the y has no stage-2 completion.
Also reports the index of the sampled-difference lattice (theory/q523_gauss.py) in the exact one, which validates the
lattice used by the Gaussian estimates.
usage: python3 theory/q523_lattice.py F q s y.json [y2.json ...]"""
import sys, os, json, math, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, pz, target_P
from q523_gauss import echelon, coords
from q523_levels import sample_P

def exact_lattice(F, q, s, Y, seed=0):
    cy = Cyclo(q, s); rng = np.random.default_rng(seed)
    alph = [[zcoords(cy.alphabet(int(y))) for y in Yk] for Yk in Y]
    Z0 = np.array([[alph[k][g][0] for g in range(F)] for k in range(4)])
    base = pz(Z0, cy.T).sum(0).ravel()
    gens = []
    for k in range(4):
        Pk0 = pz(Z0[k], cy.T).ravel()
        single = {}
        for g in range(F):
            for t in range(1, len(alph[k][g])):
                Z = Z0[k].copy(); Z[g] = alph[k][g][t]
                d = pz(Z, cy.T).ravel() - Pk0; single[(g, t)] = d; gens.append(d)
        keys = list(single)
        for (g1, t1), (g2, t2) in it.combinations(keys, 2):
            if g1 == g2: continue
            Z = Z0[k].copy(); Z[g1] = alph[k][g1][t1]; Z[g2] = alph[k][g2][t2]
            x = pz(Z, cy.T).ravel() - Pk0 - single[(g1, t1)] - single[(g2, t2)]
            if x.any(): gens.append(x)
    B, piv = echelon(gens)
    return cy, base, B, piv

def covol_log2(B, piv):
    return sum(math.log2(abs(B[i][p])) for i, p in enumerate(piv)) if B else 0.0

def main(F, q, s, Y):
    v = F * q + 1
    cy, base, B, piv = exact_lattice(F, q, s, Y)
    T = target_P(F, v, s).ravel()
    ct = coords(B, piv, (T - base).tolist())
    # sampled lattice (as in the Gaussian estimator) and its index in the exact lattice
    P = sample_P(cy, Y, 401, np.random.default_rng(1)).reshape(401, -1)
    Bs, ps = echelon((P[1:] - P[0]).tolist())
    same_rank = len(Bs) == len(B)
    idx = round(covol_log2(Bs, ps) - covol_log2(B, piv), 3) if same_rank and ps == piv else None
    return dict(v=v, rank=len(B), target_in_lattice=ct is not None, sampled_rank=len(Bs), sampled_index_log2=idx)

if __name__ == "__main__":
    F, q, s = map(int, sys.argv[1:4])
    for f in sys.argv[4:]:
        d = json.load(open(f)); Y = d['y'] if isinstance(d, dict) else d
        if Y is None: continue
        r = main(F, q, s, Y); r['file'] = os.path.basename(f); print(json.dumps(r), flush=True)
