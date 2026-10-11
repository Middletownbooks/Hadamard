r"""S/D-split cost model with per-pair S-alphabet restrictions [E] (report/paley_523_answer.md, section 3).

Pairs of Z_q-cosets (g, g + F2) of Z_F = Z_2 x Z_F2 carry S = Z(g) + Z(g + F2) and D = Z(g) - Z(g + F2) (eta-coords).
Given y, a pair has a finite set of (c, c') combinations; grouping by S gives the distinct S-values s with
multiplicity m(s) = number of D-options (the D-freedom given S).  Sampling model used by the joins:
    S uniform over the ALLOWED distinct S-values of each pair (the S-lists hold distinct configurations),
    D uniform over the m(s) options given S.
Alphabet restrictions ("min_mult"): only S-values with m(s) >= min_mult are allowed (m = 1 values are the pairs with
c = c', which leave no D-freedom).
Reported (Gaussian, lattice-exact): C_S (S-level, characters trivial on 9), C_Z3 and C_Z1 (compressions), C_mid =
C_S - C_Z3 (the order-F2 S-step), C(D|S), h_D = E sum log2 m, and the cost exponent of the S/D route with one big-fibre
Z3 parent:  cost = (C_mid + C(D|S) - h_D) / 2, with the fibre size it needs, (C_mid + C(D|S) - h_D) / 4 per block.
usage: python3 theory/q523_sdmodel.py F q s y.json [--N 40000] [--min-mult 1]"""
import sys, os, json, math, itertools as it
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, pz, target_P
from q523_levels import level_matrix
from q523_gauss import gauss_bits

def pair_tables(cy, Y, F, min_mult=1):
    """per block, per pair g < F2: list of (S, [D options]) with D signed by (-1)^g, restricted to m >= min_mult"""
    F2 = F // 2; out = []
    for k in range(4):
        alph = [zcoords(cy.alphabet(int(y))) for y in Y[k]]
        blk = []
        for g in range(F2):
            sg = 1 if g % 2 == 0 else -1
            byS = {}
            for a in alph[g]:
                for b in alph[g + F2]:
                    S = tuple((a + b).tolist()); D = tuple((sg * (a - b)).tolist())
                    byS.setdefault(S, []).append(D)
            blk.append([(S, Ds) for S, Ds in byS.items() if len(Ds) >= min_mult])
        out.append(blk)
    return out

def sample(tables, F, s, N, rng):
    """returns Z (N, 4, F, s) built from (S, D) samples: Z(g) = (S + D')/2, Z(g+F2) = (S - D')/2 with D' = sg D"""
    F2 = F // 2; Z = np.zeros((N, 4, F, s), np.int64); hD = np.zeros(N)
    for k in range(4):
        for g in range(F2):
            sg = 1 if g % 2 == 0 else -1
            L = tables[k][g]; i = rng.integers(len(L), size=N)
            Ss = np.array([S for S, _ in L]); m = np.array([len(Ds) for _, Ds in L])
            j = (rng.random(N) * m[i]).astype(np.int64)
            Dsel = np.array([L[a][1][b] for a, b in zip(i, j)])
            S = Ss[i]; Dp = sg * Dsel
            Z[:, k, g] = (S + Dp) // 2; Z[:, k, g + F2] = (S - Dp) // 2
            hD += np.log2(m[i])
    return Z, hD

def main(F, q, s, Y, N=40000, min_mult=1, seed=1):
    cy = Cyclo(q, s); v = F * q + 1; rng = np.random.default_rng(seed)
    tab = pair_tables(cy, Y, F, min_mult)
    Z, hD = sample(tab, F, s, N, rng)
    P = np.zeros((N, F, s), np.int64)
    for c0 in range(0, N, 4000): P[c0:c0 + 4000] = pz(Z[c0:c0 + 4000], cy.T).sum(1)
    P = P.reshape(N, -1); T = target_P(F, v, s).ravel()
    gal = [d for d in range(1, s + 1) if s % d == 0]
    divs = [d for d in range(1, F + 1) if F % d == 0]; F2 = F // 2
    Sd = [d for d in divs if F2 % d == 0]
    def C(ds):
        M = level_matrix(F, s, ds, gal); return gauss_bits(P @ M.T, M @ T)["bits"]
    res = dict(v=v, min_mult=min_mult, N=N)
    res["C_all"] = C(divs); res["C_S"] = C(Sd); res["C_Z1"] = C([1])
    if F2 % 3 == 0: res["C_Z3"] = C([1, 3])
    else: res["C_Z3"] = res["C_Z1"]
    res["C_mid"] = round(res["C_S"] - res["C_Z3"], 2); res["C(D|S)"] = round(res["C_all"] - res["C_S"], 2)
    res["h_D"] = round(float(hD.mean()), 2)
    res["H_S_lists"] = round(sum(math.log2(len(tab[k][g])) for k in range(4) for g in range(F2)), 2)
    res["yield"] = round(res["H_S_lists"] + res["h_D"] - res["C_all"], 2)
    res["D_deficit"] = round(res["C(D|S)"] - res["h_D"], 2)
    res["cost_log2"] = round((res["C_mid"] + max(0.0, res["D_deficit"])) / 2, 2)
    res["fibre_needed_per_block_log2"] = round((res["C_mid"] + max(0.0, res["D_deficit"])) / 4, 2)
    return res

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "40000", "--min-mult": "1", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y["y"] if isinstance(Y, dict) else Y
    print(json.dumps(main(F, q, s, Y, int(opt["--N"]), int(opt["--min-mult"]), int(opt["--seed"]))))
