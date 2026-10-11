r"""Hierarchical (level-by-level) stage-2 cost model for the generalised-Paley slice on Z_F x Z_q, F = 2 * F2 [E + exact].

Levels are nested sets of Z_F-character orders (theory/q523_levels.py).  For each level the per-block data that
determines its equations is exact:
    {1}           block total (Z1 compression)          {1, 3}        Z3 compression (F2 divisible by 3)
    S = {d | F2}  S(g) = Z(g) + Z(g + F2) on Z_F2       S + {2}       (S, D-total),  D-total = sum_g (-1)^g Z(g)
    S + {2, 6}    (S, D3), D3 = Z3 compression of D(g) = (-1)^g (Z(g) - Z(g + F2))      all     the configuration
Per-block entropies are computed exactly where the data factor over cosets (compressions, S, (S, D3) by rows of three
pairs) and by Monte Carlo over S for H(D-total | S) (exact conditional convolution per sample).  Constraint bits per
level come from the lattice-exact Gaussian (theory/q523_levels.py).

Cost model (log2, Gaussian): level i adds per-block fibre entropy dh_i (summed over blocks: dH_i) and constraint dc_i,
yield y_i = dH_i - dc_i.  Working back from one final solution, N_{i-1} = max(1, N_i / 2^{y_i}) parents are needed;
each parent costs a 2+2 meet in the middle over its fibres: 2^{(dh_0 + dh_1)} pairs per side when every child is
needed (full fibres), or 2^{(dc_i + log N_i)/2} with balanced sublists when one parent suffices.
usage: python3 theory/q523_chain.py F q s y.json [--N 40000] [--mc 2000]"""
import sys, os, json, math, itertools as it
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, target_P
from q523_levels import sample_P, level_matrix, dist_sum, H
from q523_gauss import gauss_bits

def per_block_data_entropies(cy, Yk, F, rng, mc=2000):
    F2 = F // 2
    alph = [[tuple(r) for r in zcoords(cy.alphabet(int(y)))] for y in Yk]
    add = lambda a, b: tuple(x + y for x, y in zip(a, b)); sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
    out = {}
    # compressions Z1, Z3 (exact)
    for d in (1, 3):
        if F2 % d and d != 1: continue
        out[f"Z{d}"] = sum(H(dist_sum([alph[g] for g in range(r, F, d)])) for r in range(d))
    # per pair (g, g + F2): joint distribution of (S, D) with D = (-1)^g (Z(g) - Z(g+F2)) (g < F2 so sign by parity)
    pairs = []
    for g in range(F2):
        sg = 1 if g % 2 == 0 else -1
        J = Counter()
        for a in alph[g]:
            for b in alph[g + F2]:
                d = sub(a, b); d = d if sg == 1 else tuple(-x for x in d)
                J[(add(a, b), d)] += 1
        n = sum(J.values()); pairs.append({k: v / n for k, v in J.items()})
    HS = 0.0; HSD = 0.0
    for P in pairs:
        cS = Counter()
        for (S, D), p in P.items(): cS[S] += p
        HS += H(cS); HSD += H(Counter(P))
    out["S"] = HS; out["full"] = HSD
    # (S, D3): rows of pairs g = r, r + 3, r + 6 (mod F2), exact
    if F2 % 3 == 0:
        HSD3 = 0.0
        for r in range(3):
            rows = [pairs[g] for g in range(r, F2, 3)]
            J = Counter()
            for combo in it.product(*[list(P.items()) for P in rows]):
                S = tuple(c[0][0] for c in combo); D3 = tuple(sum(x) for x in zip(*[c[0][1] for c in combo]))
                p = 1.0
                for c in combo: p *= c[1]
                J[(S, D3)] += p
            HSD3 += H(J)
        out["S+D3"] = HSD3
    # (S, D-total): Monte Carlo over S, exact conditional distribution of D-total given S
    cond = []
    for P in pairs:
        byS = {}
        for (S, D), p in P.items(): byS.setdefault(S, Counter())[D] += p
        Ss = list(byS); w = np.array([sum(byS[S].values()) for S in Ss])
        cond.append((Ss, w / w.sum(), [Counter({D: p / sum(byS[S].values()) for D, p in byS[S].items()}) for S in Ss]))
    hs = []
    for _ in range(mc):
        dist = Counter({tuple([0] * cy.s): 1.0})
        for Ss, w, cd in cond:
            i = rng.choice(len(Ss), p=w); nxt = Counter()
            for a, pa in dist.items():
                for b, pb in cd[i].items(): nxt[add(a, b)] += pa * pb
            dist = nxt
        hs.append(H(dist))
    out["S+D1"] = HS + float(np.mean(hs))
    return {k: round(v, 3) for k, v in out.items()}

CHAINS = {
    "SD": [("Z1", [1]), ("Z3", [1, 3]), ("S", [1, 3, 9]), ("S+D1", [1, 2, 3, 9]), ("S+D3", [1, 2, 3, 6, 9]), ("full", [1, 2, 3, 6, 9, 18])],
    "SD_noZ3": [("Z1", [1]), ("S", [1, 3, 9]), ("S+D1", [1, 2, 3, 9]), ("S+D3", [1, 2, 3, 6, 9]), ("full", [1, 2, 3, 6, 9, 18])],
    "SD_direct": [("Z1", [1]), ("Z3", [1, 3]), ("S", [1, 3, 9]), ("full", [1, 2, 3, 6, 9, 18])],
}

def chain_cost(levels, blocks_h, C):
    """levels: list of names; blocks_h[k][name] per-block entropies; C[name] constraint bits.  Returns per-level rows
    and the log2 total cost (sum of pair-enumeration work)."""
    rows = []; prevH = [0.0] * 4; prevC = 0.0
    for nm in levels:
        dh = [blocks_h[k][nm] - prevH[k] for k in range(4)]; dc = C[nm] - prevC
        rows.append(dict(level=nm, dh=[round(x, 2) for x in dh], dH=round(sum(dh), 2), dc=round(dc, 2), y=round(sum(dh) - dc, 2)))
        prevH = [blocks_h[k][nm] for k in range(4)]; prevC = C[nm]
    N = 0.0; total = []
    for r in reversed(rows):
        need_parents = max(0.0, N - r["y"])
        dhs = sorted(r["dh"], reverse=True); full_side = dhs[0] + dhs[3]   # pair the largest with the smallest
        full_side = max(dhs[0] + dhs[3], dhs[1] + dhs[2])
        if need_parents > 0:
            per = full_side; cost = need_parents + per
        else:
            per = min(full_side, (r["dc"] + N) / 2); cost = per
        r.update(parents_needed_log2=round(need_parents, 2), cost_log2=round(cost, 2))
        total.append(cost); N = need_parents
    tot = math.log2(sum(2 ** c for c in total))
    return rows, round(tot, 2)

def main(F, q, s, Y, N=40000, mc=2000, seed=1):
    cy = Cyclo(q, s); v = F * q + 1; rng = np.random.default_rng(seed)
    bh = [per_block_data_entropies(cy, Y[k], F, rng, mc) for k in range(4)]
    P = sample_P(cy, Y, N, rng).reshape(N, -1); T = target_P(F, v, s).ravel()
    names = {}
    for ch in CHAINS.values():
        for nm, ds in ch: names[nm] = ds
    C = {}
    for nm, ds in names.items():
        M = level_matrix(F, s, ds, [d for d in range(1, s + 1) if s % d == 0])
        C[nm] = gauss_bits(P @ M.T, M @ T)["bits"]
    out = dict(v=v, block_entropies=bh, C=C, chains={})
    for cn, ch in CHAINS.items():
        rows, tot = chain_cost([nm for nm, _ in ch], bh, C)
        out["chains"][cn] = dict(rows=rows, total_log2=tot)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "40000", "--mc": "2000", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y["y"] if isinstance(Y, dict) else Y
    print(json.dumps(main(F, q, s, Y, int(opt["--N"]), int(opt["--mc"]), int(opt["--seed"]))))
