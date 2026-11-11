#!/usr/bin/env python3
"""Border for a Q-level solution by 'fix E, enumerate Z' (see fable_border5.py for the setting).
All candidate rows = 4 blocks, each a +-1 function with Walsh support in T (30^4 for the 599 slice).
For a trial basis E (28 orthogonal candidate rows) the admissible Z-rows are the candidates z with
E M^T z^T in {+-4k}^b; a border is an orthogonal b-subset of them.  Trial bases: tensor bases (rows of a 4 x 4
Hadamard matrix times an orthogonal basis of block functions) and random orthogonal bases grown greedily.  [M]
Usage: fable_border6.py solutions.json index s ntrials seed
"""
import sys, json, time, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2

def load(path, idx):
    d = json.load(open(path)); T = d['T']; K_ = 16; a = 4
    sol = d['solutions'][idx]; cellorb = d['cell_orbits']
    G = B2.Grp2(a, 1); ch = G.chars()
    pos = [G.Qel.index(tuple((q >> (a - 1 - i)) & 1 for i in range(a))) for q in range(K_)]
    ys = []
    for c in sol['cells']:
        y = np.zeros(K_, dtype=np.int64)
        for j, ob in enumerate(cellorb):
            for q in ob: y[pos[q]] = c[j]
        ys.append(y)
    M = FB.gs_array(G, np.array(ys))
    E_ = [sum(int(ch[m] @ y) ** 2 for y in ys) for m in range(K_)]
    low = [m for m in range(K_) if E_[m] == 4 * len(T)]
    funcs = []
    for bits in range(1 << K_):
        f = np.array([1 - 2 * ((bits >> q) & 1) for q in range(K_)], dtype=np.int8)
        if all(int(ch[m] @ f) == 0 for m in range(K_) if m not in low): funcs.append(f)
    return d, sol, M, np.array(funcs), len(T)

def orth_clique(Gm, need, order, maxnodes=3_000_000):
    n = Gm.shape[0]; nodes = [0]
    def rec(cset, chosen):
        if len(chosen) == need: return chosen
        idx = np.nonzero(cset)[0]
        if len(idx) + len(chosen) < need: return None
        nodes[0] += 1
        if nodes[0] > maxnodes: return None
        for x in idx:
            m2 = cset & Gm[x]; m2[:x + 1] = False
            r = rec(m2, chosen + [x])
            if r is not None: return r
            if nodes[0] > maxnodes: return None
        return None
    return rec(np.ones(n, dtype=bool), []), nodes[0]

def main():
    path = sys.argv[1]; idx = int(sys.argv[2]); s_ = int(sys.argv[3]); ntr = int(sys.argv[4]); seed = int(sys.argv[5])
    d, sol, M, F, t = load(path, idx)
    K_ = 16; b = 4 * t; n = 4 * (K_ * s_ + t); no = len(F)
    rng = np.random.default_rng(seed)
    cand_idx = np.array(list(itertools.product(range(no), repeat=4)), dtype=np.int16)
    C = np.concatenate([F[cand_idx[:, a]] for a in range(4)], axis=1).astype(np.int8)       # (no^4, 64)
    print(f'{no} block functions, {len(C)} candidate rows', flush=True)
    MT = M.T.astype(np.int32)
    H4s = [FB.sylv(2)]
    # orthogonal bases of block functions inside F (7 mutually orthogonal functions): characters and random ones
    FG = (F.astype(np.int32) @ F.T.astype(np.int32)) == 0
    def block_bases(k):
        out = []
        for _ in range(k):
            order = rng.permutation(no); chosen = []
            for x in order:
                if all(FG[x, y] for y in chosen): chosen.append(int(x))
                if len(chosen) == t: break
            if len(chosen) == t: out.append(chosen)
        return out
    trials = []
    for bb in block_bases(max(1, ntr // 2)):
        for H4 in H4s:
            rows = [np.concatenate([H4[i, a] * F[m] for a in range(4)]) for i in range(4) for m in bb]
            trials.append(('tensor', np.array(rows)))
    for _ in range(ntr - len(trials)):
        # random greedy orthogonal basis from candidates
        chosen = []
        pool = rng.permutation(len(C))[:200000]
        for x in pool:
            v = C[x].astype(np.int32)
            if all(int(v @ c) == 0 for c in chosen): chosen.append(v)
            if len(chosen) == b: break
        if len(chosen) == b: trials.append(('greedy', np.array(chosen)))
    best = 0
    for ti, (kind_, E) in enumerate(trials):
        G_ = (E.astype(np.int32) @ MT)                           # (b, 64)
        vals = C.astype(np.int32) @ G_.T                         # (cands, b)
        adm = np.nonzero(np.all(np.abs(vals) == 4 * K_, axis=1))[0]
        best = max(best, len(adm))
        msg = f'trial {ti} ({kind_}): {len(adm)} admissible Z-rows'
        if len(adm) >= b:
            A = C[adm].astype(np.int32); Gm = (A @ A.T) == 0
            cl, nodes = orth_clique(Gm, b, None)
            msg += f'; orthogonal {b}-subset: {"FOUND" if cl else "none"} ({nodes} nodes)'
            if cl:
                Z = A[cl]; Km = (E.astype(np.int64) @ M.T @ Z.T) // (4 * K_)
                okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
                Fm = -Z.T
                top = Km @ Km.T + s_ * (E @ E.T); cross = Km @ Fm.T + E @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
                good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * K_, dtype=int))
                print(msg + f'; K Hadamard {okK}; compressed matrix Gram = {n} I: {good}', flush=True)
                json.dump(dict(T=d['T'], s=s_, cells=sol['cells'], K=Km.tolist(), E=E.tolist(), Z=Z.tolist()),
                          open(f'logs/fable_fresh/q599_border_{idx}.json', 'w'))
                return
        print(msg, flush=True)
    print('no border found in these trials (not a zero); most admissible rows in a trial:', best)
if __name__ == '__main__':
    main()
