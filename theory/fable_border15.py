#!/usr/bin/env python3
"""sigma-invariant borders for the 2396 slice, decided by clique search instead of CP-SAT.  [M, exact]

Same question as theory/fable_border10.py (rows of E and of Z are unions of orbits of sigma, the element of order 3
that defines the slice; NFIX fixed rows and (28 - NFIX)/3 orbits of size 3 on each side; K commutes with the two
permutation actions, so both sides have the same NFIX), where CP-SAT ended UNKNOWN for NFIX = 1, 4, 10.

Method (as in fable_border13/14.py).  Every sigma-invariant border has a sigma-fixed row (28 = 1 mod 3).  Take a
fixed row x1 of E.  All of Z lies in compat(x1) = {y : K(x1, y) = +-1}.  Z is a set of n3 orbit representatives and
NFIX fixed rows, pairwise orthogonal (orbit rows: orthogonal to all three images).  The program peels compat(x1) to
its (n3 - 1)-core for the orbit rows, enumerates the cliques, completes them with fixed rows, and for every complete
Z looks for E among the rows compatible with all of Z.  Exhaustive when it ends below NODECAP.
Rows of size-3 orbits use only block functions that are not sigma-fixed [P, report 7.1].
Control: 'singer' (Z_2^2, spike compressions, sigma the 3-cycle; NFIX = 0 or 3).
Usage: fable_border15.py q599 IDX NFIX [nodecap] | singer NFIX
"""
import sys, json, time, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border6 as B6, fable_tcensus2 as C2, fable_clique as FC

def main():
    if sys.argv[1] == 'q599':
        idx = int(sys.argv[2]); nfix = int(sys.argv[3]); nodecap = int(sys.argv[4]) if len(sys.argv) > 4 else 50_000_000
        d, sol, M, F, t = B6.load('logs/fable_fresh/q599_solutions.json', idx)
        k = 16; s_ = 37; T = d['T']; name = f'q599_case{idx}'
        img = [x for x in C2.gl_odd_classes() if (x[0], x[1]) == (3, 1)][0][2]
        M = np.array(M, dtype=np.int64); F = np.array(F, dtype=np.int64)
    else:
        nfix = int(sys.argv[2]); nodecap = 50_000_000
        k = 4; s_ = 7; G = B2.Grp2(2, 1); ch = G.chars(); T = [1, 2, 3]; t = 3; name = 'singer'
        ys = [np.array(y, dtype=np.int64) for y in ([1, 1, 1, 1], [3, 1, 1, 1], [3, 1, 1, 1], [3, 1, 1, 1])]
        M = FB.gs_array(G, np.array(ys)); img = [0, 2, 3, 1]
        F = np.array([f for f in (np.array([1 - 2 * ((bits >> q) & 1) for q in range(k)]) for bits in range(1 << k))
                      if all(int(ch[mm] @ f) == 0 for mm in range(k) if mm not in T)], dtype=np.int64)
        sol = dict(cells=[y.tolist() for y in ys])
    b = 4 * t; n = 4 * (k * s_ + t); assert (b - nfix) % 3 == 0
    n3 = (b - nfix) // 3
    Ps = np.zeros((k, k), dtype=np.int64)
    for q in range(k): Ps[q, img[q]] = 1
    P = np.kron(np.eye(4, dtype=np.int64), Ps)
    assert np.array_equal(P @ M, M @ P), 'M does not commute with sigma'
    no = len(F); key = {tuple(f): i for i, f in enumerate(F)}
    pi = np.array([key[tuple(F[o] @ Ps.T)] for o in range(no)]); neg = np.array([key[tuple(-f)] for f in F])
    fixed = [o for o in range(no) if pi[o] == o]; nonfix = [o for o in range(no) if pi[o] != o]
    print(f'{name}: b = {b}, NFIX = {nfix}, orbits of size 3: {n3}; {no} block functions, {len(fixed)} sigma-fixed', flush=True)
    code = lambda R: ((R[:, 0] * no + R[:, 1]) * no + R[:, 2]) * no + R[:, 3]
    vec = lambda R: np.concatenate([F[R[:, c]] for c in range(4)], axis=1).astype(np.float32)
    # orbit rows: all four blocks non-fixed; representative = least code among the six images (+-, sigma^j)
    R0 = np.array(list(itertools.product(nonfix, repeat=4)), dtype=np.int64) if nonfix else np.zeros((0, 4), dtype=np.int64)
    if len(R0):
        ims = [R0, pi[R0], pi[pi[R0]]]; ims += [neg[x] for x in ims]
        cs = np.stack([code(x) for x in ims], axis=1)
        R0 = R0[cs[:, 0] == cs.min(axis=1)]
    OV = [vec(R0), vec(pi[R0]), vec(pi[pi[R0]])] if len(R0) else [np.zeros((0, 4 * k), dtype=np.float32)] * 3
    # each orbit row must be orthogonal to its images
    if len(R0):
        ok = (np.einsum('ij,ij->i', OV[0], OV[1]) == 0)
        R0 = R0[ok]; OV = [v[ok] for v in OV]
    Rf = np.array(list(itertools.product(fixed, repeat=4)), dtype=np.int64).reshape(-1, 4)
    cf = np.stack([code(Rf), code(neg[Rf])], axis=1); Rf = Rf[cf[:, 0] == cf.min(axis=1)]
    FV = vec(Rf)
    NO, NF = len(R0), len(Rf)
    print(f'orbit representatives {NO}, fixed rows (mod sign) {NF}', flush=True)
    MT = M.T.astype(np.float32); Mf = M.astype(np.float32); fk = float(4 * k)
    # K(x, y) = x M^T y^T / 4k
    def comp_from_E(xv):          # Z-candidates compatible with an E-row vector xv (and its images, if an orbit row)
        w = xv @ MT
        co = np.ones(NO, dtype=bool)
        for j in range(3): co &= (np.abs(OV[j] @ w) == fk)
        cfx = np.abs(FV @ w) == fk
        return co, cfx
    def comp_from_Z(yv):          # E-candidates compatible with a Z-row vector yv
        w = MT @ yv
        co = np.ones(NO, dtype=bool)
        for j in range(3): co &= (np.abs(OV[j] @ w) == fk)
        cfx = np.abs(FV @ w) == fk
        return co, cfx
    nodes = [0]; t0 = time.time(); found = []; exhausted = True; nz = 0; trunc = [False]
    def orbit_cliques(idx, need):
        if need == 0: return [[]]
        V = [OV[j][idx] for j in range(3)]
        A = np.ones((len(idx), len(idx)), dtype=bool)
        for j in range(3): A &= ((V[0] @ V[j].T) == 0)
        np.fill_diagonal(A, False)
        alive = np.ones(len(idx), dtype=bool)
        while True:
            bad = alive & ((A & alive[None, :]).sum(axis=1) < need - 1)
            if not bad.any(): break
            alive &= ~bad
        cand = [c for c in range(len(idx)) if alive[c]]
        if len(cand) < need: return []
        out, nn, tr = FC.kcliques(A, need, cand=cand, maxout=200000, maxnodes=max(1, nodecap - nodes[0]))
        nodes[0] += nn
        if tr: trunc[0] = True
        return [[int(idx[c]) for c in cl] for cl in out]
    def fixed_cliques(idx, need, must=None):
        if need == 0: return [[]]
        V = FV[idx]; A = (V @ V.T) == 0; np.fill_diagonal(A, False)
        if must is None:
            out, nn, tr = FC.kcliques(A, need, maxout=200000, maxnodes=max(1, nodecap - nodes[0]))
        else:
            m = int(np.nonzero(idx == must)[0][0])
            cand = [c for c in range(len(idx)) if A[m, c]]
            out, nn, tr = FC.kcliques(A, need - 1, cand=cand, maxout=200000, maxnodes=max(1, nodecap - nodes[0]))
            out = [[m] + cl for cl in out]
        nodes[0] += nn
        if tr: trunc[0] = True
        return [[int(idx[c]) for c in cl] for cl in out]
    def complete_sides(co, cfx, need3, needf, must_fixed=None):
        """all (orbit list, fixed list) with pairwise orthogonality inside the candidate sets.
        The fixed rows are chosen first (they cut the orbit candidates down), then the orbit cliques."""
        res = []
        if int(co.sum()) < need3: return res
        ifx = np.nonzero(cfx)[0]
        if len(ifx) < needf: return res
        if needf and must_fixed is not None and not cfx[must_fixed]: return res
        for fc in (fixed_cliques(ifx, needf, must_fixed) if needf else [[]]):
            co2 = co.copy()
            for f in fc: co2 &= ((OV[0] @ FV[f]) == 0)
            io = np.nonzero(co2)[0]
            if len(io) < need3: continue
            for oc in orbit_cliques(io, need3):
                res.append((oc, fc))
            if nodes[0] > nodecap or trunc[0]: break
        return res
    firsts = range(NF) if nfix > 0 else range(NO)
    ck = f'logs/fable_fresh/ckpt_b15_{name}_f{nfix}.json'; start = 0; lastp = time.time()
    try:
        cj = json.load(open(ck)); start = cj['done']; nz = cj['nz']; nodes[0] = cj['nodes']
        print(f'resuming after {start} first rows (checkpoint)', flush=True)
    except Exception: pass
    for cnt, x1 in enumerate(firsts):
        if cnt < start: continue
        if nodes[0] > nodecap or trunc[0]: exhausted = False; break
        if cnt > start and (time.time() - lastp > 120 or cnt % 20 == 0):
            json.dump(dict(done=cnt, nz=nz, nodes=nodes[0]), open(ck, 'w')); lastp = time.time()
            print(f'  first row {cnt}/{len(firsts)}: complete Z so far {nz}, nodes {nodes[0]}, {time.time()-t0:.0f}s', flush=True)
        xv = FV[x1] if nfix > 0 else OV[0][x1]
        co, cfx = comp_from_E(xv)
        for Zo, Zf in complete_sides(co, cfx, n3, nfix):
            nz += 1
            eo = np.ones(NO, dtype=bool); ef = np.ones(NF, dtype=bool)
            for o in Zo:
                for j in range(3):
                    a_, b_ = comp_from_Z(OV[j][o]); eo &= a_; ef &= b_
            for f in Zf:
                a_, b_ = comp_from_Z(FV[f]); eo &= a_; ef &= b_
            if nfix > 0:
                if not ef[x1]: continue
                cand = complete_sides(eo, ef, n3, nfix, must_fixed=x1)
            else:
                if not eo[x1]: continue
                cand = [c for c in complete_sides(eo, ef, n3, nfix) if x1 in c[0]]
            if cand:
                found.append((cand[0], (Zo, Zf))); break
        if found: break
    tag = f'sigma-invariant border ({name}, NFIX {nfix})'
    if trunc[0]: exhausted = False
    if not found:
        print(f'{tag}: ' + ('NONE (search exhaustive)' if exhausted and nodes[0] <= nodecap else f'UNDECIDED (a cap was hit: nodes {nodecap} or 200000 cliques in one list)')
              + f'; complete Z found {nz}, nodes {nodes[0]}, {time.time()-t0:.0f}s', flush=True)
        return
    (Eo, Ef), (Zo, Zf) = found[0]
    def expand(Lo, Lf):
        rows = []
        for o in Lo:
            r = R0[o]; rows += [r, pi[r], pi[pi[r]]]
        rows += [Rf[f] for f in Lf]
        return np.array([np.concatenate([F[r[c]] for c in range(4)]) for r in rows], dtype=np.int64)
    Em = expand(Eo, Ef); Zm = expand(Zo, Zf)
    Km = (Em @ M.T @ Zm.T) // (4 * k)
    okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=np.int64)) \
        and np.array_equal(Em @ M.T @ Zm.T, 4 * k * Km)
    Fm = -Zm.T
    top = Km @ Km.T + s_ * (Em @ Em.T); cross = Km @ Fm.T + Em @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
    good = np.array_equal(top, n * np.eye(b, dtype=np.int64)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * k, dtype=np.int64))
    print(f'{tag}: FOUND (nodes {nodes[0]}, {time.time()-t0:.0f}s); K Hadamard {okK}; compressed matrix of order {b + 4*k}: Gram = {n} I is {good}', flush=True)
    json.dump(dict(name=name, nfix=nfix, T=T, s=s_, cells=sol['cells'], M=M.tolist(), K=Km.tolist(), E=Em.tolist(), Z=Zm.tolist(),
                   checks=dict(K_hadamard=bool(okK), compressed_gram=bool(good))),
              open(f'logs/fable_fresh/sigma15_border_{name}_f{nfix}.json', 'w'))

if __name__ == '__main__':
    main()
