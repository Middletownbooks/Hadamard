#!/usr/bin/env python3
"""Exact search: sigma-invariant SINGLE-CHARACTER borders for the 2396 slice (case 0: Walsh values constant on T, so
the compressed core acts on the border space as Gamma (x) I_7).

A single-character row is (label l_c in T, sign s_c) for the blocks c = 0..3.  With characters orthonormal (norm 16):
   <a, b> / 16            = sum_c [l_c = l'_c] s_c s'_c                    (must be 0 for distinct rows of one side)
   K(e, z) = <e, M^T z>/64 = (1/4) sum_{c,c'} Gamma[c][c'] [l_c = l'_{c'}] s_c s'_{c'}   (must be +-1)
sigma acts on labels (T = {0} u {1,2,3} u {4,8,12}, orbits of the dual action) and keeps signs.  A sigma-invariant
side = NFIX fixed rows (all labels 0) + NORB orbits {r, sigma r, sigma^2 r}.  Representatives are rows modulo sigma and
global sign.  Needed: on each side, representatives pairwise orbit-orthogonal and self-orbit-orthogonal; every
E-representative compatible with every Z-row.  Exhaustive depth-first search with candidate-set intersection.  [M]
Usage: fable_border11.py nfix [case_index] [max_nodes]
"""
import sys, json, itertools, time
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border6 as B6, fable_tcensus2 as C2

def main():
    nfix = int(sys.argv[1]); idx = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    maxnodes = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10 ** 9
    d, sol, M, F, t = B6.load('logs/fable_fresh/q599_solutions.json', idx)
    k = 16; b = 4 * t; s_ = 37; n = 4 * (k * s_ + t); T = d['T']
    norb = (b - nfix) // 3; assert 3 * norb + nfix == b
    G = B2.Grp2(4, 1); ch = G.chars()
    def vec(c, m):
        v = np.zeros(4 * k, dtype=np.int64); v[c * k:(c + 1) * k] = ch[m]; return v
    Gam = [np.array([[int(vec(c, m) @ M.T @ vec(c2, m)) // k for c2 in range(4)] for c in range(4)]) for m in T]
    assert all(np.array_equal(Gam[0], g) for g in Gam), 'Walsh values not constant on T: this case needs the general model'
    Gm = Gam[0]
    sig = [x for x in C2.gl_odd_classes() if (x[0], x[1]) == (3, 1)][0]
    dual = sig[3]                                   # sigma on characters
    lab = {m: i for i, m in enumerate(T)}
    sl = [lab[dual[m]] for m in T]                  # sigma on label indices
    assert sorted(sl) == list(range(t))
    # all rows
    rows = [(L, S) for L in itertools.product(range(t), repeat=4) for S in itertools.product((1, -1), repeat=4)]
    rid = {r: i for i, r in enumerate(rows)}
    sg = lambda r: (tuple(sl[l] for l in r[0]), r[1])
    neg = lambda r: (r[0], tuple(-x for x in r[1]))
    reps = []; seen = set(); fixed = []
    for r in rows:
        if r in seen: continue
        orb = {r, sg(r), sg(sg(r))}; orbn = {neg(x) for x in orb}
        seen |= orb | orbn
        (fixed if len(orb) == 1 else reps).append(min(orb | orbn))
    print(f'case {idx}: Gamma {Gm.tolist()}; sigma on labels {sl}; {len(reps)} orbit representatives, {len(fixed)} fixed rows (mod sign)', flush=True)
    def arr(rs):
        return np.array([r[0] for r in rs]), np.array([r[1] for r in rs])
    def ip(A, B):    # <a, b>/16 for all pairs
        (La, Sa), (Lb, Sb) = A, B
        return ((La[:, None, :] == Lb[None, :, :]) * (Sa[:, None, :] * Sb[None, :, :])).sum(2)
    def kk(A, B):    # 4 K(a, b)
        (La, Sa), (Lb, Sb) = A, B
        out = np.zeros((len(La), len(Lb)), dtype=np.int64)
        for c in range(4):
            for c2 in range(4):
                out += Gm[c, c2] * (La[:, None, c] == Lb[None, :, c2]) * (Sa[:, None, c] * Sb[None, :, c2])
        return out
    R0 = arr(reps); R1 = arr([sg(r) for r in reps]); R2 = arr([sg(sg(r)) for r in reps]); FX = arr(fixed)
    nr = len(reps); nf = len(fixed)
    # self-orbit orthogonality of representatives
    selfok = (np.diagonal(ip(R0, R1)) == 0) & (np.diagonal(ip(R0, R2)) == 0)
    # orbit-orthogonality between representatives (a against b, sigma b, sigma^2 b)
    oo = (ip(R0, R0) == 0) & (ip(R0, R1) == 0) & (ip(R0, R2) == 0)
    np.fill_diagonal(oo, False)
    ofx = (ip(R0, FX) == 0)                                   # rep vs fixed row
    ffx = (ip(FX, FX) == 0); np.fill_diagonal(ffx, False)     # fixed vs fixed
    # compatibility: |K| = 1 means kk = +-4; e-rep vs z-rep needs all three images; with fixed rows
    cc = (np.abs(kk(R0, R0)) == 4) & (np.abs(kk(R0, R1)) == 4) & (np.abs(kk(R0, R2)) == 4)
    cfx = (np.abs(kk(R0, FX)) == 4)                           # e-rep vs z-fixed (sigma-images give the same values)
    fxc = (np.abs(kk(FX, R0)) == 4) & (np.abs(kk(FX, R1)) == 4) & (np.abs(kk(FX, R2)) == 4)   # e-fixed vs z-rep orbit
    ffc = (np.abs(kk(FX, FX)) == 4)
    good = selfok
    print(f'self-orthogonal reps {int(good.sum())}; orbit-orthogonal rep pairs {int(oo[np.ix_(good, good)].sum()) // 2}; '
          f'compatible rep pairs {int(cc[np.ix_(good, good)].sum())}', flush=True)
    t0 = time.time(); nodes = [0]; found = []
    gidx = np.nonzero(good)[0]
    # choose the fixed rows first (few), then E reps with running Z-candidate sets
    for EF in itertools.combinations(range(nf), nfix):
        if nfix > 1 and not all(ffx[i, j] for i, j in itertools.combinations(EF, 2)): continue
        for ZF in itertools.combinations(range(nf), nfix):
            if nfix > 1 and not all(ffx[i, j] for i, j in itertools.combinations(ZF, 2)): continue
            if not all(ffc[i, j] for i in EF for j in ZF): continue
            ecand = good.copy(); zcand = good.copy()
            for i in EF: ecand &= ofx[:, i]; zcand &= fxc[i]
            for j in ZF: zcand &= ofx[:, j]; ecand &= cfx[:, j]
            def zclique(zc, need, chosen):
                if need == 0: return chosen
                ids = np.nonzero(zc)[0]
                if len(ids) < need: return None
                for x in ids:
                    nodes[0] += 1
                    m2 = zc & oo[x]; m2[:x + 1] = False
                    r = zclique(m2, need - 1, chosen + [x])
                    if r is not None: return r
                return None
            def eclique(ec, zc, need, chosen):
                if nodes[0] > maxnodes: return None
                if zc.sum() < norb: return None
                if need == 0:
                    zs = zclique(zc, norb, [])
                    return (chosen, zs) if zs is not None else None
                ids = np.nonzero(ec)[0]
                if len(ids) < need: return None
                for x in ids:
                    nodes[0] += 1
                    m2 = ec & oo[x]; m2[:x + 1] = False
                    r = eclique(m2, zc & cc[x], need - 1, chosen + [x])
                    if r is not None: return r
                    if nodes[0] > maxnodes: return None
                return None
            res = eclique(ecand, zcand, norb, [])
            if res is not None:
                found.append((EF, ZF, res)); break
            if nodes[0] > maxnodes: break
        if found or nodes[0] > maxnodes: break
    status = 'FOUND' if found else ('node cap reached: UNDECIDED' if nodes[0] > maxnodes else 'NONE (exhaustive)')
    print(f'sigma-invariant single-character border, NFIX {nfix}, case {idx}: {status}; {nodes[0]} nodes, {time.time()-t0:.0f}s', flush=True)
    if not found: return
    EF, ZF, (Er, Zr) = found[0]
    def mat(rep_ids, fx_ids):
        rs = []
        for i in rep_ids:
            r = reps[i]; rs += [r, sg(r), sg(sg(r))]
        rs += [fixed[j] for j in fx_ids]
        X = np.zeros((len(rs), 4 * k), dtype=np.int64)
        for a, (L, S) in enumerate(rs):
            for c in range(4): X[a, c * k:(c + 1) * k] = S[c] * ch[T[L[c]]]
        return X
    Em = mat(Er, EF); Zm = mat(Zr, ZF)
    Km = (Em @ M.T @ Zm.T) // (4 * k)
    okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
    Fm = -Zm.T
    top = Km @ Km.T + s_ * (Em @ Em.T); cross = Km @ Fm.T + Em @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
    gd = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * k, dtype=int))
    print(f'  K Hadamard {okK}; compressed matrix Gram = {n} I: {gd}', flush=True)
    json.dump(dict(case=idx, nfix=nfix, T=T, s=s_, cells=sol['cells'], M=M.tolist(), K=Km.tolist(), E=Em.tolist(),
                   Z=Zm.tolist(), checks=dict(K_hadamard=bool(okK), compressed_gram=bool(gd))),
              open(f'logs/fable_fresh/q599_sigma_border_exact_f{nfix}_case{idx}.json', 'w'))

if __name__ == '__main__':
    main()
