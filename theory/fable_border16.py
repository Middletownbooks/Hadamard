#!/usr/bin/env python3
"""Borders with a symmetry of order 3, for T = two 2-flats of characters of Z_2^4 meeting in the trivial character
and a constant Walsh quadruple (4, 2, 2, 2) on T (the 2396 slice of report 7.1, in normal form).  [M, exact]

On the border space the compressed core acts as Gamma (x) I with Gamma = GS(4, 2, 2, 2) (any signed permutation of
(4, 2, 2, 2) is equivalent to it: automorphisms of the quaternions, unit multiplications, transposition), so the
problem is posed here with the synthetic M = Gamma (x) I_16: K(x, y) = x M^T y^T / 64 for admissible rows x, y.

Symmetry group of this problem (apart from E <-> Z dualities): (Q8 x| Z_3) x Aut_Q, where Q8 x| Z_3 is the group of
24 signed block permutations commuting with Gamma (right multiplications by +-1, +-i, +-j, +-k and the rotation
rho: i -> j -> k) and Aut_Q the 1152 affine maps of F_2^4 keeping T.  Subgroups of order 3 up to conjugacy:
   name   block part   Q part
   a      1            3-cycle on one 2-flat only (4 fixed points on Q)
   d      1            3-cycles on both flats (fixed-point free; the sigma of the slice; decided in fable_border15.py)
   r      rho          1
   ra     rho          as a
   rd     rho          as d
This program decides, for one of them and every admissible number NFIX of fixed rows (NFIX = 1 mod 3), whether a
border invariant under it exists: a fixed row x1 of E, the fixed rows of Z inside compat(x1), then a coupled clique
search for the orbit rows (theory/fable_clique.py).  Exhaustive when no cap is reported.
Usage: fable_border16.py a|d|r|ra|rd [nodecap]
"""
import sys, json, time, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_border2 as B2, fable_clique as FC

def main():
    which = sys.argv[1]; nodecap = int(sys.argv[2]) if len(sys.argv) > 2 else 10**12
    k = 16; T = [0, 1, 2, 3, 4, 8, 12]; t = 7; b = 28
    G = B2.Grp2(4, 1); ch = G.chars()
    # admissible block functions
    F = []
    for bits in range(1 << k):
        f = np.array([1 - 2 * ((bits >> q) & 1) for q in range(k)], dtype=np.int64)
        if all(int(ch[m] @ f) == 0 for m in range(k) if m not in T): F.append(f)
    F = np.array(F); no = len(F); assert no == 30
    key = {tuple(f): i for i, f in enumerate(F)}; neg = np.array([key[tuple(-f)] for f in F])
    Gam = np.array([[4, 2, 2, 2], [-2, 4, 2, -2], [-2, -2, 4, 2], [-2, 2, -2, 4]], dtype=np.int64)
    assert np.array_equal(Gam @ Gam.T, 28 * np.eye(4, dtype=np.int64))
    M = np.kron(Gam, np.eye(k, dtype=np.int64))
    # ---- block part: the rotation rho among the signed permutations commuting with Gamma
    cen = []
    for perm in itertools.permutations(range(4)):
        for sg in itertools.product((1, -1), repeat=4):
            P = np.zeros((4, 4), dtype=np.int64)
            for c in range(4): P[c, perm[c]] = sg[c]
            if np.array_equal(P @ Gam, Gam @ P): cen.append((perm, sg, P))
    assert len(cen) == 24, len(cen)
    rho = next((p, sg) for p, sg, P in cen if p == (0, 2, 3, 1) and np.array_equal(np.linalg.matrix_power(P, 3), np.eye(4, dtype=np.int64)))
    ident = ((0, 1, 2, 3), (1, 1, 1, 1))
    # ---- Q part: linear maps of order 3 keeping the admissible functions (elements of F_2^4 as integers, the pairing
    # of characters and points is the parity of m & q, as in fable_border14.py)
    def lin_maps():
        def rec(imgs, span):
            if len(imgs) == 4: yield list(imgs); return
            for v in range(1, k):
                if v in span: continue
                yield from rec(imgs + [v], span | {v ^ w for w in span})
        yield from rec([], {0})
    ord3 = []
    for imgs in lin_maps():
        lin = [0] * k
        for q in range(k):
            r = 0
            for i in range(4):
                if (q >> i) & 1: r ^= imgs[i]
            lin[q] = r
        if lin == list(range(k)): continue
        if [lin[lin[lin[q]]] for q in range(k)] != list(range(k)): continue
        if any(tuple(f[lin]) not in key for f in F): continue
        ord3.append(lin)
    assert len(ord3) == 8, len(ord3)
    fixpts = lambda lin: sum(1 for q in range(k) if lin[q] == q)
    qa = next(l for l in ord3 if fixpts(l) == 4); qd = next(l for l in ord3 if fixpts(l) == 1)
    blk, qm = {'a': (ident, qa), 'd': (ident, qd), 'r': (rho, list(range(k))), 'ra': (rho, qa), 'rd': (rho, qd)}[which]
    (bp, bs) = blk
    fi = np.array([key[tuple(f[qm])] for f in F])                       # f -> f o qm on block functions
    def g(R):                                                           # the symmetry on rows (tuples of 4 block indices)
        out = np.empty_like(R)
        for c in range(4):
            col = fi[R[:, bp[c]]]; out[:, c] = col if bs[c] == 1 else neg[col]
        return out
    vec = lambda R: np.concatenate([F[R[:, c]] for c in range(4)], axis=1).astype(np.float32)
    code = lambda R: ((R[:, 0] * no + R[:, 1]) * no + R[:, 2]) * no + R[:, 3]
    rows = np.array(list(itertools.product(range(no), repeat=4)), dtype=np.int64)
    MT = M.T.astype(np.float32); fk = float(4 * k)
    # g is a symmetry: K(xg, yg) = K(x, y) and <xg, yg> = <x, y> on random rows, and g^3 = 1
    rng = np.random.default_rng(1); sm = rows[rng.choice(len(rows), 3000, replace=False)]
    V0 = vec(sm); V1 = vec(g(sm))
    assert np.array_equal(V0 @ MT @ V0.T, V1 @ MT @ V1.T) and np.array_equal(V0 @ V0.T, V1 @ V1.T)
    assert np.array_equal(g(g(g(sm))), sm)
    r1 = g(rows); r2 = g(r1); c0 = code(rows); c1 = code(r1); c2 = code(r2)
    isfix = c1 == c0
    Rf = rows[isfix]; cf = np.stack([code(Rf), code(neg[Rf])], axis=1); Rf = Rf[cf[:, 0] == cf.min(axis=1)]
    nf_ = ~isfix
    allc = np.stack([c0, c1, c2, code(neg[rows]), code(neg[r1]), code(neg[r2])], axis=1)
    isrep = nf_ & (c0 == allc.min(axis=1))
    R0 = rows[isrep]
    OV = [vec(R0), vec(g(R0)), vec(g(g(R0)))]
    ok = np.einsum('ij,ij->i', OV[0], OV[1]) == 0
    R0 = R0[ok]; OV = [v[ok] for v in OV]
    FV = vec(Rf); NO, NF = len(R0), len(Rf)
    dfix = int(np.linalg.matrix_rank(FV.astype(np.float64))) if NF else 0
    print(f'symmetry {which}: block part {bp} {bs}, Q part with {fixpts(qm)} fixed points; orbit representatives {NO}, '
          f'fixed rows (mod sign) {NF}, dimension of their span {dfix}', flush=True)
    def comp_orb(w):
        co = np.ones(NO, dtype=bool)
        for j in range(3): co &= (np.abs(OV[j] @ w) == fk)
        return co
    nodes = [0]; trunc = [False]
    def fixed_cliques(idx, need, must=None):
        if need == 0: return [[]]
        V = FV[idx]; A = (V @ V.T) == 0; np.fill_diagonal(A, False)
        if must is None: out, nn, tr = FC.kcliques(A, need, maxout=500000, maxnodes=max(1, nodecap - nodes[0]))
        else:
            m = int(np.nonzero(idx == must)[0][0])
            out, nn, tr = FC.kcliques(A, need - 1, cand=[c for c in range(len(idx)) if A[m, c]], maxout=500000, maxnodes=max(1, nodecap - nodes[0]))
            out = [[m] + cl for cl in out]
        nodes[0] += nn
        if tr: trunc[0] = True
        return [[int(idx[c]) for c in cl] for cl in out]
    def orth_graph(ix):
        A = np.ones((len(ix), len(ix)), dtype=bool)
        for j in range(3): A &= ((OV[0][ix] @ OV[j][ix].T) == 0)
        np.fill_diagonal(A, False); return A
    results = {}
    for nfix in range(1, min(dfix, 28) + 1, 3):
        n3 = (b - nfix) // 3; t0 = time.time(); nodes[0] = 0; trunc[0] = False; found = None; nzc = 0
        for x1 in range(NF):
            xv = FV[x1]; wE = xv @ MT
            co = comp_orb(wE); cfx = np.abs(FV @ wE) == fk                 # Z-candidates compatible with x1
            eo_base = (OV[0] @ xv) == 0; ef_base = (FV @ xv) == 0           # E-candidates orthogonal to x1
            ifx = np.nonzero(cfx)[0]
            if len(ifx) < nfix: continue
            for fc in fixed_cliques(ifx, nfix):
                eo = eo_base.copy(); ef = ef_base.copy()
                for f in fc:
                    wZ = MT @ FV[f]; eo &= comp_orb(wZ); ef &= (np.abs(FV @ wZ) == fk)
                if int(ef.sum()) < nfix - 1: continue
                zo = co.copy()
                for f in fc: zo &= ((OV[0] @ FV[f]) == 0)
                if n3 == 0:
                    cand = [(fc, [])]
                else:
                    if int(eo.sum()) < n3 or int(zo.sum()) < n3: continue
                    iz = np.nonzero(zo)[0]; ie = np.nonzero(eo)[0]
                    AZ = orth_graph(iz); AE = orth_graph(ie)
                    CM = np.ones((len(iz), len(ie)), dtype=bool); XM = OV[0][ie] @ MT
                    for j in range(3): CM &= (np.abs(OV[j][iz] @ XM.T) == fk)
                    az = np.ones(len(iz), dtype=bool); ae = np.ones(len(ie), dtype=bool)
                    while True:
                        az2 = az & ((AZ & az[None, :]).sum(axis=1) >= n3 - 1) & ((CM & ae[None, :]).sum(axis=1) >= n3)
                        ae2 = ae & ((AE & ae[None, :]).sum(axis=1) >= n3 - 1) & ((CM & az2[:, None]).sum(axis=0) >= n3)
                        if az2.sum() == az.sum() and ae2.sum() == ae.sum(): break
                        az, ae = az2, ae2
                        if az.sum() < n3 or ae.sum() < n3: break
                    if az.sum() < n3 or ae.sum() < n3: continue
                    got, _, nn, tr, zc = FC.coupled(AZ, CM, AE, n3, n3, candz=np.nonzero(az)[0].tolist(), cande=np.nonzero(ae)[0].tolist(),
                                                    maxnodes=max(1, nodecap - nodes[0]), collect=100000)
                    nodes[0] += nn
                    if tr: trunc[0] = True
                    cand = [(fc, [int(iz[c]) for c in cl]) for cl in got]
                for Zf, Zo in cand:                                      # complete Z: now the whole of E
                    nzc += 1
                    eo2 = eo_base.copy(); ef2 = ef_base.copy()
                    for f in Zf:
                        wZ = MT @ FV[f]; eo2 &= comp_orb(wZ); ef2 &= (np.abs(FV @ wZ) == fk)
                    for o in Zo:
                        for j in range(3):
                            wZ = MT @ OV[j][o]; eo2 &= comp_orb(wZ); ef2 &= (np.abs(FV @ wZ) == fk)
                    ie2 = np.nonzero(eo2)[0]; if2 = np.nonzero(ef2)[0]
                    if len(ie2) < n3 or len(if2) < nfix - 1: continue
                    # fixed rows of E: x1 and nfix - 1 more, mutually orthogonal; then n3 orbit rows orthogonal to them
                    for fe in fixed_cliques(if2, nfix - 1):
                        e3 = eo2.copy()
                        for f in fe: e3 &= ((OV[0] @ FV[f]) == 0)
                        i3 = np.nonzero(e3)[0]
                        if len(i3) < n3: continue
                        cl, nn, tr = FC.kcliques(orth_graph(i3), n3, maxout=1, maxnodes=max(1, nodecap - nodes[0]))
                        nodes[0] += nn
                        if cl: found = ((Zf, Zo), ([x1] + fe, [int(i3[c]) for c in cl[0]])); break
                        if tr and not cl: trunc[0] = True
                    if found: break
                if found or trunc[0] or nodes[0] > nodecap: break
            if found or trunc[0] or nodes[0] > nodecap: break
        if found:
            (Zf, Zo), (Ef, Eo) = found
            def expand(Lf, Lo):
                rr = [Rf[f] for f in Lf]
                for o in Lo:
                    r = R0[o:o + 1]; rr += [r[0], g(r)[0], g(g(r))[0]]
                return np.array([np.concatenate([F[r[c]] for c in range(4)]) for r in rr], dtype=np.int64)
            Em = expand(Ef, Eo); Zm = expand(Zf, Zo); Km = (Em @ M.T @ Zm.T) // (4 * k)
            okk = bool(np.array_equal(Em @ M.T @ Zm.T, 4 * k * Km) and np.array_equal(np.abs(Km), np.ones_like(Km))
                       and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=np.int64))
                       and np.array_equal(Em @ Em.T, 4 * k * np.eye(b, dtype=np.int64)) and np.array_equal(Zm @ Zm.T, 4 * k * np.eye(b, dtype=np.int64)))
            verdict = f'FOUND (K Hadamard and E, Z orthogonal: {okk})'
            json.dump(dict(symmetry=which, nfix=nfix, E=Em.tolist(), Z=Zm.tolist(), K=Km.tolist(), ok=okk), open(f'logs/fable_fresh/b16_border_{which}_f{nfix}.json', 'w'))
        elif trunc[0] or nodes[0] > nodecap: verdict = 'UNDECIDED (a cap was hit)'
        else: verdict = 'NONE (search exhaustive)'
        results[nfix] = verdict
        print(f'order-3 symmetry {which}, NFIX {nfix} ({n3} orbits): {verdict}; complete Z examined {nzc}, nodes {nodes[0]}, {time.time()-t0:.0f}s', flush=True)
    print(f'summary {which}: ' + '; '.join(f'NFIX {a_}: {v.split(" (")[0]}' for a_, v in results.items()), flush=True)

if __name__ == '__main__':
    main()
