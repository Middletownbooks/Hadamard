r"""Planted test of the Wagner S/D pipeline at full size (report/paley_523_answer.md, section 4) [M].

A random configuration given y is planted; the stage-2 target becomes its P.  The candidate rows of every block are
forced to contain the planted row sums and the planted row configurations; the middle value is the planted left Z3
partial sum, and the planted parent pairs are joined together with other parent pairs of the same middle value (so the
join runs at its production shape).  Pass = the multi-parent join returns the planted S-solution, the compiled D-check
returns the planted D among its solutions, and the assembled configuration reproduces the planted target exactly
(2-D FFT).  The search rate itself is calibrated separately (theory/q523_calib.py, theory/q523_dfine_calib.py).
usage: python3 theory/q523_planted.py F q s y.json [--trials 3] [--logN 30]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_wagner as W
from q523_sdsolve import SDInst, planted_target, seq_sig, next_prime, assemble, verify
from q523_dcheck_nb import DCheckNB
from q523_lib import zcoords, pz, expand_x, paf2

def planted_trial(F, q, s, Y, seed, logN=30.0, K=24, nrow=48):
    rng = np.random.default_rng(seed)
    inst0 = SDInst(F, q, s, Y); Pt, AC = planted_target(inst0, rng); inst = SDInst(F, q, s, Y, T=Pt)
    Z = zcoords(AC); F2 = inst.F2; T = inst.cy.T
    pidx = [np.array([next(i for i, sv in enumerate(inst.pairs[k][g]['S'].tolist()) if sv == (Z[k, g] + Z[k, g + F2]).tolist())
                      for g in range(F2)]) for k in range(4)]
    Dpl = [np.array([(1 if g % 2 == 0 else -1) * (Z[k, g] - Z[k, g + F2]) for g in range(F2)]) for k in range(4)]
    rows = []
    for k in range(4):
        rk = W.block_candidates(inst, k, 1.0, nrow, K, rng)
        for r in range(3):
            pr = tuple(int(x) for x in pidx[k][inst.rows[r]])                       # planted row configuration
            rs = tuple(sum(inst.pairs[k][g]['S'][pidx[k][g]] for g in inst.rows[r]).tolist())
            hit = [j for j, (key, idx, lm) in enumerate(rk[r]) if key == rs]
            if hit:
                key, idx, lm = rk[r][hit[0]]
                if not any(tuple(x) == pr for x in idx.tolist()): idx = np.vstack([idx, np.array(pr)[None, :]])
                rk[r][hit[0]] = (key, idx, lm)
            else:
                rk[r].append((rs, np.array(pr)[None, :], 0.0))
        rows.append(rk)
    C = [W.cand_table(inst, rows[k]) for k in range(4)]
    # planted candidate index per block
    pc = []
    for k in range(4):
        rs = [tuple(sum(inst.pairs[k][g]['S'][pidx[k][g]] for g in inst.rows[r]).tolist()) for r in range(3)]
        j = [next(t for t, e in enumerate(rows[k][r]) if e[0] == rs[r]) for r in range(3)]
        pc.append(int(np.nonzero((C[k][0] == np.array(j)[None, :]).all(1))[0][0]))
    c = C[0][1][pc[0]] + C[1][1][pc[1]]
    w = np.random.default_rng(12345).integers(1, 2 ** 62, C[0][1].shape[1], dtype=np.int64).astype(np.uint64)
    h1, h3 = W.SigHash(C[1][1], w), W.SigHash(C[3][1], w)
    k0 = (C[0][1].astype(np.uint64) * w[None, :]).sum(1); k2 = (C[2][1].astype(np.uint64) * w[None, :]).sum(1)
    L = W.parent_pairs(C[0][1], k0, h1, c); R = W.parent_pairs(C[2][1], k2, h3, inst.tsig3() - c)
    assert (pc[0], pc[1]) in L and (pc[2], pc[3]) in R
    def trim(Pp, la, lb, must):
        Pp = sorted(Pp, key=lambda ij: -(la[ij[0]] + lb[ij[1]])); keep = [must]; tot = 2.0 ** (la[must[0]] + lb[must[1]])
        for ij in Pp:
            if ij == must: continue
            if math.log2(tot) >= logN: break
            keep.append(ij); tot += 2.0 ** (la[ij[0]] + lb[ij[1]])
        return keep, math.log2(tot)
    Lk, wL = trim(L, C[0][2], C[1][2], (pc[0], pc[1])); Rk, wR = trim(R, C[2][2], C[3][2], (pc[2], pc[3]))
    tS = inst.tsigS(); D = tS.shape[0]
    wts = rng.integers(1, 2 ** 63, D, dtype=np.int64).astype(np.uint64) | np.uint64(1); rho = rng.integers(1, 1 << 30, D, dtype=np.int64)
    tbits = 20; P = next_prime(max(2, 2 ** wL / ((1 << tbits) * 0.6)))
    need = [sorted({p[0] for p in Lk}), sorted({p[1] for p in Lk}), sorted({p[0] for p in Rk}), sorted({p[1] for p in Rk})]
    sides = []; fid = []
    for k in range(4):
        fib = []
        for ci in need[k]:
            pid, S = W.materialise(inst, k, rows[k], C[k][0][ci], 10 ** 9)
            fib.append((pid, seq_sig(S, T, inst.nl)))
        sides.append(W.Side(fib, wts, rho, P)); fid.append({x: i for i, x in enumerate(need[k])})
    lpA = np.array([fid[0][i] for i, j in Lk]); lpB = np.array([fid[1][j] for i, j in Lk])
    rpC = np.array([fid[2][i] for i, j in Rk]); rpD = np.array([fid[3][j] for i, j in Rk])
    kT = (tS.astype(np.int64).astype(np.uint64) * wts).sum(); rT = int(((tS % P) * (rho % P) % P).sum() % P)
    A, B, Cc, Dd = sides; t0 = time.time()
    hits, nh, ov, pr = W.mjoin(A.key, A.res, A.off, B.key, B.rst, B.off, Cc.key, Cc.res, Cc.off, Dd.key, Dd.rst, Dd.off,
                               lpA.astype(np.int64), lpB.astype(np.int64), rpC.astype(np.int64), rpD.astype(np.int64), P, rT, kT, tbits, 1 << 22, 0, P)
    found_S = False; found_core = False; n_s = 0
    dchk = DCheckNB(inst)
    for a, b, cc, dd in hits:
        p4 = [A.pidx_all[a].astype(np.int64), B.pidx_all[b].astype(np.int64), Cc.pidx_all[cc].astype(np.int64), Dd.pidx_all[dd].astype(np.int64)]
        S4 = [W.S_of(inst, k, p4[k])[0] for k in range(4)]
        if not (sum(seq_sig(x[None], T, inst.nl)[0] for x in S4) == tS).all(): continue
        n_s += 1
        if all((p4[k] == pidx[k]).all() for k in range(4)):
            found_S = True
            sols, info = dchk.check(p4)
            planted_D = any(all((np.asarray(D4[k]) == Dpl[k]).all() for k in range(4)) for D4 in sols)
            for D4 in sols:
                Zs, ACs = assemble(inst, S4, D4)
                Q = paf2(expand_x(ACs, inst.cy))
                R = np.stack([Q[:, inst.cy.reps[l]] for l in range(s)], 1)
                ok = bool((pz(zcoords(ACs), inst.cy.T).sum(0) == Pt).all() and (R - Q[:, [0]] == Pt).all())
                found_core |= ok
    return dict(seed=seed, side_log2=(round(wL, 2), round(wR, 2)), parents=(len(Lk), len(Rk)), s_solutions=n_s,
                planted_S_found=found_S, planted_D_found=bool(found_S and planted_D), assembled_ok=found_core,
                secs=round(time.time() - t0, 1))

if __name__ == "__main__":
    a = sys.argv[1:]; tr = 3; logN = 30.0
    if '--trials' in a: i = a.index('--trials'); tr = int(a[i + 1]); a = a[:i] + a[i + 2:]
    if '--logN' in a: i = a.index('--logN'); logN = float(a[i + 1]); a = a[:i] + a[i + 2:]
    F, q, s = map(int, a[:3]); Y = json.load(open(a[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    for t in range(tr): print(json.dumps(planted_trial(F, q, s, Y, 1000 + t, logN)), flush=True)
