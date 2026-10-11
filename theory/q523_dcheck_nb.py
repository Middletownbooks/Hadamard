r"""Fully compiled D-check (numba): one call per S-solution, no Python in the inner loops
(report/paley_523_answer.md, section 3).  Same algorithm and output as theory/q523_dcheck.DChecker:
D-options per block given S -> configurations -> D3 groups (row sums of D) -> 2+2 join of the D3 values on the
Z3-compressed D-level signature -> per D3-solution 2+2 join of the group members on the full D-level signature;
solutions re-checked exactly.  `python3 theory/q523_dcheck_nb.py test` compares with DChecker on planted cases."""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_sdsolve import seq_sig

@numba.njit(cache=True)
def _sig_into(X, T, nl, out, cc):
    """PAF signature of one sequence X (n, s) at lags 0..nl-1 into out (nl * s,) (no allocation)"""
    n, s = X.shape
    for h in range(nl):
        for i in range(s):
            for j in range(s): cc[i, j] = 0
        for g in range(n):
            g2 = (g + h) % n
            for i in range(s):
                xi = X[g2, i]
                if xi == 0: continue
                for j in range(s): cc[i, j] += xi * X[g, j]
        for l in range(s):
            acc = 0
            for i in range(s):
                for j in range(s): acc += T[i, j, l] * cc[i, j]
            out[h * s + l] = acc

@numba.njit(cache=True)
def _block(pidx, dval, dstart, dcnt, rows, T, w3, wd, maxn):
    """one block: configs (n, G, s), D3 groups (order, gstart), D3 keys and signatures"""
    G = pidx.shape[0]; s = dval.shape[1]
    m = np.zeros(G, np.int64); n = 1
    for g in range(G):
        m[g] = dcnt[g, pidx[g]]; n *= m[g]
    status = 0
    if n > maxn:
        status = -1; n = 0
    cfg = np.zeros((n, G, s), np.int64)
    for l in range(n):
        rem = l
        for g in range(G - 1, -1, -1):
            t = rem % m[g]; rem //= m[g]
            base = dstart[g, pidx[g]] + t
            for j in range(s): cfg[l, g, j] = dval[base, j]
    R = rows.shape[0]; RL = rows.shape[1]
    D3 = np.zeros((n, R, s), np.int64)
    kd = np.zeros(n, np.uint64)
    for l in range(n):
        for r in range(R):
            for t in range(RL):
                g = rows[r, t]
                for j in range(s): D3[l, r, j] += cfg[l, g, j]
        acc = np.uint64(0)
        for r in range(R):
            for j in range(s): acc += np.uint64(D3[l, r, j]) * wd[r * s + j]
        kd[l] = acc
    order = np.argsort(kd)
    gstart = np.zeros(n + 1, np.int64); gfirst = np.zeros(n, np.int64); ng = 0
    for i in range(n):
        if i == 0 or kd[order[i]] != kd[order[i - 1]]:
            gstart[ng] = i; gfirst[ng] = order[i]; ng += 1
    gstart[ng] = n
    k3 = np.zeros(ng, np.uint64); sig3 = np.zeros((ng, 2 * s), np.int64)
    buf = np.zeros(2 * s, np.int64); cc = np.zeros((s, s), np.int64)
    for q in range(ng):
        _sig_into(D3[gfirst[q]], T, 2, buf, cc)
        acc = np.uint64(0)
        for t in range(2 * s):
            sig3[q, t] = buf[t]; acc += np.uint64(buf[t]) * w3[t]
        k3[q] = acc
    # distinct signatures (D3 and -D3 share one: the D-options of a block are closed under global negation)
    so = np.argsort(k3)
    ustart = np.zeros(ng + 1, np.int64); ufirst = np.zeros(ng, np.int64); nu = 0
    for i in range(ng):
        if i == 0 or k3[so[i]] != k3[so[i - 1]]:
            ustart[nu] = i; ufirst[nu] = so[i]; nu += 1
    ustart[nu] = ng
    uk = np.zeros(nu, np.uint64)
    for u in range(nu): uk[u] = k3[ufirst[u]]
    return status, cfg, order, gstart[:ng + 1], k3, sig3, uk, so, ustart[:nu + 1]

@numba.njit(cache=True)
def _fullsig(cfg, order, gstart, need, T, nl, wf):
    """full D-level signatures and keys for members of the groups flagged in need (others stay zero)"""
    n = cfg.shape[0]; s = cfg.shape[2]
    sigf = np.zeros((n, nl * s), np.int64); kf = np.zeros(n, np.uint64)
    buf = np.zeros(nl * s, np.int64); cc = np.zeros((s, s), np.int64)
    for q in range(need.shape[0]):
        if not need[q]: continue
        for i in range(gstart[q], gstart[q + 1]):
            x = order[i]
            _sig_into(cfg[x], T, nl, buf, cc)
            acc = np.uint64(0)
            for t in range(nl * s):
                sigf[x, t] = buf[t]; acc += np.uint64(buf[t]) * wf[t]
            kf[x] = acc
    return sigf, kf

@numba.njit(cache=True)
def radix_join4(kA, kB, kC, kD, kT, rbits, cap):
    """all (a, b, c, d) with kA[a] + kB[b] + kC[c] + kD[d] == kT (uint64).  Left keys and right needs are radix-
    partitioned by their top rbits bits, and each partition is joined with a small cache-resident hash table."""
    nA = kA.shape[0]; nB = kB.shape[0]; nC = kC.shape[0]; nD = kD.shape[0]
    if rbits < 1: rbits = 1
    P = 1 << rbits; sh = np.uint64(64 - rbits)
    cL = np.zeros(P + 1, np.int64); cR = np.zeros(P + 1, np.int64)
    for a in range(nA):
        for b in range(nB): cL[int((kA[a] + kB[b]) >> sh) + 1] += 1
    for c in range(nC):
        for d in range(nD): cR[int((kT - kC[c] - kD[d]) >> sh) + 1] += 1
    for i in range(P):
        cL[i + 1] += cL[i]; cR[i + 1] += cR[i]
    Lk = np.empty(nA * nB, np.uint64); Lv = np.empty(nA * nB, np.int64)
    Rk = np.empty(nC * nD, np.uint64); Rv = np.empty(nC * nD, np.int64)
    pos = cL[:P].copy()
    for a in range(nA):
        for b in range(nB):
            key = kA[a] + kB[b]; p = int(key >> sh); i = pos[p]; Lk[i] = key; Lv[i] = a * nB + b; pos[p] = i + 1
    pos = cR[:P].copy()
    for c in range(nC):
        for d in range(nD):
            key = kT - kC[c] - kD[d]; p = int(key >> sh); i = pos[p]; Rk[i] = key; Rv[i] = c * nD + d; pos[p] = i + 1
    mx = 1
    for p in range(P):
        if cL[p + 1] - cL[p] > mx: mx = cL[p + 1] - cL[p]
    size = 1
    while size < 2 * mx: size <<= 1
    mask = size - 1; umask = np.uint64(mask)
    tk = np.zeros(size, np.uint64); tv = np.zeros(size, np.int64); tst = np.zeros(size, np.int64)
    out = np.zeros((min(cap, 4096), 4), np.int64); n = 0
    for p in range(P):
        st = p + 1
        for i in range(cL[p], cL[p + 1]):
            slot = np.int64(Lk[i] & umask)
            while tst[slot] == st: slot = (slot + 1) & mask
            tst[slot] = st; tk[slot] = Lk[i]; tv[slot] = Lv[i]
        for i in range(cR[p], cR[p + 1]):
            slot = np.int64(Rk[i] & umask)
            while tst[slot] == st:
                if tk[slot] == Rk[i]:
                    if n < out.shape[0]:
                        v = tv[slot]; w = Rv[i]
                        out[n, 0] = v // nB; out[n, 1] = v % nB; out[n, 2] = w // nD; out[n, 3] = w % nD
                    n += 1
                slot = (slot + 1) & mask
    return out[:min(n, out.shape[0])], n

@numba.njit(cache=True)
def dcheck_nb(pidx4, dval, dstart, dcnt, rows, T, nl, w3, wf, wd, k3T, kfT, T3, TF, maxn, maxsol,
              tab, tst, stamp):
    """returns (n_d3_solutions, n_solutions, solutions (maxsol, 4, G, s), status).  tab (size, 2) uint64 interleaved
    key/value table and tst stamps are reused across calls (stamp distinct per call); the table actually used is the
    smallest power of two >= 2 * (left pairs), so it stays cache-resident.  Blocks are paired largest-with-smallest by
    number of D3 values (the condition is a sum over blocks)."""
    G = pidx4.shape[1]; s = dval.shape[2]
    sols = np.zeros((maxsol, 4, G, s), np.int64); ns = 0
    bl = [_block(pidx4[k], dval[k], dstart[k], dcnt[k], rows, T, w3, wd, maxn) for k in range(4)]
    for k in range(4):
        if bl[k][0] < 0: return 0, 0, sols, -1
    nq_ = np.array([bl[k][6].shape[0] for k in range(4)])
    ordk = np.argsort(nq_)                       # ascending
    # pairing: (largest, smallest) on the right (probe side), the middle two on the left (table side) or vice versa
    A = ordk[1]; B = ordk[2]; C = ordk[0]; D = ordk[3]
    if nq_[A] * nq_[B] > nq_[C] * nq_[D]:
        A, B, C, D = C, D, A, B
    kA = bl[A][6]; kB = bl[B][6]; kC = bl[C][6]; kD = bl[D][6]
    nA = kA.shape[0]; nB = kB.shape[0]
    nL = nA * nB; rbits = 0
    while (nL >> rbits) > 2048 and rbits < 16: rbits += 1
    cap = 1 << 16
    hq, nh = radix_join4(kA, kB, kC, kD, k3T, rbits, cap)
    q4 = np.zeros((cap, 4), np.int64); nq = 0
    ua = np.zeros(4, np.int64)
    for i in range(hq.shape[0]):
        ua[A] = hq[i, 0]; ua[B] = hq[i, 1]; ua[C] = hq[i, 2]; ua[D] = hq[i, 3]
        # expand each distinct signature to its member D3 groups
        for g0 in range(bl[0][8][ua[0]], bl[0][8][ua[0] + 1]):
            a = bl[0][7][g0]
            for g1 in range(bl[1][8][ua[1]], bl[1][8][ua[1] + 1]):
                b = bl[1][7][g1]
                for g2 in range(bl[2][8][ua[2]], bl[2][8][ua[2] + 1]):
                    c = bl[2][7][g2]
                    for g3 in range(bl[3][8][ua[3]], bl[3][8][ua[3] + 1]):
                        d = bl[3][7][g3]
                        ok = True
                        for t in range(T3.shape[0]):
                            if bl[0][5][a, t] + bl[1][5][b, t] + bl[2][5][c, t] + bl[3][5][d, t] != T3[t]:
                                ok = False; break
                        if ok and nq < cap:
                            q4[nq, 0] = a; q4[nq, 1] = b; q4[nq, 2] = c; q4[nq, 3] = d; nq += 1
    nd3 = nq
    if nq > cap: nq = cap
    if nq == 0: return nd3, 0, sols, 0
    b0 = bl[0]; b1 = bl[1]; b2 = bl[2]; b3 = bl[3]
    nd0 = np.zeros(b0[3].shape[0] - 1, np.bool_); nd1 = np.zeros(b1[3].shape[0] - 1, np.bool_)
    nd2 = np.zeros(b2[3].shape[0] - 1, np.bool_); nd3_ = np.zeros(b3[3].shape[0] - 1, np.bool_)
    for i in range(nq):
        nd0[q4[i, 0]] = True; nd1[q4[i, 1]] = True; nd2[q4[i, 2]] = True; nd3_[q4[i, 3]] = True
    f0 = _fullsig(b0[1], b0[2], b0[3], nd0, T, nl, wf); f1 = _fullsig(b1[1], b1[2], b1[3], nd1, T, nl, wf)
    f2 = _fullsig(b2[1], b2[2], b2[3], nd2, T, nl, wf); f3 = _fullsig(b3[1], b3[2], b3[3], nd3_, T, nl, wf)
    for i in range(nq):
        a = q4[i, 0]; b = q4[i, 1]; c = q4[i, 2]; d = q4[i, 3]
        for ia in range(b0[3][a], b0[3][a + 1]):
            xa = b0[2][ia]
            for ib in range(b1[3][b], b1[3][b + 1]):
                xb = b1[2][ib]; kab = f0[1][xa] + f1[1][xb]
                for ic in range(b2[3][c], b2[3][c + 1]):
                    xc = b2[2][ic]
                    for idd in range(b3[3][d], b3[3][d + 1]):
                        xd = b3[2][idd]
                        if kab + f2[1][xc] + f3[1][xd] != kfT: continue
                        good = True
                        for t in range(TF.shape[0]):
                            if f0[0][xa, t] + f1[0][xb, t] + f2[0][xc, t] + f3[0][xd, t] != TF[t]:
                                good = False; break
                        if good:
                            if ns < maxsol:
                                sols[ns, 0] = b0[1][xa]; sols[ns, 1] = b1[1][xb]
                                sols[ns, 2] = b2[1][xc]; sols[ns, 3] = b3[1][xd]
                            ns += 1
    return nd3, ns, sols[:min(ns, maxsol)], 0

class DCheckNB:
    def __init__(self, inst, seed=7, maxn=1 << 20):
        self.inst = inst; s = inst.s; F2 = inst.F2
        mx = max(len(inst.pairs[k][g]['D']) for k in range(4) for g in range(F2))
        tot = [sum(sum(len(d) for d in inst.pairs[k][g]['D']) for g in range(F2)) for k in range(4)]
        self.dval = np.zeros((4, max(tot), s), np.int64); self.dstart = np.zeros((4, F2, mx), np.int64)
        self.dcnt = np.zeros((4, F2, mx), np.int64)
        for k in range(4):
            off = 0
            for g in range(F2):
                for i, d in enumerate(inst.pairs[k][g]['D']):
                    self.dstart[k, g, i] = off; self.dcnt[k, g, i] = len(d); self.dval[k, off:off + len(d)] = d; off += len(d)
        self.rows = np.array(inst.rows, np.int64)
        rng = np.random.default_rng(seed)
        self.w3 = rng.integers(1, 2 ** 63, 2 * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.wf = rng.integers(1, 2 ** 63, inst.nl * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.wd = rng.integers(1, 2 ** 63, 3 * s, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        self.T3 = inst.tsigD3().astype(np.int64); self.TF = inst.tsigD().astype(np.int64)
        self.k3T = (self.T3.astype(np.uint64) * self.w3).sum(); self.kfT = (self.TF.astype(np.uint64) * self.wf).sum()
        self.maxn = maxn; self.tbits = 20; self._alloc()

    def _alloc(self):
        size = 1 << self.tbits
        self.tab = np.zeros((size, 2), np.uint64); self.tst = np.zeros(size, np.int64); self.stamp = 0

    def check(self, pidx4, maxsol=64):
        p4 = np.array([np.asarray(p, np.int64) for p in pidx4])
        while True:
            self.stamp += 1
            nd3, ns, sols, st = dcheck_nb(p4, self.dval, self.dstart, self.dcnt, self.rows, self.inst.cy.T, self.inst.nl,
                                           self.w3, self.wf, self.wd, self.k3T, self.kfT, self.T3, self.TF, self.maxn, maxsol,
                                           self.tab, self.tst, self.stamp)
            if st != -2: break
            self.tbits += 1; self._alloc()
        return [list(x) for x in sols], dict(d3_solutions=int(nd3), status=int(st))

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    from q523_sdsolve import SDInst, planted_target
    from q523_dcheck import DChecker
    from q523_lib import zcoords
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    res = {}
    for (F, q, yf) in [(6, 13, 'logs/q523/stage1_v79_s1.json'), (18, 29, 'logs/cyclo/stage1_cpsat_523.json')]:
        Y = json.load(open(os.path.join(ROOT, yf)))['y']; inst0 = SDInst(F, q, 4, Y); rng = np.random.default_rng(3)
        agree = 0; n = 0; t_ref = 0.0; t_nb = 0.0
        for t in range(12):
            Pt, AC = planted_target(inst0, rng); inst = SDInst(F, q, 4, Y, T=Pt); Z = zcoords(AC)
            pidx4 = [np.array([next(i for i, sv in enumerate(inst.pairs[k][g]['S'].tolist()) if sv == (Z[k, g] + Z[k, g + inst.F2]).tolist())
                               for g in range(inst.F2)]) for k in range(4)]
            t0 = time.time(); a, ia = DChecker(inst).check(pidx4); t_ref += time.time() - t0
            dn = DCheckNB(inst); t0 = time.time(); b, ib = dn.check(pidx4); t_nb += time.time() - t0
            ka = sorted(tuple(np.concatenate([np.ravel(x) for x in s_]).tolist()) for s_ in a)
            kb = sorted(tuple(np.concatenate([np.ravel(x) for x in s_]).tolist()) for s_ in b)
            n += 1; agree += int(ka == kb and ia['d3_solutions'] == ib['d3_solutions'] and len(kb) >= 1)
        res[f"F{F}_q{q}"] = dict(agree=agree, total=n, secs_ref=round(t_ref, 2), secs_nb=round(t_nb, 2))
    print(json.dumps(res))
