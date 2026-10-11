r"""Row-factored compiled D-check (report/paley_523_answer.md, section 3.4) [M: identical output to
theory/q523_dcheck_nb.py, `python3 theory/q523_dcheck2.py test`].

Given an S-solution, the D-configurations of a block are a product over its 9 pairs; the three Z3 rows
{0,3,6}, {1,4,7}, {2,5,8} split this product into three independent row products.  The row sums D3 of a configuration
depend on its row-configurations only, so the D3 groups of the block are exactly the triples of distinct row sums,
with members = product of the per-row member lists.  This replaces the enumeration and sort of all configurations
(theory/q523_dcheck_nb.py) by three enumerations of at most 216 row-configurations; everything else (signature
dedup, bitset 2+2 join on the Z3-level signature, exact expansion, full-level check of the members) is unchanged.
All work buffers are preallocated once per instance."""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_dcheck_nb import _sig_into, bitset_join4, radix_join4

@numba.njit(cache=True)
def _rows(pidx, dval, dstart, dcnt, rows, rc_opt, rc_sum, u_val, u_cnt, m_start, m_list, nrc):
    """per row: row-configurations (option indices, sums), distinct sums and their member lists (CSR)"""
    s = dval.shape[1]
    for r in range(3):
        g0 = rows[r, 0]; g1 = rows[r, 1]; g2 = rows[r, 2]
        m0 = dcnt[g0, pidx[g0]]; m1 = dcnt[g1, pidx[g1]]; m2 = dcnt[g2, pidx[g2]]
        n = 0
        for a in range(m0):
            ba = dstart[g0, pidx[g0]] + a
            for b in range(m1):
                bb = dstart[g1, pidx[g1]] + b
                for c in range(m2):
                    bc = dstart[g2, pidx[g2]] + c
                    rc_opt[r, n, 0] = a; rc_opt[r, n, 1] = b; rc_opt[r, n, 2] = c
                    for j in range(s): rc_sum[r, n, j] = dval[ba, j] + dval[bb, j] + dval[bc, j]
                    n += 1
        nrc[r] = n
        # distinct sums: sort the packed sums (each coordinate in [-64, 63]) and member lists (CSR)
        pk = np.empty(n, np.int64)
        for i in range(n):
            v = 0
            for j in range(s): v = (v << 7) | (rc_sum[r, i, j] + 64)
            pk[i] = v
        o = np.argsort(pk)
        nu = 0
        for ii in range(n):
            i = o[ii]
            if ii == 0 or pk[i] != pk[o[ii - 1]]:
                for j in range(s): u_val[r, nu, j] = rc_sum[r, i, j]
                m_start[r, nu] = ii; nu += 1
            m_list[r, ii] = i
        m_start[r, nu] = n
        u_cnt[r] = nu

@numba.njit(cache=True)
def _tcontract(x, y, T, out, o):
    """out[o:o+s] = T-contraction of x (x) y, i.e. sum_ij T[i, j, l] x_i y_j"""
    s = x.shape[0]
    for l in range(s):
        acc = 0
        for i in range(s):
            xi = x[i]
            if xi == 0: continue
            for j in range(s): acc += T[i, j, l] * xi * y[j]
        out[o + l] = acc

@numba.njit(cache=True)
def _triples(u_val, u_cnt, T, w3, maxn, k3, sig3, selfv, cr01, cr12, cr20):
    """Z3-level signature (lags 0, 1 of the row-sum sequence) of every triple of distinct row sums, from per-row self
    terms and per-row-pair cross terms: lag 0 = sum_r T(X_r X_r), lag 1 = T(X_1 X_0) + T(X_2 X_1) + T(X_0 X_2)"""
    U0 = u_cnt[0]; U1 = u_cnt[1]; U2 = u_cnt[2]; nt = U0 * U1 * U2
    if nt > maxn: return -1
    s = u_val.shape[2]
    for r in range(3):
        for u in range(u_cnt[r]): _tcontract(u_val[r, u], u_val[r, u], T, selfv[r, u], 0)
    for a in range(U0):
        for b in range(U1): _tcontract(u_val[1, b], u_val[0, a], T, cr01[a, b], 0)
    for b in range(U1):
        for c in range(U2): _tcontract(u_val[2, c], u_val[1, b], T, cr12[b, c], 0)
    for c in range(U2):
        for a in range(U0): _tcontract(u_val[0, a], u_val[2, c], T, cr20[c, a], 0)
    t = 0
    for a in range(U0):
        for b in range(U1):
            for c in range(U2):
                acc = np.uint64(0)
                for l in range(s):
                    v0 = selfv[0, a, l] + selfv[1, b, l] + selfv[2, c, l]
                    v1 = cr01[a, b, l] + cr12[b, c, l] + cr20[c, a, l]
                    sig3[t, l] = v0; sig3[t, s + l] = v1
                    acc += np.uint64(v0) * w3[l] + np.uint64(v1) * w3[s + l]
                k3[t] = acc; t += 1
    return nt

@numba.njit(cache=True)
def _dedup(k3, nt, so, ustart, uk):
    """distinct keys among k3[:nt]: so = order, ustart CSR into so, uk distinct keys; returns the count"""
    o = np.argsort(k3[:nt])
    for i in range(nt): so[i] = o[i]
    nu = 0
    for i in range(nt):
        if i == 0 or k3[so[i]] != k3[so[i - 1]]:
            ustart[nu] = i; uk[nu] = k3[so[i]]; nu += 1
    ustart[nu] = nt
    return nu

@numba.njit(cache=True)
def _member_keys(k, t, u_cnt, m_start, m_list, rc_opt, rows, pidx, dval, dstart, T, nl, wf, D, buf, cc,
                 mkey, msig, mrc, moff, nm):
    """full-level keys/signatures of the members of triple t of block k, appended at nm; returns new nm"""
    U1 = u_cnt[k, 1]; U2 = u_cnt[k, 2]
    a = t // (U1 * U2); b = (t // U2) % U1; c = t % U2
    s = dval.shape[1]; L = nl * s
    for i0 in range(m_start[k, 0, a], m_start[k, 0, a + 1]):
        x0 = m_list[k, 0, i0]
        for i1 in range(m_start[k, 1, b], m_start[k, 1, b + 1]):
            x1 = m_list[k, 1, i1]
            for i2 in range(m_start[k, 2, c], m_start[k, 2, c + 1]):
                x2 = m_list[k, 2, i2]
                xs = (x0, x1, x2)
                for r in range(3):
                    for q in range(3):
                        g = rows[r, q]; opt = rc_opt[k, r, xs[r], q]; base = dstart[g, pidx[g]] + opt
                        for j in range(s): D[g, j] = dval[base, j]
                _sig_into(D, T, nl, buf, cc)
                acc = np.uint64(0)
                for q in range(L):
                    msig[nm, q] = buf[q]; acc += np.uint64(buf[q]) * wf[q]
                mkey[nm] = acc; mrc[nm, 0] = x0; mrc[nm, 1] = x1; mrc[nm, 2] = x2; nm += 1
    moff[0] = nm
    return nm

@numba.njit(cache=True)
def _row_self(k, r, x, rows, rc_opt, pidx, dval, dstart, T, out, Dr):
    """lags 0 and 3 of the Z9 signature restricted to row r (positions r, r+3, r+6) for row-configuration x of block k:
    out[0:s] = T(sum_g D_g D_g), out[s:2s] = T(sum_g D_{g+3} D_g) (g + 3 stays in the row)"""
    s = dval.shape[1]
    for q in range(3):
        g = rows[r, q]; opt = rc_opt[k, r, x, q]; base = dstart[g, pidx[g]] + opt
        for j in range(s): Dr[q, j] = dval[base, j]
    for l in range(2 * s): out[l] = 0
    for q in range(3):
        q3 = (q + 1) % 3                                   # position g + 3 within the row (r, r+3, r+6 cyclic)
        for l in range(s):
            a0 = 0; a3 = 0
            for i in range(s):
                for j in range(s):
                    t = T[i, j, l]
                    if t != 0:
                        a0 += t * Dr[q, i] * Dr[q, j]; a3 += t * Dr[q3, i] * Dr[q, j]
            out[l] += a0; out[s + l] += a3

@numba.njit(cache=True)
def dcheck2(pidx4, dval, dstart, dcnt, rows, T, nl, w3, wf, k3T, kfT, T3, TF, maxn, maxsol,
            rc_opt, rc_sum, u_val, u_cnt, m_start, m_list, nrc, k3, sig3, so, ustart, uk, nuk,
            bs, htk, htv, mkey, msig, mrc, sols, selfv, cr01, cr12, cr20, q4, stamp_rs, rs_val, stage=9):
    """returns (n_d3_solutions, n_solutions, status); solutions written to sols[:min(n, maxsol)].
    stage < 9 stops early (profiling only)."""
    s = dval.shape[2]; G = pidx4.shape[1]
    X = np.zeros((3, s), np.int64); buf = np.zeros(nl * s, np.int64); cc = np.zeros((s, s), np.int64)
    for k in range(4):
        _rows(pidx4[k], dval[k], dstart[k], dcnt[k], rows, rc_opt[k], rc_sum[k], u_val[k], u_cnt[k], m_start[k],
              m_list[k], nrc[k])
        nt = _triples(u_val[k], u_cnt[k], T, w3, maxn, k3[k], sig3[k], selfv, cr01, cr12, cr20)
        if nt < 0: return 0, 0, -1
        nuk[k] = _dedup(k3[k], nt, so[k], ustart[k], uk[k])
    if stage == 1: return 0, 0, 0
    # 2+2 join on distinct Z3-level keys (pairing as in theory/q523_dcheck_nb.py)
    nq_ = np.array([nuk[0], nuk[1], nuk[2], nuk[3]])
    ordk = np.argsort(nq_)
    A = ordk[1]; B = ordk[2]; C = ordk[0]; Dd = ordk[3]
    if nq_[A] * nq_[B] > nq_[C] * nq_[Dd]:
        A, B, C, Dd = C, Dd, A, B
    kA = uk[A, :nuk[A]]; kB = uk[B, :nuk[B]]; kC = uk[C, :nuk[C]]; kD = uk[Dd, :nuk[Dd]]
    nL = nuk[A] * nuk[B]; cap = 1 << 16
    bb = 16
    while (1 << bb) < 32 * nL and bb < 24: bb += 1
    hq, nh, ov = bitset_join4(kA, kB, kC, kD, k3T, bs, bb, htk, htv, cap)
    if ov:
        rbits = 0
        while (nL >> rbits) > 2048 and rbits < 16: rbits += 1
        hq, nh = radix_join4(kA, kB, kC, kD, k3T, rbits, cap)
    if stage == 2: return 0, 0, 0
    # expansion to triples with the exact Z3-level check
    nq = 0; ua = np.zeros(4, np.int64)
    for i in range(hq.shape[0]):
        ua[A] = hq[i, 0]; ua[B] = hq[i, 1]; ua[C] = hq[i, 2]; ua[Dd] = hq[i, 3]
        for g0 in range(ustart[0, ua[0]], ustart[0, ua[0] + 1]):
            t0 = so[0, g0]
            for g1 in range(ustart[1, ua[1]], ustart[1, ua[1] + 1]):
                t1 = so[1, g1]
                for g2 in range(ustart[2, ua[2]], ustart[2, ua[2] + 1]):
                    t2 = so[2, g2]
                    for g3 in range(ustart[3, ua[3]], ustart[3, ua[3] + 1]):
                        t3 = so[3, g3]
                        ok = True
                        for q in range(T3.shape[0]):
                            if sig3[0, t0, q] + sig3[1, t1, q] + sig3[2, t2, q] + sig3[3, t3, q] != T3[q]:
                                ok = False; break
                        if ok:
                            if nq < cap:
                                q4[nq, 0] = t0; q4[nq, 1] = t1; q4[nq, 2] = t2; q4[nq, 3] = t3
                            nq += 1
    nd3 = nq
    if nq > cap: nq = cap
    if nq == 0 or stage == 3: return nd3, 0, 0
    # full level.  Lags 0 and 3 of the Z9 signature are sums of per-row terms (_row_self, cached per row-configuration);
    # every candidate is first checked on them, and only survivors get the full signature (exact check of all lags).
    D = np.zeros((G, s), np.int64); Dr = np.zeros((3, s), np.int64); tmp = np.zeros(2 * s, np.int64)
    ns = 0; nmax = mkey.shape[1]
    st = np.zeros(5, np.int64)
    for k in range(4):
        for r in range(3):
            for x in range(nrc[k, r]): stamp_rs[k, r, x] = 0
    T03 = np.zeros(2 * s, np.int64)
    for l in range(s): T03[l] = TF[l]; T03[s + l] = TF[3 * s + l]
    for i in range(nq):
        nm = 0
        for k in range(4):
            st[k] = nm
            t = q4[i, k]; U1 = u_cnt[k, 1]; U2 = u_cnt[k, 2]
            a = t // (U1 * U2); b = (t // U2) % U1; c = t % U2
            cnt = (m_start[k, 0, a + 1] - m_start[k, 0, a]) * (m_start[k, 1, b + 1] - m_start[k, 1, b]) * \
                  (m_start[k, 2, c + 1] - m_start[k, 2, c])
            if nm + cnt > nmax: return nd3, ns, -3
            us = (a, b, c)
            for i0 in range(m_start[k, 0, a], m_start[k, 0, a + 1]):
                x0 = m_list[k, 0, i0]
                for i1 in range(m_start[k, 1, b], m_start[k, 1, b + 1]):
                    x1 = m_list[k, 1, i1]
                    for i2 in range(m_start[k, 2, c], m_start[k, 2, c + 1]):
                        x2 = m_list[k, 2, i2]
                        xs = (x0, x1, x2)
                        for l in range(2 * s): msig[0, nm, l] = 0
                        for r in range(3):
                            x = xs[r]
                            if stamp_rs[k, r, x] == 0:
                                _row_self(k, r, x, rows, rc_opt, pidx4[k], dval[k], dstart[k], T, tmp, Dr)
                                for l in range(2 * s): rs_val[k, r, x, l] = tmp[l]
                                stamp_rs[k, r, x] = 1
                            for l in range(2 * s): msig[0, nm, l] += rs_val[k, r, x, l]
                        mrc[0, nm, 0] = x0; mrc[0, nm, 1] = x1; mrc[0, nm, 2] = x2; nm += 1
        st[4] = nm
        for ia in range(st[0], st[1]):
            for ib in range(st[1], st[2]):
                for ic in range(st[2], st[3]):
                    for idd in range(st[3], st[4]):
                        good = True
                        for l in range(2 * s):
                            if msig[0, ia, l] + msig[0, ib, l] + msig[0, ic, l] + msig[0, idd, l] != T03[l]:
                                good = False; break
                        if not good: continue
                        # survivor: full signatures of the four members, exact check of every lag
                        mem = (ia, ib, ic, idd); tot = np.zeros(TF.shape[0], np.int64)
                        for k in range(4):
                            x = mem[k]
                            for r in range(3):
                                for qq in range(3):
                                    g = rows[r, qq]; opt = rc_opt[k, r, mrc[0, x, r], qq]
                                    base = dstart[k, g, pidx4[k, g]] + opt
                                    for j in range(s): D[g, j] = dval[k, base, j]
                            _sig_into(D, T, nl, buf, cc)
                            for q in range(TF.shape[0]): tot[q] += buf[q]
                            if ns < maxsol:
                                for g in range(G):
                                    for j in range(s): sols[ns, k, g, j] = D[g, j]
                        full = True
                        for q in range(TF.shape[0]):
                            if tot[q] != TF[q]:
                                full = False; break
                        if full: ns += 1
    return nd3, ns, 0

class DCheck2:
    """drop-in replacement for theory/q523_dcheck_nb.DCheckNB (same check() interface and output)"""
    def __init__(self, inst, seed=7, maxn=1 << 18):
        from q523_dcheck_nb import DCheckNB
        self.ref = DCheckNB(inst, seed)                              # same tables and hash weights
        r = self.ref; self.inst = inst; s = inst.s; G = inst.F2
        self.maxn = maxn; R = 216
        self.rc_opt = np.zeros((4, 3, R, 3), np.int64); self.rc_sum = np.zeros((4, 3, R, s), np.int64)
        self.u_val = np.zeros((4, 3, R, s), np.int64); self.u_cnt = np.zeros((4, 3), np.int64)
        self.m_start = np.zeros((4, 3, R + 1), np.int64); self.m_list = np.zeros((4, 3, R), np.int64)
        self.nrc = np.zeros((4, 3), np.int64)
        self.k3 = np.zeros((4, maxn), np.uint64); self.sig3 = np.zeros((4, maxn, 2 * s), np.int64)
        self.so = np.zeros((4, maxn), np.int64); self.ustart = np.zeros((4, maxn + 1), np.int64)
        self.uk = np.zeros((4, maxn), np.uint64); self.nuk = np.zeros(4, np.int64)
        nmax = 1 << 17
        self.mkey = np.zeros((1, nmax), np.uint64); self.msig = np.zeros((1, nmax, inst.nl * s), np.int64)
        self.mrc = np.zeros((1, nmax, 3), np.int64)
        self.rows = np.array(inst.rows, np.int64)
        self.selfv = np.zeros((3, R, s), np.int64); self.cr01 = np.zeros((R, R, s), np.int64)
        self.cr12 = np.zeros((R, R, s), np.int64); self.cr20 = np.zeros((R, R, s), np.int64)
        self.q4 = np.zeros((1 << 16, 4), np.int64); self.sols = np.zeros((64, 4, inst.F2, s), np.int64)
        self.stamp_rs = np.zeros((4, 3, R), np.int64); self.rs_val = np.zeros((4, 3, R, 2 * s), np.int64)

    def check(self, pidx4, maxsol=64):
        r = self.ref
        p4 = np.array([np.asarray(p, np.int64) for p in pidx4])
        sols = self.sols if maxsol == 64 else np.zeros((maxsol, 4, self.inst.F2, self.inst.s), np.int64)
        nd3, ns, st = dcheck2(p4, r.dval, r.dstart, r.dcnt, self.rows, self.inst.cy.T, self.inst.nl, r.w3, r.wf, r.k3T, r.kfT,
                              r.T3, r.TF, self.maxn, maxsol, self.rc_opt, self.rc_sum, self.u_val, self.u_cnt, self.m_start,
                              self.m_list, self.nrc, self.k3, self.sig3, self.so, self.ustart, self.uk, self.nuk, r.bs,
                              r.htk, r.htv, self.mkey, self.msig, self.mrc, sols, self.selfv, self.cr01, self.cr12,
                              self.cr20, self.q4, self.stamp_rs, self.rs_val)
        if st < 0:                                                   # oversized instance: reference implementation
            return r.check(pidx4, maxsol)
        return [list(x.copy()) for x in sols[:min(ns, maxsol)]], dict(d3_solutions=int(nd3), status=int(st))

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    from q523_sdsolve import SDInst, planted_target
    from q523_dcheck_nb import DCheckNB
    from q523_lib import zcoords
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    res = {}
    for (F, q, yf) in [(18, 13, 'logs/q523/y_v235_a.json'), (18, 29, 'logs/cyclo/stage1_cpsat_523.json')]:
        Y = json.load(open(os.path.join(ROOT, yf)))['y']; inst0 = SDInst(F, q, 4, Y); rng = np.random.default_rng(3)
        agree = 0; n = 0
        for t in range(12):
            Pt, AC = planted_target(inst0, rng); inst = SDInst(F, q, 4, Y, T=Pt); Z = zcoords(AC)
            pidx4 = [np.array([next(i for i, sv in enumerate(inst.pairs[k][g]['S'].tolist()) if sv == (Z[k, g] + Z[k, g + inst.F2]).tolist())
                               for g in range(inst.F2)]) for k in range(4)]
            a, ia = DCheckNB(inst).check(pidx4); b, ib = DCheck2(inst).check(pidx4)
            ka = sorted(tuple(np.concatenate([np.ravel(x) for x in s_]).tolist()) for s_ in a)
            kb = sorted(tuple(np.concatenate([np.ravel(x) for x in s_]).tolist()) for s_ in b)
            n += 1; agree += int(ka == kb and ia['d3_solutions'] == ib['d3_solutions'] and len(kb) >= 1)
        res[f"planted_F{F}_q{q}"] = dict(agree=agree, total=n)
    if len(sys.argv) > 2:                                            # saved real S-solutions: identical results, timing
        Y = json.load(open(os.path.join(ROOT, sys.argv[3])))['y']; inst = SDInst(18, int(sys.argv[4]), 4, Y)
        P4 = np.load(sys.argv[2]); a_ = DCheckNB(inst, 1); b_ = DCheck2(inst, 1)
        a_.check(P4[0]); b_.check(P4[0]); agree = 0
        t0 = time.time(); RA = [a_.check(p) for p in P4]; t1 = time.time(); RB = [b_.check(p) for p in P4]; t2 = time.time()
        for (sa, ia), (sb, ib) in zip(RA, RB): agree += int(len(sa) == len(sb) and ia['d3_solutions'] == ib['d3_solutions'])
        res["real"] = dict(agree=agree, total=len(P4), ms_nb=round(1e3 * (t1 - t0) / len(P4), 2), ms_v2=round(1e3 * (t2 - t1) / len(P4), 2))
    print(json.dumps(res))
