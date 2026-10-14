r"""How far does the Gaussian tail of the D-level hold?  (report/paley_523_answer.md, section 3.4) [M against E].

theory/q523_dtail_check.py compares exact hit rates with the D-level Gaussian on lag subsets up to 17 bits, by
enumerating every D-candidate.  Here the production D-check (theory/q523_dcheck2.py, same code) also counts, among the
candidates that pass its exact lag-{0, 3} prefilter, those that match the target on lags {0, 1, 3}, {0, 2, 3} and
{0, 1, 2, 3} (about 24, 24 and 37 bits).  Run on a few hundred thousand production S-solutions this tests the Gaussian
at about 24 bits, 7 more than direct enumeration can.
usage: python3 theory/q523_dtail_prod.py F q s y.json P4.npy [--n 300000] [--off 0]"""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_dcheck_nb import _sig_into, bitset_join4, radix_join4
from q523_dcheck2 import _rows, _triples, _dedup, _row_self, _member_keys, _tcontract, DCheck2

@numba.njit(cache=True)
def dcount2(pidx4, dval, dstart, dcnt, rows, T, nl, w3, wf, k3T, kfT, T3, TF, maxn, maxsol,
            rc_opt, rc_sum, u_val, u_cnt, m_start, m_list, nrc, k3, sig3, so, ustart, uk, nuk,
            bs, htk, htv, mkey, msig, mrc, sols, selfv, cr01, cr12, cr20, q4, stamp_rs, rs_val, pcnt, stage=9):
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
                        # counts: pcnt[0] lag-{0,3} survivors; pcnt[1] {0,1,3}; pcnt[2] {0,2,3}; pcnt[3] {0,1,2,3}; pcnt[4] all lags
                        pcnt[0] += 1
                        m1 = True; m2 = True; m4 = True
                        for l in range(s):
                            if tot[s + l] != TF[s + l]: m1 = False
                            if tot[2 * s + l] != TF[2 * s + l]: m2 = False
                        for q in range(TF.shape[0]):
                            if tot[q] != TF[q]: m4 = False
                        if m1: pcnt[1] += 1
                        if m2: pcnt[2] += 1
                        if m1 and m2: pcnt[3] += 1
                        if m4:
                            pcnt[4] += 1; ns += 1
    return nd3, ns, 0


if __name__ == "__main__":
    from q523_sdsolve import SDInst
    a = sys.argv[1:]; opt = {"--n": "300000", "--off": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y); P4 = np.load(pos[4], mmap_mode='r')
    off = int(opt["--off"]); n = min(int(opt["--n"]), len(P4) - off)
    d = DCheck2(inst, 1); r = d.ref; cnt = np.zeros(5, np.int64); t0 = time.time(); nd3 = 0; ncalls = 0; fallback = 0
    for j in range(off, off + n):
        p4 = np.asarray(P4[j], np.int64)
        res = dcount2(p4, r.dval, r.dstart, r.dcnt, d.rows, inst.cy.T, inst.nl, r.w3, r.wf, r.k3T, r.kfT, r.T3, r.TF, d.maxn, 64,
                      d.rc_opt, d.rc_sum, d.u_val, d.u_cnt, d.m_start, d.m_list, d.nrc, d.k3, d.sig3, d.so, d.ustart, d.uk,
                      d.nuk, r.bs, r.htk, r.htv, d.mkey, d.msig, d.mrc, d.sols, d.selfv, d.cr01, d.cr12, d.cr20, d.q4,
                      d.stamp_rs, d.rs_val, cnt)
        nd3 += res[0]; ncalls += 1; fallback += int(res[2] < 0)
    names = ["lags 0,3 (prefilter survivors)", "lags 0,1,3", "lags 0,2,3", "lags 0,1,2,3", "all lags"]
    print(json.dumps(dict(s_solutions=n, d3_solutions=int(nd3), fallback=fallback, secs=round(time.time() - t0, 1),
                          counts={k: int(v) for k, v in zip(names, cnt)})), flush=True)
