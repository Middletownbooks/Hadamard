r"""Gate 1, second method, of report/prompt_blockswap_863.md: an exhaustive decision of stage 1 of the
block-exchanging ansatz (A, B, C, tau C) at v = 863 that shares no code with theory/blockswap_stage1_exact.py.

Credit.  Ansatz and cost formula: the OpenAI model (report/openai_deficit_answer.md, codex/order-3452-compression,
5b26dae).  863 set-up and certificate: OpenAI Codex.  Bordered construction: Kumar's form of Theorem 7 of EKLO.

Differences from the primary method (deliberate):
  * formulation: the 20 H-level equations are replaced by their Q-character form [P].  For psi in Q^ write
    x_psi(u) = sum_q psi(q) x(q, u) (a sequence on Z5) and P(s) = (sum s_u^2, sum s_u s_{u+1}, sum s_u s_{u+2}).
    Then sum_a psi(a) PAF_H(a, b) = PAF_{Z5}(x_psi)(b), and under the ansatz (tau C)_psi(u) = C_{psi o tau}(-u),
    A_01 = A_10 reversed.  So the H-level is exactly the 9 equations
        P(A_00) + P(B_00) + 2 P(C_00) = (3452, 0, 0)
        P(A_11) + P(B_11) + 2 P(C_11) = (2764, -688, -688)
        P(A_01) + P(B_01) + P(C_01) + P(C_10) = (2764, -688, -688);
  * keys: exact bit-packed integers (13 bits per value, offset 4096; values 0..3 in k1, 4..7 in k2, value 8 compared
    directly), no hashing; the C side is lexsorted and searched by binary search;
  * enumeration: A and B by a recursive per-orbit search with exact reachability pruning of the remaining row sums
    and energy; C rows by a depth-first search with min/max pruning, C itself by a compiled loop over row pairs;
  * every solution is re-checked on the explicit length-20 arrays by the 20-lag definition (in this file).
usage: python3 theory/blockswap_stage1_check2.py run FIRST LAST OUT.jsonl
       python3 theory/blockswap_stage1_check2.py planted IN_PLANTED.jsonl OUT.jsonl   (targets of the primary planted log)
"""
import sys, os, json, time
import numpy as np
import numba as nb

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
LET = [1, -1, 13, -13, 15, -15, 27, -27, 29, -29, 41, -41, 43, -43]
LETS = sorted(LET)
CHI = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]   # psi_j(q), q = 2 q1 + q2, j same coding
# row of j = 1 (psi = 01): (-1)^(q2); j = 2 (psi = 10): (-1)^(q1); j = 3 (psi = 11): (-1)^(q1 + q2)
OFF = 4096


def cell(q, u):
    return 5 * q + u


def tau_cell(q, u):
    q1, q2 = q >> 1, q & 1
    return 5 * (2 * q2 + q1) + (-u) % 5


# ------------------------------------------------------------------ the 9-value signature
def zpaf(s):
    return (sum(s[u] * s[u] for u in range(5)), sum(s[u] * s[(u + 1) % 5] for u in range(5)),
            sum(s[u] * s[(u + 2) % 5] for u in range(5)))


def twist(y, j):
    return [sum(CHI[j][q] * y[cell(q, u)] for q in range(4)) for u in range(5)]


def sig_fixed(y):
    return zpaf(twist(y, 0)) + zpaf(twist(y, 3)) + zpaf(twist(y, 1))


def sig_free(y):
    p0, p3, p1, p2 = zpaf(twist(y, 0)), zpaf(twist(y, 3)), zpaf(twist(y, 1)), zpaf(twist(y, 2))
    return tuple(2 * v for v in p0) + tuple(2 * v for v in p3) + tuple(a + b for a, b in zip(p1, p2))


TARGET9 = (3452, 0, 0, 2764, -688, -688, 2764, -688, -688)


def pack(vals):
    k1 = 0
    for v in vals[0:4]:
        k1 = (k1 << 13) | (v + OFF)
    k2 = 0
    for v in vals[4:8]:
        k2 = (k2 << 13) | (v + OFF)
    return k1, k2, vals[8]


# ------------------------------------------------------------------ A, B: recursive per-orbit enumeration
def orbits_fixed():
    """tau-orbits of the 20 cells, each with its row q (or rows) and multiplicity"""
    seen, orbs = set(), []
    for q in range(4):
        for u in range(5):
            c = cell(q, u)
            if c in seen:
                continue
            o = sorted({c, tau_cell(q, u)})
            seen.update(o)
            orbs.append(o)
    return orbs


ORB = orbits_fixed()


def enum_fixed(z, cap):
    """every tau-invariant y with row sums z and excess sum_h (y_h^2 - 1) <= cap"""
    # the row of an orbit: orbits are inside one row (q = 00 or 11) or straddle rows 01/10 with one cell each
    rows_of = [[c // 5 for c in o] for o in ORB]
    # remaining-sum bounds per row: each cell of row q not yet assigned contributes in [-43, 43] and is odd
    out = []
    y = [0] * 20
    n = len(ORB)
    # pre-compute, for suffixes, the number of remaining cells per row
    rem = [[0] * 4 for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        rem[i] = list(rem[i + 1])
        for q in rows_of[i]:
            rem[i][q] += 1
    sums = [0] * 4

    def rec(i, e):
        if i == n:
            if sums == list(z):
                out.append(list(y))
            return
        for a in LET:
            de = len(ORB[i]) * (a * a - 1)
            if e + de > cap:
                continue
            for c in ORB[i]:
                y[c] = a
                sums[c // 5] += a
            ok = True
            for q in range(4):
                r = z[q] - sums[q]
                k = rem[i + 1][q]
                if abs(r) > 43 * k or (r - k) % 2:
                    ok = False
                    break
                if k == 0 and r != 0:
                    ok = False
                    break
            if ok:
                rec(i + 1, e + de)
            for c in ORB[i]:
                sums[c // 5] -= a
                y[c] = 0
    rec(0, 0)
    return out


# ------------------------------------------------------------------ C rows: depth-first with pruning
def rows_dfs(s, cap):
    """all 5-tuples over the alphabet with sum s and sum (y^2 - 1) <= cap"""
    out = []
    cur = []

    def rec(k, r, e):
        if k == 5:
            if r == 0:
                out.append(tuple(cur))
            return
        left = 5 - k - 1
        for a in LETS:
            ne = e + a * a - 1
            if ne > cap:
                continue
            rr = r - a
            if abs(rr) > 43 * left or (left == 0 and rr != 0) or (rr - left) % 2:
                continue
            cur.append(a)
            rec(k + 1, rr, ne)
            cur.pop()
    rec(0, s, 0)
    return out


def row_min_excess(s):
    r = rows_dfs(s, 4 * 1848 + 1)
    return min(sum(a * a - 1 for a in t) for t in r) if r else None


def fixed_min(z):
    # q = 00 / 11 rows: f0 + 2 p1 + 2 p2; q = 01 row (mirrored, counted twice)
    best = 0
    for q in (0, 3):
        m = None
        for p1 in LET:
            for p2 in LET:
                f0 = z[q] - 2 * p1 - 2 * p2
                if f0 in LET:
                    e = f0 * f0 - 1 + 2 * (p1 * p1 - 1) + 2 * (p2 * p2 - 1)
                    m = e if m is None else min(m, e)
        if m is None:
            return None
        best += m
    m01 = row_min_excess(z[1])
    if m01 is None or z[1] != z[2]:
        return None
    return best + 2 * m01


def free_min(z):
    t = 0
    for s in z:
        m = row_min_excess(s)
        if m is None:
            return None
        t += m
    return 2 * t


# ------------------------------------------------------------------ C signatures by compiled loop over row pairs
@nb.njit(cache=True)
def _zp(s0, s1, s2, s3, s4, lag):
    s = (s0, s1, s2, s3, s4)
    t = 0
    for u in range(5):
        t += s[u] * s[(u + lag) % 5]
    return t


@nb.njit(cache=True)
def _c_side(R0, R1, R2, R3, E0, E1, E2, E3, cap, K1, K2, V8, IX, fill):
    """C = rows (R0, R1, R2, R3) for q = 00, 01, 10, 11; signature 2P(C00), 2P(C11), P(C01) + P(C10)"""
    n = 0
    tw = np.zeros((4, 5), np.int64)
    for a in range(R0.shape[0]):
        ea = E0[a]
        if ea > cap:
            continue
        for b in range(R1.shape[0]):
            eb = ea + E1[b]
            if eb > cap:
                continue
            for c in range(R2.shape[0]):
                ec = eb + E2[c]
                if ec > cap:
                    continue
                for d in range(R3.shape[0]):
                    if ec + E3[d] > cap:
                        continue
                    if fill:
                        for u in range(5):
                            x0 = R0[a, u]; x1 = R1[b, u]; x2 = R2[c, u]; x3 = R3[d, u]
                            tw[0, u] = x0 + x1 + x2 + x3
                            tw[1, u] = x0 - x1 + x2 - x3     # psi = 01: (-1)^q2
                            tw[2, u] = x0 + x1 - x2 - x3     # psi = 10: (-1)^q1
                            tw[3, u] = x0 - x1 - x2 + x3     # psi = 11
                        v = np.zeros(9, np.int64)
                        for lag in range(3):
                            p0 = 0; p1 = 0; p2 = 0; p3 = 0
                            for u in range(5):
                                p0 += tw[0, u] * tw[0, (u + lag) % 5]
                                p1 += tw[1, u] * tw[1, (u + lag) % 5]
                                p2 += tw[2, u] * tw[2, (u + lag) % 5]
                                p3 += tw[3, u] * tw[3, (u + lag) % 5]
                            v[lag] = 2 * p0
                            v[3 + lag] = 2 * p3
                            v[6 + lag] = p1 + p2
                        inr = True
                        for i in range(9):
                            if v[i] < -4096 or v[i] >= 4096:
                                inr = False
                        k1 = np.int64(0)
                        for i in range(4):
                            k1 = (k1 << 13) | (v[i] + 4096)
                        k2 = np.int64(0)
                        for i in range(4, 8):
                            k2 = (k2 << 13) | (v[i] + 4096)
                        if not inr:
                            # [P] no solution has a C-signature value outside +-3452 (lag-0 values are >= 0 and sum
                            # to the target; |P(b)| <= P(0) by Cauchy-Schwarz), so such C never match: sentinel -1
                            k1 = np.int64(-1); k2 = np.int64(-1)
                        K1[n] = k1; K2[n] = k2; V8[n] = v[8]
                        IX[n, 0] = a; IX[n, 1] = b; IX[n, 2] = c; IX[n, 3] = d
                    n += 1
    return n


@nb.njit(cache=True)
def _probe(SA, SB, T, K1s, K2s, V8s, out):
    """for every (a, b): need = T - SA[a] - SB[b] (9 values); find all C with equal (k1, k2, v8)"""
    n = 0
    N = K1s.shape[0]
    for a in range(SA.shape[0]):
        for b in range(SB.shape[0]):
            ok = True
            k1 = np.int64(0)
            for i in range(4):
                w = T[i] - SA[a, i] - SB[b, i]
                if w < -4096 or w >= 4096:
                    ok = False
                    break
                k1 = (k1 << 13) | (w + 4096)
            if not ok:
                continue
            k2 = np.int64(0)
            for i in range(4, 8):
                w = T[i] - SA[a, i] - SB[b, i]
                if w < -4096 or w >= 4096:
                    ok = False
                    break
                k2 = (k2 << 13) | (w + 4096)
            if not ok:
                continue
            v8 = T[8] - SA[a, 8] - SB[b, 8]
            lo = 0; hi = N
            while lo < hi:
                mid = (lo + hi) >> 1
                if K1s[mid] < k1 or (K1s[mid] == k1 and K2s[mid] < k2):
                    lo = mid + 1
                else:
                    hi = mid
            p = lo
            while p < N and K1s[p] == k1 and K2s[p] == k2:
                if V8s[p] == v8:
                    if n < out.shape[0]:
                        out[n, 0] = a; out[n, 1] = b; out[n, 2] = p
                    n += 1
                p += 1
    return n


def check20(ys, tgt20):
    """the 20-lag definition on explicit arrays (A, B, C, tau C)"""
    for a in range(4):
        for b in range(5):
            t = 0
            for y in ys:
                for q in range(4):
                    for u in range(5):
                        t += y[cell(q, u)] * y[cell(q ^ a, (u + b) % 5)]
            if t != tgt20[5 * a + b]:
                return False
    return True


def tgt9_from20(t20):
    """sum_a psi(a) T(a, b) for psi = 00, 11, 01 and b = 0, 1, 2"""
    out = []
    for j in (0, 3, 1):
        for b in range(3):
            out.append(sum(CHI[j][a] * t20[5 * a + b] for a in range(4)))
    return tuple(out)


TGT20 = [2936 if g == 0 else (-516 if g < 5 else 172) for g in range(20)]


def decide(zA, zB, zC, tgt20=None, etot=2856):
    tgt20 = TGT20 if tgt20 is None else list(tgt20)
    T9 = np.array(tgt9_from20(tgt20), np.int64)
    t0 = time.time()
    mA, mB, mC = fixed_min(zA), fixed_min(zB), free_min(zC)
    r = dict(z=[list(zA), list(zB), list(zC)], mins=[mA, mB, mC])
    if None in (mA, mB, mC) or mA + mB + mC > etot:
        r.update(status='EMPTY (Q sums or energy)', sizes=None, nsol=0, solutions=[])
        return r
    LA = enum_fixed(zA, etot - mB - mC)
    LB = enum_fixed(zB, etot - mA - mC)
    capC = etot - mA - mB
    rows = [rows_dfs(s, capC // 2) for s in zC]
    R = [np.array(x, np.int64).reshape(-1, 5) for x in rows]
    E = [2 * (Rq * Rq - 1).sum(1) for Rq in R]
    n = _c_side(*R, *E, capC, np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.int64),
                np.zeros((0, 4), np.int64), False)
    K1 = np.zeros(n, np.int64); K2 = np.zeros(n, np.int64); V8 = np.zeros(n, np.int64); IX = np.zeros((n, 4), np.int64)
    _c_side(*R, *E, capC, K1, K2, V8, IX, True)
    o = np.lexsort((K2, K1))
    K1, K2, V8, IX = K1[o], K2[o], V8[o], IX[o]
    SA = np.array([sig_fixed(y) for y in LA], np.int64).reshape(-1, 9)
    SB = np.array([sig_fixed(y) for y in LB], np.int64).reshape(-1, 9)
    t1 = time.time()
    out = np.zeros((10 ** 6, 3), np.int64)
    nh = _probe(SA, SB, T9, K1, K2, V8, out)
    sols = []
    for a, b, p in out[:min(nh, len(out))]:
        C = [0] * 20
        for q in range(4):
            for u in range(5):
                C[cell(q, u)] = int(R[q][IX[p, q], u])
        tC = [C[tau_cell(c // 5, c % 5)] for c in range(20)]
        ok = check20([LA[a], LB[b], C, tC], tgt20) and sig_free(C) == tuple(int(x) for x in (T9 - SA[a] - SB[b]))
        sols.append(dict(A=LA[a], B=LB[b], C=C, check20=ok))
    r.update(status='SOLUTIONS' if sols else 'EMPTY (exhaustive join)', sizes=[len(LA), len(LB), int(n)],
             nsol=len(sols), nverified=sum(s['check20'] for s in sols), secs_lists=round(t1 - t0, 1),
             secs=round(time.time() - t0, 1), solutions=sols)
    return r


def survivors():
    d = json.load(open(os.path.join(ROOT, 'logs', 'blockswap', 'gate0_m5.json')))
    return [t['z'] for t in d['tuples'] if t['survives']]


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        first, last, outp = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        Z = survivors()
        for i in range(first, min(last, len(Z))):
            zA, zB, zC = [tuple(z) for z in Z[i][:3]]
            r = decide(zA, zB, zC)
            r['index'] = i
            s = json.dumps(r)
            print(s[:400], flush=True)
            with open(outp, 'a') as f:
                f.write(s + '\n')
    elif sys.argv[1] == 'planted':
        for line in open(sys.argv[2]):
            p = json.loads(line)
            t20 = p['target']
            zA, zB, zC = [tuple(z) for z in p['z'][:3]]
            r = decide(zA, zB, zC, tgt20=t20, etot=t20[0] - 80)
            pl = p['plant']
            found = any(s['A'] == pl['A'] and s['B'] == pl['B'] and s['C'] == pl['C'] for s in r['solutions'])
            rec = dict(trial=p['trial'], seed=p['seed'], sizes=r['sizes'], nsol=r['nsol'], nverified=r.get('nverified'),
                       primary_hits=p['hits'], planted_recovered=found, secs=r.get('secs'))
            s = json.dumps(rec)
            print(s, flush=True)
            with open(sys.argv[3], 'a') as f:
                f.write(s + '\n')
