r"""Gate 1 (primary method) of report/prompt_blockswap_863.md: exact stage 1 (the H-level, H = Z2^2 x Z5) of the
block-exchanging ansatz (x0, x1, x2, x3) = (A, B, C, tau C), tau A = A, tau B = B, at v = 863.

Credit.  Ansatz and cost formula: the OpenAI model (report/openai_deficit_answer.md, codex/order-3452-compression,
5b26dae).  863 set-up and certificate: OpenAI Codex.  Bordered construction: Kumar's form of Theorem 7 of EKLO.
The method follows theory/h863_s14_stage1_exact.py (exhaustive energy-capped lists, exact join), but this file is
new code and shares no code with that script or with theory/openai_deficit_*.py.

For a Q tuple (zA, zB, zC) (cell sums z(q), q = 2 q1 + q2 in 00, 01, 10, 11):
  * A: every tau-invariant y on the 20 cells h = (q, u) (11 orbit values in the alphabet +-{1,13,15,27,29,41,43})
    with Q row sums zA and excess energy sum_h (y^2 - 1) <= 2856 - min(B) - min(C); same for B;
  * C: every y on the 20 cells with row sums zC and 2 * excess <= 2856 - min(A) - min(B) (C and tau C both count);
  * join A + B against C on the exact 20-lag H-level vector, where C contributes PAF_C + PAF_C o tau (tau C is not
    an independent list).  Keys are linear 64-bit hashes of the integer vectors (equal vectors have equal keys, so
    nothing is lost); every key match is re-checked on the integer vectors inside the join, and every hit is
    re-checked again in pure Python on the four explicit length-20 arrays (A, B, C, tau C).
Every solution of the ansatz has, after block negations and A <-> B, a Q tuple among the 48 of Gate 0; the 8 that
fail the joint energy gate cannot hold a solution.  So the 40 survivors decide the ansatz at the H-level.

usage: python3 theory/blockswap_stage1_exact.py run FIRST LAST OUT.jsonl     (survivor indices, Gate 0 order)
       python3 theory/blockswap_stage1_exact.py planted N SEED OUT.jsonl
       python3 theory/blockswap_stage1_exact.py planted_real N SEED OUT.jsonl INDEX   (plants from survivor INDEX)
"""
import sys, os, json, time, itertools as it
import numpy as np
import numba as nb

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
ALPH = np.array(sorted(s * v for v in (1, 13, 15, 27, 29, 41, 43) for s in (1, -1)), np.int64)
ASET = set(ALPH.tolist())
M = 5
NC = 20
SWQ = (0, 2, 1, 3)
TAU = np.array([SWQ[h // M] * M + (-(h % M)) % M for h in range(NC)], np.int64)        # cell tau(h)
# SH[g, h] = index of cell h + g  (g = lag (a, b), cell (q, u) -> (q ^ a, u + b))
SH = np.array([[((h // M) ^ (g // M)) * M + (h % M + g % M) % M for h in range(NC)] for g in range(NC)], np.int64)
TGT = np.array([2936 if g == 0 else (-516 if g // M == 0 else 172) for g in range(NC)], np.int64)
EXC = 2936 - 80


# ------------------------------------------------------------------ rows
def rows_fixed(s):
    """(f0, p1, p2) with f0 + 2 p1 + 2 p2 = s; cells u = 0, {1, 4}, {2, 3}"""
    out = []
    for p1 in ALPH:
        for p2 in ALPH:
            f0 = s - 2 * int(p1) - 2 * int(p2)
            if f0 in ASET:
                out.append((f0, int(p1), int(p1), int(p2), int(p2)))     # order of cells u = 0, 1, 4, 2, 3 -> fix below
    res = []
    for f0, p1, _, p2, _ in out:
        res.append((f0, p1, p2, p2, p1))                                 # u = 0, 1, 2, 3, 4
    return res


_R5 = {}
def rows_free(s):
    if s not in _R5:
        out = []
        for m4 in it.product(ALPH.tolist(), repeat=4):
            last = s - sum(m4)
            if last in ASET:
                out.append(m4 + (last,))
        _R5[s] = out
    return _R5[s]


def excess(row):
    return sum(v * v - 1 for v in row)


def fixed_block_rows(z):
    """rows of a tau-invariant block: q = 00 row, q = 01 row (q = 10 row is y(10, u) = y(01, -u)), q = 11 row"""
    if z[1] != z[2]:
        raise ValueError('tau-invariant block needs z01 = z10')
    r00 = rows_fixed(z[0]); r01 = rows_free(z[1]); r11 = rows_fixed(z[3])
    return r00, r01, r11


def expand_fixed(r00, r01, r11):
    y = [0] * NC
    for u in range(M):
        y[0 * M + u] = r00[u]
        y[1 * M + u] = r01[u]
        y[2 * M + u] = r01[(-u) % M]
        y[3 * M + u] = r11[u]
    return y


def fixed_list(z, cap):
    r00, r01, r11 = fixed_block_rows(z)
    e00 = [excess(r) for r in r00]
    e01 = [2 * excess(r) for r in r01]           # the q = 01 row and its mirror q = 10
    e11 = [excess(r) for r in r11]
    o01 = np.argsort(e01, kind='stable')
    Y = []
    for i, r in enumerate(r00):
        for k, w in enumerate(r11):
            rest = cap - e00[i] - e11[k]
            if rest < 0:
                continue
            for j in o01:
                if e01[j] > rest:
                    break
                Y.append(expand_fixed(r, r01[j], w))
    return np.array(Y, np.int64).reshape(-1, NC)


def min_fixed(z):
    r00, r01, r11 = fixed_block_rows(z)
    if not (r00 and r01 and r11):
        return None
    return min(map(excess, r00)) + 2 * min(map(excess, r01)) + min(map(excess, r11))


def min_free(z):
    rs = [rows_free(s) for s in z]
    if not all(rs):
        return None
    return 2 * sum(min(map(excess, r)) for r in rs)


# ------------------------------------------------------------------ signatures and keys
def sig_block(Y):
    """H-level PAF vectors of the rows of Y (n x 20): S[n, g] = sum_h Y[n, h] Y[n, h + g]"""
    Y = np.asarray(Y, np.int64)
    return np.stack([(Y * Y[:, SH[g]]).sum(1) for g in range(NC)], 1)


TAULAG = TAU   # tau acts on lags exactly as on cells (it is an automorphism of H)


@nb.njit(cache=True)
def _enum_C(R0, R1, R2, R3, E0, E1, E2, E3, cap, SH, TAULAG, RK, out_key, out_idx, out_e):
    """all (i0, i1, i2, i3) with E0[i0] + ... <= cap (E already doubled, rows sorted by E); key of PAF_C + PAF_C o tau"""
    n = 0
    y = np.zeros(20, np.int64)
    wk = np.zeros(20, np.int64)
    for g in range(20):
        wk[g] = RK[g] + RK[TAULAG[g]]           # key(PAF o tau) = sum_g RK[g] PAF(tau g) = sum_g RK[tau g] PAF(g)
    for i0 in range(R0.shape[0]):
        e0 = E0[i0]
        if e0 > cap:
            break
        for i1 in range(R1.shape[0]):
            e1 = e0 + E1[i1]
            if e1 > cap:
                break
            for i2 in range(R2.shape[0]):
                e2 = e1 + E2[i2]
                if e2 > cap:
                    break
                for i3 in range(R3.shape[0]):
                    e3 = e2 + E3[i3]
                    if e3 > cap:
                        break
                    for u in range(5):
                        y[u] = R0[i0, u]; y[5 + u] = R1[i1, u]; y[10 + u] = R2[i2, u]; y[15 + u] = R3[i3, u]
                    k = np.int64(0)
                    for g in range(20):
                        p = np.int64(0)
                        for h in range(20):
                            p += y[h] * y[SH[g, h]]
                        k += wk[g] * p
                    if n < out_key.shape[0]:
                        out_key[n] = k
                        out_idx[n, 0] = i0; out_idx[n, 1] = i1; out_idx[n, 2] = i2; out_idx[n, 3] = i3
                        out_e[n] = e3
                    n += 1
    return n


@nb.njit(cache=True)
def _build_table(ukeys, bits):
    size = 1 << bits
    mask = size - 1
    tab = -np.ones(size, np.int64)
    for i in range(ukeys.shape[0]):
        h = (ukeys[i] * np.int64(-7046029254386353131)) >> np.int64(64 - bits)
        h &= mask
        while tab[h] >= 0:
            h = (h + 1) & mask
        tab[h] = i
    return tab


@nb.njit(cache=True)
def _join(kA, kB, tk, ukeys, tab, bits, out):
    size = 1 << bits
    mask = size - 1
    n = 0
    for a in range(kA.shape[0]):
        base = tk - kA[a]
        for b in range(kB.shape[0]):
            need = base - kB[b]
            h = (need * np.int64(-7046029254386353131)) >> np.int64(64 - bits)
            h &= mask
            while tab[h] >= 0:
                if ukeys[tab[h]] == need:
                    if n < out.shape[0]:
                        out[n, 0] = a; out[n, 1] = b; out[n, 2] = tab[h]
                    n += 1
                    break
                h = (h + 1) & mask
    return n


def recheck_direct(A, B, C, zA, zB, zC, tgt):
    """pure-Python re-check on the explicit arrays (A, B, C, tau C): all 20 lags, alphabet, tau-invariance, Q sums"""
    tC = [C[int(TAU[h])] for h in range(NC)]
    X = [list(map(int, A)), list(map(int, B)), list(map(int, C)), tC]
    for g in range(NC):
        tot = 0
        for x in X:
            for q in range(4):
                for u in range(M):
                    a, b = divmod(g, M)
                    tot += x[q * M + u] * x[(q ^ a) * M + (u + b) % M]
        if tot != int(tgt[g]):
            return False
    if any(v not in ASET for x in X for v in x):
        return False
    if any(A[h] != A[int(TAU[h])] or B[h] != B[int(TAU[h])] for h in range(NC)):
        return False
    for x, z in ((X[0], zA), (X[1], zB), (X[2], zC), (X[3], (zC[0], zC[2], zC[1], zC[3]))):
        if [sum(x[q * M:(q + 1) * M]) for q in range(4)] != list(z):
            return False
    return True


def decide(zA, zB, zC, tgt=None, etot=None, seed=5, max_hits=10 ** 6):
    tgt = TGT if tgt is None else np.asarray(tgt, np.int64)
    etot = EXC if etot is None else etot
    t0 = time.time()
    mA, mB, mC = min_fixed(zA), min_fixed(zB), min_free(zC)
    res = dict(z=[list(zA), list(zB), list(zC), [zC[0], zC[2], zC[1], zC[3]]], mins=[mA, mB, mC])
    if None in (mA, mB, mC) or mA + mB + mC > etot:
        res.update(status='EMPTY (Q sums or energy)', sizes=None, hits=0, solutions=[])
        return res
    LA = fixed_list(zA, etot - mB - mC)
    LB = fixed_list(zB, etot - mA - mC)
    capC = etot - mA - mB
    rows = [rows_free(s) for s in zC]
    Rs, Es = [], []
    for r in rows:
        e = np.array([2 * excess(x) for x in r], np.int64)
        o = np.argsort(e, kind='stable')
        Rs.append(np.array(r, np.int64)[o]); Es.append(e[o])
    rng = np.random.default_rng(seed)
    RK = rng.integers(1, 2 ** 62, size=NC, dtype=np.int64) | 1
    # count first (cheap pass with zero-length buffers), then fill
    nC = _enum_C(*Rs, *Es, capC, SH, TAULAG, RK, np.zeros(0, np.int64), np.zeros((0, 4), np.int32), np.zeros(0, np.int64))
    kC = np.zeros(nC, np.int64); iC = np.zeros((nC, 4), np.int32); eC = np.zeros(nC, np.int64)
    _enum_C(*Rs, *Es, capC, SH, TAULAG, RK, kC, iC, eC)
    t_enum = time.time() - t0
    SA, SB = sig_block(LA), sig_block(LB)
    kA = (SA * RK).sum(1); kB = (SB * RK).sum(1)              # int64 wraparound: linear in the vector
    tk = np.int64((tgt * RK).sum())
    order = np.argsort(kC, kind='stable')
    kCs = kC[order]
    ukeys, start = np.unique(kCs, return_index=True)
    bits = max(4, int(np.ceil(np.log2(max(2, len(ukeys))))) + 1)
    tab = _build_table(ukeys, bits)
    out = np.zeros((max_hits, 3), np.int64)
    nh = _join(kA, kB, tk, ukeys, tab, bits, out)
    t_join = time.time() - t0 - t_enum
    sols = []
    nverified = 0
    stop = np.append(start[1:], len(kCs))
    for a, b, u in out[:min(nh, max_hits)]:
        for pos in range(start[u], stop[u]):
            c = order[pos]
            Cy = np.concatenate([Rs[q][iC[c, q]] for q in range(4)])
            sv = SA[a] + SB[b] + sig_block(Cy[None, :])[0] + sig_block(Cy[TAU][None, :])[0]
            if not (sv == tgt).all():
                continue                                        # key collision: discarded by the integer check
            ok = recheck_direct(LA[a].tolist(), LB[b].tolist(), Cy.tolist(), zA, zB, zC, tgt)
            nverified += ok
            sols.append(dict(A=LA[a].tolist(), B=LB[b].tolist(), C=Cy.tolist(), recheck=ok))
    res.update(status='SOLUTIONS' if sols else 'EMPTY (exhaustive join)', sizes=[len(LA), len(LB), int(nC)],
               distinct_C_sigs=int(len(ukeys)), key_matches=int(nh), hits=len(sols), verified=nverified,
               secs_enum=round(t_enum, 1), secs_join=round(t_join, 1), secs=round(time.time() - t0, 1),
               solutions=sols)
    return res


def survivors():
    d = json.load(open(os.path.join(ROOT, 'logs', 'blockswap', 'gate0_m5.json')))
    return [t['z'] for t in d['tuples'] if t['survives']]


# ------------------------------------------------------------------ planted test
def random_ansatz(rng, nbig=(2, 2, 2)):
    """random tau-invariant A, B and free C over the alphabet: mostly +-1, a few letters 13/15 (random positions)"""
    def letter():
        return int(rng.choice([13, 15])) * int(rng.choice([-1, 1]))
    out = []
    for blk in range(3):
        y = [int(v) for v in rng.choice([-1, 1], size=NC)]
        pos = rng.choice(NC, size=nbig[blk], replace=False)
        for p in pos:
            y[p] = letter()
        if blk < 2:                                          # symmetrise: y(h) := y(tau h) on each orbit
            for h in range(NC):
                if TAU[h] < h:
                    y[h] = y[int(TAU[h])]
        out.append(y)
    return out


def planted(n, seed, outp):
    rng = np.random.default_rng(seed)
    for i in range(n):
        A, B, C = random_ansatz(rng, nbig=tuple(int(x) for x in rng.integers(1, 4, size=3)))
        tC = [C[int(TAU[h])] for h in range(NC)]
        Y = np.array([A, B, C, tC], np.int64)
        tgt = sig_block(Y).sum(0)
        etot = int(tgt[0]) - 80
        zs = [tuple(int(sum(y[q * M:(q + 1) * M])) for q in range(4)) for y in (A, B, C)]
        r = decide(zs[0], zs[1], zs[2], tgt=tgt, etot=etot, seed=seed + i)
        found = any(s['A'] == A and s['B'] == B and s['C'] == C for s in r['solutions'])
        rec = dict(trial=i, seed=seed, z=r['z'], sizes=r['sizes'], hits=r['hits'], verified=r.get('verified'),
                   planted_recovered=found, secs=r.get('secs'), target=tgt.tolist(), plant=dict(A=A, B=B, C=C))
        s = json.dumps(rec)
        print(s[:300], flush=True)
        with open(outp, 'a') as f:
            f.write(s + '\n')


def planted_real(n, seed, outp, index=0):
    """plant (A, B, C) drawn uniformly from the real lists of survivor `index` (C uniform among the energy-capped C
    list, then A, B uniform among the pairs completing the excess 2856 exactly), set the target to its 20-lag sum, and
    decide with exactly the real lists (the default energy 2856 holds because the plant's lag 0 is 2936)."""
    rng = np.random.default_rng(seed)
    zA, zB, zC, _ = [tuple(z) for z in survivors()[index]]
    mA, mB, mC = min_fixed(zA), min_fixed(zB), min_free(zC)
    LA = fixed_list(zA, EXC - mB - mC); LB = fixed_list(zB, EXC - mA - mC)
    eA = (LA * LA - 1).sum(1); eB = (LB * LB - 1).sum(1)
    capC = EXC - mA - mB
    rows = [rows_free(s_) for s_ in zC]
    Rs, Es = [], []
    for r in rows:
        e = np.array([2 * excess(x) for x in r], np.int64)
        o = np.argsort(e, kind='stable')
        Rs.append(np.array(r, np.int64)[o]); Es.append(e[o])
    RK = np.ones(NC, np.int64)
    nC = _enum_C(*Rs, *Es, capC, SH, TAULAG, RK, np.zeros(0, np.int64), np.zeros((0, 4), np.int32), np.zeros(0, np.int64))
    kC = np.zeros(nC, np.int64); iC = np.zeros((nC, 4), np.int32); eCl = np.zeros(nC, np.int64)
    _enum_C(*Rs, *Es, capC, SH, TAULAG, RK, kC, iC, eCl)
    del kC
    done = 0
    while done < n:
        c = int(rng.integers(nC))                     # uniform in the energy-capped C list
        Cy = np.concatenate([Rs[q][iC[c, q]] for q in range(4)]).tolist()
        ec = int(eCl[c])
        ia = rng.integers(len(LA))
        jb = np.nonzero(eB == EXC - ec - eA[ia])[0]
        if len(jb) == 0:
            continue
        ib = int(rng.choice(jb))
        A, B = LA[ia].tolist(), LB[ib].tolist()
        tC = [Cy[int(TAU[h])] for h in range(NC)]
        tgt = sig_block(np.array([A, B, Cy, tC], np.int64)).sum(0)
        assert tgt[0] == 2936
        r = decide(zA, zB, zC, tgt=tgt, seed=seed + done)
        found = any(s_['A'] == A and s_['B'] == B and s_['C'] == Cy for s_ in r['solutions'])
        rec = dict(trial=done, seed=seed, index=index, z=r['z'], sizes=r['sizes'], hits=r['hits'],
                   verified=r.get('verified'), planted_recovered=found, secs=r.get('secs'), target=tgt.tolist(),
                   plant=dict(A=A, B=B, C=Cy))
        s_ = json.dumps(rec)
        print(s_[:300], flush=True)
        with open(outp, 'a') as f:
            f.write(s_ + '\n')
        done += 1


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'run':
        first, last, outp = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        Z = survivors()
        for n in range(first, min(last, len(Z))):
            zA, zB, zC, _ = [tuple(z) for z in Z[n]]
            r = decide(zA, zB, zC)
            r['index'] = n
            s = json.dumps(r)
            print(s[:400], flush=True)
            with open(outp, 'a') as f:
                f.write(s + '\n')
    elif cmd == 'planted':
        planted(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif cmd == 'planted_real':
        planted_real(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], index=int(sys.argv[5]))
