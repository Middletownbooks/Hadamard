"""Exact Williamson search on Z_n, n = p*q, gcd(p,q)=1 (pre-registered in williamson/PREREGISTRATION.md).

Stage A: symmetric compressions onto Z_p and Z_q (PSD bound), joined into complementary quadruples.
Stage B: combine Z_p x Z_q quadruples with the same row sums; reduce mod (Z_p)*/+-1 x (Z_q)*/+-1 (acts independently
         by CRT, so orbit reps of pairs = pairs of orbit reps).
Stage C: per block, every symmetric +-1 p x q array with the given margins (exact MITM over columns), PSD-filtered.
Stage D: meet in the middle (A,B) vs (C,D) on the PAF sum plus the product condition a_j b_j = eps c_j d_j (j != 0).
Nothing is sampled or capped. Usage: wil_search.py n p q [--no-mult] [--out FILE]"""
import sys, json, time, itertools, math, numpy as np
from wil_common import expand, verify, rowsum_cases, epsilon

T0 = time.time(); LAST = [T0]
import os
JOIN = os.environ.get("WIL_JOIN", "nb")
CACHE_LIMIT = float(os.environ.get("WIL_CACHE_GB", "2")) * 2**30; CB = [0, 0]   # bytes, clears
if JOIN == "nb":
    from wil_join_nb import Joiner
    JOINER = Joiner()
def log(**k):
    k["t"] = round(time.time() - T0, 1); print(json.dumps(k), flush=True)

# ---------------- stage A ----------------
def comp_seqs(n, m):
    """all symmetric compressed halves on Z_m with entries of parity n/m, |.|<=n/m, passing PSD<=4n"""
    k = n // m; hm = (m + 1) // 2; vals = list(range(-k, k + 1, 2))
    H = np.array(list(itertools.product(vals, repeat=hm)), dtype=np.int64)
    r = np.arange(1, hm)
    Xh = H[:, [0]] + 2 * H[:, 1:] @ np.cos(2 * np.pi * np.outer(r, r) / m) if hm > 1 else H[:, :1] * 0
    ok = np.all(Xh ** 2 <= 4 * n + 1e-6, 1)
    H, Xh = H[ok], Xh[ok]
    S = H[:, 0] + 2 * H[:, 1:].sum(1)
    X = expand(H, m)
    P = np.stack([np.sum(X * np.roll(X, -s, 1), 1) for s in range(1, hm)], 1) if hm > 1 else np.zeros((len(H), 0), np.int64)
    return H, S, P, Xh

def pair_keys(PA, PB, XA, XB, n, W, sign):
    """all pairs (i,j) with pair-PSD <= 4n; returns i, j, hashed PAF-sum key (sign * (PA_i+PB_j))"""
    out_i, out_j, out_k = [], [], []
    step = max(1, 4_000_000 // max(1, len(PB)))
    for a in range(0, len(PA), step):
        ps = XA[a:a + step, None, :] ** 2 + XB[None, :, :] ** 2
        ii, jj = np.nonzero(np.all(ps <= 4 * n + 1e-6, 2))
        ii += a
        key = (sign * (PA[ii] + PB[jj])).astype(np.uint64) @ W if PA.shape[1] else np.zeros(len(ii), np.uint64)
        out_i.append(ii); out_j.append(jj); out_k.append(key)
    cat = lambda L, dt: np.concatenate(L).astype(dt) if L else np.zeros(0, dt)
    return cat(out_i, np.int64), cat(out_j, np.int64), cat(out_k, np.uint64)

def match(kab, kcd):
    """all (x,y) with kab[x]==kcd[y]"""
    o = np.argsort(kcd, kind="stable"); ks = kcd[o]
    lo = np.searchsorted(ks, kab, "left"); hi = np.searchsorted(ks, kab, "right")
    cnt = hi - lo; x = np.repeat(np.arange(len(kab)), cnt)
    if len(x) == 0: return x, x
    off = np.arange(len(x)) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    return x, o[np.repeat(lo, cnt) + off]

def comp_quads(n, m, case, rng):
    H, S, P, Xh = comp_seqs(n, m)
    W = rng.integers(1, 2**63, P.shape[1], dtype=np.uint64) | np.uint64(1)
    L = [np.where(S == s)[0] for s in case]
    ai, bj, kab = pair_keys(P[L[0]], P[L[1]], Xh[L[0]], Xh[L[1]], n, W, 1)
    ci, dj, kcd = pair_keys(P[L[2]], P[L[3]], Xh[L[2]], Xh[L[3]], n, W, -1)
    x, y = match(kab, kcd)
    Q = np.stack([L[0][ai[x]], L[1][bj[x]], L[2][ci[y]], L[3][dj[y]]], 1) if len(x) else np.zeros((0, 4), np.int64)
    # exact check: compressed PAF sum = 0 off zero, all four PSDs sum to 4n
    good = np.all(P[Q].sum(1) == 0, 1) & np.all(np.abs((Xh[Q] ** 2).sum(1) - 4 * n) < 1e-6, 1)
    return H, H[Q[good]]                     # (all halves, quadruples as halves [N,4,hm])

def canon(Qh, m):
    """lexicographic min over u in (Z_m)*/+-1 of the quadruple (block order fixed)"""
    us = [u for u in range(1, (m + 1) // 2) if math.gcd(u, m) == 1]
    X = expand(Qh, m)                         # [N,4,m]
    best = None
    for u in us:
        perm = (u * np.arange(m)) % m          # y[i] = x[u i]
        Y = X[..., perm][..., :(m + 1) // 2].reshape(len(X), -1)
        best = Y if best is None else np.where(lex_less(Y, best)[:, None], Y, best)
    return best.reshape(Qh.shape)

def lex_less(A, B):
    d = A != B; first = np.argmax(d, 1); has = d.any(1)
    return has & (A[np.arange(len(A)), first] < B[np.arange(len(A)), first])

def uniq_rows(Q):
    if len(Q) == 0: return Q
    return np.unique(Q.reshape(len(Q), -1), axis=0).reshape((-1,) + Q.shape[1:])

# ---------------- stage C ----------------
def lifts(n, p, q, R, C):
    """all symmetric +-1 p x q arrays (rows Z_p, cols Z_q) with half row sums R and half col sums C.
    Returns full sequences on Z_n (int8, index i <-> (i mod p, i mod q))."""
    hp, hq = (p + 1) // 2, (q + 1) // 2
    cols = []                                   # per free column: (options [k, p] as +-1 on rows, contrib [k, hp])
    # column 0: M[0,0], M[a,0]=M[-a,0] a=1..hp-1
    o = np.array(list(itertools.product([1, -1], repeat=hp)), dtype=np.int64)
    s0 = o[:, 0] + 2 * o[:, 1:].sum(1)
    o = o[s0 == C[0]]
    full0 = np.zeros((len(o), p), np.int64); full0[:, 0] = o[:, 0]
    for a in range(1, hp): full0[:, a] = o[:, a]; full0[:, p - a] = o[:, a]
    cols.append((full0, o.copy()))              # contrib: row0 += M00, row a += M[a,0]
    for b in range(1, hq):
        o = np.array(list(itertools.product([1, -1], repeat=p)), dtype=np.int64)
        o = o[o.sum(1) == C[b]]
        con = np.zeros((len(o), hp), np.int64); con[:, 0] = 2 * o[:, 0]
        for a in range(1, hp): con[:, a] = o[:, a] + o[:, p - a]
        cols.append((o, con))
    if any(len(c[0]) == 0 for c in cols): return np.zeros((0, n), np.int8)
    half = len(cols) // 2
    def side(cs):
        idx = np.array(list(itertools.product(*[range(len(c[0])) for c in cs])), dtype=np.int64).reshape(-1, len(cs))
        tot = sum(cs[t][1][idx[:, t]] for t in range(len(cs))) if cs else np.zeros((1, hp), np.int64)
        return idx, tot
    i1, t1 = side(cols[:half]); i2, t2 = side(cols[half:])
    base = np.int64(4 * q + 1); w = base ** np.arange(hp)
    k1 = (t1 + 2 * q) @ w; k2 = (np.asarray(R) - t2 + 2 * q) @ w
    x, y = match(k1, k2)
    if len(x) == 0: return np.zeros((0, n), np.int8)
    M = np.zeros((len(x), p, q), np.int8)
    for t, (opts, _) in enumerate(cols):
        ch = i1[x, t] if t < half else i2[y, t - half]
        M[:, :, t] = opts[ch]
        if t: M[:, :, q - t] = opts[ch][:, (-np.arange(p)) % p]   # M[a,-b] = M[-a,b]
    i = np.arange(n); X = M[:, i % p, i % q]
    assert np.all(X.sum(1) == np.asarray(R)[0] + 2 * np.asarray(R)[1:].sum())
    return X

def lift_table(n, X):
    """PSD filter, then (half PAF lags 1..(n-1)/2, -1 mask on positions 1..(n-1)/2)"""
    F = np.fft.rfft(X.astype(np.float64), axis=1)
    ok = np.all(np.abs(F[:, 1:]) ** 2 <= 4 * n + 1e-6, 1)
    X = X[ok]; hn = (n + 1) // 2
    P = np.rint(np.fft.irfft(np.abs(np.fft.rfft(X.astype(np.float64), axis=1)) ** 2, n=n, axis=1)).astype(np.int64)[:, 1:hn]
    bits = (X[:, 1:hn] < 0).astype(np.uint64) @ (np.uint64(1) << np.arange(hn - 1, dtype=np.uint64))
    F = np.ascontiguousarray(np.abs(np.fft.rfft(X.astype(np.float64), axis=1)[:, 1:]) ** 2)
    return X, P, bits, F

# ---------------- stage D ----------------
def join4(n, T, eps, W, K1):
    """T: 4 lift tables. returns list of solution quadruples (full sequences)."""
    (XA, PA, MA, FA), (XB, PB, MB, FB), (XC, PC, MC, FC), (XD, PD, MD, FD) = T
    hn = (n + 1) // 2
    E = np.uint64((1 << (hn - 1)) - 1) if eps == -1 else np.uint64(0)
    def keys(P1, M1, P2, M2, sign, flip):
        i = np.repeat(np.arange(len(P1)), len(P2)); j = np.tile(np.arange(len(P2)), len(P1))
        k = ((M1[i] ^ M2[j] ^ flip) * K1) + (sign * (P1[i] + P2[j])).astype(np.uint64) @ W
        return i, j, k
    sols = []; work = len(XA) * len(XB) + len(XC) * len(XD)
    if min(len(XA), len(XB), len(XC), len(XD)) == 0: return sols, work
    if JOIN == "nb":
        h = [(P.astype(np.uint64) @ W) for P in (PA, PB, PC, PD)]
        F = (FA, FB, FC, FD)
        for a, b, c, d in JOINER.run(h[0], MA, F[0], h[1], MB, F[1], h[2], MC, F[2], h[3], MD, F[3], int(E), int(K1), 4.0 * n + 1e-6):
            Q = np.stack([XA[a], XB[b], XC[c], XD[d]]).astype(np.int64)
            if verify(Q): sols.append(Q)
        return sols, work
    ci, dj, kcd = keys(PC, MC, PD, MD, -1, E)
    o = np.argsort(kcd); kcd = kcd[o]; ci, dj = ci[o], dj[o]
    step = max(1, 20_000_000 // max(1, len(XB)))
    for a in range(0, len(XA), step):
        ai, bj, kab = keys(PA[a:a + step], MA[a:a + step], PB, MB, 1, np.uint64(0))
        lo = np.searchsorted(kcd, kab, "left"); hi = np.searchsorted(kcd, kab, "right")
        for t in np.nonzero(hi > lo)[0]:
            for u in range(lo[t], hi[t]):
                Q = np.stack([XA[a + ai[t]], XB[bj[t]], XC[ci[u]], XD[dj[u]]]).astype(np.int64)
                if verify(Q): sols.append(Q)
    return sols, work

# ---------------- driver ----------------
def main():
    n, p, q = map(int, sys.argv[1:4]); assert p * q == n and math.gcd(p, q) == 1 and p < q
    mult = "--no-mult" not in sys.argv; cost_only = "--cost-only" in sys.argv
    global SHI, SHK
    SHI, SHK = map(int, sys.argv[sys.argv.index("--shard") + 1].split("/")) if "--shard" in sys.argv else (0, 1)
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else f"logs/wil_n{n}{'' if mult else '_nomult'}{'_cost' if cost_only else ''}{('_c' + sys.argv[sys.argv.index('--case') + 1]) if '--case' in sys.argv else ''}{f'_s{SHI}of{SHK}' if SHK > 1 else ''}.json"
    rng = np.random.default_rng(12345)
    hn = (n + 1) // 2
    W = rng.integers(1, 2**63, hn - 1, dtype=np.uint64) | np.uint64(1)
    K1 = np.uint64(rng.integers(1, 2**63, dtype=np.uint64) | 1)
    global CKPT
    CKPT = out + ".ckpt" if "--no-ckpt" not in sys.argv else os.devnull
    rec = {"n": n, "p": p, "q": q, "mult_reduction": mult, "cost_only": cost_only, "shard": [SHI, SHK], "cases": [], "solutions": [], "caps_hit": False}
    lift_cache = {}
    CASES = rowsum_cases(n)
    if "--case" in sys.argv: CASES = [CASES[int(sys.argv[sys.argv.index("--case") + 1])]]
    for case in CASES:
        eps = epsilon(n, case)
        _, Qp = comp_quads(n, p, case, rng); _, Qq = comp_quads(n, q, case, rng)
        nqp, nqq = len(Qp), len(Qq)
        if mult: Qp, Qq = uniq_rows(canon(Qp, p)), uniq_rows(canon(Qq, q))
        cr = {"case": case, "eps": eps, "Qp": nqp, "Qq": nqq, "Qp_reps": len(Qp), "Qq_reps": len(Qq),
              "combos": len(Qp) * len(Qq), "combos_this_shard": sum(1 for t in range(len(Qp)) if t % SHK == SHI) * len(Qq), "liftable": 0, "work": 0, "ledger": 0, "sols": 0}
        log(stage="A", **{k: cr[k] for k in ("case", "Qp", "Qq", "Qp_reps", "Qq_reps")})
        ck = {}                                   # checkpoint: finished Z_p reps of this case (survives restarts)
        if os.path.exists(CKPT):
            for line in open(CKPT):
                try: d = json.loads(line)
                except ValueError: continue          # a torn last line from a kill is ignored (that rep reruns)
                if tuple(d["case"]) == tuple(case): ck[d["ip"]] = d
        for ip, qp in enumerate(Qp):
            if ip % SHK != SHI: continue
            if ip in ck:
                d = ck[ip]
                for k in ("ledger", "liftable", "work", "sols"): cr[k] += d[k]
                rec["solutions"] += d["solutions"]; cr["resumed_reps"] = cr.get("resumed_reps", 0) + 1
                continue
            before = {k: cr[k] for k in ("ledger", "liftable", "work", "sols")}; nsol0 = len(rec["solutions"])
            for iq, qq in enumerate(Qq):
                T = []
                for blk in range(4):
                    key = (qp[blk].tobytes(), qq[blk].tobytes())
                    if key not in lift_cache:
                        lt = lift_table(n, lifts(n, p, q, qp[blk], qq[blk]))
                        lift_cache[key] = (np.empty(len(lt[0]), np.int8), None, None, None) if cost_only else lt
                        CB[0] += sum(a.nbytes for a in lift_cache[key] if a is not None)
                        if CB[0] > CACHE_LIMIT:            # bounded memory; results unaffected (lifts are recomputed)
                            lift_cache.clear(); CB[0] = 0; CB[1] += 1
                            lift_cache[key] = lt if not cost_only else (np.empty(len(lt[0]), np.int8), None, None, None)
                    T.append(lift_cache[key])
                cr["ledger"] += 1
                if time.time() - LAST[0] > 60:
                    LAST[0] = time.time()
                    log(stage="progress", case=case, ledger=cr["ledger"], of=cr["combos"], work=cr["work"], sols=cr["sols"], cache=len(lift_cache))
                if min(len(t[0]) for t in T) == 0: continue
                cr["liftable"] += 1
                if cost_only:
                    L = [len(t[0]) for t in T]; sols, work = [], L[0] * L[1] + L[2] * L[3]
                    cr["max_pairs"] = max(cr.get("max_pairs", 0), L[0] * L[1], L[2] * L[3])
                else:
                    sols, work = join4(n, T, eps, W, K1)
                cr["work"] += work; cr["sols"] += len(sols)
                for Q in sols: rec["solutions"].append({"case": case, "halves": Q[:, :hn].tolist()})
            d = {"case": list(case), "ip": ip, **{k: cr[k] - before[k] for k in before}, "solutions": rec["solutions"][nsol0:]}
            with open(CKPT, "a") as fh: fh.write(json.dumps(d) + "\n"); fh.flush(); os.fsync(fh.fileno())
        log(stage="case_done", **cr)
        rec["cases"].append(cr)
    rec["N_found"] = len(rec["solutions"]); rec["N_expanded"] = orbit_count(n, [s["halves"] for s in rec["solutions"]]); rec["secs"] = round(time.time() - T0, 1)
    rec["lift_cache_entries"] = len(lift_cache); rec["cache_clears"] = CB[1]
    json.dump(rec, open(out, "w"))
    log(stage="done", N_found=rec["N_found"], N_expanded=rec["N_expanded"], out=out)


def orbit_count(n, sols_halves):
    """number of distinct quadruples in the union of (Z_n)*-orbits of the given solutions (block order fixed)"""
    hn = (n + 1) // 2; seen = set()
    for h in sols_halves:
        X = expand(np.array(h), n)
        for u in range(1, n):
            if math.gcd(u, n) != 1: continue
            seen.add(X[:, (u * np.arange(n)) % n][:, :hn].tobytes())
    return len(seen)

if __name__ == "__main__":
    main()
