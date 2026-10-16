"""2-list search for a b = 0 Goethals-Seidel core of prime order v in the 3-cycle-twisted multiplier slice
(scripts/gs_twist.py, g = 3-cycle on (B, C, D), A fixed):
    U = <u> of order d (3 | d);  A is U-invariant;  B is U3-invariant (U3 = <u^3>);
    C(x) = B(u^-1 x),  D(x) = B(u^-2 x).
ΣPAF = 4v·δ  <=>  for every class c of <U, -1> on Z_v^*:  PAF_A(c) + PAF_B(c) + PAF_B(u^-1 c) + PAF_B(u^-2 c) = 0.
Row sums (s_A, s_B, s_B, s_B) with s_A^2 + 3 s_B^2 = 4v.

A-table: every A with the required row sum and A(0) fixed (complete list). B: uniform random samples with the required
row sum and B(0) fixed. Each B is matched against the A-table by a hashed key of its class values; every hit is checked
exactly on all v lags, the full GS array is built and H·Hᵀ = 4v·I is verified in integer arithmetic.
Restart-safe: sampling is memoryless; each run appends JSON lines to OUT and writes every solution to OUT.sol.

Usage: b0_twist3_search.py v d sA sB OUT [--seed S] [--batch B] [--maxb N] [--exhaustB]
"""
import os, sys, json, math, time, itertools, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tmat_twist import gen_of_order


def setup(v, d):
    u = gen_of_order(v, d)
    def orbits(gen, order):
        seen, out = {0}, []
        for x in range(1, v):
            if x in seen: continue
            o = sorted({x * pow(gen, j, v) % v for j in range(order)}); seen.update(o); out.append(o)
        return out
    OA = orbits(u, d)                          # U-orbits (A)
    OB = orbits(pow(u, 3, v), d // 3)          # U3-orbits (B)
    seen, reps = set(), []
    for x in range(1, v):
        if x in seen: continue
        seen.update({(s * x * pow(u, j, v)) % v for j in range(d) for s in (1, -1)}); reps.append(x)
    ui = pow(u, -1, v)
    lagsB = np.array([[c, c * ui % v, c * ui * ui % v] for c in reps])   # lags whose PAF_B values sum into class c
    return u, OA, OB, np.array(reps), lagsB


def zero_value(s, size, m):
    """(x0, k): x0 = value at 0, k = sum of the m orbit signs, with x0 + size*k = s."""
    out = []
    for x0 in (1, -1):
        if (s - x0) % size: continue
        k = (s - x0) // size
        if abs(k) <= m and (k - m) % 2 == 0: out.append((x0, k))
    return out


def paf_all(X):
    F = np.fft.rfft(X, axis=1)
    return np.rint(np.fft.irfft(F.real ** 2 + F.imag ** 2, n=X.shape[1], axis=1)).astype(np.int64)


def build_core(v, u, A, B):
    ui = pow(u, -1, v)
    C = np.array([B[(ui * x) % v] for x in range(v)]); D = np.array([B[(ui * ui * x) % v] for x in range(v)])
    return A, B, C, D


def verify_full(v, A, B, C, D):
    X = np.stack([A, B, C, D]).astype(np.int64)
    P = np.zeros(v, np.int64)
    for k in range(4):
        for s in range(v): P[s] += int(np.dot(X[k], np.roll(X[k], -s)))
    if not (P[0] == 4 * v and (P[1:] == 0).all()): return False
    # Goethals-Seidel array with the back-identity R; check H Hᵀ = 4v I exactly
    def circ(a): return np.stack([np.roll(a, i) for i in range(v)])
    Ac, Bc, Cc, Dc = (circ(a) for a in (A, B, C, D))
    R = np.fliplr(np.eye(v, dtype=np.int64))
    H = np.block([[Ac, Bc @ R, Cc @ R, Dc @ R],
                  [-Bc @ R, Ac, Dc.T @ R, -Cc.T @ R],
                  [-Cc @ R, -Dc.T @ R, Ac, Bc.T @ R],
                  [-Dc @ R, Cc.T @ R, -Bc.T @ R, Ac]])
    G = H @ H.T
    return bool((G == 4 * v * np.eye(4 * v, dtype=np.int64)).all())


def main():
    v, d, sA, sB, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    arg = lambda f, dflt, ty=int: ty(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else dflt
    seed, batch, maxb = arg("--seed", 1), arg("--batch", 1 << 14), arg("--maxb", 10 ** 9)
    exhaustB = "--exhaustB" in sys.argv
    u, OA, OB, reps, lagsB = setup(v, d)
    mA, mB, C = len(OA), len(OB), len(reps)
    rng = np.random.default_rng(seed)
    W = rng.integers(1, 2 ** 40, C).astype(np.int64)
    t0 = time.time()
    log = open(OUT, "a")
    def emit(**kw):
        kw["secs"] = round(time.time() - t0, 1); s = json.dumps(kw); print(s, flush=True); log.write(s + "\n"); log.flush()
    zA, zB = zero_value(sA, d, mA), zero_value(sB, d // 3, mB)
    emit(stage="setup", v=v, d=d, u=u, orbitsA=mA, orbitsB=mB, classes=C, sA=sA, sB=sB, zeroA=zA, zeroB=zB, seed=seed)
    if not zA or not zB:
        emit(stage="infeasible"); return
    # A-table (complete)
    a0, kA = zA[0]; npos = (mA + kA) // 2
    combs = np.array(list(itertools.combinations(range(mA), npos)), dtype=np.int64)
    SA = -np.ones((len(combs), mA), np.int8); np.put_along_axis(SA, combs, 1, axis=1)
    XA = np.zeros((len(SA), v), np.float64); XA[:, 0] = a0
    for j, o in enumerate(OA): XA[:, o] = SA[:, [j]]
    PA = paf_all(XA)[:, reps]                                   # |A| x C
    keyA = PA @ W
    order = np.argsort(keyA, kind="stable"); keyA_s = keyA[order]
    emit(stage="A-table", size=len(SA), log2=round(math.log2(len(SA)), 2))
    b0, kB = zB[0]; nposB = (mB + kB) // 2
    idxB = [np.array(o) for o in OB]
    # Fast path: T_B(c) = beta^T Q_c beta, beta = (B(0), b_1..b_mB), Q_c = sum over the 3 lags of E^T S_lag E (sym.)
    Q = None
    if "--fft" not in sys.argv:
        E = np.zeros((v, mB + 1), np.float64); E[0, 0] = 1
        for j, o in enumerate(OB): E[o, j + 1] = 1
        Qs = []
        for c in range(C):
            M = np.zeros((mB + 1, mB + 1))
            for lag in lagsB[c]:
                M += E.T @ np.roll(E, -int(lag), axis=0)          # (E^T S E)[i,j] = #{x in i: x+lag in j}
            Qs.append((M + M.T) / 2)
        Q = np.concatenate(Qs, axis=1).astype(np.float32)          # (mB+1) x (C*(mB+1))
        emit(stage="fast-path", Q_shape=list(Q.shape))
    nsol, done = 0, 0
    if exhaustB:
        allc = itertools.combinations(range(mB), nposB)
    for it in range(maxb):
        if exhaustB:
            chunk = list(itertools.islice(allc, batch))
            if not chunk: break
            SB = -np.ones((len(chunk), mB), np.int8); np.put_along_axis(SB, np.array(chunk), 1, axis=1)
        else:
            # uniform random subsets of size nposB
            r = rng.random((batch, mB)); pos = np.argpartition(r, nposB - 1, axis=1)[:, :nposB]
            SB = -np.ones((batch, mB), np.int8); np.put_along_axis(SB, pos, 1, axis=1)
        if Q is not None:                                            # orbit-level quadratic forms (fast path)
            beta = np.concatenate([np.full((len(SB), 1), b0, np.float32), SB.astype(np.float32)], axis=1)
            Y = (beta @ Q).reshape(len(SB), C, mB + 1)
            TB = np.rint(np.einsum("ncj,nj->nc", Y, beta)).astype(np.int64)
            XB = None
        else:
            XB = np.zeros((len(SB), v), np.float64); XB[:, 0] = b0
            for j, o in enumerate(idxB): XB[:, o] = SB[:, [j]]
            PB = paf_all(XB)
            TB = PB[:, lagsB[:, 0]] + PB[:, lagsB[:, 1]] + PB[:, lagsB[:, 2]]   # batch x C
        key = -(TB @ W)
        pos_ = np.searchsorted(keyA_s, key); pos_[pos_ >= len(keyA_s)] = 0
        hit = np.flatnonzero(keyA_s[pos_] == key)
        for h in hit:
            lo = np.searchsorted(keyA_s, key[h], "left"); hi = np.searchsorted(keyA_s, key[h], "right")
            for ia in order[lo:hi]:
                if not ((PA[ia] + TB[h]) == 0).all(): continue
                A = XA[ia].astype(np.int64)
                if XB is None:
                    B = np.zeros(v, np.int64); B[0] = b0
                    for j, o in enumerate(idxB): B[o] = int(SB[h, j])
                else:
                    B = XB[h].astype(np.int64)
                A, B, Cc, D = build_core(v, u, A, B)
                ok = verify_full(v, A, B, Cc, D)
                nsol += 1
                rec = dict(stage="SOLUTION", verified=ok, n=nsol, A="".join("+" if x > 0 else "-" for x in A),
                           B="".join("+" if x > 0 else "-" for x in B))
                emit(**rec)
                with open(OUT + ".sol", "a") as f:
                    f.write(json.dumps(dict(v=v, d=d, u=u, verified=ok, A=rec["A"], B=rec["B"],
                                            C="".join("+" if x > 0 else "-" for x in Cc),
                                            D="".join("+" if x > 0 else "-" for x in D))) + "\n")
        done += len(SB)
        if it % 50 == 0 or exhaustB:
            emit(stage="progress", B_tested=done, log2_B=round(math.log2(done), 2), solutions=nsol)
    emit(stage="end", B_tested=done, solutions=nsol)


if __name__ == "__main__":
    main()
