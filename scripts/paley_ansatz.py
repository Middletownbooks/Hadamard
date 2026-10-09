"""Paley-coset ansatz for GS-type cores over G = G' x Z_p (cyclic G' = Z_n', gcd(n', p) = 1, so G = Z_{n'p}).

Forensics (report/alpoge_methods_forensics.md, verified independently): the published v=359 core is EXACTLY of this form
on its Z_11 factor (116 cosets Legendre-twisted, 12 constant), v=419 approximately.

Ansatz, per block k and coset g in G':   x(g, 0) = alpha_g,   x(g, z) = beta_g * chi(z)  (twisted)  or  gamma_g (constant)
for z != 0, chi = Legendre symbol mod p.  8 labels per coset: alpha in +-1, type in {twisted, constant}, sign in +-1.
Write x(g, .) = alpha_g*delta + u_g*chi + w_g*eps  (eps = indicator of Z_p^*; u = 0 or w = 0 by type).
PAF at lag (h, s) depends on s only through s = 0, chi(s) = +1, chi(s) = -1, so the core condition reduces EXACTLY to
3 equations per h in G':   R_k(h, c) = sum_g sum_{f1,f2} A_k[g,f1] A_k[g+h,f2] C[f1,f2,c],  c in {0, QR, NQR},
with C the cross-correlation table of (delta, chi, eps) on Z_p (computed numerically).
Targets: b=0 GS core (sum PAF = 0 off 0); b=1 core (sum PAF = -4 off 0, m = v - 1).
Search: exhaustive 4-list meet-in-the-middle over all 8^n' labelings per block (exact; feasible for n' <= ~5).
usage: paley_ansatz.py validate | paley_ansatz.py search v b n' p"""
import sys, json, itertools, time
import numpy as np

def legendre(p):
    chi = np.zeros(p, dtype=np.int64)
    for z in range(1, p): chi[z] = 1 if pow(z, (p - 1) // 2, p) == 1 else -1
    return chi

def corr_table(p):
    chi = legendre(p); delta = np.zeros(p, np.int64); delta[0] = 1; eps = 1 - delta
    F = [delta, chi, eps]
    nqr = next(z for z in range(1, p) if chi[z] == -1)
    C = np.zeros((3, 3, 3), np.int64)
    for i, j in itertools.product(range(3), repeat=2):
        for ci, s in enumerate([0, 1, nqr]):
            C[i, j, ci] = int(np.dot(F[i], np.roll(F[j], -s)))      # sum_z f_i(z) f_j(z + s)
    return C, nqr

LABELS = [(a, t, s) for a in (1, -1) for t in (0, 1) for s in (1, -1)]   # alpha, type (1 = twisted), sign
def label_coeffs():
    A = np.zeros((8, 3), np.int64)
    for i, (a, t, s) in enumerate(LABELS):
        A[i] = [a, s if t else 0, 0 if t else s]
    return A
LC = label_coeffs()

def all_blocks(n):
    idx = np.arange(8 ** n, dtype=np.int64)
    return np.stack([(idx // 8 ** j) % 8 for j in range(n)], 1)             # labels per coset

def reduced_paf(lab, n, C):
    """lab: (N, n) labels -> R: (N, n, 3) = PAF at (h, 0), (h, QR), (h, NQR) for h in Z_n'"""
    A = LC[lab]                                                              # (N, n, 3)
    R = np.zeros((lab.shape[0], n, 3), np.int64)
    for h in range(n):
        B = np.roll(A, -h, axis=1)                                           # B[:, g] = A[:, g + h]
        R[:, h, :] = np.einsum('ngi,ngj,ijc->nc', A, B, C)
    return R

def full_sequence(lab_row, n, p):
    chi = legendre(p); m = n * p; x = np.zeros(m, np.int64)
    for i in range(m):
        g, z = i % n, i % p
        a, t, s = LABELS[lab_row[g]]
        x[i] = a if z == 0 else (s * chi[z] if t else s)
    return x

def paf(x):
    F = np.fft.fft(x); return np.rint(np.fft.ifft(np.abs(F) ** 2).real).astype(np.int64)

def validate():
    rng = np.random.default_rng(0); ok = True
    for n, p in [(3, 7), (4, 11), (2, 13), (5, 3), (6, 7)]:
        C, nqr = corr_table(p); lab = rng.integers(0, 8, (20, n))
        R = reduced_paf(lab, n, C)
        for r in range(20):
            x = full_sequence(lab[r], n, p); P = paf(x)
            for h in range(n):
                for ci, s in enumerate([0, 1, nqr]):
                    i = next(i for i in range(n * p) if i % n == h and i % p == s)
                    ok &= P[i] == R[r, h, ci]
                # every s != 0 with the same chi must give the same PAF (the reduction is exact)
                for s in range(1, p):
                    i = next(i for i in range(n * p) if i % n == h and i % p == s)
                    ok &= P[i] == R[r, h, 1 if legendre(p)[s] == 1 else 2]
    print(json.dumps({"validate_reduced_paf_vs_fft": bool(ok)}))

def search(v, b, n, p, max_sols=5):
    m = n * p
    assert np.gcd(n, p) == 1 and m == v - b and b in (0, 1)
    C, nqr = corr_table(p)
    tgt = np.full((n, 3), -4 * b, np.int64); tgt[0, 0] = 4 * m               # sum_k PAF_k target (lag 0 = 4m)
    lab = all_blocks(n); R = reduced_paf(lab, n, C).reshape(len(lab), -1)
    T = tgt.reshape(-1)
    # 4-list MITM: pair sums S12 = R_i + R_j, need S12 + S34 = T. Hash pair sums.
    t0 = time.time(); N = len(R)
    keys = {}
    Rb = [tuple(r) for r in R]
    # dictionary of pair sums -> one witness (i, j)
    pair = {}
    for i in range(N):
        ri = R[i]
        S = R + ri
        for j in range(i, N):
            pair.setdefault(S[j].tobytes(), (i, j))
    sols = []
    for key, (i, j) in pair.items():
        need = (T - np.frombuffer(key, np.int64)).tobytes()
        if need in pair:
            k, l = pair[need]
            blocks = [full_sequence(lab[q], n, p) for q in (i, j, k, l)]
            tot = sum(paf(x) for x in blocks)
            if (tot[1:] == -4 * b).all():
                sols.append([x.tolist() for x in blocks])
                if len(sols) >= max_sols: break
    return {"v": v, "b": b, "n": n, "p": p, "configs_per_block": N, "distinct_pair_sums": len(pair),
            "solutions": len(sols), "secs": round(time.time() - t0, 1), "examples": sols[:1]}

if __name__ == "__main__" and sys.argv[1] in ("validate", "search"):
    if sys.argv[1] == "validate": validate()
    elif sys.argv[1] == "search":
        v, b, n, p = map(int, sys.argv[2:6]); res = search(v, b, n, p)
        print(json.dumps({k: res[k] for k in res if k != "examples"}))
        if res["solutions"]:
            json.dump(res, open(f"logs/paley_ansatz_v{v}_b{b}_n{n}_p{p}.json", "w"))


# ---------------- Gaussian estimate (calibrated against exact counts) ----------------
def exact_count(v, b, n, p):
    """number of ordered 4-tuples of labelings (exact, n small)"""
    C, nqr = corr_table(p); m = n * p
    tgt = np.full((n, 3), -4 * b, np.int64); tgt[0, 0] = 4 * m; T = tgt.reshape(-1)
    R = reduced_paf(all_blocks(n), n, C).reshape(8 ** n, -1)
    from collections import Counter
    pairs = Counter()
    for i in range(len(R)):
        S = R + R[i]
        for row in S: pairs[row.tobytes()] += 1
    tot = 0
    for key, cnt in pairs.items():
        need = (T - np.frombuffer(key, np.int64)).tobytes()
        tot += cnt * pairs.get(need, 0)
    return tot

def rowsum_feasible(v, b, n, p):
    """row sum of a block = sum_g alpha_g + (p-1) * sum_g w_g (twisted cosets sum to 0 on Z_p^*); exact DP over cosets"""
    per = {a + (p - 1) * w for a, t, s_ in LABELS for w in [0 if t else s_]}
    S = {0}
    for _ in range(n): S = {x + y for x in S for y in per}
    if b == 0:
        import itertools as it
        return any(a*a + b_*b_ + c*c + d*d == 4 * v for a, b_, c, d in it.product(sorted(S), repeat=4))
    return (2 in S or -2 in S) and 0 in S                                 # b=1 row sums (+-2, 0, 0, 0)

def gauss_estimate(v, b, n, p, N=200000, seed=1):
    if not rowsum_feasible(v, b, n, p): return {"v": v, "b": b, "n": n, "p": p, "feasible": False, "why": "row sums"}
    C, nqr = corr_table(p); m = n * p
    tgt = np.full((n, 3), -4 * b, np.int64); tgt[0, 0] = 4 * m; T = tgt.reshape(-1)
    rng = np.random.default_rng(seed)
    S = sum(reduced_paf(rng.integers(0, 8, (N, n)), n, C).reshape(N, -1) for _ in range(4))
    keep = [i for i in range(S.shape[1]) if i != 0]                       # drop lag (0,0) (always 4m)
    S = S[:, keep]; Tk = T[keep]
    # remove exactly duplicated coordinates (PAF(h,c) = PAF(-h,c')), checking target consistency
    cols, seen = [], {}
    for j in range(S.shape[1]):
        key = S[:2000, j].tobytes()
        if key in seen:
            if Tk[j] != Tk[seen[key]]: return {"feasible": False, "why": "duplicate coordinate with different target"}
            continue
        seen[key] = j; cols.append(j)
    S = S[:, cols]; Tk = Tk[cols]
    mu = S.mean(0); Cv = np.atleast_2d(np.cov(S.T))
    w, V = np.linalg.eigh(Cv); keepw = w > 1e-6 * w.max()
    d = S - S[0]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j])))) for j in range(S.shape[1])])
    if not all((Tk[j] - S[0, j]) % steps[j] == 0 for j in range(len(steps))): return {"feasible": False, "why": "target off lattice"}
    pr = V.T @ (Tk - mu)
    if (~keepw).any() and np.abs(pr[~keepw]).max() > 1e-6 * max(1, np.abs(pr).max()): return {"feasible": False, "why": "null direction"}
    import math
    log2P = (np.log(steps).sum() - 0.5 * keepw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[keepw]).sum()
             - 0.5 * (pr[keepw] ** 2 / w[keepw]).sum()) / math.log(2)
    nbits = -log2P
    return {"v": v, "b": b, "n": n, "p": p, "free_bits": 12 * n, "constraint_bits": round(nbits, 1),
            "log2E": round(12 * n - nbits, 1), "eqs": int(keepw.sum())}

if __name__ == "__main__" and sys.argv[1] in ("count", "estimate"):
    v, b, n, p = map(int, sys.argv[2:6])
    out = gauss_estimate(v, b, n, p)
    if sys.argv[1] == "count":
        import math; c = exact_count(v, b, n, p); out["exact_count"] = c; out["exact_log2"] = round(math.log2(c), 1) if c else None
    print(json.dumps(out), flush=True)


# ---------------- row-sum-conditioned estimate (as nbits.py does) ----------------
def block_rowsum(lab, p):
    a = np.array([LABELS[i][0] for i in range(8)]); w = np.array([0 if LABELS[i][1] else LABELS[i][2] for i in range(8)])
    return a[lab].sum(1) + (p - 1) * w[lab].sum(1)

def count_rowsum(n, p):
    """exact number of labelings per block row sum (DP over cosets)"""
    from collections import Counter
    per = Counter(LABELS[i][0] + (p - 1) * (0 if LABELS[i][1] else LABELS[i][2]) for i in range(8))
    D = Counter({0: 1})
    for _ in range(n):
        E = Counter()
        for x, cx in D.items():
            for y, cy in per.items(): E[x + y] += cx * cy
        D = E
    return D

def sample_rowsum(n, p, a, N, rng, tries=200):
    out = []
    for _ in range(tries):
        lab = rng.integers(0, 8, (max(N, 20000), n)); lab = lab[block_rowsum(lab, p) == a]
        out.append(lab)
        if sum(len(x) for x in out) >= N: break
    lab = np.concatenate(out)[:N]
    return lab if len(lab) >= min(N, 2000) else None

def gauss_estimate2(v, b, n, p, N=60000, seed=1, maxprof=8):
    import math, itertools as it
    C, nqr = corr_table(p); m = n * p; rng = np.random.default_rng(seed)
    tgt = np.full((n, 3), -4 * b, np.int64); tgt[0, 0] = 4 * m; T = tgt.reshape(-1)
    cnt = count_rowsum(n, p); A = sorted(a for a in cnt if a >= 0)
    if b == 0:
        profs = [q for q in it.combinations_with_replacement(A, 4) if sum(x * x for x in q) == 4 * v]
    else:
        profs = [(2, 0, 0, 0)] if (2 in cnt and 0 in cnt) else []
    if not profs: return {"v": v, "b": b, "n": n, "p": p, "feasible": False, "why": "row sums"}
    profs = sorted(profs, key=lambda q: -sum(math.log2(cnt[x]) for x in q))[:maxprof]
    stats = {}
    def st(a):
        if a not in stats:
            lab = sample_rowsum(n, p, a, N, rng)
            if lab is None: stats[a] = None
            else:
                R = reduced_paf(lab, n, C).reshape(len(lab), -1); stats[a] = (R.mean(0), np.cov(R.T), R)
        return stats[a]
    best = None; tot = 0.0
    for q in profs:
        S = [st(a) for a in q]
        if any(s is None for s in S): continue
        perms = len(set(it.permutations(q))); signs = 2 ** sum(1 for x in q if x != 0)
        mu = sum(s[0] for s in S); Cv = sum(s[1] for s in S)
        keep = [j for j in range(len(T)) if j != 0]
        # drop exactly duplicated coordinates (PAF(h,c) = PAF(-h,c'))
        Rcat = sum(s[2][:2000] for s in S)
        cols, seen = [], {}
        for j in keep:
            key = Rcat[:, j].tobytes()
            if key in seen: continue
            seen[key] = j; cols.append(j)
        mu_ = mu[cols]; C_ = Cv[np.ix_(cols, cols)]; T_ = T[cols]
        w, V = np.linalg.eigh(C_); kw = w > 1e-6 * w.max(); pr = V.T @ (T_ - mu_)
        if (~kw).any() and np.abs(pr[~kw]).max() > 1e-6 * max(1, np.abs(pr).max()): continue
        d = Rcat[:, cols] - Rcat[0, cols]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j])))) for j in range(len(cols))])
        lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum() - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
        lE = sum(math.log2(cnt[x]) for x in q) + lp + math.log2(perms * signs)
        tot += 2 ** lE
        if best is None or lE > best[1]: best = (list(q), round(lE, 1), round(-lp, 1))
    if best is None: return {"v": v, "b": b, "n": n, "p": p, "feasible": False, "why": "no admissible profile"}
    return {"v": v, "b": b, "n": n, "p": p, "free_bits": 12 * n, "log2E_total": round(math.log2(tot), 1) if tot > 0 else None,
            "best_profile": best[0], "best_profile_log2E": best[1], "best_profile_constraint_bits": best[2]}

if __name__ == "__main__" and sys.argv[1] in ("estimate2", "count2"):
    v, b, n, p = map(int, sys.argv[2:6]); out = gauss_estimate2(v, b, n, p)
    if sys.argv[1] == "count2":
        import math; c = exact_count(v, b, n, p); out["exact_log2"] = round(math.log2(c), 1) if c else None
    print(json.dumps(out), flush=True)


# ---------------- exact screen: row sums + within-coset (h = 0) equations ----------------
def exact_h0_screen(v, b, n, p):
    """Per-coset additive invariants: row sum r(L), within-coset PAF e(L, QR), e(L, NQR) (lag (0, s), s != 0).
    Block value = sum over its n cosets; DP over cosets gives the reachable set per block; then require a
    row-sum profile and sum_k e_k = target at (0, QR) and (0, NQR).  Exact necessary condition."""
    import itertools as it
    C, nqr = corr_table(p)
    r = [LABELS[L][0] + (p - 1) * (0 if LABELS[L][1] else LABELS[L][2]) for L in range(8)]
    e = [[int(sum(LC[L, i] * LC[L, j] * C[i, j, c] for i in range(3) for j in range(3))) for c in (1, 2)] for L in range(8)]
    S = {(0, 0, 0)}
    for _ in range(n):
        S = {(x[0] + r[L], x[1] + e[L][0], x[2] + e[L][1]) for x in S for L in range(8)}
    tgtQ = tgtN = -4 * b
    by_r = {}
    for rs, eq, en in S: by_r.setdefault(rs, set()).add((eq, en))
    if b == 0:
        rows = sorted(by_r)
        profs = [q for q in it.product([x for x in rows if x >= 0], repeat=4) if sum(x * x for x in q) == 4 * v and list(q) == sorted(q)]
        cand = [[s * x for s, x in zip(signs, q)] for q in profs for signs in it.product((1, -1), repeat=4)]
    else:
        cand = [[a, 0, 0, 0] for a in (2, -2)]
    for q in cand:
        sets = [by_r.get(a) for a in q]
        if any(s is None for s in sets): continue
        s12 = {(x[0] + y[0], x[1] + y[1]) for x in sets[0] for y in sets[1]}
        s34 = {(x[0] + y[0], x[1] + y[1]) for x in sets[2] for y in sets[3]}
        if any((tgtQ - x[0], tgtN - x[1]) in s34 for x in s12):
            return {"v": v, "b": b, "n": n, "p": p, "h0_screen": True, "profile": [int(x) for x in q]}
    return {"v": v, "b": b, "n": n, "p": p, "h0_screen": False}

if __name__ == "__main__" and sys.argv[1] == "h0":
    v, b, n, p = map(int, sys.argv[2:6]); print(json.dumps(exact_h0_screen(v, b, n, p)))


# ---------------- fast exact h0 screen (aggregate form, 2-D boolean DP + convolution) ----------------
def exact_h0_screen_fast(v, b, n, p):
    """Same condition as exact_h0_screen, computed via aggregates.
    Per coset: twisted (alpha, beta): r = alpha, s = -1, d = alpha*beta*[p = 1 mod 4]
               constant (alpha, gamma): r = alpha + (p-1)*gamma, s = p - 2 + 2*alpha*gamma, d = 0
    (s = within-coset PAF averaged over QR/NQR, d = their half-difference / 2).  Target: sum_k s_k = -4b, sum_k d_k = 0,
    row sums: profile with sum a^2 = 4v (b=0) or (+-2, 0, 0, 0) (b=1)."""
    import itertools as it
    from scipy.signal import fftconvolve
    p1 = (p % 4 == 1)
    R = n * p; S0 = n; Smax = n * p                       # r in [-R, R], s in [-n, n*p]
    # constant part: G[Nc] = bool array over (r, s) for Nc constant cosets
    const_steps = [(a + (p - 1) * g, p - 2 + 2 * a * g) for a in (1, -1) for g in (1, -1)]
    shape = (2 * R + 1, Smax + S0 + 1)
    G = np.zeros(shape, bool); G[R, S0] = True            # zero cosets: r=0, s=0
    consts = [G]
    for _ in range(n):
        H = np.zeros(shape, bool)
        for dr, ds in const_steps:
            src = consts[-1]
            rs = slice(max(0, -dr), shape[0] - max(0, dr)); rd = slice(max(0, dr), shape[0] - max(0, -dr))
            ss = slice(max(0, -ds), shape[1] - max(0, ds)); sd = slice(max(0, ds), shape[1] - max(0, -ds))
            H[rd, sd] |= src[rs, ss]
        consts.append(H)
    # block reachable set as dict r -> bool array over (s, d)  (d in [-n, n])
    D0 = n
    block = {}
    for Nc in range(n + 1):
        Nt = n - Nc
        C = consts[Nc]
        # twisted part: r_t, d with parity Nt, |r_t| <= Nt, |d| <= Nt (d = 0 if not p1); s contribution -Nt
        rts = range(-Nt, Nt + 1, 2)
        ds = list(range(-Nt, Nt + 1, 2)) if p1 else [0]
        if not p1 and Nt % 2 == 1 and False: pass
        for rt in rts:
            rows = np.nonzero(C.any(1))[0]
            for ri in rows:
                rc = ri - R; rtot = rc + rt
                svec = C[ri]                                  # over s index (s = idx - S0)
                arr = block.setdefault(rtot, np.zeros((Smax + S0 + 1, 2 * D0 + 1), bool))
                sidx = np.nonzero(svec)[0] - Nt               # shift by -Nt
                sidx = sidx[sidx >= 0]
                for d in ds:
                    if (not p1) or True: arr[sidx, d + D0] = True
    def prof_list():
        rows = sorted(block)
        if b == 0:
            pos = [x for x in rows if x >= 0 and -x in block]
            return [list(q) for q in it.combinations_with_replacement(pos, 4) if sum(x * x for x in q) == 4 * v]
        return [[2, 0, 0, 0]] if 2 in block and 0 in block else []
    tS = -4 * b
    for q in prof_list():
        conv = None
        for a in q:
            A = block[a].astype(float)
            conv = A if conv is None else (fftconvolve(conv, A) > 0.5).astype(float)
        # conv index: s_total + 4*S0 ; d_total + 4*D0
        si, di = tS + 4 * S0, 4 * D0
        if 0 <= si < conv.shape[0] and conv[si, di] > 0.5:
            return {"v": v, "b": b, "n": n, "p": p, "h0_screen": True, "profile": [int(x) for x in q]}
    return {"v": v, "b": b, "n": n, "p": p, "h0_screen": False}

if __name__ == "__main__" and sys.argv[1] == "h0fast":
    v, b, n, p = map(int, sys.argv[2:6]); print(json.dumps(exact_h0_screen_fast(v, b, n, p)))


# ---------------- multiplier-symmetric slices: labels constant on orbits of g -> u*g on G' ----------------
# u in (Z/m)^*, u mod p a quadratic residue (so chi(u z) = chi(z) and the ansatz is preserved):
#   x(u g, u z) = x(g, z)  <=>  label(u g mod n) = label(g).
def mult_orbits(n, un):
    seen, orbs = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * un) % n
        seen.update(o); orbs.append(o)
    return orbs

def expand(lab_orb, orbs, n):
    lab = np.zeros((lab_orb.shape[0], n), np.int64)
    for j, o in enumerate(orbs): lab[:, o] = lab_orb[:, [j]]
    return lab

def count_rowsum_orb(orbs, p):
    from collections import Counter
    D = Counter({0: 1})
    for o in orbs:
        per = Counter(len(o) * (LABELS[i][0] + (p - 1) * (0 if LABELS[i][1] else LABELS[i][2])) for i in range(8))
        E = Counter()
        for x, cx in D.items():
            for y, cy in per.items(): E[x + y] += cx * cy
        D = E
    return D

def gauss_estimate_sym(v, b, n, p, un, N=60000, seed=1):
    import math, itertools as it
    orbs = mult_orbits(n, un); k = len(orbs)
    C, nqr = corr_table(p); m = n * p; rng = np.random.default_rng(seed)
    tgt = np.full((n, 3), -4 * b, np.int64); tgt[0, 0] = 4 * m; T = tgt.reshape(-1)
    cnt = count_rowsum_orb(orbs, p)
    if b == 1:
        profs = [(2, 0, 0, 0)] if (2 in cnt and 0 in cnt) else []
    else:
        A = sorted(a for a in cnt if a >= 0)
        profs = [q for q in it.combinations_with_replacement(A, 4) if sum(x * x for x in q) == 4 * v]
    if not profs: return {"v": v, "n": n, "p": p, "u": un, "orbits": k, "feasible": False, "why": "row sums"}
    sizes = np.array([len(o) for o in orbs])
    rsum = lambda lo: (np.array([LABELS[i][0] for i in range(8)])[lo] * sizes).sum(1) + (p - 1) * (np.array([0 if LABELS[i][1] else LABELS[i][2] for i in range(8)])[lo] * sizes).sum(1)
    def samp(a):
        out = []
        for _ in range(400):
            lo = rng.integers(0, 8, (40000, k)); lo = lo[rsum(lo) == a]; out.append(lo)
            if sum(len(x) for x in out) >= N: break
        lo = np.concatenate(out)[:N]
        return None if len(lo) < 1000 else reduced_paf(expand(lo, orbs, n), n, C).reshape(len(lo), -1)
    best = None; tot = 0.0; cache = {}
    for q in profs[:8]:
        Rs = [cache.setdefault(a, samp(a)) for a in q]
        if any(r is None for r in Rs): continue
        mu = sum(r.mean(0) for r in Rs); Cv = sum(np.cov(r.T) for r in Rs)
        Rc = sum(r[:1000] for r in Rs)
        cols, seen = [], {}
        for j in range(1, len(T)):
            key = Rc[:, j].tobytes()
            if key in seen: continue
            seen[key] = j; cols.append(j)
        mu_ = mu[cols]; C_ = Cv[np.ix_(cols, cols)]; T_ = T[cols]
        w, V = np.linalg.eigh(C_); kw = w > 1e-6 * w.max(); pr = V.T @ (T_ - mu_)
        if (~kw).any() and np.abs(pr[~kw]).max() > 1e-6 * max(1, np.abs(pr).max()): continue
        d = Rc[:, cols] - Rc[0, cols]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j])))) for j in range(len(cols))])
        lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum() - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
        perms = len(set(it.permutations(q))); signs = 2 ** sum(1 for x in q if x != 0)
        lE = sum(math.log2(cnt[x]) for x in q) + lp + math.log2(perms * signs)
        tot += 2 ** lE
        if best is None or lE > best[1]: best = (list(q), round(lE, 1), round(-lp, 1), int(kw.sum()))
    if best is None: return {"v": v, "n": n, "p": p, "u": un, "orbits": k, "feasible": False, "why": "no admissible profile"}
    return {"v": v, "b": b, "n": n, "p": p, "u_mod_n": un, "orbits": k, "bits_per_block": round(math.log2(8) * k, 1),
            "log2E": round(math.log2(tot), 1), "constraint_bits": best[2], "eqs": best[3], "profile": best[0]}

def sym_options(n, p):
    """multipliers u mod n' (units, order d > 1) that lift to u in (Z/m)^* with u mod p a QR: always possible (take u = 1 mod p)."""
    from math import gcd
    out = {}
    for un in range(2, n):
        if gcd(un, n) != 1: continue
        orbs = mult_orbits(n, un); out.setdefault(len(orbs), un)
    return sorted(out.items())                                    # (number of orbits, a representative u)

if __name__ == "__main__" and sys.argv[1] == "estsym":
    v, b, n, p = map(int, sys.argv[2:6])
    for k, un in sym_options(n, p):
        if k < 4: continue
        print(json.dumps(gauss_estimate_sym(v, b, n, p, un)), flush=True)
