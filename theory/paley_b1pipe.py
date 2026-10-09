r"""Exact two-stage pipeline for b = 1 Paley-coset cores over G' = Z_n (p = 3 mod 4), completed to H(4v) by the standard
bordered Goethals-Seidel array (row sums (+-2, 0, 0, 0)) and verified.  Positive control for report section 7.

Reductions (report sections 2, 4):
  * norm identity -> (K, N-, N+);   Lemma 4 -> s(h) = sum_k PAF(w_k)(h) in S(h) (small), sum_k t_k^2 = K + sum s;
  * row sums r = (+-2, 0, 0, 0): |r_k - (p-1) t_k| <= n  (t_k = sum w_k);
  * per-block translations and reflections preserve every PAF and X (X_k is preserved by translation; reflection is
    combined with u -> -u, the NQR multiplier), so each block's gadget is normalised up to translation/reflection;
  * stage 1: per gadget, all alpha with the right sums, exact 2+2 join on sum_k PAF(y_k) = T_J;
  * Lemma 6 (S2 - S1 at the real characters): the u row sums rho_k = sum u_k satisfy sum_k (r_k - p t_k)^2 + p sum rho_k^2
    = 4v with rho_k = |T_k| (mod 2); at psi = (-1)^g the same with (y^, w^)  -- filters gadgets and stage-1 solutions;
  * stage 2: u on twisted positions, exact 2+2 join on (sum PAF(beta) + p PAF(u), X) = (4v delta, 0).
usage: python3 theory/paley_b1pipe.py v p n [max_cores]
       python3 theory/paley_b1pipe.py verify data/core_v*_b1_paley.txt"""
import sys, os, json, time, itertools as it, math
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
import paley_ls as LSm
from paley_ansatz import LABELS, legendre
from paley_wscreen import KN_options, S_allowed

def pafc(x):                                  # cyclic PAF, rows
    x = np.atleast_2d(x); n = x.shape[1]
    return np.stack([(x * np.roll(x, -h, axis=1)).sum(1) for h in range(n)], 1)

def ccx(f, g):                                # CC(f,g)(h) = sum f(x) g(x+h), rows of g
    g = np.atleast_2d(g); n = g.shape[1]
    return np.stack([(f[None, :] * np.roll(g, -h, axis=1)).sum(1) for h in range(n)], 1)

def join2x2(P, T, seed=5):
    """all (i, j) with P0[a0] + P1[a1] + P2[a2] + P3[a3] == T, i = a0*|P1| + a1, j = a2*|P3| + a3 (linear uint64 hash,
    sort + searchsorted, exact re-check)"""
    wts = np.random.default_rng(seed).integers(1, 2 ** 62, P[0].shape[1], dtype=np.int64).astype(np.uint64)
    key = lambda M: (M.astype(np.uint64) * wts[None, :]).sum(1)
    k = [key(M) for M in P]
    K01 = (k[0][:, None] + k[1][None, :]).reshape(-1); order = np.argsort(K01, kind='stable'); Ks = K01[order]
    kT = (np.asarray(T).astype(np.uint64) * wts).sum()
    n1, n3 = len(P[1]), len(P[3])
    step = max(1, (1 << 24) // max(n3, 1))                 # chunk the (2,3) side: at most ~2^24 pair keys at a time
    for c0 in range(0, len(P[2]), step):
        need = (kT - k[2][c0:c0 + step, None] - k[3][None, :]).reshape(-1)
        lo = np.searchsorted(Ks, need, 'left'); hi = np.searchsorted(Ks, need, 'right')
        for jj in np.nonzero(hi > lo)[0]:
            a2, a3 = divmod(int(jj), n3); a2 += c0
            for i in order[lo[jj]:hi[jj]]:
                a0, a1 = divmod(int(i), n1)
                if (P[0][a0] + P[1][a1] + P[2][a2] + P[3][a3] == T).all(): yield int(i), a2 * n3 + a3

def block_shapes(n, kmax):
    """signed point sets (positions, w-signs) up to translation and reflection, |set| <= kmax"""
    out = {0: [((), ())]}
    for k in range(1, kmax + 1):
        seen = set(); lst = []
        for rest in it.combinations(range(1, n), k - 1):
            pos = (0,) + rest
            for sg in it.product((1, -1), repeat=k):
                if sg[0] != 1 and False: pass
                w = np.zeros(n, np.int64); w[list(pos)] = sg
                cands = []
                for rf in (False, True):
                    ww = w[(-np.arange(n)) % n] if rf else w
                    for tr in range(n): cands.append(tuple(np.roll(ww, -tr)))
                key = min(cands)
                if key in seen: continue
                seen.add(key); lst.append(key)
        out[k] = [(tuple(np.nonzero(np.array(kk))[0]), tuple(np.array(kk)[np.nonzero(np.array(kk))[0]])) for kk in lst]
    return out

def sq4(Q, sizes):
    """exists rho in Z^4 with sum rho_k^2 = Q, rho_k = sizes_k (mod 2), |rho_k| <= sizes_k  (Lemma 6)"""
    if Q < 0: return False
    rng = [[r for r in range(sz % 2, sz + 1, 2) if r * r <= Q] for sz in sizes]
    for r0 in rng[0]:
        for r1 in rng[1]:
            for r2 in rng[2]:
                q = Q - r0 * r0 - r1 * r1 - r2 * r2
                if q < 0: break
                r3 = math.isqrt(q)
                if r3 * r3 == q and r3 % 2 == sizes[3] % 2 and r3 <= sizes[3]: return True
    return False

def lemma6_real_chars(A, W, p, n):
    """S2 - S1 at the real character psi(g) = (-1)^g (n even):  sum_k u^_k(psi)^2 = 2 sum_k y^_k w^_k - p sum_k w^_k^2,
    u^_k(psi) = |T_k| (mod 2), |u^_k(psi)| <= |T_k| (T_k = twisted cosets of block k)"""
    if n % 2: return True
    sg = 1 - 2 * (np.arange(n) % 2)
    Y = A + (p - 1) * W; yh = (Y * sg).sum(1); wh = (W * sg).sum(1)
    R = int(2 * (yh * wh).sum() - p * (wh * wh).sum())
    return sq4(R, [int((W[k] == 0).sum()) for k in range(4)])

def run(v, p, n, max_cores=1, max_stage1=10 ** 9, seed=0, record=None, lemma6=True, tag=''):
    b = 1; m = n * p; assert m == v - 1
    t0 = time.time(); TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m + (p - 1) * (-4)
    KN = KN_options(v, b, p, n)
    log = open(os.path.join(ROOT, 'theory', 'logs', f'paley_b1pipe_v{v}{tag}.log'), 'w')
    def L(**kw): log.write(json.dumps(kw) + '\n'); log.flush()
    L(v=v, p=p, n=n, KN=KN)
    found = []; n_gadgets = 0; n_s1 = 0; n_s2_calls = 0; n_l6 = 0; n_l6y = 0
    for (K, Nm, Np) in KN:
        Soff = S_allowed(-4 * p, p, n, K)
        shapes = block_shapes(n, min(K, 4))
        # t_k range from row sums: r in (2,0,0,0) (block 0 carries +2; others 0 -- block order fixed WLOG up to perm)
        rs = [2, 0, 0, 0]
        def tok(r, t): return abs(r - (p - 1) * t) <= n and (r - (p - 1) * t - n) % 2 == 0
        # enumerate size partitions of K over 4 blocks
        for sizes in it.product(range(0, min(K, 4) + 1), repeat=4):
            if sum(sizes) != K: continue
            for shp in it.product(*[shapes[sz] for sz in sizes]):
                W = np.zeros((4, n), np.int64)
                for k, (pos, sg) in enumerate(shp): W[k, list(pos)] = sg
                ts = W.sum(1)
                if not all(tok(rs[k], int(ts[k])) for k in range(4)): continue
                s = pafc(W).sum(0)
                if not all(int(s[h]) in Soff for h in range(1, n)): continue
                if int((ts ** 2).sum()) != K + int(s[1:].sum()): continue
                if lemma6:                                  # Lemma 6 at the trivial character: u row sums
                    num = 4 * v - sum((rs[k] - p * int(ts[k])) ** 2 for k in range(4))
                    if num % p or not sq4(num // p, [n - int((W[k] != 0).sum()) for k in range(4)]): n_l6 += 1; continue
                # J/J' assignment: choose which gadget positions are J' (alpha = -w): exactly Nm of them
                gpos = [(k, g) for k in range(4) for g in np.nonzero(W[k])[0]]
                for jm in it.combinations(range(len(gpos)), Nm):
                    n_gadgets += 1
                    Jp = np.zeros((4, n), bool)
                    for i in jm: Jp[gpos[i]] = True                          # J' positions
                    # stage 1: alpha lists per block with sum alpha_k = r_k - (p-1) t_k
                    lists = []
                    for k in range(4):
                        fixed = {int(g): (-int(W[k, g]) if Jp[k, g] else int(W[k, g])) for g in np.nonzero(W[k])[0]}
                        free = [g for g in range(n) if g not in fixed]
                        need = rs[k] - (p - 1) * int(ts[k]) - sum(fixed.values())
                        nm_ = (len(free) - need)
                        if nm_ % 2 or abs(need) > len(free): lists = None; break
                        A = []
                        for minus in it.combinations(free, nm_ // 2):
                            a = np.ones(n, np.int64)
                            for g, val in fixed.items(): a[g] = val
                            a[list(minus)] = -1; A.append(a)
                        A = np.array(A)
                        lists.append((A, pafc(A + (p - 1) * W[k][None, :])))
                    if lists is None: continue
                    if n_gadgets % 200 == 0: L(progress=True, gadgets=n_gadgets, K=K, stage1=n_s1, secs=round(time.time() - t0, 1))
                    for i, j in join2x2([l[1] for l in lists], TJ):
                            a0, a1 = divmod(i, len(lists[1][0])); a2, a3 = divmod(j, len(lists[3][0]))
                            Aall = np.array([lists[0][0][a0], lists[1][0][a1], lists[2][0][a2], lists[3][0][a3]])
                            n_s1 += 1
                            if lemma6 and not lemma6_real_chars(Aall, W, p, n):
                                n_l6y += 1
                                if record is not None: record.append((Aall.copy(), W.copy(), -1))
                                continue
                            Us = stage2(Aall, W, p, v, n); n_s2_calls += 1
                            if record is not None: record.append((Aall.copy(), W.copy(), len(Us)))
                            for U in Us:
                                rec = finish(Aall, U, W, v, p, n)
                                found.append(rec); L(CORE=True, K=K, Nminus=Nm, **{k_: rec[k_] for k_ in rec if k_ not in ('X',)}, secs=round(time.time() - t0, 1))
                                if len(found) >= max_cores:
                                    L(summary=True, gadgets=n_gadgets, stage1=n_s1, cores=len(found), lemma6_gadgets_cut=n_l6, lemma6_y_cut=n_l6y, secs=round(time.time() - t0, 1)); log.close()
                                    return found
                            if n_s1 % 50 == 0: L(progress=True, gadgets=n_gadgets, stage1=n_s1, secs=round(time.time() - t0, 1))
                            if n_s1 >= max_stage1: break
    L(summary=True, gadgets=n_gadgets, stage1=n_s1, cores=len(found), lemma6_gadgets_cut=n_l6, lemma6_y_cut=n_l6y, secs=round(time.time() - t0, 1)); log.close()
    return found

def stage2(A, W, p, v, n, max_sols=1):
    beta = A - W; T = np.concatenate([np.eye(1, n, 0, dtype=np.int64)[0] * 4 * v, np.zeros(n, np.int64)])
    lists = []
    for k in range(4):
        tw = np.nonzero(W[k] == 0)[0]; mm = len(tw); idx = np.arange(2 ** mm, dtype=np.int64)
        Ub = np.zeros((2 ** mm, n), np.int64); Ub[:, tw] = 1 - 2 * ((idx[:, None] >> np.arange(mm)) & 1)
        sig = np.concatenate([pafc(beta[k])[0][None, :] + p * pafc(Ub), ccx(beta[k], Ub) - np.stack([(Ub * np.roll(beta[k], -h)[None, :]).sum(1) for h in range(n)], 1)], 1)
        lists.append((Ub, sig))
    wts = np.random.default_rng(3).integers(1, 2 ** 62, 2 * n, dtype=np.int64).astype(np.uint64)
    key = lambda M: (M.astype(np.uint64) * wts[None, :]).sum(1)
    K01 = (key(lists[0][1])[:, None] + key(lists[1][1])[None, :]).reshape(-1); order = np.argsort(K01); Ks = K01[order]
    need = ((T.astype(np.uint64) * wts).sum() - key(lists[2][1])[:, None] - key(lists[3][1])[None, :]).reshape(-1)
    lo = np.searchsorted(Ks, need, 'left'); hi = np.searchsorted(Ks, need, 'right'); sols = []
    for jj in np.nonzero(hi > lo)[0]:
        c_, d = divmod(int(jj), len(lists[3][0]))
        for h in order[lo[jj]:hi[jj]]:
            a, b_ = divmod(int(h), len(lists[1][0]))
            if (lists[0][1][a] + lists[1][1][b_] + lists[2][1][c_] + lists[3][1][d] == T).all():
                sols.append(np.array([lists[0][0][a], lists[1][0][b_], lists[2][0][c_], lists[3][0][d]]))
                if len(sols) >= max_sols: return sols
    return sols

def finish(A, U, W, v, p, n):
    """expand to G = Z_n x Z_p (any n: element (g, z) <-> index g*p + z; cyclic or not), verify the core by a 2-D FFT,
    build the bordered GS array over G and verify H H^T = N I"""
    chi = legendre(p); m = n * p
    X = np.zeros((4, n, p), np.int64)
    for k in range(4):
        for g in range(n):
            for z in range(p):
                X[k, g, z] = A[k, g] if z == 0 else (U[k, g] * chi[z] if W[k, g] == 0 else W[k, g])
    P = np.stack([np.rint(np.fft.ifft2(np.abs(np.fft.fft2(x)) ** 2).real).astype(np.int64) for x in X]).sum(0).reshape(-1)
    core_ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
    rows = X.reshape(4, -1).sum(1).tolist()
    Hm = bordered_gs(X.reshape(4, -1), n, p)
    N = 4 * v
    had = bool(Hm is not None and np.array_equal(Hm @ Hm.T, N * np.eye(N, dtype=np.int64)))
    if had:
        grp = f'G = Z{n} x Z{p}' + (f' (cyclic, = Z{m})' if math.gcd(n, p) == 1 else ' (not cyclic)')
        with open(os.path.join(ROOT, 'data', f'core_v{v}_b1_paley.txt'), 'w') as fh:
            fh.write(f'v={v} N={N} b=1 m={m} n={n} p={p} {grp}; entry (g, z) at index g*{p} + z; Paley-coset ansatz, theory/paley_b1pipe.py\n')
            for nm, row in zip('ABCD', X.reshape(4, -1)): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')
    return dict(core_ok=core_ok, rowsums=rows, hadamard=had, N=N, X=X.reshape(4, -1))

def group_tables(n, p):
    """index g*p + z for (g, z) in Z_n x Z_p: SUB[a, b] = index of b - a, NEG[a] = index of -a"""
    g = np.arange(n * p) // p; z = np.arange(n * p) % p
    SUB = ((g[None, :] - g[:, None]) % n) * p + (z[None, :] - z[:, None]) % p
    NEG = ((-g) % n) * p + (-z) % p
    return SUB, NEG

def bordered_gs(X, n, p):
    """standard bordered Goethals-Seidel (b = 1), order 4m + 4, over the abelian group G = Z_n x Z_p.  Core = GS array of
    group-developed matrices M[a, b] = x(b - a) with R[a, b] = [a + b = 0] (valid for any abelian group):
       [[A, BR, CR, DR], [-BR, A, D^T R, -C^T R], [-CR, -D^T R, A, B^T R], [-DR, C^T R, -B^T R, A]].
    Border rows = [Cn | S (x) 1^T], border columns = [Cn ; T^T (x) 1]: orthogonality forces S S^T = T T^T = 4I and
    Cn = -S Sigma^T T / 4, Sigma[I, J] = row sum of core sub-block (I, J); we search the 768^2 pairs of H(4)'s."""
    m = X.shape[1]
    SUB, NEG = group_tables(n, p)
    R = np.zeros((m, m), np.int64); R[np.arange(m), NEG] = 1
    A, B, C, D = [x[SUB] for x in X]
    core = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    Sig = np.array([[core[I * m, J * m:(J + 1) * m].sum() for J in range(4)] for I in range(4)])
    H4s = [np.array(M).reshape(4, 4) for M in it.product((1, -1), repeat=16)]
    H4s = [M for M in H4s if np.array_equal(M @ M.T, 4 * np.eye(4, dtype=np.int64))]
    N = 4 * m + 4
    for S in H4s:
        for T in H4s:
            num = -(S @ Sig.T @ T)
            if (num % 4).any(): continue
            Cn = num // 4
            if not (np.abs(Cn) == 1).all(): continue
            Hm = np.zeros((N, N), np.int64); Hm[:4, :4] = Cn
            Hm[:4, 4:] = np.repeat(S, m, axis=1); Hm[4:, :4] = np.repeat(T.T, m, axis=0); Hm[4:, 4:] = core
            if np.array_equal(Hm @ Hm.T, N * np.eye(N, dtype=np.int64)): return Hm
    return None

def verify(path):
    """independent check of a data/core_v*_b1_paley.txt file: 2-D PAF sum = 4m at 0 and -4 elsewhere, row sums
    (+-2, 0, 0, 0) up to order, and the bordered GS matrix is Hadamard"""
    lines = open(path).read().split('\n')
    hd = dict(t.split('=') for t in lines[0].split(';')[0].split() if '=' in t and t.split('=')[0] in ('v', 'N', 'm', 'n', 'p'))
    v, N, m, n, p = (int(hd[k]) for k in ('v', 'N', 'm', 'n', 'p'))
    X = np.array([[1 if c == '+' else -1 for c in l.split('=')[1]] for l in lines[1:5]], np.int64)
    assert X.shape == (4, m) and m == n * p and N == 4 * v == 4 * m + 4
    P = np.stack([np.rint(np.fft.ifft2(np.abs(np.fft.fft2(x.reshape(n, p))) ** 2).real).astype(np.int64) for x in X]).sum(0).reshape(-1)
    core_ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
    rows = sorted(np.abs(X.sum(1)).tolist())
    Hm = bordered_gs(X, n, p)
    had = bool(Hm is not None and np.array_equal(Hm @ Hm.T, N * np.eye(N, dtype=np.int64)))
    return dict(file=os.path.basename(path), v=v, N=N, n=n, p=p, core_ok=core_ok, abs_rowsums=rows, hadamard=had)

if __name__ == "__main__":
    if sys.argv[1] == 'verify':
        for f in sys.argv[2:]: print(json.dumps(verify(f)))
        sys.exit()
    v, p, n = map(int, sys.argv[1:4]); mc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    res = run(v, p, n, mc)
    print(json.dumps([{k: r[k] for k in r if k != 'X'} for r in res]))
