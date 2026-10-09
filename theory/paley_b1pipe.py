r"""Exact two-stage pipeline for b = 1 Paley-coset cores over G' = Z_n (p = 3 mod 4), completed to H(4v) by the standard
bordered Goethals-Seidel array (row sums (+-2, 0, 0, 0)) and verified.  Positive control for report section 7.

Reductions (report sections 2, 4):
  * norm identity -> (K, N-, N+);   Lemma 4 -> s(h) = sum_k PAF(w_k)(h) in S(h) (small), sum_k t_k^2 = K + sum s;
  * row sums r = (+-2, 0, 0, 0): |r_k - (p-1) t_k| <= n  (t_k = sum w_k);
  * per-block translations and reflections preserve every PAF and X (X_k is preserved by translation; reflection is
    combined with u -> -u, the NQR multiplier), so each block's gadget is normalised up to translation/reflection;
  * stage 1: per gadget, all alpha with the right sums, exact 2+2 join on sum_k PAF(y_k) = T_J;
  * stage 2: u on twisted positions, exact 2+2 join on (sum PAF(beta) + p PAF(u), X) = (4v delta, 0).
usage: python3 theory/paley_b1pipe.py v p n [max_cores]"""
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

def run(v, p, n, max_cores=1, max_stage1=10 ** 9, seed=0):
    b = 1; m = n * p; assert m == v - 1
    t0 = time.time(); TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m + (p - 1) * (-4)
    KN = KN_options(v, b, p, n)
    log = open(os.path.join(ROOT, 'theory', 'logs', f'paley_b1pipe_v{v}.log'), 'w')
    def L(**kw): log.write(json.dumps(kw) + '\n'); log.flush()
    L(v=v, p=p, n=n, KN=KN)
    found = []; n_gadgets = 0; n_s1 = 0; n_s2_calls = 0
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
                    S01 = (lists[0][1][:, None, :] + lists[1][1][None, :, :]).reshape(-1, n)
                    H = defaultdict(list)
                    for i, r in enumerate(S01): H[r.tobytes()].append(i)
                    S23 = (TJ[None, None, :] - lists[2][1][:, None, :] - lists[3][1][None, :, :]).reshape(-1, n)
                    for j, r in enumerate(S23):
                        for i in H.get(r.tobytes(), []):
                            a0, a1 = divmod(i, len(lists[1][0])); a2, a3 = divmod(j, len(lists[3][0]))
                            Aall = np.array([lists[0][0][a0], lists[1][0][a1], lists[2][0][a2], lists[3][0][a3]])
                            n_s1 += 1
                            Us = stage2(Aall, W, p, v, n); n_s2_calls += 1
                            for U in Us:
                                rec = finish(Aall, U, W, v, p, n)
                                found.append(rec); L(CORE=True, K=K, Nminus=Nm, **{k_: rec[k_] for k_ in rec if k_ not in ('X',)}, secs=round(time.time() - t0, 1))
                                if len(found) >= max_cores:
                                    L(summary=True, gadgets=n_gadgets, stage1=n_s1, cores=len(found), secs=round(time.time() - t0, 1)); log.close()
                                    return found
                            if n_s1 % 50 == 0: L(progress=True, gadgets=n_gadgets, stage1=n_s1, secs=round(time.time() - t0, 1))
                            if n_s1 >= max_stage1: break
    L(summary=True, gadgets=n_gadgets, stage1=n_s1, cores=len(found), secs=round(time.time() - t0, 1)); log.close()
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
    """expand to G = Z_n x Z_p (CRT order: index i <-> (i mod n, i mod p) needs gcd(n,p)=1), verify core, build bordered GS"""
    chi = legendre(p); m = n * p
    X = np.zeros((4, m), np.int64)
    for k in range(4):
        for i in range(m):
            g, z = i % n, i % p
            X[k, i] = A[k, g] if z == 0 else (U[k, g] * chi[z] if W[k, g] == 0 else W[k, g])
    P = np.stack([np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X]).sum(0)
    core_ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
    rows = X.sum(1).tolist()
    Hm = bordered_gs(X)
    N = 4 * v
    had = bool(Hm is not None and np.array_equal(Hm @ Hm.T, N * np.eye(N, dtype=np.int64)))
    if had:
        with open(os.path.join(ROOT, 'data', f'core_v{v}_b1_paley.txt'), 'w') as fh:
            fh.write(f'v={v} N={N} b=1 m={m} cyclic Z{m} = Z{n} x Z{p} (CRT), Paley-coset ansatz p={p}, theory/paley_b1pipe.py\n')
            for nm, row in zip('ABCD', X): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')
    return dict(core_ok=core_ok, rowsums=rows, hadamard=had, N=N, X=X)

def bordered_gs(X):
    """standard bordered Goethals-Seidel (b = 1), order 4m + 4.  Core = GS array of circulants (R = back-identity):
       [[A, BR, CR, DR], [-BR, A, D^T R, -C^T R], [-CR, -D^T R, A, B^T R], [-DR, C^T R, -B^T R, A]].
    Border rows = [Cn | S (x) 1^T], border columns = [Cn ; T^T (x) 1]: orthogonality forces S S^T = T T^T = 4I and
    Cn = -S Sigma^T T / 4, Sigma[I, J] = row sum of core sub-block (I, J); we search the 768^2 pairs of H(4)'s."""
    m = X.shape[1]
    circ = lambda x: np.stack([np.roll(x, i) for i in range(m)])
    R = np.eye(m, dtype=np.int64)[::-1]
    A, B, C, D = [circ(x) for x in X]
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

if __name__ == "__main__":
    v, p, n = map(int, sys.argv[1:4]); mc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    res = run(v, p, n, mc)
    print(json.dumps([{k: r[k] for k in r if k != 'X'} for r in res]))
