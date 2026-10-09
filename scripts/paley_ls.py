"""Tabu local search on the REDUCED Paley-coset equations (3 per h' in G'), bordered b >= 3 or b in {0, 1}.

State: label (alpha, type, sign) per block k and coset g' in G' (paley_ansatz.LABELS); x(g', 0) = alpha,
x(g', z) = sign*chi(z) (twisted) or sign (constant) for z != 0.  The reduced PAF R[h', c] (c = lag class 0, QR, NQR)
must equal the target: 4m at (0, 0); -4b on H' x {0, QR, NQR} \ (0, 0); +4 elsewhere (b >= 3: H' of index c in G');
b in {0, 1}: -4b everywhere off (0, 0).  Objective = sum over ALL lags of G (weights 1, (p-1)/2, (p-1)/2) of
(R - T)^2.  Moves: change one label; all 7*n'*4 moves scored exactly per step (vectorised delta).
A zero-cost state is expanded to four +-1 sequences on G' x Z_p and verified by FFT.

usage: paley_ls.py v b p mods... [--Hgen e1,e2,...] [--seed s] [--steps N] [--restarts R] [--sa 1 [--T0 t]]
  --sa 1: simulated annealing (numba Metropolis, geometric cooling T0 -> 0.5 over N proposals) instead of tabu
  Hgen: generators of H' as mixed-radix element coordinates "a.b.c" (default for b >= 3: the index-c subgroup c*G'
  when G' is cyclic).  Example (359 control): paley_ls.py 359 7 11 2 2 8 --Hgen 0.0.2"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_ansatz import legendre, corr_table, LABELS, LC
from paley_bge3 import AG

def subgroup(G, gens):
    S = {0}; frontier = [0]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = int(G.add[x, g])
            if y not in S: S.add(y); frontier.append(y)
    return S

class Prob:
    def __init__(s, v, b, p, mods, Hgen=None):
        s.v, s.b, s.p = v, b, p; s.G = AG(mods); n = s.n = s.G.n; s.m = v - b
        assert n * p == s.m, (n, p, s.m)
        s.C, _ = corr_table(p); s.C = s.C.astype(np.int64)
        T = np.full((n, 3), -4 * b, np.int64)
        if b >= 3:
            c = b + 1
            if Hgen is None:
                assert len(mods) == 1; H = {(c * j) % n for j in range(n // c)}
            else:
                H = subgroup(s.G, [int(s.G.idx(np.array(g))) for g in Hgen])
            assert len(H) * c == n, (len(H), c, n)
            T[:] = 4; T[sorted(H)] = -4 * b
            s.H = sorted(H)
        T[0, 0] = 4 * s.m; s.T = T
        s.w = np.array([1, (p - 1) // 2, (p - 1) // 2], np.int64)
        s.add = s.G.add; s.sub = np.array([[int(s.G.add[g, s.neg(h)]) for h in range(n)] for g in range(n)])  # g - h

    def neg(s, h):
        return int(s.G.idx(-s.G.co[h]))

    def block_R(s, A):
        """A: (n, 3) -> (n, 3) reduced PAF of one block"""
        R = np.zeros((s.n, 3), np.int64)
        for h in range(s.n): R[h] = np.einsum('gi,gj,ijc->c', A, A[s.add[:, h]], s.C)
        return R

    def deltas(s, A):
        """dR[g, L, h, c] for setting coset g to label L (one block)"""
        D = LC[None, :, :] - A[:, None, :]                                    # (n, 8, 3)
        M1 = A[s.add]                                                          # M1[g, h] = A[g + h]
        M2 = A[s.sub]                                                          # M2[g, h] = A[g - h]
        dR = np.einsum('gli,ghj,ijc->glhc', D, M1, s.C) + np.einsum('ghi,glj,ijc->glhc', M2, D, s.C)
        dR[:, :, 0, :] += np.einsum('gli,glj,ijc->glc', D, D, s.C)
        return dR

def full_seq(P, lab_k):
    chi = legendre(P.p); n, p = P.n, P.p
    x = np.zeros((n, p), np.int64)
    for g in range(n):
        a, t, sg = LABELS[lab_k[g]]
        x[g, 0] = a; x[g, 1:] = sg * chi[1:] if t else sg
    return x                                                                   # shape (n', p); group G' x Z_p

def verify(P, lab):
    mods = tuple(P.G.mods) + (P.p,)
    X = [full_seq(P, lab[k]).reshape(mods) for k in range(4)]
    S = sum(np.rint(np.fft.ifftn(np.abs(np.fft.fftn(x)) ** 2).real).astype(np.int64) for x in X).reshape(P.n, P.p)
    Tf = np.zeros((P.n, P.p), np.int64); chi = legendre(P.p)
    for z in range(P.p): Tf[:, z] = P.T[:, 0] if z == 0 else P.T[:, 1 if chi[z] == 1 else 2]
    return bool((S == Tf).all()), X

import numba

@numba.njit(cache=True)
def _best_move(A, E, add, sub, C, w, LCa, tl, it, best, rnd, allow):
    """A: (4, n, 3) components, E: (n, 3) residual.  Returns (newcost_delta, k, g, L) of the best non-tabu move
    (aspiration if it beats best).  Ties broken by the random key rnd (4, n, 8)."""
    K, n = A.shape[0], A.shape[1]
    bestd = 1 << 62; bk = -1; bg = -1; bl = -1; bkey = 2.0
    dR = np.zeros((n, 3), np.int64)
    for k in range(K):
        for g in range(n):
            for L in range(8):
                if not allow[k, g, L]: continue
                d0 = LCa[L, 0] - A[k, g, 0]; d1 = LCa[L, 1] - A[k, g, 1]; d2 = LCa[L, 2] - A[k, g, 2]
                if d0 == 0 and d1 == 0 and d2 == 0: continue
                tot = 0
                for h in range(n):
                    a = add[g, h]; b_ = sub[g, h]
                    for c in range(3):
                        val = 0
                        for j in range(3):
                            val += (d0 * C[0, j, c] + d1 * C[1, j, c] + d2 * C[2, j, c]) * A[k, a, j]
                            val += A[k, b_, j] * (d0 * C[j, 0, c] + d1 * C[j, 1, c] + d2 * C[j, 2, c])
                        if h == 0:
                            val += d0 * (d0 * C[0, 0, c] + d1 * C[0, 1, c] + d2 * C[0, 2, c]) \
                                 + d1 * (d0 * C[1, 0, c] + d1 * C[1, 1, c] + d2 * C[1, 2, c]) \
                                 + d2 * (d0 * C[2, 0, c] + d1 * C[2, 1, c] + d2 * C[2, 2, c])
                        tot += w[c] * (2 * E[h, c] * val + val * val)
                if tl[k, g] > it and tot + best[1] >= best[0]: continue          # tabu unless it beats the best
                key = rnd[k, g, L]
                if tot < bestd or (tot == bestd and key < bkey):
                    bestd = tot; bk = k; bg = g; bl = L; bkey = key
    return bestd, bk, bg, bl

def search(P, seed=0, steps=200000, tabu=None, log_every=20000, lab0=None, allow=None):
    """lab0: starting labels (default random); allow: (4, n', 8) bool mask of permitted labels per coset"""
    rng = np.random.default_rng(seed); n = P.n
    tabu = tabu or max(5, n // 2)
    lab = rng.integers(0, 8, (4, n)) if lab0 is None else np.array(lab0)
    allow = np.ones((4, n, 8), np.bool_) if allow is None else allow
    A = LC[lab].astype(np.int64)
    E = sum(P.block_R(A[k]) for k in range(4)) - P.T
    cost = int(((E ** 2) * P.w).sum()); best = cost
    tl = np.zeros((4, n), np.int64) - 1
    LCa = LC.astype(np.int64); C = P.C.astype(np.int64); w = P.w.astype(np.int64)
    add = P.add.astype(np.int64); sub = P.sub.astype(np.int64)
    t0 = time.time()
    for it in range(steps):
        dc, k, g, L = _best_move(A, E, add, sub, C, w, LCa, tl, it, np.array([best, cost], np.int64), rng.random((4, n, 8)), allow)
        if k < 0: tl[:] = -1; continue
        Rold = P.block_R(A[k]); lab[k, g] = L; A[k, g] = LCa[L]
        E = E + P.block_R(A[k]) - Rold; cost += int(dc)
        tl[k, g] = it + tabu + rng.integers(0, tabu)
        if cost < best: best = cost
        if cost == 0:
            ok, X = verify(P, lab)
            return dict(found=True, verified=ok, steps=it + 1, secs=round(time.time() - t0, 1), lab=lab.tolist())
        if log_every and (it + 1) % log_every == 0:
            assert cost == int(((E ** 2) * P.w).sum())
            print(json.dumps(dict(it=it + 1, cost=cost, best=best, secs=round(time.time() - t0, 1))), flush=True)
    return dict(found=False, best=best, steps=steps, secs=round(time.time() - t0, 1))

@numba.njit(cache=True)
def _anneal(A, lab, E, add, sub, C, w, LCa, steps, T0, T1, seed, allow):
    """Metropolis simulated annealing on the exact reduced objective sum_{h,c} w_c E[h,c]^2; single-label moves.
    A (4, n, 3) components and lab (4, n) labels are updated in place; returns the final cost."""
    np.random.seed(seed)
    K, n = lab.shape
    cost = 0
    for h in range(n):
        for c in range(3): cost += w[c] * E[h, c] * E[h, c]
    dv = np.zeros((n, 3), np.int64)
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        k = np.random.randint(K); g = np.random.randint(n); L = np.random.randint(8)
        if L == lab[k, g] or not allow[k, g, L]: continue
        d0 = LCa[L, 0] - A[k, g, 0]; d1 = LCa[L, 1] - A[k, g, 1]; d2 = LCa[L, 2] - A[k, g, 2]
        tot = 0
        for h in range(n):
            a = add[g, h]; b_ = sub[g, h]
            for c in range(3):
                val = 0
                for j in range(3):
                    val += (d0 * C[0, j, c] + d1 * C[1, j, c] + d2 * C[2, j, c]) * A[k, a, j]
                    val += A[k, b_, j] * (d0 * C[j, 0, c] + d1 * C[j, 1, c] + d2 * C[j, 2, c])
                if h == 0:
                    val += d0 * (d0 * C[0, 0, c] + d1 * C[0, 1, c] + d2 * C[0, 2, c]) \
                         + d1 * (d0 * C[1, 0, c] + d1 * C[1, 1, c] + d2 * C[1, 2, c]) \
                         + d2 * (d0 * C[2, 0, c] + d1 * C[2, 1, c] + d2 * C[2, 2, c])
                dv[h, c] = val
                tot += w[c] * (2 * E[h, c] * val + val * val)
        if tot <= 0 or np.random.random() < np.exp(-tot / T):
            A[k, g, 0] += d0; A[k, g, 1] += d1; A[k, g, 2] += d2; lab[k, g] = L
            for h in range(n):
                for c in range(3): E[h, c] += dv[h, c]
            cost += tot
            if cost == 0: return 0
    return cost

def anneal(P, seed=0, steps=10**8, T0=None, T1=0.5, lab0=None, allow=None):
    rng = np.random.default_rng(seed); n = P.n
    lab = rng.integers(0, 8, (4, n)) if lab0 is None else np.array(lab0)
    allow = np.ones((4, n, 8), np.bool_) if allow is None else allow
    A = LC[lab].astype(np.int64)
    E = (sum(P.block_R(A[k]) for k in range(4)) - P.T).astype(np.int64)
    T0 = T0 or float(P.p * P.p * 8)
    t0 = time.time()
    c = _anneal(A, lab, E, P.add.astype(np.int64), P.sub.astype(np.int64), P.C.astype(np.int64), P.w.astype(np.int64),
                LC.astype(np.int64), steps, T0, T1, int(rng.integers(1 << 30)), allow)
    out = dict(found=c == 0, final_cost=int(c), steps=steps, secs=round(time.time() - t0, 1))
    if c == 0:
        ok, X = verify(P, lab); out.update(verified=ok, lab=lab.tolist())
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    opt = {"--seed": "0", "--steps": "200000", "--restarts": "1", "--Hgen": None, "--log": "20000", "--sa": "0", "--T0": "0"}
    pos = []
    i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, b, p = pos[:3]; mods = tuple(pos[3:])
    Hgen = None if opt["--Hgen"] is None else [tuple(map(int, e.split('.'))) for e in opt["--Hgen"].split(',')]
    P = Prob(v, b, p, mods, Hgen)
    for r in range(int(opt["--restarts"])):
        if int(opt["--sa"]):
            res = anneal(P, seed=int(opt["--seed"]) + r, steps=int(float(opt["--steps"])), T0=float(opt["--T0"]) or None)
        else:
            res = search(P, seed=int(opt["--seed"]) + r, steps=int(opt["--steps"]), log_every=int(opt["--log"]))
        out = {k: res[k] for k in res if k != "lab"}
        print(json.dumps(dict(v=v, b=b, p=p, mods=list(mods), seed=int(opt["--seed"]) + r, **out)), flush=True)
        if res.get("found") and res.get("verified"):
            fn = f"logs/paley_ls_v{v}_b{b}_p{p}_{'x'.join(map(str, mods))}_s{int(opt['--seed']) + r}.json"
            json.dump(dict(v=v, b=b, p=p, mods=list(mods), H=getattr(P, "H", None), lab=res["lab"]), open(fn, "w"))
            break
