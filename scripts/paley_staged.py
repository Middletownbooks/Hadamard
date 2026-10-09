"""Two-stage search for Paley-coset cores (p = 3 mod 4), exact decomposition of the reduced equations.

Per coset g' of Z_p: x = alpha*delta + u*chi + w*eps (u = +-1 twisted / 0 constant; w = +-1 constant / 0 twisted).
With R0 = PAF at (h, 0) and Y, X the symmetric / antisymmetric parts at (h, QR), (h, NQR) (R_QR = Y + X, R_NQR = Y - X):
   R0 = PAF_a + (p-1)(PAF_u + PAF_w)
   Y  = -PAF_u + (p-2) PAF_w + CC(a, w) + CC(w, a)
   X  = CC(a - w, u) - CC(u, a - w)                          (CC(f, g)(h) = sum_g f(g) g(g + h); sums over blocks)
and R0 + (p-1) Y = PAF(y) with y = a + (p-1) w = the Z_p-coset sums.  Hence, exactly:
   stage 1  (compression onto G'):   sum_k PAF(y_k) = T0 + (p-1) Ts     y in {+-1 (twisted), +-(p-2), +-p (constant)}
   stage 2  (twisted signs u | y):   sum_k PAF(u_k) = (p-2) PAF_w + CC(a,w) + CC(w,a) - Ts,   X = 0
Stage 1 is a small integer-valued GS-type problem on G'; stage 2 is a +-1 GS-type problem with a linear side condition.
Stage 2 is run with paley_ls's exact reduced objective restricted to twisted-sign moves.

usage: paley_staged.py v b p mods... [--Hgen a.b.c,...] [--seed s] [--s1steps N] [--s2steps N] [--tries T] [--s2tries R]"""
import sys, os, json, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_ansatz import LABELS, LC
import paley_ls as LSm

@numba.njit(cache=True)
def _s1_move(Y, E, add, sub, vals, tl, it, best, cost, rnd):
    K, n = Y.shape
    bestd = 1 << 62; bk = -1; bg = -1; bv = -1; bkey = 2.0
    for k in range(K):
        for g in range(n):
            for vi in range(vals.shape[0]):
                d = vals[vi] - Y[k, g]
                if d == 0: continue
                tot = 0
                for h in range(n):
                    val = d * (Y[k, add[g, h]] + Y[k, sub[g, h]])
                    if h == 0: val += d * d
                    tot += 2 * E[h] * val + val * val
                if tl[k, g] > it and cost + tot >= best: continue
                if tot < bestd or (tot == bestd and rnd[k, g, vi] < bkey):
                    bestd = tot; bk = k; bg = g; bv = vi; bkey = rnd[k, g, vi]
    return bestd, bk, bg, bv

def paf_g(P, f):
    return np.array([(f * f[P.add[:, h]]).sum() for h in range(P.n)])

def cc_g(P, f, g):
    return np.array([(f * g[P.add[:, h]]).sum() for h in range(P.n)])

@numba.njit(cache=True)
def _anneal_pos(Y, E, add, sub, V, nv, steps, T0, T1, seed, mask):
    """as _s1_anneal, but position (k, g) may only take the values V[k, g, :nv[k, g]]"""
    np.random.seed(seed)
    K, n = Y.shape
    cost = 0
    for h in range(n):
        if mask[h]: cost += E[h] * E[h]
    dv = np.zeros(n, np.int64)
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        k = np.random.randint(K); g = np.random.randint(n)
        if nv[k, g] < 2: continue
        d = V[k, g, np.random.randint(nv[k, g])] - Y[k, g]
        if d == 0: continue
        tot = 0
        for h in range(n):
            val = d * (Y[k, add[g, h]] + Y[k, sub[g, h]])
            if h == 0: val += d * d
            dv[h] = val
            if mask[h]: tot += 2 * E[h] * val + val * val
        if tot <= 0 or np.random.random() < np.exp(-tot / T):
            Y[k, g] += d
            for h in range(n): E[h] += dv[h]
            cost += tot
            if cost == 0: return 0
    return cost

@numba.njit(cache=True)
def _s1_anneal(Y, E, add, sub, vals, steps, T0, T1, seed, mask):
    """Metropolis on the cost sum_{h: mask[h]} E[h]^2 (E = sum_k PAF(y_k) - target); geometric cooling;
    returns final cost (E is kept up to date on every h, masked or not)"""
    np.random.seed(seed)
    K, n = Y.shape
    cost = 0
    for h in range(n):
        if mask[h]: cost += E[h] * E[h]
    dv = np.zeros(n, np.int64)
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        k = np.random.randint(K); g = np.random.randint(n); vi = np.random.randint(vals.shape[0])
        d = vals[vi] - Y[k, g]
        if d == 0: continue
        tot = 0
        for h in range(n):
            val = d * (Y[k, add[g, h]] + Y[k, sub[g, h]])
            if h == 0: val += d * d
            dv[h] = val
            if mask[h]: tot += 2 * E[h] * val + val * val
        # note: E[h] and E[-h] are the same equation; tot double-counts symmetrically, which only rescales T
        if tot <= 0 or np.random.random() < np.exp(-tot / T):
            Y[k, g] = vals[vi]
            for h in range(n): E[h] += dv[h]
            cost += tot
            if cost == 0: return 0
    return cost

def stage1(P, rng, steps, T0=None, T1=1.0):
    n, p = P.n, P.p
    T1v = P.T[:, 0] + (p - 1) * P.T[:, 1]
    vals = np.array([1, -1, p - 2, -(p - 2), p, -p], np.int64)
    Y = rng.choice(np.array([1, -1]), (4, n)).astype(np.int64)
    E = (sum(paf_g(P, Y[k]) for k in range(4)) - T1v).astype(np.int64)
    T0 = T0 or float(4 * p * p * n)
    c = _s1_anneal(Y, E, P.add.astype(np.int64), P.sub.astype(np.int64), vals, steps, T0, T1, int(rng.integers(1 << 30)),
                    np.ones(n, np.bool_))
    if c == 0:
        assert ((sum(paf_g(P, Y[k]) for k in range(4)) - T1v) == 0).all()
        return Y, steps
    return None, steps

@numba.njit(cache=True)
def _u_anneal(U, Q, Xr, beta, add, sub, steps, T0, T1, seed):
    """Stage 2 as a +-1 problem on the twisted cosets (U = 0 on constant cosets, never moved):
    Q = sum_k PAF(u_k) - target_u,  Xr = sum_k [CC(beta_k, u_k) - CC(u_k, beta_k)],  beta = alpha - w.
    Cost = sum_h Q[h]^2 + Xr[h]^2 (proportional to the exact reduced objective once stage 1 holds)."""
    np.random.seed(seed)
    K, n = U.shape
    cost = 0
    for h in range(n): cost += Q[h] * Q[h] + Xr[h] * Xr[h]
    dq = np.zeros(n, np.int64); dx = np.zeros(n, np.int64)
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        k = np.random.randint(K); g = np.random.randint(n)
        if U[k, g] == 0: continue
        d = -2 * U[k, g]
        tot = 0
        for h in range(n):
            a = add[g, h]; b_ = sub[g, h]
            q = d * (U[k, a] + U[k, b_])
            if h == 0: q += d * d
            x = d * (beta[k, b_] - beta[k, a])
            dq[h] = q; dx[h] = x
            tot += 2 * Q[h] * q + q * q + 2 * Xr[h] * x + x * x
        if tot <= 0 or np.random.random() < np.exp(-tot / T):
            U[k, g] += d
            for h in range(n):
                Q[h] += dq[h]; Xr[h] += dx[h]
            cost += tot
            if cost == 0: return 0
    return cost

def stage2_setup(P, Y):
    """from stage-1 y: alpha, w, beta, target_u (= (p-2) PAF_w + CC(a,w) + CC(w,a) - Ts)"""
    p = P.p; W = np.where(np.abs(Y) > 1, np.sign(Y), 0).astype(np.int64); A = (Y - (p - 1) * W).astype(np.int64)
    tu = -P.T[:, 1].astype(np.int64) + sum((p - 2) * paf_g(P, W[k]) + cc_g(P, A[k], W[k]) + cc_g(P, W[k], A[k]) for k in range(4))
    return A, W, (A - W).astype(np.int64), tu

def stage2(P, Y, rng, steps, T0=200.0, T1=0.5):
    A, W, beta, tu = stage2_setup(P, Y)
    U = np.where(W == 0, rng.choice(np.array([1, -1]), Y.shape), 0).astype(np.int64)
    Q = (sum(paf_g(P, U[k]) for k in range(4)) - tu).astype(np.int64)
    Xr = sum(cc_g(P, beta[k], U[k]) - cc_g(P, U[k], beta[k]) for k in range(4)).astype(np.int64)
    c = _u_anneal(U, Q, Xr, beta, P.add.astype(np.int64), P.sub.astype(np.int64), steps, T0, T1, int(rng.integers(1 << 30)))
    if c != 0: return None, c
    lab = np.zeros(Y.shape, np.int64)
    for (k, g), y in np.ndenumerate(Y):
        lab[k, g] = LABELS.index((int(A[k, g]), 1, int(U[k, g]))) if W[k, g] == 0 else LABELS.index((int(A[k, g]), 0, int(W[k, g])))
    return lab, 0

def labels_from_y(Y, p, rng):
    """y -> (alpha, type, w); twisted signs random.  Returns labels and the move mask (twisted-sign flips only)."""
    lab = np.zeros(Y.shape, np.int64); allow = np.zeros(Y.shape + (8,), np.bool_)
    for (k, g), y in np.ndenumerate(Y):
        if abs(y) == 1:
            a = int(y); s = int(rng.choice([1, -1]))
            lab[k, g] = LABELS.index((a, 1, s))
            allow[k, g, LABELS.index((a, 1, 1))] = allow[k, g, LABELS.index((a, 1, -1))] = True
        else:
            w = 1 if y > 0 else -1; a = int(y - (p - 1) * w)
            lab[k, g] = LABELS.index((a, 0, w)); allow[k, g, lab[k, g]] = True
    return lab, allow

def run(v, b, p, mods, Hgen, seed, s1steps, s2steps, tries, s2tries):
    assert p % 4 == 3
    P = LSm.Prob(v, b, p, mods, Hgen); rng = np.random.default_rng(seed); t0 = time.time()
    for t in range(tries):
        Y, st1 = stage1(P, rng, s1steps)
        if Y is None:
            print(json.dumps(dict(v=v, try_=t, stage1="fail", steps=st1, secs=round(time.time() - t0, 1))), flush=True); continue
        K = int((np.abs(Y) > 1).sum())
        for r in range(s2tries):
            lab, allow = labels_from_y(Y, p, rng)
            res = LSm.search(P, seed=int(rng.integers(1 << 30)), steps=s2steps, log_every=0, lab0=lab, allow=allow)
            rec = dict(v=v, b=b, p=p, mods=list(mods), try_=t, s2try=r, stage1_steps=st1, constants=K,
                       **{k_: res[k_] for k_ in res if k_ != "lab"}, total_secs=round(time.time() - t0, 1))
            print(json.dumps(rec), flush=True)
            if res.get("found") and res.get("verified"):
                fn = f"logs/paley_staged_v{v}_b{b}_p{p}_{'x'.join(map(str, mods))}_s{seed}.json"
                json.dump(dict(v=v, b=b, p=p, mods=list(mods), H=getattr(P, "H", None), lab=res["lab"]), open(fn, "w"))
                return True
    return False

if __name__ == "__main__":
    a = sys.argv[1:]
    opt = {"--seed": "0", "--s1steps": "20000000", "--s2steps": "100000", "--tries": "20", "--s2tries": "3", "--Hgen": None}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, b, p = pos[:3]; mods = tuple(pos[3:])
    Hgen = None if opt["--Hgen"] is None else [tuple(map(int, e.split('.'))) for e in opt["--Hgen"].split(',')]
    run(v, b, p, mods, Hgen, int(opt["--seed"]), int(opt["--s1steps"]), int(opt["--s2steps"]), int(opt["--tries"]),
        int(opt["--s2tries"]))
