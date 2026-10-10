"""Joint simulated annealing inside a Paley(p) + multiplier slice (b = 1): no staging.

Variables: one label (alpha, type, sign) per block and per orbit of g -> mu g on G' = Z_n (4 x #orbits labels, 3 bits
each).  Every block is then mu-invariant on G' (and Paley-twisted on Z_p), so the reduced PAF R[h, c] (c = lag class
0 / QR / NQR on Z_p) is constant on <mu>-orbits of lags h: we keep R only at one representative per lag orbit and
weight it by the orbit size, which gives exactly the full objective  sum_{h, c} w_c (R - T)^2  (w = 1, (p-1)/2, (p-1)/2;
T = 4m at (0, 0), -4 elsewhere).
Move: change one block's label on one whole orbit; the change of R at the representative lags is computed from the
pairs touching that orbit only.  Metropolis with geometric cooling, restarts.  A zero-cost state is expanded and checked
by FFT on G' x Z_p, then completed to a bordered GS Hadamard matrix (theory/paley_b1pipe.finish).
usage: slice_sa.py v p n mu [--tabu 1] [--steps 2e8] [--restarts 10] [--T0 t] [--T1 0.5] [--seed s]"""
import sys, os, json, time, math
import numpy as np
import numba
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
from paley_ansatz import LABELS, LC, corr_table

def mu_orbits(n, mu):
    seen, out = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * mu) % n
        seen.update(o); out.append(sorted(o))
    return out

class Slice:
    def __init__(s, v, p, n, mu):
        s.v, s.p, s.n, s.mu = v, p, n, mu; s.m = n * p; assert s.m == v - 1
        s.orbs = mu_orbits(n, mu); s.r = len(s.orbs)
        s.oid = np.zeros(n, np.int64)
        for j, o in enumerate(s.orbs): s.oid[o] = j
        s.ostart = np.cumsum([0] + [len(o) for o in s.orbs]).astype(np.int64)
        s.oelems = np.array([g for o in s.orbs for g in o], np.int64)
        lag_orbs = mu_orbits(n, mu)                       # lag orbits = same partition (mu acts on lags the same way)
        s.H = np.array([o[0] for o in lag_orbs], np.int64); s.Hw = np.array([len(o) for o in lag_orbs], np.int64)
        C, _ = corr_table(p); s.C = C.astype(np.int64)
        s.w = np.array([1, (p - 1) // 2, (p - 1) // 2], np.int64)
        T = np.full((len(s.H), 3), -4, np.int64); T[0, 0] = 4 * s.m      # H[0] = 0 (orbit of 0 listed first)
        assert s.H[0] == 0
        s.T = T

    def R_all(s, A):
        """A (4, n, 3) -> E-part: total reduced PAF at representative lags (len(H), 3)"""
        n = s.n; R = np.zeros((len(s.H), 3), np.int64)
        for k in range(4):
            for i, h in enumerate(s.H):
                B = A[k][(np.arange(n) + h) % n]
                R[i] += np.einsum('gi,gj,ijc->c', A[k], B, s.C)
        return R

@numba.njit(cache=True)
def _anneal(lab, A, E, LCa, C, w, H, Hw, oelems, ostart, oid, n, steps, T0, T1, seed):
    np.random.seed(seed)
    K = lab.shape[0]; r = lab.shape[1]; nh = H.shape[0]
    cost = 0
    for i in range(nh):
        for c in range(3): cost += Hw[i] * w[c] * E[i, c] * E[i, c]
    dR = np.zeros((nh, 3), np.int64)
    best = cost
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        k = np.random.randint(K); j = np.random.randint(r); L = np.random.randint(8)
        if L == lab[k, j]: continue
        a0 = LCa[L, 0]; a1 = LCa[L, 1]; a2 = LCa[L, 2]
        for i in range(nh):
            for c in range(3): dR[i, c] = 0
        # pairs (x, x+h) touching orbit j: x in O (any x+h), or x+h in O with x not in O
        for t in range(ostart[j], ostart[j + 1]):
            x = oelems[t]
            o0 = A[k, x, 0]; o1 = A[k, x, 1]; o2 = A[k, x, 2]
            for i in range(nh):
                y = (x + H[i]) % n
                if oid[y] == j:
                    q0 = a0; q1 = a1; q2 = a2                 # new value at y too
                else:
                    q0 = A[k, y, 0]; q1 = A[k, y, 1]; q2 = A[k, y, 2]
                p0 = A[k, y, 0]; p1 = A[k, y, 1]; p2 = A[k, y, 2]
                for c in range(3):
                    newv = (a0 * (C[0, 0, c] * q0 + C[0, 1, c] * q1 + C[0, 2, c] * q2)
                            + a1 * (C[1, 0, c] * q0 + C[1, 1, c] * q1 + C[1, 2, c] * q2)
                            + a2 * (C[2, 0, c] * q0 + C[2, 1, c] * q1 + C[2, 2, c] * q2))
                    oldv = (o0 * (C[0, 0, c] * p0 + C[0, 1, c] * p1 + C[0, 2, c] * p2)
                            + o1 * (C[1, 0, c] * p0 + C[1, 1, c] * p1 + C[1, 2, c] * p2)
                            + o2 * (C[2, 0, c] * p0 + C[2, 1, c] * p1 + C[2, 2, c] * p2))
                    dR[i, c] += newv - oldv
                # pair (x - h, x) with x - h outside O
                z = (x - H[i]) % n
                if oid[z] != j:
                    b0 = A[k, z, 0]; b1 = A[k, z, 1]; b2 = A[k, z, 2]
                    for c in range(3):
                        newv = (b0 * (C[0, 0, c] * a0 + C[0, 1, c] * a1 + C[0, 2, c] * a2)
                                + b1 * (C[1, 0, c] * a0 + C[1, 1, c] * a1 + C[1, 2, c] * a2)
                                + b2 * (C[2, 0, c] * a0 + C[2, 1, c] * a1 + C[2, 2, c] * a2))
                        oldv = (b0 * (C[0, 0, c] * o0 + C[0, 1, c] * o1 + C[0, 2, c] * o2)
                                + b1 * (C[1, 0, c] * o0 + C[1, 1, c] * o1 + C[1, 2, c] * o2)
                                + b2 * (C[2, 0, c] * o0 + C[2, 1, c] * o1 + C[2, 2, c] * o2))
                        dR[i, c] += newv - oldv
        tot = 0
        for i in range(nh):
            for c in range(3):
                tot += Hw[i] * w[c] * (2 * E[i, c] * dR[i, c] + dR[i, c] * dR[i, c])
        if tot <= 0 or np.random.random() < np.exp(-tot / T):
            lab[k, j] = L
            for t in range(ostart[j], ostart[j + 1]):
                x = oelems[t]; A[k, x, 0] = a0; A[k, x, 1] = a1; A[k, x, 2] = a2
            for i in range(nh):
                for c in range(3): E[i, c] += dR[i, c]
            cost += tot
            if cost < best: best = cost
            if cost == 0: return 0, best, it
    return cost, best, steps

@numba.njit(cache=True)
def _delta(k, j, L, A, LCa, C, H, oelems, ostart, oid, n, dR):
    nh = H.shape[0]
    a0 = LCa[L, 0]; a1 = LCa[L, 1]; a2 = LCa[L, 2]
    for i in range(nh):
        for c in range(3): dR[i, c] = 0
    for t in range(ostart[j], ostart[j + 1]):
        x = oelems[t]
        o0 = A[k, x, 0]; o1 = A[k, x, 1]; o2 = A[k, x, 2]
        for i in range(nh):
            y = (x + H[i]) % n
            p0 = A[k, y, 0]; p1 = A[k, y, 1]; p2 = A[k, y, 2]
            if oid[y] == j:
                q0 = a0; q1 = a1; q2 = a2
            else:
                q0 = p0; q1 = p1; q2 = p2
            for c in range(3):
                dR[i, c] += (a0 * (C[0, 0, c] * q0 + C[0, 1, c] * q1 + C[0, 2, c] * q2)
                             + a1 * (C[1, 0, c] * q0 + C[1, 1, c] * q1 + C[1, 2, c] * q2)
                             + a2 * (C[2, 0, c] * q0 + C[2, 1, c] * q1 + C[2, 2, c] * q2)) \
                            - (o0 * (C[0, 0, c] * p0 + C[0, 1, c] * p1 + C[0, 2, c] * p2)
                               + o1 * (C[1, 0, c] * p0 + C[1, 1, c] * p1 + C[1, 2, c] * p2)
                               + o2 * (C[2, 0, c] * p0 + C[2, 1, c] * p1 + C[2, 2, c] * p2))
            z = (x - H[i]) % n
            if oid[z] != j:
                b0 = A[k, z, 0]; b1 = A[k, z, 1]; b2 = A[k, z, 2]
                for c in range(3):
                    dR[i, c] += (b0 * (C[0, 0, c] * (a0 - o0) + C[0, 1, c] * (a1 - o1) + C[0, 2, c] * (a2 - o2))
                                 + b1 * (C[1, 0, c] * (a0 - o0) + C[1, 1, c] * (a1 - o1) + C[1, 2, c] * (a2 - o2))
                                 + b2 * (C[2, 0, c] * (a0 - o0) + C[2, 1, c] * (a1 - o1) + C[2, 2, c] * (a2 - o2)))

@numba.njit(cache=True)
def _tabu(lab, A, E, LCa, C, w, H, Hw, oelems, ostart, oid, n, steps, tenure, seed):
    """best-move tabu search (aspiration: a tabu move is allowed if it beats the best cost); returns (cost, best, it)"""
    np.random.seed(seed)
    K = lab.shape[0]; r = lab.shape[1]; nh = H.shape[0]
    cost = 0
    for i in range(nh):
        for c in range(3): cost += Hw[i] * w[c] * E[i, c] * E[i, c]
    best = cost
    tl = np.zeros((K, r), np.int64)
    dR = np.zeros((nh, 3), np.int64); bdR = np.zeros((nh, 3), np.int64)
    for it in range(steps):
        bt = 1 << 62; bk = -1; bj = -1; bL = -1; bkey = 2.0
        for k in range(K):
            for j in range(r):
                for L in range(8):
                    if L == lab[k, j]: continue
                    _delta(k, j, L, A, LCa, C, H, oelems, ostart, oid, n, dR)
                    tot = 0
                    for i in range(nh):
                        for c in range(3): tot += Hw[i] * w[c] * (2 * E[i, c] * dR[i, c] + dR[i, c] * dR[i, c])
                    if tl[k, j] > it and cost + tot >= best: continue
                    key = np.random.random()
                    if tot < bt or (tot == bt and key < bkey):
                        bt = tot; bk = k; bj = j; bL = L; bkey = key
                        for i in range(nh):
                            for c in range(3): bdR[i, c] = dR[i, c]
        if bk < 0:
            tl[:, :] = 0; continue
        lab[bk, bj] = bL
        for t in range(ostart[bj], ostart[bj + 1]):
            x = oelems[t]; A[bk, x, 0] = LCa[bL, 0]; A[bk, x, 1] = LCa[bL, 1]; A[bk, x, 2] = LCa[bL, 2]
        for i in range(nh):
            for c in range(3): E[i, c] += bdR[i, c]
        cost += bt
        tl[bk, bj] = it + tenure + np.random.randint(tenure + 1)
        if cost < best: best = cost
        if cost == 0: return 0, best, it
    return cost, best, steps

def expand(S, lab):
    A = np.zeros((4, S.n, 3), np.int64)
    for k in range(4):
        for j, o in enumerate(S.orbs): A[k, o] = LC[lab[k, j]]
    return A

TABU = 0
def run(v, p, n, mu, steps, restarts, T0, T1, seed):
    S = Slice(v, p, n, mu); rng = np.random.default_rng(seed)
    LCa = LC.astype(np.int64)
    for r in range(restarts):
        lab = rng.integers(0, 8, (4, S.r)).astype(np.int64)
        A = expand(S, lab); E = S.R_all(A) - S.T
        t0 = time.time()
        if TABU:
            c, best, it = _tabu(lab, A, E, LCa, S.C, S.w, S.H, S.Hw, S.oelems, S.ostart, S.oid, n, steps,
                                max(3, S.r // 3), int(rng.integers(1 << 30)))
        else:
            c, best, it = _anneal(lab, A, E, LCa, S.C, S.w, S.H, S.Hw, S.oelems, S.ostart, S.oid, n, steps,
                                  T0 or float(8 * p * p), T1, int(rng.integers(1 << 30)))
        rec = dict(v=v, p=p, n=n, mu=mu, restart=r, final=int(c), best=int(best), steps=int(it), secs=round(time.time() - t0, 1))
        if c == 0:
            assert (S.R_all(A) == S.T).all()                       # exact recheck of the reduced equations
            import paley_b1pipe as BP
            Al, Ul, Wl = A[:, :, 0], A[:, :, 1], A[:, :, 2]
            fin = BP.finish(Al, Ul, Wl, v, p, n)
            rec.update(core_ok=fin['core_ok'], hadamard=fin['hadamard'])
            print(json.dumps(rec), flush=True)
            if fin['hadamard']:
                print(json.dumps(dict(FOUND=True, v=v)), flush=True); return True
        else:
            print(json.dumps(rec), flush=True)
    return False

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--tabu": "0", "--steps": "2e8", "--restarts": "10", "--T0": "0", "--T1": "0.5", "--seed": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, p, n, mu = pos
    TABU = int(opt['--tabu'])
    run(v, p, n, mu, int(float(opt["--steps"])), int(opt["--restarts"]), float(opt["--T0"]), float(opt["--T1"]), int(opt["--seed"]))
