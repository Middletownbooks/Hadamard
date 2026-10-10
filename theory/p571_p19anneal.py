r"""Heuristic stage-1 search (simulated annealing, numba) for the Paley(p) ansatz on Z_n at a FIXED gadget class (W, eps):
find A (A = eps W on gadgets, +-1 on twisted cosets, row sums r_k - (p-1) t_k) with
    sum_k [PAF(A_k)(h) + (p-1) X_k(h)] = -4p - (p-1)^2 s(h),   h = 1..n/2.
Moves swap a +1 and a -1 of one block (row sums stay exact).  Cost = sum_h (L(h) - target(h))^2.
Not a proof of anything when it fails; used (i) as a calibration on the published 419 class (p = 11, n = 38), where
stage-1 solutions exist, and (ii) as a probe of the 571 gadget classes that CP-SAT could not refute (report section 7).
usage: python3 theory/p571_p19anneal.py control419 [secs]
       python3 theory/p571_p19anneal.py classes.pkl --only idxfile [--secs 20]"""
import sys, os, json, time, pickle
import numpy as np
import numba

@numba.njit(cache=True)
def _lags(A, W, n, p, H):
    L = np.zeros(H + 1, np.int64)
    for h in range(1, H + 1):
        tot = 0
        for k in range(4):
            for g in range(n):
                tot += A[k, g] * A[k, (g + h) % n]
                if W[k, g] != 0:
                    tot += (p - 1) * W[k, g] * (A[k, (g + h) % n] + A[k, (g - h) % n])
        L[h] = tot
    return L

@numba.njit(cache=True)
def _anneal(A, W, free_pos, nfree, n, p, H, target, iters, T0, T1, seed):
    np.random.seed(seed)
    L = _lags(A, W, n, p, H)
    cost = 0
    for h in range(1, H + 1): cost += (L[h] - target[h]) ** 2
    best = cost; bestA = A.copy()
    # wx[k, h, i] = sum of W(g) over gadgets g with g + h = i or g - h = i
    wx = np.zeros((4, H + 1, n), np.int64)
    for k in range(4):
        for g in range(n):
            if W[k, g] != 0:
                for h in range(1, H + 1):
                    wx[k, h, (g + h) % n] += W[k, g]; wx[k, h, (g - h) % n] += W[k, g]
    dL = np.zeros(H + 1, np.int64)
    for it in range(iters):
        T = T0 * (T1 / T0) ** (it / iters)
        k = np.random.randint(4)
        if nfree[k] < 2: continue
        i = free_pos[k, np.random.randint(nfree[k])]; j = free_pos[k, np.random.randint(nfree[k])]
        if A[k, i] == A[k, j]: continue
        ai = A[k, i]; aj = A[k, j]
        dc = 0
        for h in range(1, H + 1):
            # PAF change: sum_g A(g)A(g+h); flipping i and j
            d = 0
            ip = (i + h) % n; im = (i - h) % n; jp = (j + h) % n; jm = (j - h) % n
            if 2 * h == n:                          # each position lies in two terms (g and g + n/2)
                d += -4 * ai * A[k, ip] - 4 * aj * A[k, jp]
                if ip == j: d += 8 * ai * aj        # both terms join i and j: unchanged
            else:
                d += -2 * ai * (A[k, ip] + A[k, im]) - 2 * aj * (A[k, jp] + A[k, jm])
                if ip == j or im == j: d += 4 * ai * aj
            d += (p - 1) * (-2 * ai * wx[k, h, i] - 2 * aj * wx[k, h, j])
            dL[h] = d
            dc += (L[h] + d - target[h]) ** 2 - (L[h] - target[h]) ** 2
        if dc <= 0 or np.random.random() < np.exp(-dc / T):
            A[k, i] = -ai; A[k, j] = -aj
            for h in range(1, H + 1): L[h] += dL[h]
            cost += dc
            if cost < best:
                best = cost; bestA = A.copy()
                if best == 0: break
    return best, bestA

def run(W, eps, n, p, R, secs=20.0, seed=0, T0=400.0, T1=2.0, iters=2_000_000):
    W = np.asarray(W, np.int64); eps = np.asarray(eps, np.int64); H = n // 2
    s = np.array([sum(int((W[k] * np.roll(W[k], -h)).sum()) for k in range(4)) for h in range(H + 1)])
    target = np.array([0] + [-4 * p - (p - 1) ** 2 * int(s[h]) for h in range(1, H + 1)], np.int64)
    t = W.sum(1); rng = np.random.default_rng(seed)
    free_pos = np.zeros((4, n), np.int64); nfree = np.zeros(4, np.int64)
    best_all = None; t0 = time.time(); restarts = 0
    while time.time() - t0 < secs:
        A = np.zeros((4, n), np.int64)
        for k in range(4):
            fr = np.nonzero(W[k] == 0)[0]; nfree[k] = len(fr); free_pos[k, :len(fr)] = fr
            A[k, W[k] != 0] = (eps[k] * W[k])[W[k] != 0]
            need = R[k] - (p - 1) * int(t[k]) - int(A[k, W[k] != 0].sum())
            npos = (len(fr) + need) // 2
            if (len(fr) + need) % 2 or npos < 0 or npos > len(fr): return dict(error='row sum infeasible')
            vals = np.array([1] * npos + [-1] * (len(fr) - npos)); rng.shuffle(vals); A[k, fr] = vals
        best, bestA = _anneal(A, W, free_pos, nfree, n, p, H, target, iters, T0, T1, int(rng.integers(1 << 30)))
        Lb = _lags(bestA, W, n, p, H); assert int(((Lb[1:] - target[1:]) ** 2).sum()) == best, "incremental cost drift"
        restarts += 1
        if best_all is None or best < best_all[0]: best_all = (int(best), bestA.copy())
        if best == 0: break
    return dict(best_cost=best_all[0], restarts=restarts, secs=round(time.time() - t0, 1), A=best_all[1])

def control419(secs=120):
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
    from paley_anatomy419 import load
    from paley_algebra import components
    v, p, n = 419, 11, 38; m_ = n * p
    X = load(); Y = []
    for x, sh in zip(X, [0, 208, 208, 208]):
        xs = np.roll(x, -sh); y = np.zeros((n, p), np.int64)
        for i in range(m_): y[i % n, i % p] = xs[i]
        Y.append(y)
    comps = [components(y, p) for y in Y]
    A = np.array([c[0] for c in comps]); W = np.array([c[2] for c in comps]); eps = np.where(W != 0, np.sign(A * W), 0)
    R = (A + (p - 1) * W).sum(1)
    r = run(W, eps, n, p, R, secs=secs)
    ok = None
    if r.get('best_cost') == 0:
        Aa = r['A']; y = Aa + (p - 1) * W
        P = sum(np.array([int((y[k] * np.roll(y[k], -h)).sum()) for h in range(n)]) for k in range(4))
        ok = bool((P[1:] == -4 * p).all()); r['differs_from_published'] = bool((Aa != A).any())
    r.pop('A'); r['s1_exact'] = ok
    print(json.dumps(dict(control419=r)), flush=True)

if __name__ == "__main__":
    if sys.argv[1] == 'control419':
        control419(float(sys.argv[2]) if len(sys.argv) > 2 else 120); sys.exit()
    a = sys.argv[1:]; opt = {"--only": "", "--secs": "20", "--log": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    classes = pickle.load(open(pos[0], 'rb'))
    keys = sorted(classes, key=lambda c: (classes[c]['K'], str(classes[c]['dist']), str(c)))
    idxs = [int(x) for x in open(opt['--only']).read().split()]
    logf = open(opt['--log'], 'a') if opt['--log'] else None
    for idx in idxs:
        info = classes[keys[idx]]
        r = run(np.array(info['W']), np.array(info['eps']), 30, 19, (2, 0, 0, 0), secs=float(opt['--secs']), seed=idx)
        Aa = r.pop('A')
        rec = dict(i=idx, dist=info['dist'], **r)
        if r.get('best_cost') == 0: rec['A'] = Aa.tolist()
        line = json.dumps(rec); print(line, flush=True)
        if logf: logf.write(line + '\n'); logf.flush()
