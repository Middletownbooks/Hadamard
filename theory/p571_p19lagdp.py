r"""Independent single-lag certificates for the <61> = Paley(19) slice of v = 571 (report section 7): for a gadget class
(W, eps) and a lag h, compute EXACTLY the set of values of
    sigma_k(h) = PAF(A_k)(h) + 18 X_k(h)
over all A_k in {+-1}^30 with the gadget values fixed (A = eps W) and the row sum sum_g A_k(g) = r_k - 18 t_k, by a
dynamic programme along the cycles of g -> g + h (PAF(h) only couples consecutive cycle elements; X_k(h) is a weighted
sum of the A-values next to gadgets; the row sum is additive).  The class is infeasible if for some h the target
-76 - 324 s(h) is not in the sumset sigma_0 + sigma_1 + sigma_2 + sigma_3.  No solver is involved.
usage: python3 theory/p571_p19lagdp.py classes.pkl [--start i] [--count c] [--log file]
       python3 theory/p571_p19lagdp.py estimate classes.pkl 'verdict-file-glob[,glob...]' [nsample]
         (e.g. 'theory/logs/p571/p19/verdicts_pass1.jsonl.gz,theory/logs/p571/p19/cpsat_s1_pass2.jsonl')
         (exact-marginal estimate log2 [N_A prod_h f_h] of the number of stage-1 solutions per class, for random
          samples of refuted and unresolved classes; f_h = exact fraction of row-sum-consistent A meeting lag h)"""
import sys, os, json, time, pickle
import numpy as np

n = 30; R = (2, 0, 0, 0)

SIG = 300; NS = 2 * SIG + 1; SUM0 = 30; NU = 2 * SUM0 + 1

def _shift(a, ds, du):
    """a[i, j] -> out[i + ds, j + du] (zero fill; values never leave the window by the bounds |sigma| <= 282)"""
    out = np.zeros_like(a)
    i0, i1 = max(0, ds), NS + min(0, ds); j0, j1 = max(0, du), NU + min(0, du)
    out[i0:i1, j0:j1] = a[i0 - ds:i1 - ds, j0 - du:j1 - du]
    return out

def _orconv(a, b):
    from scipy.signal import fftconvolve
    c = fftconvolve(a.astype(np.float64), b.astype(np.float64))
    return c[SIG:SIG + NS, SUM0:SUM0 + NU] > 0.5

_cache = {}
def block_sigma(W, eps, k, h):
    """set of sigma_k(h) = PAF(A_k)(h) + 18 X_k(h) values over A_k with the row sum r_k - 18 t_k (exact DP)"""
    key = (k == 0, tuple(int(x) for x in W[k]), tuple(int(x) for x in eps[k]), h)
    if key in _cache: return _cache[key]
    Aval = [int(eps[k][g] * W[k][g]) for g in range(n)]
    wx = [0] * n
    for g in range(n):
        if W[k][g]:
            wx[(g + h) % n] += int(W[k][g]); wx[(g - h) % n] += int(W[k][g])
    if h == 15: cycles = [[st, st + 15] for st in range(15)]; mult = 2
    else:
        c = int(np.gcd(h, n)); L = n // c; mult = 1
        cycles = [[(st + i * h) % n for i in range(L)] for st in range(c)]
    opts = lambda q: (Aval[q],) if Aval[q] != 0 else (1, -1)
    total = None
    for cyc in cycles:
        C = np.zeros((NS, NU), bool)
        for a0 in opts(cyc[0]):
            D = {1: np.zeros((NS, NU), bool), -1: np.zeros((NS, NU), bool)}
            D[a0][SIG + 18 * wx[cyc[0]] * a0, SUM0 + a0] = True
            for q in cyc[1:]:
                ND = {1: np.zeros((NS, NU), bool), -1: np.zeros((NS, NU), bool)}
                for a2 in opts(q):
                    for a1 in (1, -1):
                        if D[a1].any(): ND[a2] |= _shift(D[a1], a1 * a2 * mult + 18 * wx[q] * a2, a2)
                D = ND
            for a1 in (1, -1):
                if D[a1].any(): C |= _shift(D[a1], 0 if h == 15 else a1 * a0, 0)
        total = C if total is None else _orconv(total, C)
    t = int(sum(W[k])); need = R[k] - 18 * t
    res = set((np.nonzero(total[:, SUM0 + need])[0] - SIG).tolist()) if abs(need) <= SUM0 else set()
    _cache[key] = res
    return res

def brute_sigma(W, eps, k, h, samples=20000, seed=0):
    """random-sample check: every sampled A_k (with the row sum) gives a sigma in block_sigma"""
    rng = np.random.default_rng(seed); S = block_sigma(W, eps, k, h)
    t = int(sum(W[k])); need = R[k] - 18 * t; fixed = W[k] != 0; free = np.nonzero(~fixed)[0]
    hit = 0
    for _ in range(samples):
        A = np.where(fixed, eps[k] * W[k], rng.choice((-1, 1), size=n))
        if A.sum() != need: continue
        hit += 1
        paf = int((A * np.roll(A, -h)).sum()); x = int((W[k] * (np.roll(A, -h) + np.roll(A, h))).sum())
        if paf + 18 * x not in S: return False, hit
    return True, hit

def class_certificate(W, eps):
    """first lag h whose target is unreachable (None if every single lag is feasible)"""
    s = [sum(int((W[k] * np.roll(W[k], -h)).sum()) for k in range(4)) for h in range(16)]
    for h in range(1, 16):
        target = -76 - 324 * s[h]
        sets = [block_sigma(W, eps, k, h) for k in range(4)]
        S01 = {a + b for a in sets[0] for b in sets[1]}
        ok = any((target - c - d) in S01 for c in sets[2] for d in sets[3])
        if not ok: return h
    return None

# ---------------------------------------------------------------- exact per-lag marginals (counts) for estimates
def block_sigma_counts(W, eps, k, h, n_=30, p=19, R_=(2, 0, 0, 0)):
    """dict sigma -> number of A_k (gadgets fixed, row sum r_k - (p-1) t_k) with PAF(A_k)(h) + (p-1) X_k(h) = sigma"""
    from scipy.signal import fftconvolve
    W = np.asarray(W); eps = np.asarray(eps)
    SG = 30 * n_ // 30 + (p - 1) * 2 * int(np.count_nonzero(W[k])) + 2; NSG = 2 * SG + 1; S0 = n_; NSU = 2 * S0 + 1
    Aval = [int(eps[k][g] * W[k][g]) for g in range(n_)]
    wx = [0] * n_
    for g in range(n_):
        if W[k][g]:
            wx[(g + h) % n_] += int(W[k][g]); wx[(g - h) % n_] += int(W[k][g])
    if 2 * h == n_: cycles = [[st, st + h] for st in range(h)]; mult = 2
    else:
        c = int(np.gcd(h, n_)); L = n_ // c; mult = 1
        cycles = [[(st + i * h) % n_ for i in range(L)] for st in range(c)]
    opts = lambda q: (Aval[q],) if Aval[q] != 0 else (1, -1)
    def sh(a, ds, du):
        out = np.zeros_like(a)
        i0, i1 = max(0, ds), NSG + min(0, ds); j0, j1 = max(0, du), NSU + min(0, du)
        out[i0:i1, j0:j1] = a[i0 - ds:i1 - ds, j0 - du:j1 - du]; return out
    total = None
    for cyc in cycles:
        C = np.zeros((NSG, NSU))
        for a0 in opts(cyc[0]):
            D = {1: np.zeros((NSG, NSU)), -1: np.zeros((NSG, NSU))}
            D[a0][SG + (p - 1) * wx[cyc[0]] * a0, S0 + a0] = 1.0
            for q in cyc[1:]:
                ND = {1: np.zeros((NSG, NSU)), -1: np.zeros((NSG, NSU))}
                for a2 in opts(q):
                    for a1 in (1, -1):
                        if D[a1].any(): ND[a2] += sh(D[a1], a1 * a2 * mult + (p - 1) * wx[q] * a2, a2)
                D = ND
            for a1 in (1, -1):
                if D[a1].any(): C += sh(D[a1], 0 if 2 * h == n_ else a1 * a0, 0)
        if total is None: total = C
        else:
            cc = fftconvolve(total, C); total = np.rint(cc[SG:SG + NSG, S0:S0 + NSU])
    t = int(W[k].sum()); need = R_[k] - (p - 1) * t
    col = total[:, S0 + need] if abs(need) <= S0 else np.zeros(NSG)
    return {int(i - SG): float(v) for i, v in enumerate(col) if v > 0.5}

def marginal_estimate(W, eps, n_=30, p=19, R_=(2, 0, 0, 0)):
    """log2 of N_A * prod_h f_h, f_h = exact fraction of row-sum-consistent A satisfying the lag-h equation
    (independent-lags estimate of the number of stage-1 solutions of the class) and the per-lag log2 f_h"""
    import math
    W = np.asarray(W); H = n_ // 2
    s = [sum(int((W[k] * np.roll(W[k], -h)).sum()) for k in range(4)) for h in range(H + 1)]
    logN = None; logs = []
    for h in range(1, H + 1):
        cnts = [block_sigma_counts(W, eps, k, h, n_, p, R_) for k in range(4)]
        if logN is None: logN = sum(math.log2(sum(c.values())) for c in cnts)
        # distribution of sigma_0 + sigma_1, then count target
        tgt = -4 * p - (p - 1) ** 2 * s[h]
        d01 = {}
        for a, ca in cnts[0].items():
            for b, cb in cnts[1].items(): d01[a + b] = d01.get(a + b, 0.0) + ca * cb
        hit = 0.0
        for c, cc in cnts[2].items():
            for d, cd in cnts[3].items(): hit += cc * cd * d01.get(tgt - c - d, 0.0)
        tot = math.prod(sum(c.values()) for c in cnts)
        logs.append(math.log2(hit / tot) if hit > 0 else float('-inf'))
    return logN + sum(logs), logs, logN


def estimate_cli(pkl, verdict_glob, nsample, seed=1):
    """exact-marginal estimates for random samples of refuted (INFEASIBLE) and unresolved (UNKNOWN) classes"""
    import glob, random, statistics
    cl = pickle.load(open(pkl, 'rb'))
    keys = sorted(cl, key=lambda c: (cl[c]['K'], str(cl[c]['dist']), str(c)))
    import gzip
    verdict = {}
    files = [f for g in verdict_glob.split(',') for f in sorted(glob.glob(g))]
    for f in files:                                   # later files override UNKNOWN verdicts of earlier ones
        for l in (gzip.open(f, 'rt') if f.endswith('.gz') else open(f)):
            r = json.loads(l)
            if r['st'] != 'UNKNOWN' or r['i'] not in verdict: verdict[r['i']] = r['st']
    random.seed(seed); res = {}
    for tag in ('INFEASIBLE', 'UNKNOWN'):
        pool = [i for i, v in verdict.items() if v == tag]
        vals = []
        for idx in random.sample(pool, min(nsample, len(pool))):
            info = cl[keys[idx]]
            est, logs, logN = marginal_estimate(np.array(info['W']), np.array(info['eps']))
            vals.append(est)
            print(json.dumps(dict(tag=tag, i=idx, dist=info['dist'], est_log2=round(est, 2), log2_NA=round(logN, 2))), flush=True)
        fin = [v for v in vals if v > -1e300]
        res[tag] = dict(pool=len(pool), sampled=len(vals), single_lag_impossible=len(vals) - len(fin),
                        median_log2=round(statistics.median(fin), 1) if fin else None, max_log2=round(max(fin), 1) if fin else None,
                        mean_E=float(sum(2.0 ** v for v in fin) / len(vals)) if vals else None)
    print(json.dumps(dict(summary=True, **res)), flush=True)

if __name__ == "__main__" and sys.argv[1] == 'estimate':
    estimate_cli(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 100); sys.exit()

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--start": "0", "--count": "1000000", "--log": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    classes = pickle.load(open(pos[0], 'rb'))
    keys = sorted(classes, key=lambda c: (classes[c]['K'], str(classes[c]['dist']), str(c)))
    st0 = int(opt['--start']); cnt = int(opt['--count']); logf = open(opt['--log'], 'a') if opt['--log'] else None
    t0 = time.time(); hist = {}
    for idx in range(st0, min(len(keys), st0 + cnt)):
        info = classes[keys[idx]]
        W = np.array(info['W']); eps = np.array(info['eps'])
        h = class_certificate(W, eps)
        hist[h] = hist.get(h, 0) + 1
        line = json.dumps(dict(i=idx, cert_lag=h))
        if logf: logf.write(line + '\n'); logf.flush()
    print(json.dumps(dict(summary=True, range=[st0, min(len(keys), st0 + cnt)], cert_lag_hist={str(k): v for k, v in sorted(hist.items(), key=lambda t: (t[0] is None, t[0] or 0))},
                          secs=round(time.time() - t0, 1))), flush=True)
