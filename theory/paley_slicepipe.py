r"""Exact two-stage pipeline inside a Paley(p) + multiplier slice, b = 1 (report section 7.5).

The labels (alpha, type, sign) on G' = Z_n are constant on the orbits of g -> mu*g.  This is the orbit slice of the
multiplier group (squares mod p) x <mu mod n>; the published 419 core lies in it for p = 11, n = 38, mu = 7
(theory/paley_anatomy419.py), and so does the repo's d = 15 slice of 573 (p = 11, n = 52, mu = 29).

Given a gadget (constant orbits with w-sign and J/J' type per block):
  stage 1: all alpha on the twisted orbits of each block with the required block sum (exact DP enumeration), exact
           2+2 join of sum_k PAF(y_k) = T_J (theory/paley_b1pipe.join2x2);
  Lemma 6 at psi = 1 and psi = (-1)^g filters each stage-1 solution;
  stage 2: all u on the twisted orbits, exact 2+2 join of (sum PAF(beta) + p PAF(u), X) = (4v delta, 0);
  completion: bordered Goethals-Seidel (theory/paley_b1pipe.finish), exact H H^T = N I.
Commands:
  given419   positive control with the PUBLISHED 419 gadget held fixed (alpha and u searched exactly)
  random v p n mu K Nm [gadgets] [seed]   sample gadgets in the slice (Lemma 4 + Lemma 6 filtered) and run both stages
  designed v p n mu K Nm [gadgets] [seed] [B]   designed gadgets: signed orbit unions with |PAF_w| <= B per block and
           sum_k PAF(w_k) in {-2..1} off 0 (cyclotomic near-difference families), then both stages
usage: python3 theory/paley_slicepipe.py given419 | random 419 11 38 7 16 11 200 0 | designed 419 11 38 7 16 11 2000 0 2"""
import sys, os, json, time, math, itertools as it
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import paley_b1pipe as BP
from paley_wscreen import S_allowed

def orbits(n, mu):
    seen, out = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * mu) % n
        seen.update(o); out.append(sorted(o))
    return out

def sign_vectors_with_sum(osz, need):
    """all +-1 vectors a (len(osz)) with sum osz*a = need, as an int8 array (exact DP enumeration)"""
    L = len(osz)
    reach = [set() for _ in range(L + 1)]; reach[L] = {0}
    for i in range(L - 1, -1, -1): reach[i] = {t + s * osz[i] for t in reach[i + 1] for s in (1, -1)}
    if need not in reach[0]: return np.zeros((0, L), np.int8)
    out = []
    def rec(i, t, acc):
        if i == L:
            if t == 0: out.append(acc.copy())
            return
        for s in (1, -1):
            if t - s * osz[i] in reach[i + 1]:
                acc.append(s); rec(i + 1, t - s * osz[i], acc); acc.pop()
    rec(0, need, [])
    return np.array(out, np.int8)

def expand(vals, orbs_sel, n, base=None):
    """vals (L, len(orbs_sel)) -> (L, n) with each orbit filled; base (n,) added"""
    M = np.zeros((vals.shape[0], n), np.int64) if base is None else np.tile(base, (vals.shape[0], 1)).astype(np.int64)
    for j, o in enumerate(orbs_sel): M[:, o] = vals[:, [j]]
    return M

def run_gadget(W, Jp, v, p, n, orbs, max_cores=1, log=None, stats=None, stage2=True, record=None):
    """W (4, n) gadget (w on constant orbits), Jp (4, n) bool J' mask.  Returns list of finished cores."""
    m = n * p; TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m - 4 * (p - 1); rs = [2, 0, 0, 0]
    lists = []
    for k in range(4):
        fixed = np.where(W[k] != 0, np.where(Jp[k], -W[k], W[k]), 0)
        tw = [o for o in orbs if W[k, o[0]] == 0]
        need = rs[k] - (p - 1) * int(W[k].sum()) - int(fixed.sum())
        S = sign_vectors_with_sum([len(o) for o in tw], need)
        if len(S) == 0:
            if stats is not None: stats['alpha_empty'] = stats.get('alpha_empty', 0) + 1
            return []
        Ak = expand(S.astype(np.int64), tw, n, base=fixed)
        lists.append((Ak, BP.pafc(Ak + (p - 1) * W[k][None, :])))
    if stats is not None:
        stats['alpha_lists'] = [len(l[0]) for l in lists]
        stats['alpha_bits_x10'] = stats.get('alpha_bits_x10', 0) + int(10 * sum(math.log2(len(l[0])) for l in lists))
        stats['joined'] = stats.get('joined', 0) + 1
    found = []
    # pair the largest alpha list with the smallest (S1 is a symmetric sum over blocks); guard memory
    ordr = sorted(range(4), key=lambda k: len(lists[k][0]))
    perm = [ordr[3], ordr[0], ordr[2], ordr[1]]
    if len(lists[perm[0]][0]) * len(lists[perm[1]][0]) > (1 << 28):
        if stats is not None: stats['too_big'] = stats.get('too_big', 0) + 1
        return []
    L_ = [lists[k] for k in perm]
    for i, j in BP.join2x2([l[1] for l in L_], TJ):
        a0, a1 = divmod(i, len(L_[1][0])); a2, a3 = divmod(j, len(L_[3][0]))
        sel = [L_[0][0][a0], L_[1][0][a1], L_[2][0][a2], L_[3][0][a3]]
        A = np.zeros((4, n), np.int64)
        for slot, k in enumerate(perm): A[k] = sel[slot]
        if stats is not None: stats['stage1'] = stats.get('stage1', 0) + 1
        if not BP.lemma6_real_chars(A, W, p, n):
            if stats is not None: stats['lemma6_cut'] = stats.get('lemma6_cut', 0) + 1
            continue
        if record is not None: record.append((A.copy(), W.copy(), Jp.copy()))
        if not stage2: continue
        for U in stage2_slice(A, W, p, v, n, orbs):
            rec = BP.finish(A, U, W, v, p, n); found.append(rec)
            if log: log(CORE=True, **{k_: rec[k_] for k_ in rec if k_ != 'X'})
            if len(found) >= max_cores: return found
    return found

def stage2_slice(A, W, p, v, n, orbs, max_sols=1):
    beta = A - W; T = np.concatenate([np.eye(1, n, 0, dtype=np.int64)[0] * 4 * v, np.zeros(n, np.int64)])
    lists = []
    for k in range(4):
        tw = [o for o in orbs if W[k, o[0]] == 0]; L = len(tw)
        idx = np.arange(2 ** L, dtype=np.int64)
        S = 1 - 2 * ((idx[:, None] >> np.arange(L)) & 1)
        Ub = expand(S, tw, n)
        sig = np.concatenate([BP.pafc(beta[k])[0][None, :] + p * BP.pafc(Ub),
                              BP.ccx(beta[k], Ub) - np.stack([(Ub * np.roll(beta[k], -h)[None, :]).sum(1) for h in range(n)], 1)], 1)
        lists.append((Ub, sig))
    sols = []
    for i, j in BP.join2x2([l[1] for l in lists], T):
        a0, a1 = divmod(i, len(lists[1][0])); a2, a3 = divmod(j, len(lists[3][0]))
        sols.append(np.array([lists[0][0][a0], lists[1][0][a1], lists[2][0][a2], lists[3][0][a3]]))
        if len(sols) >= max_sols: break
    return sols

def published419_gadget():
    sys.path.insert(0, os.path.join(ROOT, 'theory'))
    from paley_anatomy419 import load
    from paley_algebra import components
    X = load(); v, p, n = 419, 11, 38; m = n * p
    W = np.zeros((4, n), np.int64); Jp = np.zeros((4, n), bool)
    for k, (x, s) in enumerate(zip(X, [0, 208, 208, 208])):
        xs = np.roll(x, -s); y = np.zeros((n, p), np.int64)
        for i in range(m): y[i % n, i % p] = xs[i]
        a, u, w, ok = components(y, p); assert ok
        W[k] = w; Jp[k] = (w != 0) & (w == -a)
    return W, Jp

def given419(max_cores=1):
    v, p, n, mu = 419, 11, 38, 7; orbs = orbits(n, mu)
    W, Jp = published419_gadget()
    t0 = time.time(); stats = {}
    found = run_gadget(W, Jp, v, p, n, orbs, max_cores=max_cores, stats=stats)
    out = dict(v=v, control='published 419 gadget fixed; alpha and u searched exactly in the Paley(11)+mult:7 slice',
               cores=len(found), hadamard=[r['hadamard'] for r in found], secs=round(time.time() - t0, 1), **stats)
    if found:
        X = found[0]['X']
        out['same_as_published_up_to_shift'] = None                   # filled below (informative only)
    print(json.dumps(out)); return out

def lemma4_ok(W, p, n, K):
    s = sum(BP.pafc(W[k])[0] for k in range(4))
    Soff = S_allowed(-4 * p, p, n, K)
    return all(int(s[h]) in Soff for h in range(1, n)) and int((W.sum(1) ** 2).sum()) == K + int(s[1:].sum())

def random_gadgets(v, p, n, mu, K, Nm, n_gadgets=200, seed=0, max_cores=1):
    """sample gadgets: constant orbits chosen at random across the 4 blocks with total size K and J' size Nm, random
    w signs; keep those passing row-sum feasibility, Lemma 4 and Lemma 6 (psi = 1); run both stages exactly."""
    rng = np.random.default_rng(seed); orbs = orbits(n, mu); r = len(orbs); osz = [len(o) for o in orbs]
    logf = open(os.path.join(ROOT, 'theory', 'logs', f'paley_slicepipe_v{v}_mu{mu}_K{K}_s{seed}.log'), 'w')
    def L(**kw): logf.write(json.dumps(kw) + '\n'); logf.flush()
    L(v=v, p=p, n=n, mu=mu, orbits=r, K=K, Nminus=Nm)
    t0 = time.time(); tried = 0; passed = 0; agg = defaultdict(int); found = []
    while passed < n_gadgets and tried < 10 ** 7:
        tried += 1
        perm = rng.permutation(4 * r); tot = 0; const = []
        for s_ in perm:
            if tot + osz[s_ % r] <= K and rng.random() < 0.5: const.append(int(s_)); tot += osz[s_ % r]
        if tot != K: continue
        jp = []; tj = 0
        for s_ in rng.permutation(const):
            if tj + osz[s_ % r] <= Nm and rng.random() < 0.7: jp.append(int(s_)); tj += osz[s_ % r]
        if tj != Nm: continue
        W = np.zeros((4, n), np.int64); Jp = np.zeros((4, n), bool)
        for s_ in const:
            k, j = divmod(s_, r); W[k, orbs[j]] = 1 if rng.random() < 0.5 else -1
        for s_ in jp:
            k, j = divmod(s_, r); Jp[k, orbs[j]] = True
        if not lemma4_ok(W, p, n, K): agg['lemma4_cut'] += 1; continue
        rs = [2, 0, 0, 0]; num = 4 * v - sum((rs[k] - p * int(W[k].sum())) ** 2 for k in range(4))
        if num % p or not BP.sq4(num // p, [n - int((W[k] != 0).sum()) for k in range(4)]): agg['lemma6_gadget_cut'] += 1; continue
        passed += 1; st = {}
        fnd = run_gadget(W, Jp, v, p, n, orbs, max_cores=max_cores, stats=st, log=L)
        for k_, val in st.items():
            if k_ != 'alpha_lists': agg[k_] += val
        if fnd:
            found += fnd
            L(gadget=passed, W=W.tolist(), Jp=Jp.astype(int).tolist(), secs=round(time.time() - t0, 1))
            break
        if passed % 10 == 0: L(progress=True, gadgets=passed, tried=tried, **agg, secs=round(time.time() - t0, 1))
    out = dict(summary=True, v=v, mu=mu, K=K, Nminus=Nm, gadgets=passed, tried=tried, cores=len(found),
               hadamard=[r_['hadamard'] for r_ in found], **agg, secs=round(time.time() - t0, 1))
    L(**out); print(json.dumps(out)); return out

def block_gadgets(n, orbs, kmax, B):
    """all w in {0,+-1}^orbits (constant on orbits) with support size <= kmax and |PAF_w(h)| <= B for every h != 0.
    Returns dict size -> list of (w (n,), t, paf (n,))."""
    r = len(orbs); osz = np.array([len(o) for o in orbs])
    out = defaultdict(list)
    E = np.zeros((r, n), np.int64)
    for j, o in enumerate(orbs): E[j, o] = 1
    for supp_mask in range(1 << r):
        bits = [j for j in range(r) if supp_mask >> j & 1]
        sz = int(osz[bits].sum()) if bits else 0
        if sz > kmax: continue
        if not bits: out[0].append((np.zeros(n, np.int64), 0, np.zeros(n, np.int64))); continue
        L = len(bits)
        S = 1 - 2 * ((np.arange(1 << (L - 1))[:, None] >> np.arange(L)) & 1)      # first sign fixed (+): w ~ -w
        Wm = S @ E[bits]
        P = BP.pafc(Wm)
        ok = (np.abs(P[:, 1:]) <= B).all(1)
        for i in np.nonzero(ok)[0]:
            for sg in (1, -1):
                out[sz].append((sg * Wm[i], int(sg * Wm[i].sum()), P[i]))
    return out

def designed(v, p, n, mu, K, Nm, max_gadgets=2000, seed=0, B=2, Stight=(-2, -1, 0, 1), kmax=None, max_cores=1):
    """designed gadgets: per block, signed orbit unions with |PAF_w| <= B off 0 (cyclotomic near-difference sets);
    4-block combinations with total size K, s(h) = sum_k PAF(w_k)(h) in Stight for all h != 0, sum t^2 = K + sum s,
    row-sum feasibility and Lemma 6 (psi = 1); random J' choice of total size Nm; then both exact stages."""
    rng = np.random.default_rng(seed); orbs = orbits(n, mu); r = len(orbs); osz = [len(o) for o in orbs]
    kmax = kmax or K
    t0 = time.time(); BG = block_gadgets(n, orbs, kmax, B)
    sizes = sorted(BG); nb = {s_: len(BG[s_]) for s_ in sizes}
    logf = open(os.path.join(ROOT, 'theory', 'logs', f'paley_slicepipe_designed_v{v}_mu{mu}_K{K}_B{B}_s{seed}.log'), 'w')
    def L(**kw): logf.write(json.dumps(kw) + '\n'); logf.flush()
    L(v=v, p=p, n=n, mu=mu, orbits=r, K=K, Nminus=Nm, B=B, Stight=list(Stight), block_gadgets=nb, secs=round(time.time() - t0, 1))
    m = n * p; rs = [2, 0, 0, 0]; found = []; agg = defaultdict(int); tried = 0; passed = 0
    comps = [c for c in it.product(sizes, repeat=4) if sum(c) == K]
    while passed < max_gadgets and tried < 5 * 10 ** 6:
        tried += 1
        c = comps[rng.integers(len(comps))]
        pick = [BG[c[k]][rng.integers(len(BG[c[k]]))] for k in range(4)]
        W = np.array([pk[0] for pk in pick]); ts = [pk[1] for pk in pick]
        s_ = sum(pk[2] for pk in pick)
        if not all(int(x) in Stight for x in s_[1:]): agg['s_cut'] += 1; continue
        if int(sum(t * t for t in ts)) != K + int(s_[1:].sum()): agg['t_cut'] += 1; continue
        num = 4 * v - sum((rs[k] - p * ts[k]) ** 2 for k in range(4))
        if num % p or not BP.sq4(num // p, [n - int((W[k] != 0).sum()) for k in range(4)]): agg['lemma6_cut'] += 1; continue
        # J' choice: constant orbits (block, orbit) with total size Nm
        cons = [(k, j) for k in range(4) for j in range(r) if W[k, orbs[j][0]] != 0]
        Jp = None
        for _ in range(200):
            sel = [x for x in cons if rng.random() < Nm / max(K, 1)]
            if sum(osz[j] for (_, j) in sel) == Nm:
                Jp = np.zeros((4, n), bool)
                for (k, j) in sel: Jp[k, orbs[j]] = True
                break
        if Jp is None: agg['jp_fail'] += 1; continue
        passed += 1; st = {}
        fnd = run_gadget(W, Jp, v, p, n, orbs, max_cores=max_cores, stats=st, log=L)
        for k_, val in st.items():
            if k_ != 'alpha_lists': agg[k_] += val
        if fnd:
            found += fnd; L(gadget=passed, W=W.tolist(), Jp=Jp.astype(int).tolist(), secs=round(time.time() - t0, 1)); break
        if passed % 20 == 0: L(progress=True, gadgets=passed, tried=tried, **agg, secs=round(time.time() - t0, 1))
    out = dict(summary=True, v=v, mu=mu, K=K, Nminus=Nm, B=B, gadgets=passed, tried=tried, cores=len(found),
               hadamard=[r_['hadamard'] for r_ in found], **agg, secs=round(time.time() - t0, 1))
    L(**out); print(json.dumps(out)); return out

def gadget_alpha_lists(W, Jp, p, n, orbs):
    """per-block alpha lists (orbit-constant, required block sums) for a gadget; None if some block is infeasible"""
    rs = [2, 0, 0, 0]; out = []
    for k in range(4):
        fixed = np.where(W[k] != 0, np.where(Jp[k], -W[k], W[k]), 0)
        tw = [o for o in orbs if W[k, o[0]] == 0]
        need = rs[k] - (p - 1) * int(W[k].sum()) - int(fixed.sum())
        S = sign_vectors_with_sum([len(o) for o in tw], need)
        if len(S) == 0: return None
        out.append(expand(S.astype(np.int64), tw, n, base=fixed))
    return out

def predicted_stage1(W, Jp, p, n, orbs, rng, N=3000):
    """Gaussian log2 E of stage-1 solutions for this gadget (alpha uniform in its lists) [E]"""
    from paley_compfix import gauss_logp
    AL = gadget_alpha_lists(W, Jp, p, n, orbs)
    if AL is None: return None
    m = n * p; TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m - 4 * (p - 1)
    bits = sum(math.log2(len(a)) for a in AL)
    S1 = sum(BP.pafc(AL[k][rng.integers(len(AL[k]), size=N)] + (p - 1) * W[k][None, :]) for k in range(4))[:, 1:n // 2 + 1]
    try:
        c1, _, _ = gauss_logp(S1, TJ[1:n // 2 + 1], 'S1')
    except AssertionError:
        return None
    return bits - c1

def guided(v, p, n, mu, K, Nm, max_runs=300, seed=0, B=2, Stight=(-2, -1, 0, 1), thresh=0.0, max_cores=1, time_cap=6000, kmax=None, stage2=True):
    """designed gadgets as in `designed`, restricted to |t_k| <= 1, alpha-feasible, and RANKED by the predicted stage-1
    yield: only gadgets with Gaussian log2 E(stage 1) > thresh run the exact stages"""
    rng = np.random.default_rng(seed); orbs = orbits(n, mu); r = len(orbs); osz = [len(o) for o in orbs]
    t0 = time.time(); BG = block_gadgets(n, orbs, kmax or K, B)
    BG = {sz: [g for g in lst if abs(g[1]) <= 1] for sz, lst in BG.items()}
    rec = [] if not stage2 else None
    BG = {sz: lst for sz, lst in BG.items() if lst}
    sizes = sorted(BG)
    logf = open(os.path.join(ROOT, 'theory', 'logs', f'paley_slicepipe_guided_v{v}_mu{mu}_K{K}_s{seed}.log'), 'w')
    def L(**kw): logf.write(json.dumps(kw) + '\n'); logf.flush()
    L(v=v, p=p, n=n, mu=mu, orbits=r, K=K, Nminus=Nm, B=B, thresh=thresh, block_gadgets={sz: len(BG[sz]) for sz in sizes})
    rs = [2, 0, 0, 0]; found = []; agg = defaultdict(int); runs = 0; tried = 0; preds = []
    comps = [c for c in it.product(sizes, repeat=4) if sum(c) == K]
    while runs < max_runs and time.time() - t0 < time_cap:
        tried += 1
        c = comps[rng.integers(len(comps))]
        pick = [BG[c[k]][rng.integers(len(BG[c[k]]))] for k in range(4)]
        W = np.array([pk[0] for pk in pick]); ts = [pk[1] for pk in pick]; s_ = sum(pk[2] for pk in pick)
        if not all(int(x) in Stight for x in s_[1:]): agg['s_cut'] += 1; continue
        if int(sum(t * t for t in ts)) != K + int(s_[1:].sum()): agg['t_cut'] += 1; continue
        num = 4 * v - sum((rs[k] - p * ts[k]) ** 2 for k in range(4))
        if num % p or not BP.sq4(num // p, [n - int((W[k] != 0).sum()) for k in range(4)]): agg['lemma6_cut'] += 1; continue
        cons = [(k, j) for k in range(4) for j in range(r) if W[k, orbs[j][0]] != 0]; Jp = None
        for _ in range(200):
            sel = [x for x in cons if rng.random() < Nm / max(K, 1)]
            if sum(osz[j] for (_, j) in sel) == Nm:
                Jp = np.zeros((4, n), bool)
                for (k, j) in sel: Jp[k, orbs[j]] = True
                break
        if Jp is None: agg['jp_fail'] += 1; continue
        pr = predicted_stage1(W, Jp, p, n, orbs, rng)
        if pr is None: agg['alpha_infeasible'] += 1; continue
        preds.append(round(pr, 1))
        if pr <= thresh: agg['low_yield'] += 1; continue
        runs += 1; st = {}
        n_before = len(rec) if rec is not None else 0
        fnd = run_gadget(W, Jp, v, p, n, orbs, max_cores=max_cores, stats=st, log=L, stage2=stage2, record=rec)
        for k_, val in st.items():
            if k_ != 'alpha_lists': agg[k_] += val
        if rec is not None:
            L(gadget_run=runs, predicted=round(pr, 1), alpha_lists=st.get('alpha_lists'), stage1_new=len(rec) - n_before,
              secs=round(time.time() - t0, 1))
            if len(rec) > n_before:
                np.savez_compressed(os.path.join(ROOT, 'theory', 'logs', f'paley_slicepipe_stage1_v{v}_mu{mu}_K{K}_s{seed}.npz'),
                                    A=np.array([x[0] for x in rec]), W=np.array([x[1] for x in rec]), Jp=np.array([x[2] for x in rec]))
        if fnd:
            found += fnd; L(gadget_run=runs, predicted=round(pr, 1), W=W.tolist(), Jp=Jp.astype(int).tolist(), secs=round(time.time() - t0, 1)); break
        if runs % 10 == 0: L(progress=True, runs=runs, tried=tried, **agg, secs=round(time.time() - t0, 1))
    q = np.percentile(preds, [50, 90, 99]).round(1).tolist() if preds else None
    out = dict(summary=True, v=v, mu=mu, K=K, Nminus=Nm, runs=runs, tried=tried, cores=len(found),
               hadamard=[r_['hadamard'] for r_ in found], predicted_yield_p50_p90_p99=q, **agg, secs=round(time.time() - t0, 1))
    L(**out); print(json.dumps(out)); return out

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == 'given419': given419()
    elif cmd == 'designed':
        a = list(map(int, sys.argv[2:]))
        v, p, n, mu, K, Nm = a[:6]; ng = a[6] if len(a) > 6 else 2000; sd = a[7] if len(a) > 7 else 0; B = a[8] if len(a) > 8 else 2
        designed(v, p, n, mu, K, Nm, ng, sd, B)
    elif cmd == 'guided':
        a = list(map(int, sys.argv[2:]))
        v, p, n, mu, K, Nm = a[:6]; ng = a[6] if len(a) > 6 else 300; sd = a[7] if len(a) > 7 else 0
        guided(v, p, n, mu, K, Nm, ng, sd)
    elif cmd == 'guided1':                      # stage 1 only (large slices): record y's, no stage-2 join
        a = list(map(int, sys.argv[2:]))
        v, p, n, mu, K, Nm = a[:6]; ng = a[6] if len(a) > 6 else 50; sd = a[7] if len(a) > 7 else 0
        km = a[8] if len(a) > 8 else 14; cap = a[9] if len(a) > 9 else 6000
        guided(v, p, n, mu, K, Nm, ng, sd, kmax=km, stage2=False, time_cap=cap)
    elif cmd == 'random':
        a = list(map(int, sys.argv[2:]))
        v, p, n, mu, K, Nm = a[:6]; ng = a[6] if len(a) > 6 else 200; sd = a[7] if len(a) > 7 else 0
        random_gadgets(v, p, n, mu, K, Nm, ng, sd)
