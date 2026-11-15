#!/usr/bin/env python3
"""Exact core search for a character-border slice of theory/fable_tcensus3.py.  [M; a cap is reported as UNDECIDED]
(Plan item 1, logs/fable_fresh/PLAN_phase3.md.)

Core: four +-1 sequences on G = Q x S (Q = Z_2^a, S = Z_s), constant on the orbits of H = Sigma x U, with
  sum_k PAF(x_k)(g) = -4 sum_{rho in T} rho(q)   for every g = (q, z) != 0        [report section 7]
(equivalently sum_k |xhat_k(chi)|^2 = 4t on T and 4p off T, n = 4p, p = 2^a s + t).

Method.
 1. Q-level: every block y (one value per Sigma-orbit of Q, values in Y = {sum over U-orbits B of +-|B|}) whose
    Walsh values satisfy yhat^2 <= target on every character orbit (and |yhat| <= sqrt(4t) on T), by a meet in the
    middle; 2 + 2 join on the energy vectors; every solution expanded to its y-quadruples (blocks up to sign).
 2. For each y-quadruple: per block all lifts (a sign for every cell (A, B), A a Sigma-orbit, B a U-orbit, with
    sum_B eps |B| = y_A for every A); PAF on one lag per orbit of <H, -1> on G \\ 0, computed as eps^T N_g eps.
 3. 2 + 2 join of the four lift lists on the full PAF vector (64-bit linear hash, exact re-check).  A hit is
    re-checked by the definition (PAF on every lag of G) and saved.
Caps: lifts per block (--maxlift), pairs per side (--maxpairs), y-quadruples (--maxquad); hitting one makes the
verdict UNDECIDED for that y-quadruple.
Usage: fable_tcore.py KEY [--maxlift N] [--maxpairs N] [--maxquad N] [--planted]     (KEY: a key of tcensus3.jsonl)
       fable_tcore.py --control112          (Z_2^3 x Z_3, T = {0,1,2,4}, no multipliers: must find cores)
"""
import sys, os, json, itertools, math, time, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fable_tcensus3 as C3

def chi(m, q): return 1 - 2 * (bin(m & q).count('1') & 1)

def build(a, label, T, s, u):
    K = 1 << a
    if label == 'trivial':
        qorbs = [[q] for q in range(K)]; corbs = [[m] for m in range(K)]; sgens = []
    else:
        g = next(x for x in C3.sigma_groups(a, 12) if x['label'] == label)
        qorbs, corbs = g['qorbs'], g['corbs']
        # generators on Q: recover from the orbit partition is not possible; rebuild the group's permutations
        sgens = sigma_perms(a, label)
    seen = [False] * s; uorbs = []
    for z in range(s):
        if seen[z]: continue
        o = []; y = z
        while not seen[y]:
            seen[y] = True; o.append(y); y = (y * u) % s
        uorbs.append(o)
    return dict(a=a, K=K, s=s, u=u, T=sorted(T), qorbs=qorbs, corbs=corbs, uorbs=uorbs, sgens=sgens)

def sigma_perms(a, label):
    """the generating permutations of the group with this label (same construction as fable_tcensus3)"""
    K = 1 << a
    kind_, rest = label.split(' ', 1)
    parts = json.loads(rest.split(']')[0] + ']')
    cs = json.loads('[' + rest.split('[', 2)[2].split(']')[0] + ']')
    jj = json.loads(rest.split('j=')[1]) if 'j=' in rest else None
    offs = [sum(parts[:i]) for i in range(len(parts))]
    blocks = [(d, o) for d, o in zip(parts, offs) if d >= 2]
    gen = {d: C3.gf_gen(d) for d, _ in blocks}
    def block_map(i, c, j=1):
        d, o = blocks[i]; h = C3.gf_pow(gen[d], ((1 << d) - 1) // c, d); hj = C3.gf_pow(h, j, d)
        return [(q & ~(((1 << d) - 1) << o)) | (C3.gf_mul((q >> o) & ((1 << d) - 1), hj, d) << o) for q in range(K)]
    if kind_ == 'prod':
        return [block_map(i, c) for i, c in enumerate(cs) if c > 1]
    nz = [i for i, c in enumerate(cs) if c > 1]
    img = list(range(K))
    for i, j in zip(nz, [1] + jj):
        m = block_map(i, cs[i], j); img = [m[x] for x in img]
    return [img]

# ------------------------------------------------------------------ Q-level with y vectors
def qlevel_solutions(B, p, t, maxy=200000):
    a, K = B['a'], B['K']; qorbs, corbs, T = B['qorbs'], B['corbs'], set(B['T'])
    N = len(qorbs)
    P = np.array([[sum(chi(co[0], q) for q in qo) for qo in qorbs] for co in corbs], dtype=np.int64)
    Torb = [i for i, co in enumerate(corbs) if co[0] in T]
    tgt = np.array([4 * t if i in Torb else 4 * p for i in range(N)], dtype=np.int64)
    Y = C3.sums_alphabet([len(o) for o in B['uorbs']])
    lowb = int(math.isqrt(4 * t)); hib = int(math.isqrt(4 * p))
    h = N // 2
    Yv = np.array(Y, dtype=np.int64)
    if len(Y) ** max(h, N - h) > 6_000_000: return 'cap: block space', None, None
    GL = np.array(list(itertools.product(range(len(Y)), repeat=h)), dtype=np.int32).reshape(-1, h)
    GR = np.array(list(itertools.product(range(len(Y)), repeat=N - h)), dtype=np.int32).reshape(-1, N - h)
    VL = Yv[GL]; VR = Yv[GR]
    WL = VL @ P[:, :h].T; WR = VR @ P[:, h:].T
    bound = np.array([lowb if i in Torb else hib for i in range(N)], dtype=np.int64)
    # keep halves that can still meet the bound on the T rows
    dR = collections.defaultdict(list)
    for idx, kk in enumerate(map(bytes, WR[:, Torb].astype(np.int64))): dR[kk].append(idx)
    keys = list(dR.keys())
    RK = np.array([np.frombuffer(k_, dtype=np.int64) for k_ in keys]).reshape(len(keys), len(Torb))
    ys = []; en = []
    CH = max(1, 2_000_000 // max(1, len(RK)))
    for c0 in range(0, len(VL), CH):
        ok = np.all(np.abs(WL[c0:c0 + CH, Torb][:, None, :] + RK[None, :, :]) <= lowb, axis=2)
        for il_, ik in zip(*np.nonzero(ok)):
            il = c0 + int(il_); lst = dR[keys[int(ik)]]
            W = WL[il][None, :] + WR[lst]
            good = np.all(np.abs(W) <= bound[None, :], axis=1)
            for jx in np.nonzero(good)[0]:
                yv = np.concatenate([VL[il], VR[lst[jx]]])
                # sign normalisation: first nonzero entry positive
                nzi = np.nonzero(yv)[0]
                if len(nzi) and yv[nzi[0]] < 0: continue
                ys.append(yv); en.append(W[jx] * W[jx])
                if len(ys) > maxy: return 'cap: y vectors', None, None
    if not ys: return 0, [], tgt
    ys = np.array(ys); en = np.array(en)
    groups = collections.defaultdict(list)
    for i, e in enumerate(map(bytes, en)): groups[e].append(i)
    ekeys = list(groups.keys()); EV = np.array([np.frombuffer(k_, dtype=np.int64) for k_ in ekeys])
    pairs = collections.defaultdict(list)
    for i in range(len(EV)):
        S2 = EV[i][None, :] + EV[i:]
        for jj in np.nonzero(np.all(S2 <= tgt[None, :], axis=1))[0]:
            pairs[S2[jj].tobytes()].append((i, i + int(jj)))
    quads = []
    for kk, lst in pairs.items():
        comp = (tgt - np.frombuffer(kk, dtype=np.int64)).tobytes()
        if comp not in pairs or kk > comp: continue
        for (i, j) in lst:
            for (k2, l2) in pairs[comp]:
                quads.append((i, j, k2, l2))
    # expand to y-quadruples (as index quadruples into ys), unordered
    out = set()
    for (i, j, k2, l2) in quads:
        for combo in itertools.product(groups[ekeys[i]], groups[ekeys[j]], groups[ekeys[k2]], groups[ekeys[l2]]):
            out.add(tuple(sorted(combo)))
            if len(out) > 2_000_000: return 'cap: y-quadruples', None, None
    return len(out), (ys, sorted(out)), tgt

# ------------------------------------------------------------------ lifts and PAF
def cell_setup(B):
    K, s = B['K'], B['s']; qorbs, uorbs = B['qorbs'], B['uorbs']
    cells = [(j, l) for j in range(len(qorbs)) for l in range(len(uorbs))]
    qo_of = {q: j for j, o in enumerate(qorbs) for q in o}; uo_of = {z: l for l, o in enumerate(uorbs) for z in o}
    nG = K * s
    cell_of = np.array([qo_of[q] * len(uorbs) + uo_of[z] for q in range(K) for z in range(s)], dtype=np.int64)
    # lag classes: orbits of <Sigma, U, -1> on G \ 0
    gens = []
    for g in B['sgens']: gens.append([g[q] * s + z for q in range(K) for z in range(s)])
    gens.append([q * s + (z * B['u']) % s for q in range(K) for z in range(s)])
    gens.append([q * s + (-z) % s for q in range(K) for z in range(s)])
    orbs = C3.group_orbits(gens, nG)
    reps = [o[0] for o in orbs if 0 not in o]
    add = lambda g, h: (((g // s) ^ (h // s)) * s) + ((g % s + h % s) % s)
    R = len(cells)
    Nm = np.zeros((len(reps), R, R), dtype=np.float32)
    for li, g in enumerate(reps):
        for h in range(nG): Nm[li, cell_of[h], cell_of[add(h, g)]] += 1
    T = B['T']
    tg = np.array([-4 * sum(chi(m, g // s) for m in T) for g in reps], dtype=np.int64)
    return dict(cells=cells, cell_of=cell_of, reps=reps, Nm=Nm, tgt=tg, nG=nG, add=add)

def lifts(B, yv, maxlift):
    """sign vectors over cells with sum_B eps |B| = y_A for every Sigma-orbit A"""
    usz = [len(o) for o in B['uorbs']]
    per = []
    for yA in yv:
        opts = [e for e in itertools.product((1, -1), repeat=len(usz)) if sum(x * z for x, z in zip(e, usz)) == yA]
        if not opts: return np.zeros((0, len(yv) * len(usz)), dtype=np.int8)
        per.append(opts)
    tot = 1
    for o in per: tot *= len(o)
    if tot > maxlift: return None
    return np.array([sum(c, ()) for c in itertools.product(*per)], dtype=np.int8)

def paf_features(X, CS):
    Nm = CS['Nm']; L, R, _ = Nm.shape
    Nflat = np.transpose(Nm, (1, 0, 2)).reshape(R, L * R)
    out = np.empty((len(X), L), dtype=np.int64)
    for c0 in range(0, len(X), 20000):
        E = X[c0:c0 + 20000].astype(np.float32)
        Yb = (E @ Nflat).reshape(len(E), L, R)
        out[c0:c0 + 20000] = np.rint(np.einsum('nlr,nr->nl', Yb, E)).astype(np.int64)
    return out

def window(Fs, tgt):
    """exact range/window filter [P]: in a solution, block k has coordinate l within
    [tgt_l - sum_{others} max, tgt_l - sum_{others} min]; iterate to a fixed point.  Returns the filtered lists
    (None if some list becomes empty or a range test fails)."""
    Fs = list(Fs)
    while True:
        if any(len(f) == 0 for f in Fs): return None
        lo = [f.min(0) for f in Fs]; hi = [f.max(0) for f in Fs]
        L = sum(lo); H = sum(hi)
        if np.any(tgt < L) or np.any(tgt > H): return None
        ch = False
        for k in range(4):
            a_ = tgt - (H - hi[k]); b_ = tgt - (L - lo[k])
            keep = np.all((Fs[k] >= a_) & (Fs[k] <= b_), axis=1)
            if not keep.all(): Fs[k] = Fs[k][keep]; ch = True
        if not ch: return Fs

def window_cols(Fs, tgt):
    """as window(), but returns (kept indices per list, filtered lists)"""
    idx = [np.arange(len(f)) for f in Fs]; Fs = list(Fs)
    while True:
        if any(len(f) == 0 for f in Fs): return None
        lo = [f.min(0) for f in Fs]; hi = [f.max(0) for f in Fs]
        L = sum(lo); H = sum(hi)
        if np.any(tgt < L) or np.any(tgt > H): return None
        ch = False
        for k in range(4):
            a_ = tgt - (H - hi[k]); b_ = tgt - (L - lo[k])
            keep = np.all((Fs[k] >= a_) & (Fs[k] <= b_), axis=1)
            if not keep.all(): Fs[k] = Fs[k][keep]; idx[k] = idx[k][keep]; ch = True
        if not ch: return idx, Fs

def join4(Fs, tgt, maxpairs, rng):
    n = [len(f) for f in Fs]
    best = min([((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))], key=lambda pr: max(n[pr[0][0]] * n[pr[0][1]], n[pr[1][0]] * n[pr[1][1]]))
    (i, j), (k, l) = best
    if n[i] * n[j] > maxpairs or n[k] * n[l] > maxpairs: return 'undecided', None
    np.seterr(over='ignore')
    w = rng.integers(-2 ** 62, 2 ** 62, size=len(tgt), dtype=np.int64)
    H = [(f * w).sum(1) for f in Fs]
    right = (H[k][:, None] + H[l][None, :]).ravel()
    order = np.argsort(right); rs = right[order]
    need = (tgt * w).sum() - (H[i][:, None] + H[j][None, :]).ravel()
    pos = np.minimum(np.searchsorted(rs, need), len(rs) - 1)
    hit = np.nonzero(rs[pos] == need)[0]
    for hh in hit:
        a_, b_ = divmod(int(hh), n[j]); c_, d_ = divmod(int(order[pos[hh]]), n[l])
        if np.array_equal(Fs[i][a_] + Fs[j][b_] + Fs[k][c_] + Fs[l][d_], tgt):
            sel = {i: a_, j: b_, k: c_, l: d_}
            return 'found', [sel[0], sel[1], sel[2], sel[3]]
    return 'none', None

def join4_stream(Fs, tgt, maxbucket, rng, maxtotal=2 * 10 ** 10):
    """2 + 2 join in buckets of the value of one coordinate c: pairs (i, j) with F_i[c] + F_j[c] = v are matched only
    against pairs (k, l) with F_k[c] + F_l[c] = tgt[c] - v.  Exact (every match in a bucket is found); memory per
    bucket at most maxbucket pairs per side, else 'undecided'."""
    n = [len(f) for f in Fs]
    (i, j), (k, l) = min([((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))],
                         key=lambda pr: n[pr[0][0]] * n[pr[0][1]] + n[pr[1][0]] * n[pr[1][1]])
    if n[i] * n[j] + n[k] * n[l] > maxtotal: return 'undecided', None
    # coordinate with the most distinct pair-sum values (finest buckets)
    c = int(np.argmax([len(np.unique(Fs[i][:, x])) * len(np.unique(Fs[k][:, x])) for x in range(Fs[0].shape[1])]))
    np.seterr(over='ignore')
    w = rng.integers(-2 ** 62, 2 ** 62, size=len(tgt), dtype=np.int64)
    H = [(f * w).sum(1) for f in Fs]
    def groups(F):
        vals, inv = np.unique(F[:, c], return_inverse=True)
        return {int(v): np.nonzero(inv == q)[0] for q, v in enumerate(vals)}
    gi, gj, gk, gl = groups(Fs[i]), groups(Fs[j]), groups(Fs[k]), groups(Fs[l])
    left_vals = sorted({a + b for a in gi for b in gj})
    rset = {a + b for a in gk for b in gl}
    th = (tgt * w).sum()
    for v in left_vals:
        u = int(tgt[c]) - v
        if u not in rset: continue
        Lp = [(gi[a], gj[v - a]) for a in gi if (v - a) in gj]
        Rp = [(gk[a], gl[u - a]) for a in gk if (u - a) in gl]
        nl = sum(len(x) * len(y) for x, y in Lp); nr = sum(len(x) * len(y) for x, y in Rp)
        if nl > maxbucket or nr > maxbucket: return 'undecided', None
        Rh = np.concatenate([(H[k][x][:, None] + H[l][y][None, :]).ravel() for x, y in Rp])
        Ri = np.concatenate([np.stack(np.meshgrid(x, y, indexing='ij'), -1).reshape(-1, 2) for x, y in Rp])
        order = np.argsort(Rh); rs = Rh[order]
        for x, y in Lp:
            need = th - (H[i][x][:, None] + H[j][y][None, :]).ravel()
            pos = np.minimum(np.searchsorted(rs, need), len(rs) - 1)
            for hh in np.nonzero(rs[pos] == need)[0]:
                a_, b_ = x[hh // len(y)], y[hh % len(y)]
                c_, d_ = Ri[order[pos[hh]]]
                if np.array_equal(Fs[i][a_] + Fs[j][b_] + Fs[k][c_] + Fs[l][d_], tgt):
                    sel = {i: int(a_), j: int(b_), k: int(c_), l: int(d_)}
                    return 'found', [sel[0], sel[1], sel[2], sel[3]]
    return 'none', None

def check_core(B, CS, Xs):
    nG = CS['nG']; s = B['s']; T = B['T']
    full = [np.array([x[CS['cell_of'][h]] for h in range(nG)], dtype=np.int64) for x in Xs]
    ok = True
    for g in range(1, nG):
        sh = np.array([CS['add'](h, g) for h in range(nG)])
        v = sum(int((f * f[sh]).sum()) for f in full)
        if v != -4 * sum(chi(m, g // s) for m in T): ok = False; break
    return ok, [f.tolist() for f in full]

def solve(B, p, t, maxlift, maxpairs, maxquad, planted=False, tag=''):
    t0 = time.time()
    cnt, sols, tgtq = qlevel_solutions(B, p, t)
    print(f'{tag} Q-level: {cnt} y-quadruples ({time.time()-t0:.0f}s)', flush=True)
    if isinstance(cnt, str) or cnt == 0: return dict(qlevel=cnt, verdict='UNDECIDED' if isinstance(cnt, str) else 'NONE (Q-level empty)')
    ys, quads = sols
    CS = cell_setup(B)
    print(f'{tag} {len(CS["reps"])} lag classes, {len(CS["cells"])} cells per block', flush=True)
    rng = np.random.default_rng(1)
    nund = 0; ndone = 0; found = None; cache = {}; NWIN[0] = 0
    for qi, quad in enumerate(quads[:maxquad]):
        Fs = []; Ls = []; bad = False
        for bi in quad:
            if bi not in cache:
                if len(cache) > CACHE: cache.clear()         # memory: keep only a few lift lists
                X = lifts(B, ys[bi], maxlift)
                cache[bi] = None if X is None else (X, paf_features(X, CS) if len(X) else np.zeros((0, len(CS['reps'])), dtype=np.int64))
            if cache[bi] is None: bad = True; break
            Ls.append(cache[bi][0]); Fs.append(cache[bi][1])
        if bad: nund += 1; continue
        if any(len(f) == 0 for f in Fs): ndone += 1; continue
        if planted and qi == 0 and not ESTIMATE:
            pick = [rng.integers(len(f)) for f in Fs]
            tp = sum(f[i] for f, i in zip(Fs, pick))
            wv = window_cols(Fs, tp)
            okw = wv is not None and all(i in set(ki.tolist()) for i, ki in zip(pick, wv[0]))
            print(f'{tag} planted control: window keeps the plant: {okw}', flush=True)
            if not okw: return dict(verdict='CONTROL FAILED (window)')
            Fp = [f[ki] for f, ki in zip(Fs, wv[0])]
            st, sel = join4(Fp, tp, maxpairs, rng)
            print(f'{tag} planted control: join {st}', flush=True)
            if st == 'none': return dict(verdict='CONTROL FAILED (join)')
        if not ESTIMATE:
            Fw = window_cols(Fs, CS['tgt'])
            if Fw is None:
                ndone += 1; NWIN[0] += 1; continue
            keep_idx, Fs = Fw
            Ls = [L_[ki] for L_, ki in zip(Ls, keep_idx)]
        if ESTIMATE and qi < 8:
            # Gaussian estimate of the number of joins hitting the target [E]: independent lag classes, lattice step
            # from the observed values of each coordinate
            lg = 0.0; ns = [len(f) for f in Fs]
            for l in range(len(CS['tgt'])):
                mu = sum(f[:, l].mean() for f in Fs); var = sum(f[:, l].var() for f in Fs)
                vals = np.concatenate([f[:, l] - f[0, l] for f in Fs]); step = int(np.gcd.reduce(np.abs(vals[vals != 0]))) if np.any(vals) else 1
                if var <= 0:
                    lg += 0.0 if abs(mu - CS['tgt'][l]) < 1e-9 else -999; continue
                lg += math.log2(step / math.sqrt(2 * math.pi * var)) - (CS['tgt'][l] - mu) ** 2 / (2 * var) / math.log(2)
            tot = sum(math.log2(n) for n in ns)
            print(f'{tag} y-quadruple {qi}: lift list sizes {ns}; log2(product) {tot:.1f}; log2 E[hits] {tot + lg:.1f} [E]', flush=True)
            if qi == 7: return dict(verdict='ESTIMATE ONLY')
            continue
        st, sel = join4(Fs, CS['tgt'], maxpairs, rng)
        if st == 'undecided' and STREAM: st, sel = join4_stream(Fs, CS['tgt'], maxpairs, rng)
        if st == 'undecided': nund += 1; continue
        ndone += 1
        if st == 'found':
            Xs = [Ls[b][sel[b]] for b in range(4)]
            ok, full = check_core(B, CS, Xs)
            found = dict(quad=[ys[b].tolist() for b in quad], core=full, check=ok)
            print(f'{tag} CORE FOUND for y-quadruple {qi}; definition check {ok}', flush=True)
            break
    nq = min(len(quads), maxquad)
    if found: verdict = 'FOUND'
    elif nund == 0 and len(quads) <= maxquad: verdict = 'NONE (exhaustive)'
    else: verdict = f'UNDECIDED ({nund} y-quadruples capped, {nq} of {len(quads)} examined)'
    return dict(qlevel=cnt, examined=ndone, killed_by_window=NWIN[0], capped=nund, verdict=verdict, found=found, secs=round(time.time() - t0, 1))

ESTIMATE = False
STREAM = False
CACHE = 12
NWIN = [0]

def main():
    global ESTIMATE, STREAM, CACHE
    a_ = sys.argv[1:]
    ESTIMATE = '--estimate' in a_
    STREAM = '--stream' in a_
    CACHE = int(a_[a_.index('--cache') + 1]) if '--cache' in a_ else 12
    get = lambda name, d: a_[a_.index(name) + 1] if name in a_ else d
    maxlift = int(float(get('--maxlift', '3e5'))); maxpairs = int(float(get('--maxpairs', '6e7')))
    maxquad = int(float(get('--maxquad', '1e5')))
    if '--control112' in a_:
        B = build(3, 'trivial', [0, 1, 2, 4], 3, 1)
        res = solve(B, 28, 4, maxlift, maxpairs, maxquad=200, planted=True, tag='control112')
        print('control112:', {k: v for k, v in res.items() if k != 'found'}, flush=True)
        if res.get('found'):
            print('  core check by definition:', res['found']['check'])
        return
    key = a_[0]; p, a, label, T, s, u, uo = json.loads(key)
    B = build(a, label, T, s, u); t = len(T)
    res = solve(B, p, t, maxlift, maxpairs, maxquad, planted='--planted' in a_, tag=f'p={p} a={a}')
    res['key'] = key
    print(json.dumps({k: v for k, v in res.items() if k != 'found'}), flush=True)
    with open('logs/fable_fresh/tcore.jsonl', 'a') as fo: fo.write(json.dumps({k: v for k, v in res.items() if k != 'found'}) + '\n')
    if res.get('found'):
        json.dump(res, open(f'logs/fable_fresh/tcore_found_{p}_{a}.json', 'w'))

if __name__ == '__main__':
    main()
