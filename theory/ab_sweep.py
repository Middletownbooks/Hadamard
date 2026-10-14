r"""Sweep of generalised-Paley (GP) type slices for the Alpoege-style builders, with exact screens [M].

Configurations (b >= 3, split G = Q x H):  v, c (even, c | v + 1), h = (v + 1)/c - 1, b = c - 1, m = c h;
  Q abelian of order c (axes: its cyclic factors), H = H' x Z_q with q a prime, q || h (Z_q a direct factor of H),
  M <= F_q^* of index s (order f = (q - 1)/s) acting on Z_q, and a product group U of per-axis scalars on G' = Q x H'
  (trivial, or one cyclic subgroup of units on one axis of H' (or of Q), or two such).  S = U x M.
b = 0 at composite v = F q: G = Z_F x Z_q, same M, U on Z_F (trivial or {+-1}).
Exact screens (all necessary conditions of the core equations; see report/alpoge_design_answer.md):
  A  Q-level:  per block, the Q-compression z (c coset sums of y) and the lag-0 energy E = sum_b y(b)^2 over all
     letter assignments (DP per Q-coset); per-character bound |z^(psi)|^2 <= 4b (psi != 1); 4-block join of
     (PAF_Q(z), E) against (T_Q, T1(0)).  For c = 1 (b = 0) this is lag 0 + rows.
  B  alpha-congruence: f0 = largest f with every letter = +-1 (mod 2f); alpha = y (mod 2 f0) in {+-1};
     sum_k PAF_{G'}(alpha_k) = T1 (mod 4 f0) at every lag, rows of alpha = rows of y (mod 2 f0) for a row quadruple
     admitted by A.  (U-invariant alpha only.)
usage: python3 theory/ab_sweep.py --v 659 [--v ...] [--c 4,6] [--out FILE] [--maxGeff 24] [--b0]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import AG, Slice, factorint, primitive_root, alpha_congruence_mod


def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield []; return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k): yield [k] + rest


def abelian_types(n):
    """list of tuples of prime-power cyclic factors"""
    per = [[[p ** e for e in lam] for lam in partitions(k)] for p, k in factorint(n).items()]
    out = []
    for combo in itertools.product(*per):
        out.append(tuple(sorted([x for part in combo for x in part], reverse=True)))
    return out if out else [tuple()]


def unit_subgroup_gens(n):
    """generators (one per cyclic subgroup) of the cyclic subgroups of (Z/n)^*, excluding the trivial one;
    n a prime power (odd) or 2^k; returns list of (order, generator)"""
    if n <= 2: return []
    units = [u for u in range(1, n) if math.gcd(u, n) == 1]
    seen = set(); out = []
    for u in units:
        if u == 1: continue
        H = set(); x = 1
        while True:
            x = x * u % n; H.add(x)
            if x == 1: break
        key = frozenset(H)
        if key in seen: continue
        seen.add(key); out.append((len(H), u))
    return sorted(out)


def mgen(q, s):
    """generator of the index-s subgroup of F_q^*"""
    g = primitive_root(q)
    return pow(g, s, q)


def configs_bge3(v, cs=(4, 6), maxU=1):
    out = []
    for c in cs:
        if c % 2 or (v + 1) % c: continue
        h = (v + 1) // c - 1; b = c - 1; m = c * h
        fh = factorint(h)
        for Qt in abelian_types(c):
            for q, e in fh.items():
                if e != 1 or q == 2: continue
                hp = h // q
                for Ht in abelian_types(hp):
                    gp_axes = list(Qt) + list(Ht)
                    mods = tuple(gp_axes) + (q,)
                    nQ = len(Qt)
                    for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                        f = (q - 1) // s
                        if f < 2: continue
                        # U menu: trivial + one subgroup on one axis (+ optionally two)
                        Umenus = [[]]
                        axis_opts = []
                        for ai, n in enumerate(gp_axes):
                            for (o, u) in unit_subgroup_gens(n):
                                axis_opts.append((ai, o, u))
                        for k in range(1, maxU + 1):
                            for combo in itertools.combinations(axis_opts, k):
                                if len({a for a, _, _ in combo}) < k: continue
                                Umenus.append(list(combo))
                        # -1 on a subset (>= 2) of the axes of order > 2 (one generator acting on several axes)
                        big = [ai for ai, n in enumerate(gp_axes) if n > 2]
                        for k in range(2, len(big) + 1):
                            for A in itertools.combinations(big, k):
                                Umenus.append([('neg', A)])
                        for U in Umenus:
                            out.append(dict(v=v, c=c, b=b, h=h, m=m, Q=Qt, H1=Ht, q=q, s=s, f=f, mods=mods, nQ=nQ,
                                            U=U))
    return out


def configs_b0(v, maxU=1):
    out = []
    fv = factorint(v)
    if len(fv) < 2 and list(fv.values()) == [1]: return out
    for q, e in fv.items():
        if e != 1: continue
        F = v // q
        if F == 1: continue
        for Ft in abelian_types(F):
            mods = tuple(Ft) + (q,)
            for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                f = (q - 1) // s
                if f < 2: continue
                Umenus = [[]]
                axis_opts = [(ai, o, u) for ai, n in enumerate(Ft) for (o, u) in unit_subgroup_gens(n)]
                for k in range(1, maxU + 1):
                    for combo in itertools.combinations(axis_opts, k):
                        if len({a for a, _, _ in combo}) < k: continue
                        Umenus.append(list(combo))
                for U in Umenus:
                    out.append(dict(v=v, c=1, b=0, h=None, m=v, Q=(), H1=Ft, q=q, s=s, f=f, mods=mods, nQ=0, U=U))
    return out


def build(cfg):
    mods = cfg['mods']; G = AG(mods); qa = len(mods) - 1
    gens = []
    g = [1] * len(mods); g[qa] = mgen(cfg['q'], cfg['s']); gens.append(tuple(g))
    for gen in cfg['U']:
        g = [1] * len(mods)
        if gen[0] == 'neg':
            for a in gen[1]: g[a] = mods[a] - 1
        else:
            a, o, u = gen; g[a] = u
        gens.append(tuple(g))
    if cfg['b'] >= 3:
        nQ = cfg['nQ']
        Hmask = (G.co[:, :nQ] == 0).all(1) if nQ else np.ones(G.n, bool)
        sl = Slice(mods, cfg['v'], cfg['b'], gens, qa, Hmask)
    else:
        sl = Slice(mods, cfg['v'], cfg['b'], gens, qa, None)
    return sl


def qcoset_map(sl, nQ):
    """Q-coset index of each element of G' (split: Q-coordinates); Q group"""
    Gp = sl.Gp
    if nQ == 0: return np.zeros(Gp.n, np.int64), AG(())
    Qg = AG(Gp.mods[:nQ])
    return Qg.idx(Gp.co[:, :nQ]), Qg


def _bits_sum(A, B, cap):
    """Minkowski sum of two sets of non-negative ints given as python-int bitsets, truncated at cap"""
    out = 0; a = A; i = 0
    mask = (1 << (cap + 1)) - 1
    while a:
        if a & 1: out |= (B << i)
        a >>= 1; i += 1
    return out & mask


def screen_A(sl, cfg, max_vec=400000):
    """exact Q-level screen (lag-0 energy + Q-compression); returns feasibility and admissible row quadruples"""
    nQ = cfg['nQ']; b = cfg['b']; v = cfg['v']
    qc, Qg = qcoset_map(sl, nQ); nq = Qg.n
    T0 = int(sl.T1[0])
    TQ = np.zeros(nq, np.int64); np.add.at(TQ, qc, sl.T1)
    # U-orbits on Q (U restricted to the Q axes)
    if nq > 1:
        Uq = sorted({tuple(u[:nQ]) for u in sl.Sp})
        Qperm = [Qg.act(u) for u in Uq]
        qorb = -np.ones(nq, np.int64); qreps = []
        for t in range(nq):
            if qorb[t] >= 0: continue
            for P in Qperm: qorb[P[t]] = len(qreps)
            qreps.append(t)
    else:
        qorb = np.zeros(1, np.int64); qreps = [0]
    nqo = len(qreps); qosz = np.bincount(qorb)
    # G'-orbits over each Q-orbit
    per = [[] for _ in range(nqo)]
    for j in range(len(sl.breps)):
        el = np.nonzero(sl.borb == j)[0]
        qs = set(qc[el].tolist()); o = {int(qorb[t]) for t in qs}
        assert len(o) == 1
        o = o.pop(); w = len(el); nper = w // len(qs)
        assert nper * len(qs) == w and len(qs) == qosz[o]
        per[o].append((w, nper, sl.letters[j]))
    minE = sum(w * min(y * y for y in L) for lst in per for (w, n_, L) in lst)
    Emax = T0 - 3 * minE
    if Emax < minE: return dict(A=False, reason="lag0 energy below minimum", minE_block=minE, T0=T0)
    # DP per Q-orbit: zval -> bitset of E
    tabs = []
    minE_o = [sum(w * min(y * y for y in L) for (w, n_, L) in lst) for lst in per]
    for o in range(nqo):
        cap = Emax - (minE - minE_o[o])
        st = {0: 1}
        for (w, n_, L) in per[o]:
            new = defaultdict(int)
            for z, bits in st.items():
                for y in L:
                    sh = w * y * y
                    if sh > cap: continue
                    nb = (bits << sh) & ((1 << (cap + 1)) - 1)
                    if nb: new[z + n_ * y] |= nb
            st = dict(new)
        tabs.append(st)
    # near-constant z vectors: |z^(psi)|^2 <= 4b (psi != 1)  =>  ||z - zbar||^2 <= 4b (nq - 1)/nq; z_q all of one parity
    keys0 = sorted(tabs[0])
    if nq > 1:
        D = math.isqrt(4 * b * (nq - 1) // nq) + 1
        win = np.arange(-2 * D, 2 * D + 1, 2)
        pats = np.array(list(itertools.product(win, repeat=nqo - 1)), np.int64).reshape(-1, nqo - 1)
        full = np.concatenate([np.zeros((len(pats), 1), np.int64), pats], 1)[:, qorb]
        P = np.abs(np.fft.fftn(full.reshape((-1,) + Qg.mods).astype(float), axes=tuple(range(1, 1 + len(Qg.mods))))) ** 2
        pats = pats[(P.reshape(len(pats), -1)[:, 1:] <= 4 * b + 1e-6).all(1)]
    else:
        pats = np.zeros((1, 0), np.int64)
    sigE = defaultdict(int); rows_of = defaultdict(set); nvec = 0
    for z0 in keys0:
        for pt in pats:
            zv = (z0,) + tuple(int(z0 + d) for d in pt)
            if any(zv[o] not in tabs[o] for o in range(1, nqo)): continue
            zfull = np.array([zv[qorb[t]] for t in range(nq)], np.int64)
            pq = Qg.paf(zfull[None, :])[0] if nq > 1 else zfull ** 2
            E = tabs[0][zv[0]]
            for o in range(1, nqo): E = _bits_sum(E, tabs[o][zv[o]], Emax)
            if not E: continue
            key = tuple(pq.tolist())
            sigE[key] |= E; rows_of[key].add(abs(int(zfull.sum())))
            nvec += 1
            if nvec > max_vec: return dict(A=None, reason="too many z vectors", nvec=nvec)
    sigs = list(sigE)
    if not sigs: return dict(A=False, reason="no block z-vector", nvec=nvec)
    if len(sigs) > 20000: return dict(A=None, reason="too many Q-signatures", nsig=len(sigs))
    Sarr = np.array(sigs, np.int64); ns = len(sigs)
    rng = np.random.default_rng(777)
    wv = rng.integers(1, 2 ** 62, size=Sarr.shape[1], dtype=np.int64)
    hs = Sarr @ wv                                         # wraps mod 2^64: fine for hashing, exact check below
    th = int((TQ.astype(np.int64) @ wv))
    I, J = np.triu_indices(ns)
    ph = hs[I] + hs[J]
    o = np.argsort(ph, kind='stable'); ps = ph[o]
    need = np.int64(th) - ph
    lo = np.searchsorted(ps, need, 'left'); hi = np.searchsorted(ps, need, 'right')
    cand = np.nonzero(hi > lo)[0]
    Ebits = [sigE[k] for k in sigs]
    def rev(x):
        r = 0; i = 0
        while x:
            if x & 1 and T0 - i >= 0: r |= 1 << (T0 - i)
            x >>= 1; i += 1
        return r
    pairE = {}
    def PE(p):
        if p not in pairE: pairE[p] = _bits_sum(Ebits[I[p]], Ebits[J[p]], T0)
        return pairE[p]
    quads = set(); feasible = False; checks = 0
    for p in cand:
        for qi in o[lo[p]:hi[p]]:
            if qi < p: continue
            if not ((Sarr[I[p]] + Sarr[J[p]] + Sarr[I[qi]] + Sarr[J[qi]]) == TQ).all(): continue
            checks += 1
            if checks > 2 * 10 ** 6: return dict(A=None, reason="E-check cap", nsig=ns)
            if PE(p) & rev(PE(qi)):
                feasible = True
                for r in itertools.product(rows_of[sigs[I[p]]], rows_of[sigs[J[p]]], rows_of[sigs[I[qi]]], rows_of[sigs[J[qi]]]):
                    if sum(x_ * x_ for x_ in r) == 4 * v: quads.add(tuple(sorted(r)))
                if len(quads) > 40: break
        if len(quads) > 40: break
        if len(pairE) > 200000: pairE.clear()
    return dict(A=feasible, nsig=len(sigs), nvec=nvec, T0=T0, TQ=TQ.tolist(), row_quads=sorted(quads)[:40])


def lag_classes(Gp, Sp_perms):
    """representatives of the <U, -1>-orbits of G' (0 first)"""
    seen = np.zeros(Gp.n, bool); reps = []
    negp = Gp.idx(-Gp.co)
    for g in range(Gp.n):
        if seen[g]: continue
        reps.append(g)
        for P in Sp_perms:
            seen[P[g]] = True; seen[negp[P[g]]] = True
    return reps


def screen_B(sl, cfg, row_quads, max_bits=24):
    """alpha-congruence screen (exact); row_quads from screen A (list of sorted |row| tuples)"""
    f0 = alpha_congruence_mod(sl.letters)
    if f0 < 2: return dict(B=None, reason="f0 < 2", f0=f0)
    mod, rmod = 4 * f0, 2 * f0
    L1 = len(sl.breps)
    if L1 > max_bits: return dict(B=None, reason=f"too many G'-orbits ({L1})", f0=f0)
    # alpha letters allowed per orbit: alpha with some letter = alpha mod 2f0
    allowed = [[a for a in (1, -1) if any((y - a) % rmod == 0 for y in L)] for L in sl.letters]
    Up = [sl.Gp.act(u) for u in sl.Sp]
    reps = lag_classes(sl.Gp, Up)
    # orbits with a single allowed alpha are fixed; enumerate the free ones only
    if any(len(al) == 0 for al in allowed): return dict(B=False, reason="an orbit admits no alpha", f0=f0)
    free = [j for j in range(L1) if len(allowed[j]) == 2]
    n = 1 << len(free)
    vals = np.empty((n, L1), np.int8)
    for j in range(L1):
        if len(allowed[j]) == 1: vals[:, j] = allowed[j][0]
    idx = np.arange(n, dtype=np.int64)
    for t, j in enumerate(free): vals[:, j] = (2 * ((idx >> t) & 1) - 1).astype(np.int8)
    Ps = []; Rs = []
    for c0 in range(0, len(vals), 1 << 15):
        A = vals[c0:c0 + (1 << 15)].astype(np.int64)[:, sl.borb]            # (chunk, |G'|)
        Ps.append(sl.Gp.paf(A)[:, reps] % mod); Rs.append(A.sum(1))
    P = np.concatenate(Ps) if Ps else np.zeros((0, len(reps)), np.int64)
    R = np.concatenate(Rs) if Rs else np.zeros(0, np.int64)
    T = sl.T1[reps] % mod
    # row residue options from admissible row quadruples: each block's row = +-r for r in the quad
    if not row_quads: return dict(B=False, reason="no row quadruple from A", f0=f0)
    # keys per block: (row mod rmod, P)
    key = np.concatenate([(R % rmod)[:, None], P], 1)
    uk, cnt = np.unique(key, axis=0, return_counts=True)
    # 2+2 join: for each row quadruple and sign pattern, need rows of alpha = +-r_k mod rmod
    left = defaultdict(list)
    for i in range(len(uk)):
        for j in range(len(uk)):
            sP = tuple(((uk[i, 1:] + uk[j, 1:]) % mod).tolist())
            left[sP].append((int(uk[i, 0]), int(uk[j, 0])))
    feasible_rows = set()
    for rq in row_quads:
        for perm in set(itertools.permutations(rq)):
            for sg in itertools.product((1, -1), repeat=4):
                feasible_rows.add(tuple((s_ * r) % rmod for s_, r in zip(sg, perm)))
    hit = None
    for sP, lst in left.items():
        need = tuple(((T - np.array(sP)) % mod).tolist())
        if need in left:
            rl = left[need]
            ls = set(lst); rs = set(rl)
            for (r1, r2) in ls:
                for (r3, r4) in rs:
                    if (r1, r2, r3, r4) in feasible_rows: hit = (sP, (r1, r2, r3, r4)); break
                if hit: break
        if hit: break
    return dict(B=hit is not None, f0=f0, n_alpha=int(len(vals)), n_keys=int(len(uk)), lag_classes=len(reps))


def run(cfg, do_B=True):
    t0 = time.time()
    sl = build(cfg)
    rec = dict(v=cfg['v'], c=cfg['c'], b=cfg['b'], Q=cfg['Q'], H1=cfg['H1'], q=cfg['q'], s=cfg['s'], f=cfg['f'],
               U=[list(g) if g[0] == 'neg' else [g[0], g[1]] for g in cfg['U']],
               Ugen=[None if g[0] == 'neg' else g[2] for g in cfg['U']], mods=cfg['mods'],
               orbits=int(sl.r), Gp=len(sl.breps), Gp_n=sl.Gp.n, sGeff=cfg['s'] * len(sl.breps),
               letters=sorted({tuple(L) for L in sl.letters}, key=len)[0][:12])
    A = screen_A(sl, cfg)
    rec.update({k: A[k] for k in A if k in ('A', 'reason', 'nsig', 'row_quads')})
    if A.get('A') and do_B:
        B = screen_B(sl, cfg, A['row_quads'])
        rec.update({('B_' + k if k != 'B' else 'B'): B[k] for k in B})
    rec['secs'] = round(time.time() - t0, 2)
    return rec


if __name__ == "__main__":
    a = sys.argv[1:]
    vs = []; cs = (4, 6); out = None; maxGeff = 10 ** 9; b0 = False; maxU = 1
    i = 0
    while i < len(a):
        if a[i] == '--v': vs.append(int(a[i + 1])); i += 2
        elif a[i] == '--c': cs = tuple(int(x) for x in a[i + 1].split(',')); i += 2
        elif a[i] == '--out': out = a[i + 1]; i += 2
        elif a[i] == '--maxGeff': maxGeff = int(a[i + 1]); i += 2
        elif a[i] == '--maxU': maxU = int(a[i + 1]); i += 2
        elif a[i] == '--b0': b0 = True; i += 1
        else: raise SystemExit("bad arg " + a[i])
    import resource
    memgb = float(os.environ.get('AB_MEMGB', '3.5'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30), int(memgb * 2 ** 30)))
    def ckey(d):
        U = [list(g) if g[0] == 'neg' else [g[0], g[1]] for g in d['U']] if d.get('U') and not isinstance(d['U'][0], list) else d.get('U')
        return json.dumps([d['v'], d['c'], list(d['Q']), list(d['H1']), d['q'], d['s'],
                           [[x if not isinstance(x, tuple) else list(x) for x in g] for g in (U or [])]])
    done = set()
    if out and os.path.exists(out):
        for l in open(out):
            try: d = json.loads(l)
            except Exception: continue
            if 'cfgkey' in d: done.add(d['cfgkey'])
            elif 'U' in d and 'Q' in d: done.add(ckey(d))
    fo = open(out, 'a') if out else None
    for v in vs:
        cfgs = configs_b0(v, maxU) if b0 else configs_bge3(v, cs, maxU)
        for cfg in cfgs:
            key = ckey(dict(cfg, U=[list(g) if g[0] == 'neg' else [g[0], g[1]] for g in cfg['U']]))
            if key in done: continue
            try:
                sl_n = None
                rec = run(cfg)
            except (Exception, MemoryError) as e:
                rec = dict(v=v, c=cfg['c'], Q=cfg['Q'], H1=cfg['H1'], q=cfg['q'], s=cfg['s'],
                           U=[list(g) if g[0] == 'neg' else [g[0], g[1]] for g in cfg['U']], error=repr(e)[:200])
            rec['cfgkey'] = key
            if rec.get('Gp', 0) * rec.get('s', 1) > maxGeff * 10: pass
            s = json.dumps(rec)
            print(s, flush=True)
            if fo: fo.write(s + "\n"); fo.flush()
