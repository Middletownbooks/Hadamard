r"""Open-order screen, step 1-2: enumerate staged GP slices at the 85 open orders below 10000 and run the exact
cheap screens [P, M].

Families (report/prompt_open_screen.md section 2), all as theory/ab_lib.Slice objects:
  b0  Goethals-Seidel core on an abelian G of order v (v odd composite; G = G' x Z_q);
  b1  Legendre quadruple / bordered GS on G of order v - 1 (G = G' x Z_q);
  k4  EKLO Theorem 7 / Kumar's bordered family with t = 3: G of order v - 3 = 4h, v = 3 (mod 4), K = H of index 4,
      split G = Q x H with Q = Z2^2 or Z4 (theory/ab_sweep.configs_bge3 with c = 4).
Slice: S = U x M, M the index-s subgroup of F_q^* on the staging axis Z_q (q a prime with q || |G|), U trivial or one
cyclic subgroup of units on one axis of G', or -1 on several axes (the menu of theory/ab_sweep.py, by EKLO/Kumar
and the repository's earlier sweeps).  Orbits per block r = n_B (1 + s), n_B = number of U-orbits on G' [P]
(M acts freely on Z_q \ 0, independently of U).  Kept: 40 <= r <= 160.

Screens (exact, cheapest first; all necessary conditions):
  A  lag-0 energy + row sums (+ Q-level for k4) of stage 1 (ab_sweep.screen_A; for b = 1 the row condition is
     sum_k r_k^2 = sum T1 = 4, passed to screen_A as v = 1, which uses v only in that test);
  B  alpha-congruence (ab_sweep.screen_B).
Extension of ab_sweep (minimal): configs_b1 (b = 1, any abelian G' of order (v - 1)/q), the orbit-count
prefilter, and restart-safe JSONL output.  Diagonal multipliers acting on several axes with one generator (as at
647) are not enumerated.
usage: python3 theory/open_screen_enum.py --list [--fam b0,b1,k4]         (count configurations)
       python3 theory/open_screen_enum.py --run OUT.jsonl [--shard i/n] [--fam ...] [--secs 300]"""
import sys, os, json, math, time, itertools, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import factorint, AG
import ab_sweep as SW

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN = json.load(open(os.path.join(ROOT, 'data/open_orders_lt10000.json')))['open']
RMIN, RMAX = 24, 160


def umenu(axes, maxU=1):
    Umenus = [[]]
    axis_opts = [(ai, o, u) for ai, n in enumerate(axes) for (o, u) in SW.unit_subgroup_gens(n)]
    for k in range(1, maxU + 1):
        for combo in itertools.combinations(axis_opts, k):
            if len({a for a, _, _ in combo}) < k: continue
            Umenus.append(list(combo))
    big = [ai for ai, n in enumerate(axes) if n > 2]
    for k in range(2, len(big) + 1):
        for A in itertools.combinations(big, k): Umenus.append([('neg', A)])
    return Umenus


def configs_b1(v, maxU=1):
    out = []; m = v - 1
    for q, e in factorint(m).items():
        if e != 1 or q == 2: continue
        F = m // q
        for Ft in SW.abelian_types(F):
            mods = tuple(Ft) + (q,)
            for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                f = (q - 1) // s
                if f < 2: continue
                for U in umenu(list(Ft), maxU):
                    out.append(dict(v=v, c=1, b=1, h=None, m=m, Q=(), H1=Ft, q=q, s=s, f=f, mods=mods, nQ=0, U=U))
    return out


def nB(cfg):
    """number of U-orbits on G' (U as in ab_sweep.build)"""
    axes = list(cfg['mods'][:-1])
    if not axes: return 1
    G = AG(axes); gens = []
    for gen in cfg['U']:
        g = [1] * len(axes)
        if gen[0] == 'neg':
            for a in gen[1]: g[a] = axes[a] - 1
        else:
            a, o, u = gen; g[a] = u
        gens.append(g)
    if not gens: return G.n
    seen = np.zeros(G.n, bool); n = 0
    perms = []
    for g in gens: perms.append(G.act(g))
    for x in range(G.n):
        if seen[x]: continue
        n += 1; st = [x]; seen[x] = True
        while st:
            y = st.pop()
            for P in perms:
                z = int(P[y])
                if not seen[z]: seen[z] = True; st.append(z)
    return n


def all_configs(fams=('b0', 'b1', 'k4')):
    out = []
    for N in OPEN:
        v = N // 4
        cs = []
        if 'b0' in fams and v % 2 == 1: cs += SW.configs_b0(v)
        if 'b1' in fams: cs += configs_b1(v)
        if 'k4' in fams and v % 4 == 3: cs += SW.configs_bge3(v, cs=(4,))
        for c in cs:
            n_b = nB(c); r = n_b * (1 + c['s'])
            if RMIN <= r <= RMAX:
                c['nB'] = n_b; c['r'] = r; c['sGeff'] = n_b * c['s']; c['fam'] = {0: 'b0', 1: 'b1', 3: 'k4'}[c['b']]
                out.append(c)
    return out


def _unique_order(n, o):
    """True if (Z/n)^* has exactly one cyclic subgroup of order o"""
    return sum(1 for (oo, _) in SW.unit_subgroup_gens(n) if oo == o) == 1


def ukey(U, axes=None):
    """U as stored in the key: ['neg', axes] or [axis, order], plus the generator when the cyclic subgroup of that
    order on that axis is not unique (2-power axes of order >= 8)"""
    out = []
    for g in U:
        if g[0] == 'neg': out.append(list(g))
        elif axes is not None and not _unique_order(axes[g[0]], g[1]): out.append([g[0], g[1], g[2]])
        else: out.append([g[0], g[1]])
    return out


def cfgkey(c):
    axes = list(c['Q']) + list(c['H1'])
    return json.dumps([c['fam'], c['v'], list(c['Q']), list(c['H1']), c['q'], c['s'],
                       [[x if not isinstance(x, tuple) else list(x) for x in g] for g in ukey(c['U'], axes)]])


def screen_B_fast(sl, cfg, row_quads, max_bits=22, max_keys=4000):
    """ab_sweep.screen_B (alpha-congruence, exact) with a vectorised 2+2 join over the distinct keys; B = None when
    the distinct-key count exceeds max_keys (a cap, not a zero)"""
    from ab_lib import alpha_congruence_mod
    f0 = alpha_congruence_mod(sl.letters)
    if f0 < 2: return dict(B=None, reason="f0 < 2", f0=f0)
    mod, rmod = 4 * f0, 2 * f0
    L1 = len(sl.breps)
    allowed = [[a for a in (1, -1) if any((y - a) % rmod == 0 for y in L)] for L in sl.letters]
    if any(len(al) == 0 for al in allowed): return dict(B=False, reason="an orbit admits no alpha", f0=f0)
    free = [j for j in range(L1) if len(allowed[j]) == 2]
    if len(free) > max_bits: return dict(B=None, reason=f"too many free alpha orbits ({len(free)})", f0=f0)
    Up = [sl.Gp.act(u) for u in sl.Sp]
    reps = SW.lag_classes(sl.Gp, Up)
    n = 1 << len(free)
    vals = np.empty((n, L1), np.int8)
    for j in range(L1):
        if len(allowed[j]) == 1: vals[:, j] = allowed[j][0]
    idx = np.arange(n, dtype=np.int64)
    for t, j in enumerate(free): vals[:, j] = (2 * ((idx >> t) & 1) - 1).astype(np.int8)
    keys = []
    for c0 in range(0, n, 1 << 14):
        A = vals[c0:c0 + (1 << 14)].astype(np.int64)[:, sl.borb]
        P = sl.Gp.paf(A)[:, reps] % mod
        keys.append(np.concatenate([(A.sum(1) % rmod)[:, None], P], 1))
    uk = np.unique(np.concatenate(keys), axis=0)
    if not row_quads: return dict(B=False, reason="no row quadruple from A", f0=f0)
    if len(uk) > max_keys: return dict(B=None, reason=f"key cap ({len(uk)} keys)", f0=f0, n_keys=int(len(uk)))
    T = sl.T1[reps] % mod
    feas = set()
    for rq in row_quads:
        for perm in set(itertools.permutations(rq)):
            for sg in itertools.product((1, -1), repeat=4):
                feas.add(tuple((s_ * r) % rmod for s_, r in zip(sg, perm)))
    nk = len(uk); d = uk.shape[1] - 1
    I, J = np.meshgrid(np.arange(nk), np.arange(nk), indexing='ij'); I = I.ravel(); J = J.ravel()
    Ps = (uk[I, 1:] + uk[J, 1:]) % mod
    base = np.array([mod ** i for i in range(d)], dtype=np.int64) if mod ** d < 2 ** 62 else None
    if base is None:
        rng = np.random.default_rng(5); base = rng.integers(1, 2 ** 40, size=d)      # hash; exact re-check below
    enc = Ps.astype(np.int64) @ base
    need = ((T[None, :] - Ps) % mod).astype(np.int64) @ base
    rr = uk[I, 0] * rmod + uk[J, 0]
    # right side: (enc, row pair) ; left queries: (need, allowed right row pair)
    rkey = enc * (rmod * rmod) + rr
    rset = np.unique(rkey)
    byleft = {}
    for (r1, r2, r3, r4) in feas: byleft.setdefault(r1 * rmod + r2, []).append(r3 * rmod + r4)
    hit = False
    for lr, rights in byleft.items():
        sel = np.nonzero(rr == lr)[0]
        if not len(sel): continue
        for r34 in rights:
            q = need[sel] * (rmod * rmod) + r34
            pos = np.searchsorted(rset, q); pos[pos >= len(rset)] = 0
            ok = np.nonzero(rset[pos] == q)[0]
            for o in ok[:50]:                            # exact re-check of the candidate
                p = sel[o]; tgt = (T - Ps[p]) % mod
                m2 = np.nonzero((rr == r34) & (Ps == tgt).all(1))[0]
                if len(m2): hit = True; break
            if hit: break
        if hit: break
    return dict(B=hit, f0=f0, n_alpha=int(n), n_keys=int(nk), lag_classes=len(reps))


class TO(Exception): pass


def _alarm(*a): raise TO()


def run_one(c, secs):
    t0 = time.time()
    rec = dict(key=cfgkey(c), fam=c['fam'], v=c['v'], N=4 * c['v'], Q=list(c['Q']), H1=list(c['H1']), q=c['q'],
               s=c['s'], f=c['f'], U=ukey(c['U'], list(c['Q']) + list(c['H1'])), Ugen=[None if g[0] == 'neg' else g[2] for g in c['U']],
               r=c['r'], nB=c['nB'], sGeff=c['sGeff'])
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(secs)
    try:
        sl = SW.build(c)
        assert sl.r == c['r'], (sl.r, c['r'])
        cA = dict(c, v=1) if c['b'] == 1 else c          # b = 1: rows sum r_k^2 = 4 (see module doc)
        A = SW.screen_A(sl, cA)
        rec.update({('A_' + k if k != 'A' else 'A'): A[k] for k in A if k in ('A', 'reason', 'nsig', 'row_quads')})
        rec['letters'] = sorted({tuple(L) for L in sl.letters}, key=len)[0][:12]
        if A.get('A'):
            B = screen_B_fast(sl, cA, A["row_quads"])
            rec.update({('B_' + k if k != 'B' else 'B'): B[k] for k in B})
        signal.alarm(0)
    except TO:
        rec['timeout'] = secs
    except (Exception, MemoryError) as e:
        signal.alarm(0); rec['error'] = repr(e)[:200]
    rec['secs'] = round(time.time() - t0, 2)
    return rec


if __name__ == "__main__":
    a = sys.argv[1:]
    fams = tuple(a[a.index('--fam') + 1].split(',')) if '--fam' in a else ('b0', 'b1', 'k4')
    cfgs = all_configs(fams)
    if '--list' in a:
        from collections import Counter
        print(len(cfgs), Counter(c['fam'] for c in cfgs))
        print(Counter((c['fam'], c['sGeff'] <= 80) for c in cfgs))
        raise SystemExit
    out = a[a.index('--run') + 1]
    secs = int(a[a.index('--secs') + 1]) if '--secs' in a else 300
    si, sn = (int(x) for x in a[a.index('--shard') + 1].split('/')) if '--shard' in a else (0, 1)
    import resource
    memgb = float(os.environ.get('OS_MEMGB', '3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30), int(memgb * 2 ** 30)))
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: done.add(json.loads(l)['key'])
            except Exception: pass
    # cheapest first: small sGeff first
    cfgs.sort(key=lambda c: (c['sGeff'], c['r']))
    fo = open(out, 'a')
    for i, c in enumerate(cfgs):
        if i % sn != si or cfgkey(c) in done: continue
        rec = run_one(c, secs)
        fo.write(json.dumps(rec) + "\n"); fo.flush()
