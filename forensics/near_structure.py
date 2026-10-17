"""Near-structure ("template + deviations") census of the published cores, with our own cores as controls.

The two kinds of hidden structure found so far in the published cores are near-structures: the half-shift arcs of
167/179 (an exact four-arc template is impossible by parity; 5 exceptions per core) and the Legendre(11) twist of
359/419 (116/128 and 48/56 pairs).  Our own methods only search exact slices.  This script looks for further
near-structures in every core, in its own group and about the common fixed point of its stabiliser U:

 (a) multiplier classes: for every class a mod U (any order), the number of U-orbits O with x(aO) = x(O), summed over
     the four blocks; a high count is near-invariance under <U, a>, a low one near-anti-invariance (a twist);
 (b) half-shift: for every translation tau of order 2, sigma(g) = x(g) x(g + tau), and its boundary (sign changes
     along every cyclic axis of length >= 3 of the core's own coordinates); a low boundary is an arc template;
 (c) b = 0 only, the multiplicative analogue: on the cyclotomic classes j = log_g(y) mod e (e = (v-1)/|U|), the
     boundary of f(j) = x(g^j) and of sigma_s(j) = f(j) f(j + s) for every shift s;
 (d) cross-block multiplier relations: agreement of x_k(aO) with x_j(O) for every ordered pair j != k and every class
     a (including the identity);
 (e) b = 0 only: the T-product x_A x_B x_C x_D = +1 and agreement with the Legendre symbol mod v.

Null for every statistic: the four blocks' orbit labels permuted independently within orbit-size classes (so U-
invariance and the label multiset are kept), P permutations; z = (observed - mean) / sd.  'z*' below is the most
extreme z of a family, and 'fam' the number of tests in it (a Bonferroni scale: |z| = 4.0 is p ~ 6e-5 per test).
usage: python3 forensics/near_structure.py [--P 400] [--only 283,487] [--json OUT]"""
import os, sys, re, glob, json, math, gzip
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from census_originals import G, natural_core, stabiliser
from gbuilder import read_seq
from load_originals import originals, load_any, header


def parse_seqs(path):
    """four +-1 sequences from a text file: lines 'A=+-+...' (or 'A: ...', or bare +- lines); '#' lines skipped"""
    op = gzip.open if path.endswith('.gz') else open
    out = []
    for l in op(path, 'rt'):
        l = l.strip()
        if not l or l.startswith('#'): continue
        mm = re.match(r'^([ABCD])?\s*[=:]?\s*([+-]{20,})$', l)
        if mm: out.append(np.array([1 if c == '+' else -1 for c in mm.group(2)], np.int64))
    assert len(out) >= 4, path
    return np.array(out[:4])


def paf_sum(X, g):
    F = np.fft.fftn(X.reshape((4,) + tuple(g.mods)), axes=tuple(range(1, len(g.mods) + 1)))
    P = np.rint(np.fft.ifftn((np.abs(F) ** 2).sum(0)).real).astype(np.int64)
    return P.reshape(-1)


def crt_cyclic(Y, mods):
    """sequence on Z_a x Z_b (entry (g, z) at g*b + z) -> cyclic Z_ab by CRT, n -> (n mod a, n mod b)"""
    a, b = mods; n = np.arange(a * b)
    return np.array([y[(n % a) * b + (n % b)] for y in Y])


CONTROLS = [  # (label, v, b, path, mods, loader)
    ('ours 487 k-tree', 487, 0, 'data/core_v487_b0_kt.txt', (487,), None),
    ('ours 883 k-tree', 883, 0, 'data/core_v883_b0_d21_0.txt', (883,), None),
    ('ours 919 twist', 919, 0, 'data/core_v919_b0_twist51_0.txt', (919,), None),
    ('ours 933 GP', 933, 0, 'pub/H3732/sequences_Z933.txt', (933,), None),
    ('Djokovic 463', 463, 0, 'data/core_v463_b0.txt', (463,), None),
    ('ours 419 Paley', 419, 1, 'data/core_v419_b1_paley_guided_K15.txt', (418,), 'crt38x11'),
    ('ours 523 quartic', 523, 1, 'data/core_v523_b1_q523_p523f_w0.txt', (522,), None),
    ('ours 573 Paley', 573, 1, 'data/core_v573_b1.txt', (572,), None),
    ('ours 647 doubly', 647, 1, 'pub/H2588/sequences_Z646.txt', (646,), None),
    ('ours 347 staged', 347, 5, 'logs/ab/prod347/core_a1_438c58a201c4084b_1.npy', None, 'npy347'),
]


def target_ok(v, b, P, g):
    """sum_k PAF_k: b = 0: 4v at 0, 0 elsewhere; b = 1: 4m at 0, -4 elsewhere; b >= 2: two-valued -4b / +4"""
    off = set(np.unique(P[1:]).tolist())
    if b == 0: return off == {0}
    if b == 1: return off == {-4}
    return off <= {-4 * b, 4}


def load_all(only=None):
    cores = []
    srcs = [('m', n) for n in originals()] + [('s', f) for f in sorted(glob.glob(os.path.join(HERE, 'originals', 'seq_v*.txt')))]
    for kind, f in srcs:
        if kind == 's':
            hdr, seqs = read_seq(f); H = None; meta = dict(re.findall(r'(\w+)=(\w+)', hdr))
        else:
            H = load_any(f); seqs = None; meta = dict(re.findall(r'(\w+)=(\w+)', header(f)))
        v, b = int(meta['v']), int(meta['b'])
        if only and v not in only: continue
        X, g = natural_core(v, b, meta['func'], H, seqs)
        cores.append(dict(label=f'Alpoge {v}', v=v, b=b, X=X, g=g, ours=False))
    for label, v, b, path, mods, how in CONTROLS:
        if only and v not in only: continue
        p = os.path.join(ROOT, path)
        if how == 'npy347':
            X = np.load(p).astype(np.int64)
            for mods_try in [(2, 3, 3, 19), (2, 3, 57), (2, 57, 3), (3, 2, 57)]:
                g = G(mods_try)
                if target_ok(v, b, paf_sum(X, g), g): break
            else: raise SystemExit('347 control: layout not recognised')
        else:
            X = parse_seqs(p)
            if how == 'crt38x11': X = crt_cyclic(X, (38, 11))
            g = G(mods)
        cores.append(dict(label=label, v=v, b=b, X=X, g=g, ours=True))
    for c in cores:
        P = paf_sum(c['X'], c['g'])
        c['paf_ok'] = bool(target_ok(c['v'], c['b'], P, c['g']))
    return cores


class Ctx:
    """a core with its stabiliser U, orbits about the fixed point (moved to 0) and the permutation null"""
    def __init__(s, core, P=400, seed=1):
        s.core = core; X, g = core['X'], core['g']
        d, c, U = stabiliser(X, g)
        sh = g.idx(g.coords + g.coords[c]); s.X = X[:, sh]; s.g = g; s.U = U; s.d = d
        s.e = int(np.lcm.reduce(g.mods)); m = g.m
        orb = -np.ones(m, np.int64); reps = []
        for y in range(m):
            if orb[y] < 0:
                o = sorted(set(int(t) for t in g.idx(np.array([(a * g.coords[y]) for a in U]))))
                orb[o] = len(reps); reps.append(o)
        s.orb = orb; s.reps = reps; s.r = len(reps)
        s.size = np.array([len(o) for o in reps]); s.rep = np.array([o[0] for o in reps])
        s.L = s.X[:, s.rep]                                                  # (4, r) orbit labels
        rng = np.random.default_rng(seed); s.Pn = P
        perm = np.empty((4, P, s.r), np.int64)
        for k in range(4):
            for p in range(P):
                lab = s.L[k].copy()
                for sz in np.unique(s.size):
                    ix = np.nonzero(s.size == sz)[0]; lab[ix] = lab[ix][rng.permutation(len(ix))]
                perm[k, p] = lab
        s.perm = perm                                                        # (4, P, r) null labels
        units = [a for a in range(1, s.e) if gcd(a, s.e) == 1]
        Uset = set(U); seen = set(); s.classes = []
        for a in units:
            key = min((a * u) % s.e for u in U)
            if key in seen: continue
            seen.add(key); s.classes.append(key)
        # orbit image under each class (orbits are permuted by multipliers since U is abelian)
        s.img = {a: s.orb[g.idx(a * g.coords[s.rep])] for a in s.classes}

    def order_mod_U(s, a):
        k, x = 1, a % s.e; Us = set(s.U)
        while x not in Us: x = (x * a) % s.e; k += 1
        return k

    def zstat(s, obs, null):
        mu, sd = float(null.mean()), float(null.std())
        return (obs - mu) / sd if sd > 0 else 0.0, mu, sd


def test_classes(c):
    """(a): sum over blocks of #{O : x(aO) = x(O), aO != O}"""
    out = []
    for a in c.classes:
        if a == 1: continue
        im = c.img[a]; mv = im != np.arange(c.r)
        if not mv.any(): continue
        obs = int(sum((c.L[k][mv] == c.L[k][im[mv]]).sum() for k in range(4)))
        null = sum((c.perm[k][:, mv] == c.perm[k][:, im[mv]]).sum(1) for k in range(4))
        z, mu, sd = c.zstat(obs, null)
        per = [int((c.L[k][mv] == c.L[k][im[mv]]).sum()) for k in range(4)]
        out.append(dict(a=int(a), ord=c.order_mod_U(a), z=round(z, 2), same=obs, null_mean=round(mu, 1),
                        n=int(4 * mv.sum()), per_block=per))
    out.sort(key=lambda t: -abs(t['z']))
    return out


def boundary(sig, g):
    """sign changes of sig along every cyclic axis of length >= 3 (sig on the group, flat index)"""
    S = sig.reshape(tuple(g.mods)); tot = 0
    for ax, n in enumerate(g.mods):
        if n >= 3: tot += int((S != np.roll(S, -1, axis=ax)).sum())
    return tot


def test_halfshift(c):
    """(b): every translation tau of order 2: boundary of sigma = x(g) x(g + tau)"""
    g = c.g; out = []
    taus = [t for t in range(1, g.m) if (2 * g.coords[t] % g.mods == 0).all()]
    for t in taus:
        sh = g.idx(g.coords + g.coords[t])
        obs = sum(boundary(c.X[k] * c.X[k][sh], g) for k in range(4))
        null = np.zeros(c.Pn)
        for p in range(min(c.Pn, 200)):
            for k in range(4):
                xk = c.perm[k, p][c.orb]; null[p] += boundary(xk * xk[sh], g)
        null = null[:min(c.Pn, 200)]
        z, mu, sd = c.zstat(obs, null)
        per = [boundary(c.X[k] * c.X[k][sh], g) for k in range(4)]
        out.append(dict(tau=[int(x) for x in g.coords[t]], z=round(z, 2), boundary=int(obs), null_mean=round(mu, 1),
                        per_block=per))
    out.sort(key=lambda t: t['z'])
    return out


def primitive_root(p):
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    for g0 in range(2, p):
        if all(pow(g0, (p - 1) // q, p) != 1 for q in fac): return g0


def test_logshift(c):
    """(c), b = 0 (prime v): f(j) = x(g^j) on Z_e; boundary of f and of sigma_s = f(j) f(j+s)"""
    v = c.g.m; g0 = primitive_root(v); e = (v - 1) // c.d
    pw = np.array([pow(g0, j, v) for j in range(e)])                   # g^j, j < e: one per class
    oi = c.orb[pw]                                                      # orbit index of class j
    F = c.L[:, oi]                                                      # (4, e)
    Fn = c.perm[:, :, oi]                                               # (4, P, e)
    res = {}
    obs = int((F != np.roll(F, -1, axis=1)).sum()); null = (Fn != np.roll(Fn, -1, axis=2)).sum((0, 2))
    z, mu, sd = c.zstat(obs, null); res['f_boundary'] = dict(z=round(z, 2), obs=obs, null_mean=round(mu, 1))
    best = []
    for s_ in range(1, e):
        S = F * np.roll(F, -s_, axis=1); Sn = Fn * np.roll(Fn, -s_, axis=2)
        ob = int((S != np.roll(S, -1, axis=1)).sum()); nl = (Sn != np.roll(Sn, -1, axis=2)).sum((0, 2))
        z, mu, sd = c.zstat(ob, nl); best.append(dict(s=s_, z=round(z, 2), obs=ob, null_mean=round(mu, 1)))
    best.sort(key=lambda t: t['z'])
    res['sigma_s_lowest'] = best[:3]; res['e'] = e; res['g'] = g0
    return res


def test_cross(c):
    """(d): agreement of x_k(aO) with x_j(O), every ordered pair j != k and every class a"""
    out = []
    for j in range(4):
        for k in range(4):
            if j == k: continue
            for a in c.classes:
                im = c.img[a]
                obs = int((c.L[k][im] == c.L[j]).sum())
                null = (c.perm[k][:, im] == c.L[j][None, :]).sum(1)
                z, mu, sd = c.zstat(obs, null)
                out.append(dict(pair='ABCD'[j] + '<-' + 'ABCD'[k], a=int(a), z=round(z, 2), agree=obs, n=c.r,
                                null_mean=round(mu, 1)))
    out.sort(key=lambda t: -abs(t['z']))
    return out


def test_b0_extra(c):
    """(e), b = 0: T-product and Legendre agreement, at the orbit level"""
    v = c.g.m; res = {}
    prod = c.L.prod(0); nprod = c.perm.prod(0)                          # (r,), (P, r)
    obs = int((prod == 1).sum()); z, mu, sd = c.zstat(obs, (nprod == 1).sum(1))
    res['T_product'] = dict(z=round(z, 2), obs=obs, r=c.r, null_mean=round(mu, 1))
    chi = np.array([0] + [1 if pow(int(y), (v - 1) // 2, v) == 1 else -1 for y in range(1, v)])
    if all(chi[o].std() == 0 for o in c.reps):                          # Legendre constant on U-orbits
        ch = chi[c.rep]; zz = []
        for k in range(4):
            obs = int((c.L[k] == ch).sum()); z, mu, sd = c.zstat(obs, (c.perm[k] == ch[None, :]).sum(1)); zz.append(round(z, 2))
        res['legendre_z_per_block'] = zz
    return res


def main():
    a = sys.argv[1:]
    P = int(a[a.index('--P') + 1]) if '--P' in a else 400
    only = set(int(x) for x in a[a.index('--only') + 1].split(',')) if '--only' in a else None
    outp = a[a.index('--json') + 1] if '--json' in a else None
    allres = []
    for core in load_all(only):
        c = Ctx(core, P=P)
        R = dict(label=core['label'], v=core['v'], b=core['b'], mods=[int(x) for x in core['g'].mods], paf_ok=core['paf_ok'],
                 U_order=c.d, orbits=c.r, classes=len(c.classes) - 1)
        R['classes_top'] = test_classes(c)[:4]
        R['cross_top'] = test_cross(c)[:4]; R['cross_family'] = 12 * len(c.classes)
        if (np.array(core['g'].mods) % 2 == 0).any(): R['halfshift'] = test_halfshift(c)
        if core['b'] == 0 and len(core['g'].mods) == 1 and all(core['g'].m % q for q in range(2, int(core['g'].m ** 0.5) + 1)):
            R['logshift'] = test_logshift(c); R['b0'] = test_b0_extra(c)
        allres.append(R)
        line = f"{R['label']:>18} v={R['v']} b={R['b']} G=Z{'xZ'.join(map(str, R['mods']))} |U|={c.d} r={c.r} paf_ok={R['paf_ok']}"
        t = R['classes_top'][0] if R['classes_top'] else None
        if t: line += f" | (a) z*={t['z']} a={t['a']} ord={t['ord']} same {t['same']}/{t['n']} (null {t['null_mean']}) fam={R['classes']}"
        x = R['cross_top'][0]; line += f" | (d) z*={x['z']} {x['pair']} a={x['a']} fam={R['cross_family']}"
        if 'halfshift' in R:
            h = R['halfshift'][0]; line += f" | (b) z*={h['z']} tau={h['tau']} bnd {h['boundary']} (null {h['null_mean']}) fam={len(R['halfshift'])}"
        if 'logshift' in R:
            L_ = R['logshift']; s0 = L_['sigma_s_lowest'][0]
            line += f" | (c) f z={L_['f_boundary']['z']}, sigma_s z*={s0['z']} s={s0['s']} fam={L_['e'] - 1}"
            line += f" | (e) T z={R['b0']['T_product']['z']} Leg z={R['b0'].get('legendre_z_per_block')}"
        print(line, flush=True)
    if outp: json.dump(allres, open(outp, 'w'), indent=1)


if __name__ == '__main__':
    main()
