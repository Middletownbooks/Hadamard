r"""Open-order screen, step 5: exact stage-1 decision of a staged GP slice by two methods, with planted controls [M].

Stage 1 of a slice (theory/ab_lib.Slice; S = U x M, U on the axes of G', M of index s on Z_q):
    y_k(b) = sum_z x_k(b, z),   sum_k PAF_{G'}(y_k) = T1,   T1(b) = sum_z T(b, z),
with y_k U-invariant and, at each b, y_k(b) = alpha + f*sigma, alpha = +-1, sigma a sum of s signs [P: the S-orbit of
(b, 0) meets the fibre of b in one point and each S-orbit of (b, z), z != 0, in f points, since U fixes the Z_q
coordinate and M acts freely on Z_q \ 0].

Method A: theory/ab_b0_stage1.stage1 (the repository's enumerator: every block within the energy bound, a 2+2 join
          on the exact PAF vector over G').
Method B (this file, shares no code with ab_lib/ab_b0_stage1): G', the U-orbits, the letters and T1 are derived from
          the definitions above; blocks are enumerated by a DFS over U-orbits with energy pruning; PAF by the
          definition (double sum over G'); the join is a dictionary of pair sums of PAF tuples split by pair energy.
Planted: a random U-invariant y-quadruple from the letters; the target T1 is replaced by its PAF sum; both methods
          must return it.  A method that hits a size cap reports UNDECIDED (never zero).
usage: python3 theory/open_screen_stage1_exact.py 'KEY' [--planted 4] [--cap 2e7] [--out FILE]"""
import sys, os, json, time, itertools, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- method B (independent)
def mk_group(axes):
    return list(itertools.product(*[range(n) for n in axes]))


def b_setup(key):
    fam, v, Q, H1, q, s, U = key
    b = {'b0': 0, 'b1': 1, 'k4': 3}[fam]
    axes = list(Q) + list(H1); f = (q - 1) // s
    els = mk_group(axes); idx = {e: i for i, e in enumerate(els)}
    # U generators: per-axis scalars (cyclic subgroup of the given order on one axis) or -1 on axes
    gens = []
    for g in U:
        mult = [1] * len(axes)
        if g[0] == 'neg':
            for a in g[1]: mult[a] = axes[a] - 1
        else:
            a, o = g; n = axes[a]
            u = next(u for u in range(2, n) if math.gcd(u, n) == 1 and _ordmod(u, n) == o and
                     _same_subgroup(u, n, o))
            mult[a] = u
        gens.append(mult)
    act = lambda m, e: tuple((x * c) % n for x, c, n in zip(e, m, axes))
    orb = [-1] * len(els); orbits = []
    for i, e in enumerate(els):
        if orb[i] >= 0: continue
        st = [e]; seen = {e}
        while st:
            x = st.pop()
            for m in gens:
                y = act(m, x)
                if y not in seen: seen.add(y); st.append(y)
        for x in seen: orb[idx[x]] = len(orbits)
        orbits.append(sorted(idx[x] for x in seen))
    sig = [s - 2 * t for t in range(s + 1)]
    letters = sorted({a + f * sg for a in (1, -1) for sg in sig})
    # T over G = G' x Z_q, summed over z: T1
    n = len(els); m = n * q
    T1 = [0] * n
    nQ = len(Q)
    for i, e in enumerate(els):
        for z in range(q):
            if b == 0: t = 4 * v if (i == 0 and z == 0) else 0
            elif b == 1: t = 4 * m if (i == 0 and z == 0) else -4
            else:
                inH = all(x == 0 for x in e[:nQ])
                t = 4 * m if (i == 0 and z == 0) else (-4 * b if inH else 4)
            T1[i] += t
    sub = lambda a, c: tuple((x - y) % nn for x, y, nn in zip(a, c, axes))
    lagmap = [[idx[tuple((x + y) % nn for x, y, nn in zip(e, h, axes))] for e in els] for h in els]
    return dict(axes=axes, els=els, orbits=orbits, orb=orb, letters=letters, T1=T1, lagmap=lagmap, n=n)


def _ordmod(u, n):
    k, x = 1, u % n
    while x != 1: x = x * u % n; k += 1
    return k


def _same_subgroup(u, n, o):
    return True                      # any element of order o generates the (unique, n a prime power or 2) cyclic subgroup
                                     # of that order when (Z/n)^* is cyclic; checked against ab_lib below


def b_paf(st, y):
    lm = st['lagmap']
    return tuple(sum(y[g] * y[lm[h][g]] for g in range(st['n'])) for h in range(st['n']))


def b_blocks(st, Emax, cap):
    orbits = st['orbits']; L = st['letters']; w = [len(o) for o in orbits]
    order = sorted(range(len(orbits)), key=lambda j: -w[j])
    out = []
    vals = [0] * len(orbits)
    def dfs(t, E):
        if len(out) > cap: raise OverflowError
        if t == len(order):
            y = [0] * st['n']
            for j, o in enumerate(orbits):
                for g in o: y[g] = vals[j]
            out.append((E, tuple(y))); return
        j = order[t]
        rest = sum(w[order[u]] for u in range(t + 1, len(order)))       # every letter has |y| >= 1
        for a in L:
            e2 = E + w[j] * a * a
            if e2 + rest > Emax: continue
            vals[j] = a; dfs(t + 1, e2)
    dfs(0, 0)
    return out


def b_solve(st, T1, cap=2 * 10 ** 7):
    minE = st['n']; T0 = T1[0]
    blocks = b_blocks(st, T0 - 3 * minE, cap)
    # normalise sign: a block and its negative have the same PAF; keep y[0] > 0
    keyed = defaultdict(list)
    for E, y in blocks:
        if y[0] < 0: continue
        keyed[b_paf(st, y)].append(y)
    keys = list(keyed)
    if len(keys) > 60000: raise OverflowError("too many distinct block PAFs for the pair join")
    pairs = defaultdict(list)
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            sm = tuple(a + b for a, b in zip(keys[i], keys[j]))
            pairs[sm].append((i, j))
    sols = set()
    for sm, lst in pairs.items():
        need = tuple(t - a for t, a in zip(T1, sm))
        if need not in pairs: continue
        for (i, j) in lst:
            for (k, l) in pairs[need]:
                for a in keyed[keys[i]]:
                    for b in keyed[keys[j]]:
                        for c in keyed[keys[k]]:
                            for d in keyed[keys[l]]:
                                sols.add(tuple(sorted((a, b, c, d))))
    for Y in sols:
        tot = [sum(b_paf(st, y)[h] for y in Y) for h in range(st['n'])]
        assert tot == list(T1)
    return sorted(sols), len(blocks)


# ---------------------------------------------------------------- method A (repository) and the driver
def a_solve(key, T1=None, cap=2 * 10 ** 7):
    import open_screen_deficit as D, ab_sweep as SW
    from ab_b0_stage1 import stage1
    sl = SW.build(D.cfg_from_key(key))
    if T1 is not None: sl.T1 = np.array(T1, np.int64)
    sols, nb = stage1(sl, cap=cap)
    return [tuple(sorted(tuple(int(x) for x in r) for r in Y)) for Y in sols], nb, sl


def planted_y(st, rng):
    """a random U-invariant quadruple from the letters with lag-0 energy at most the real T1(0) (so the planted
    enumeration has the size of the real one): start from +-1 letters, then raise random orbits to random larger
    letters while the energy stays within the bound"""
    T0 = st['T1'][0]; w = [len(o) for o in st['orbits']]; L = st['letters']
    vals = [[int(rng.choice([-1, 1])) for _ in w] for _ in range(4)]
    E = sum(w) * 4
    for t in range(200):
        k = int(rng.integers(4)); j = int(rng.integers(len(w))); a = int(rng.choice(L))
        e2 = E - w[j] * vals[k][j] ** 2 + w[j] * a * a
        if e2 <= T0: vals[k][j] = a; E = e2
    Y = []
    for k in range(4):
        y = [0] * st['n']
        for j, o in enumerate(st['orbits']):
            for g in o: y[g] = vals[k][j]
        Y.append(tuple(y))
    return Y


def canon(Y):
    """quadruple up to block order and block sign"""
    return tuple(sorted(tuple(y) if y[0] > 0 else tuple(-x for x in y) for y in Y))


def main():
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
    a = sys.argv[1:]
    key = json.loads(a[0]); npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    cap = int(float(a[a.index('--cap') + 1])) if '--cap' in a else 2 * 10 ** 7
    outp = a[a.index('--out') + 1] if '--out' in a else None
    rec = dict(key=key)
    st = b_setup(key)
    # consistency of the independent set-up with ab_lib (orbit partition of G', letters, T1)
    t0 = time.time()
    try:
        solsA, nbA, sl = a_solve(key, cap=cap); rec['A'] = dict(solutions=len(solsA), blocks=nbA, secs=round(time.time() - t0, 1))
    except Exception as e:
        rec['A'] = dict(undecided=repr(e)[:120]); sl = None
    if sl is not None:
        rec['setup_agrees'] = (sorted(sorted(o) for o in st['orbits']) == sorted(sorted(np.nonzero(sl.borb == j)[0].tolist()) for j in range(len(sl.breps)))
                               and set(st['letters']) == {x for L in sl.letters for x in L} and list(sl.T1) == st['T1'])
    t0 = time.time()
    try:
        solsB, nbB = b_solve(st, st['T1'], cap); rec['B'] = dict(solutions=len(solsB), blocks=nbB, secs=round(time.time() - t0, 1))
    except (OverflowError, MemoryError) as e:
        rec['B'] = dict(undecided=repr(e)[:120]); solsB = None
    if 'solutions' in rec['A'] and solsB is not None:
        rec['agree'] = {canon(Y) for Y in solsA} == {canon(Y) for Y in solsB}
    rec['solutions'] = [list(map(list, Y)) for Y in (solsB or [])[:20]]
    # planted controls
    rng = np.random.default_rng(11); pl = []
    for t in range(npl):
        Y = planted_y(st, rng)
        Tp = [sum(b_paf(st, y)[h] for y in Y) for h in range(st['n'])]
        r = dict(T0=Tp[0])
        try:
            sB, _ = b_solve(st, Tp, cap); r['B'] = canon(Y) in {canon(Z) for Z in sB}
        except (OverflowError, MemoryError): r['B'] = None
        try:
            sA, _, _ = a_solve(key, Tp, cap); r['A'] = canon(Y) in {canon(Z) for Z in sA}
        except Exception: r['A'] = None
        pl.append(r)
    rec['planted'] = pl
    print(json.dumps(rec), flush=True)
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
