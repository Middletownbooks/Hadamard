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
Method C (this file): CP-SAT on the integer model (ortools); INFEASIBLE decides, UNKNOWN does not.
--strong adds to C the redundant row and real-character equations and the block symmetry breaking, and to B the
per-block row and real-character pruning [P].
usage: python3 theory/open_screen_stage1_exact.py 'KEY' [--methods A,B,C] [--strong] [--planted 4] [--cap 2e7] [--csecs 600]
       [--out FILE] [--aslimit GB]"""
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
        elif len(g) == 3:                        # [axis, order, generator]: the generator is given (keys on 2-power
            a, o, u = g; n = axes[a]               # axes, where cyclic unit subgroups of one order are not unique)
            assert math.gcd(u, n) == 1 and _ordmod(u, n) == o, (g, n)
            mult[a] = u
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
    return dict(axes=axes, els=els, orbits=orbits, orb=orb, letters=letters, T1=T1, lagmap=lagmap, n=n,
                fam=fam, nQ=len(Q), q=q)


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


def nf_set(st):
    """k4 with Q of order 4 and |K| = m odd: the allowed Q-compressions of a block by the compression normal form
    (OpenAI model, report/openai_4p3_review.md section 1, [P]): +-(a, a, a, a), a odd, a^2 <= m, or
    +-(t 1 + 2 delta_s), t odd, t(t + 1) <= m, any position s; returns (Q-index of every element of G', list of
    allowed z vectors, number of Q elements) or None when the theorem does not apply"""
    nQ = st.get('nQ', 0)
    if st.get('fam') != 'k4' or nQ == 0: return None
    Qax = st['axes'][:nQ]; Qn = int(np.prod(Qax)); m = st['n'] * st['q'] // Qn
    if Qn != 4 or m % 2 == 0: return None
    qidx = [int(np.ravel_multi_index(e[:nQ], Qax)) for e in st['els']]
    S = set()
    for a in range(1, math.isqrt(m) + 1, 2):
        for sg in (1, -1): S.add(tuple([sg * a] * 4))
    for t in range(-math.isqrt(m) - 2, math.isqrt(m) + 2):
        if t % 2 == 0 or t * (t + 1) > m: continue
        for pos in range(4):
            z = [t] * 4; z[pos] += 2
            for sg in (1, -1): S.add(tuple(sg * x for x in z))
    return qidx, sorted(S), 4


def b_blocks(st, Emax, cap, T1=None, prune=False):
    """all U-invariant blocks with lag-0 energy <= Emax; with prune, also r^2 <= sum T1 and, for every real character
    chi constant on the U-orbits, (sum_g chi(g) y(g))^2 <= sum_h chi(h) T1(h) (each of the four nonnegative squares
    is at most their sum, which stage 1 fixes [P])"""
    orbits = st['orbits']; L = st['letters']; w = [len(o) for o in orbits]
    order = sorted(range(len(orbits)), key=lambda j: -w[j])
    out = []
    vals = [0] * len(orbits)
    T1 = st['T1'] if T1 is None else T1; n = st['n']
    chis = [c for c in real_characters(st) if all(c[g] == c[o[0]] for o in orbits for g in o)] if prune else []
    Rb = math.isqrt(int(sum(T1))) if prune else None
    Ab = [math.isqrt(max(0, int(sum(c[h] * T1[h] for h in range(n))))) for c in chis]
    W = [[c[o[0]] * len(o) for o in orbits] for c in chis]
    big = max(abs(x) for x in L)
    restw = [sum(w[order[u]] for u in range(t + 1, len(order))) for t in range(len(order))]
    nf = nf_set(st) if prune else None
    if nf:
        qidx, NF, nq_ = nf; NFa = np.array(NF, np.int64)
        cq = [[sum(1 for g in orbits[j] if qidx[g] == q) for q in range(nq_)] for j in range(len(orbits))]
        capq = [[big * sum(cq[order[u]][q] for u in range(t + 1, len(order))) for q in range(nq_)] for t in range(len(order))]
    zc = [0] * (nf[2] if nf else 0)
    def nf_ok(t):
        if not nf: return True
        lo = np.array([zc[q] - capq[t][q] for q in range(nf[2])]); hi = np.array([zc[q] + capq[t][q] for q in range(nf[2])])
        return bool(((NFa >= lo) & (NFa <= hi)).all(1).any())
    def dfs(t, E, r, A):
        if len(out) > cap: raise OverflowError
        if t == len(order):
            if prune and (r * r > Rb * Rb or any(a * a > b * b for a, b in zip(A, Ab))): return
            if nf and tuple(zc) not in set(NF): return
            y = [0] * st['n']
            for j, o in enumerate(orbits):
                for g in o: y[g] = vals[j]
            out.append((E, tuple(y))); return
        j = order[t]; rest = restw[t]                                    # every letter has |y| >= 1
        for a in L:
            e2 = E + w[j] * a * a
            if e2 + rest > Emax: continue
            r2 = r + w[j] * a; A2 = [x + Wc[j] * a for x, Wc in zip(A, W)]
            if prune and (abs(r2) - big * rest > Rb or any(abs(x) - big * rest > b for x, b in zip(A2, Ab))): continue
            if nf:
                for q in range(nf[2]): zc[q] += cq[j][q] * a
                ok = nf_ok(t)
                if ok: vals[j] = a; dfs(t + 1, e2, r2, A2)
                for q in range(nf[2]): zc[q] -= cq[j][q] * a
                continue
            vals[j] = a; dfs(t + 1, e2, r2, A2)
    dfs(0, 0, 0, [0] * len(chis))
    return out


def b_solve(st, T1, cap=2 * 10 ** 7, prune=False):
    minE = st['n']; T0 = T1[0]
    blocks = b_blocks(st, T0 - 3 * minE, cap, T1, prune)
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


def b_solve_np(st, T1, cap=2 * 10 ** 7, prune=True, maxpairs=6 * 10 ** 7, blocks=None):
    """method B with a vectorised, bucketed join (same blocks as b_solve).  Distinct block PAF vectors are grouped by
    |row sum| and by the invariants e = PAF(0) (lag-0 energy) and a_chi = sum_h chi(h) PAF(h) (= A_chi^2) for every
    real character chi constant on the U-orbits.  A quadruple solves stage 1 only if these invariants of its two pairs
    add up to those of T1, so for every row multiset (r_a, r_b | r_c, r_d) with sum r^2 = sum T1 the pairs are formed
    group pair by group pair and joined bucket by bucket on the pair-sum invariants [P]; inside a bucket the right
    pair sums are hashed with a random 64-bit linear hash (exact under wrap-around), sorted, the left pairs streamed
    against them, and every hash match re-checked on the exact PAF vectors."""
    minE = st['n']; T0 = T1[0]; n = st['n']
    if blocks is None: blocks = b_blocks(st, T0 - 3 * minE, cap, T1, prune)
    by = defaultdict(lambda: defaultdict(list))
    for E, y in blocks:
        if y[0] < 0: continue
        by[abs(sum(y))][b_paf(st, y)].append(y)
    chis = [c for c in real_characters(st) if all(c[g] == c[o[0]] for o in st['orbits'] for g in o)]
    nf = nf_set(st)
    if nf:                                                       # Q-level PAF (sum over the cosets of H'): linear
        qidx = nf[0]
        inv = lambda P: (P[0],) + tuple(sum(c[h] * P[h] for h in range(n)) for c in chis) + \
            tuple(sum(P[h] for h in range(n) if qidx[h] == q) for q in range(nf[2]))
    else:
        inv = lambda P: (P[0],) + tuple(sum(c[h] * P[h] for h in range(n)) for c in chis)
    tinv = inv(T1)
    R = int(sum(T1)); rows = sorted(by)
    rng = np.random.default_rng(12345)
    hv = rng.integers(-2 ** 62, 2 ** 62, size=n, dtype=np.int64)
    K = {r: list(by[r]) for r in rows}
    H = {r: (np.array(K[r], np.int64) * hv).sum(1) if K[r] else np.zeros(0, np.int64) for r in rows}
    grp = {r: defaultdict(list) for r in rows}
    for r in rows:
        for i_, P in enumerate(K[r]): grp[r][inv(P)].append(i_)
    for r in rows:
        for g in grp[r]: grp[r][g] = np.array(grp[r][g], np.int64)
    tH = np.int64((np.array(T1, np.int64) * hv).sum())
    def pair_buckets(ra, rb):
        out = defaultdict(list); ga = list(grp[ra]); gb = list(grp[rb])
        for x in ga:
            for y in gb:
                if ra == rb and x > y: continue
                out[tuple(p + q for p, q in zip(x, y))].append((x, y))
        return out
    def pairs(ra, rb, lst):
        I = []; J = []
        for x, y in lst:
            a, b = grp[ra][x], grp[rb][y]
            A_, B_ = np.meshgrid(a, b, indexing='ij'); A_ = A_.ravel(); B_ = B_.ravel()
            if ra == rb and x == y: m = A_ <= B_; A_ = A_[m]; B_ = B_[m]
            I.append(A_); J.append(B_)
        return (np.concatenate(I), np.concatenate(J)) if I else (np.zeros(0, np.int64), np.zeros(0, np.int64))
    sols = set(); stats = []
    combos = {tuple(sorted(c)) for c in itertools.product(rows, repeat=4) if sum(x * x for x in c) == R}
    for c in sorted(combos):
        ra, rb, rc, rd = c
        BL = pair_buckets(ra, rb); BR = pair_buckets(rc, rd); nL = nR = 0
        for key, lstR in BR.items():
            need_key = tuple(t - k for t, k in zip(tinv, key))
            if need_key not in BL: continue
            I, J = pairs(rc, rd, lstR)
            if len(I) > maxpairs: raise OverflowError(f"bucket of {len(I)} right pairs")
            right = H[rc][I] + H[rd][J]; o = np.argsort(right, kind='stable'); rs = right[o]; nR += len(rs)
            A_, B_ = pairs(ra, rb, BL[need_key])
            for c0 in range(0, len(A_), 1 << 22):
                a_, b_ = A_[c0:c0 + (1 << 22)], B_[c0:c0 + (1 << 22)]
                need = tH - (H[ra][a_] + H[rb][b_]); nL += len(need)
                pos = np.searchsorted(rs, need); pos[pos >= len(rs)] = len(rs) - 1
                for h in np.nonzero(rs[pos] == need)[0]:
                    lo = np.searchsorted(rs, need[h], 'left'); hi = np.searchsorted(rs, need[h], 'right')
                    for q in o[lo:hi]:
                        ka, kb, kc, kd = K[ra][a_[h]], K[rb][b_[h]], K[rc][I[q]], K[rd][J[q]]
                        if [x + y_ + z + t for x, y_, z, t in zip(ka, kb, kc, kd)] != list(T1): continue
                        for ya in by[ra][ka]:
                            for yb in by[rb][kb]:
                                for yc in by[rc][kc]:
                                    for yd in by[rd][kd]: sols.add(tuple(sorted((ya, yb, yc, yd))))
        stats.append(dict(rows=c, left_pairs=int(nL), right_pairs=int(nR)))
    for Y in sols:
        assert [sum(b_paf(st, y)[h] for y in Y) for h in range(n)] == list(T1)
    return sorted(sols), len(blocks), stats


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
    # start from the letters of least modulus (+-1 for even s; for odd s, sigma is odd and +-1 is not a letter --
    # fixed 2026-11-06: the plants of odd-s slices were not letter vectors; for even s the draws are unchanged)
    m0 = min(abs(x) for x in L); low = [x for x in (-m0, m0) if x in L]
    vals = [[int(rng.choice(low)) for _ in w] for _ in range(4)]
    E = sum(w) * 4 * m0 * m0
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


def real_characters(st):
    """the nontrivial characters of G' with values +-1 (products of e -> (-1)^e on the even-order axes)"""
    axes = st['axes']; els = st['els']; out = []
    ev = [a for a, n in enumerate(axes) if n % 2 == 0]
    for mask in itertools.product((0, 1), repeat=len(ev)):
        if not any(mask): continue
        out.append([(-1) ** sum(e[a] for m, a in zip(mask, ev) if m) for e in els])
    return out


def c_solve(st, T1, secs=600, workers=4, seed=0, strong=False):
    """method C: CP-SAT on the integer stage-1 model (y per block and U-orbit in the letter domain, lag products by
    AddMultiplicationEquality, one equation per lag); returns 'INFEASIBLE', 'FEASIBLE' (with a verified quadruple) or
    'UNKNOWN' (undecided).  Symmetry breaking: block 0 cell 0 > 0 is not imposed (sound and simple)."""
    from ortools.sat.python import cp_model
    M = cp_model.CpModel(); orbits = st['orbits']; L = st['letters']; n = st['n']; lm = st['lagmap']
    dom = cp_model.Domain.FromValues(L)
    Y = [[M.NewIntVarFromDomain(dom, f"y{k}_{j}") for j in range(len(orbits))] for k in range(4)]
    cell = st['orb']
    bound = max(abs(x) for x in L) ** 2
    prod = {}
    def pv(k, a, b):
        a, b = min(a, b), max(a, b)
        if (k, a, b) not in prod:
            z = M.NewIntVar(-bound, bound, f"p{k}_{a}_{b}")
            if a == b: M.AddMultiplicationEquality(z, [Y[k][a], Y[k][a]])
            else: M.AddMultiplicationEquality(z, [Y[k][a], Y[k][b]])
            prod[(k, a, b)] = z
        return prod[(k, a, b)]
    for h in range(n):
        terms = [pv(k, cell[g], cell[lm[h][g]]) for k in range(4) for g in range(n)]
        M.Add(sum(terms) == T1[h])
    if strong:
        # redundant exact consequences [P]: rows r_k (trivial character) and every real (order-2) character chi of G'
        # constant on U-orbits: sum_k (sum_g chi(g) y_k(g))^2 = sum_h chi(h) T1(h).  Symmetry: negating a block or
        # permuting blocks maps solutions to solutions, so impose r_k >= 0 and r_0 >= r_1 >= r_2 >= r_3.
        w = [len(o) for o in orbits]; R = int(sum(T1))
        rk = [M.NewIntVar(0, 10 ** 6, f"r{k}") for k in range(4)]
        for k in range(4): M.Add(rk[k] == sum(w[j] * Y[k][j] for j in range(len(orbits))))
        for k in range(3): M.Add(rk[k] >= rk[k + 1])
        sq = [M.NewIntVar(0, 10 ** 9, f"rs{k}") for k in range(4)]
        for k in range(4): M.AddMultiplicationEquality(sq[k], [rk[k], rk[k]])
        M.Add(sum(sq) == R)
        for chi in real_characters(st):
            vals = [chi[o[0]] for o in orbits]
            if any(chi[g] != vals[j] for j, o in enumerate(orbits) for g in o): continue
            tgt = int(sum(chi[h] * T1[h] for h in range(n)))
            Ak = [M.NewIntVar(-10 ** 6, 10 ** 6, f"A{k}") for k in range(4)]
            Aq = [M.NewIntVar(0, 10 ** 9, f"Aq{k}") for k in range(4)]
            for k in range(4):
                M.Add(Ak[k] == sum(vals[j] * w[j] * Y[k][j] for j in range(len(orbits))))
                M.AddMultiplicationEquality(Aq[k], [Ak[k], Ak[k]])
            M.Add(sum(Aq) == tgt)
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = secs
    sol.parameters.num_search_workers = workers; sol.parameters.random_seed = seed
    r = sol.Solve(M)
    if r == cp_model.INFEASIBLE: return 'INFEASIBLE', None
    if r in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        Yv = []
        for k in range(4):
            y = [0] * n
            for j, o in enumerate(orbits):
                for g in o: y[g] = sol.Value(Y[k][j])
            Yv.append(tuple(y))
        assert [sum(b_paf(st, y)[h] for y in Yv) for h in range(n)] == list(T1)
        return 'FEASIBLE', Yv
    return 'UNKNOWN', None


def canon(Y):
    """quadruple up to block order and block sign"""
    return tuple(sorted(tuple(y) if y[0] > 0 else tuple(-x for x in y) for y in Y))


def main():
    import resource
    a = sys.argv[1:]
    gb = int(a[a.index('--aslimit') + 1]) if '--aslimit' in a else 6      # address-space cap in GB (default 6)
    resource.setrlimit(resource.RLIMIT_AS, (gb << 30, gb << 30))
    key = json.loads(a[0]); npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    cap = int(float(a[a.index('--cap') + 1])) if '--cap' in a else 2 * 10 ** 7
    outp = a[a.index('--out') + 1] if '--out' in a else None
    rec = dict(key=key)
    st = b_setup(key)
    # consistency of the independent set-up with ab_lib (orbit partition of G', letters, T1)
    methods = a[a.index('--methods') + 1].split(',') if '--methods' in a else ['A', 'B']
    strong = '--strong' in a
    csecs = int(a[a.index('--csecs') + 1]) if '--csecs' in a else 600
    t0 = time.time()
    if 'C' in methods:
        r_, Yc = c_solve(st, st['T1'], csecs, strong=strong); rec['C'] = dict(status=r_, secs=round(time.time() - t0, 1))
        if Yc: rec['C']['solution'] = [list(y) for y in Yc]
    t0 = time.time()
    try:
        if 'A' not in methods: raise RuntimeError('method A not run')
        solsA, nbA, sl = a_solve(key, cap=cap); rec['A'] = dict(solutions=len(solsA), blocks=nbA, secs=round(time.time() - t0, 1))
    except Exception as e:
        rec['A'] = dict(undecided=repr(e)[:120]); sl = None
    if sl is not None:
        rec['setup_agrees'] = (sorted(sorted(o) for o in st['orbits']) == sorted(sorted(np.nonzero(sl.borb == j)[0].tolist()) for j in range(len(sl.breps)))
                               and set(st['letters']) == {x for L in sl.letters for x in L} and list(sl.T1) == st['T1'])
    t0 = time.time()
    try:
        if 'B' not in methods: raise OverflowError('method B not run')
        if strong:
            solsB, nbB, stB = b_solve_np(st, st['T1'], cap)
            rec['B'] = dict(solutions=len(solsB), blocks=nbB, join=stB, secs=round(time.time() - t0, 1))
        else:
            solsB, nbB = b_solve(st, st['T1'], cap); rec['B'] = dict(solutions=len(solsB), blocks=nbB, secs=round(time.time() - t0, 1))
    except (OverflowError, MemoryError) as e:
        rec['B'] = dict(undecided=repr(e)[:120]); solsB = None
    if 'solutions' in rec['A'] and solsB is not None:
        rec['agree'] = {canon(Y) for Y in solsA} == {canon(Y) for Y in solsB}
    rec['solutions'] = [list(map(list, Y)) for Y in (solsB or [])[:20]]
    if outp: open(outp, 'a').write(json.dumps(dict(rec, phase='decision')) + "\n")
    csecs = min(csecs, int(a[a.index('--plsecs') + 1]) if '--plsecs' in a else 600)
    # planted controls
    rng = np.random.default_rng(11); pl = []
    for t in range(npl):
        # planted quadruples drawn with sum_k r_k^2 as small as found in 4000 draws (the real target has
        # sum T1 = sum r_k^2, e.g. 4 for b = 1), so the pruned enumeration is of the real problem's size
        best = None
        if strong and nf_set(st) and 'B' in methods:
            # k4 normal form: plant four real enumerated blocks (one constant Q-compression, three singletons)
            if t == 0:
                BL = b_blocks(st, st['T1'][0] - 3 * st['n'], cap, st['T1'], True)
                qidx_ = nf_set(st)[0]
                zof = lambda y: tuple(sum(y[g] for g in range(st['n']) if qidx_[g] == q) for q in range(4))
                const = [y for E, y in BL if len(set(zof(y))) == 1]; single = [y for E, y in BL if len(set(zof(y))) > 1]
            Y = [const[int(rng.integers(len(const)))]] + [single[int(rng.integers(len(single)))] for _ in range(3)]
            best = (0, Y)
        for _ in range(4000 if (strong and best is None) else (1 if best is None else 0)):
            Yc = planted_y(st, rng); R_ = sum(sum(y) ** 2 for y in Yc)
            if best is None or abs(R_ - sum(st['T1'])) < best[0]: best = (abs(R_ - sum(st['T1'])), Yc)
            if best[0] == 0: break
        Y = best[1]
        Tp = [sum(b_paf(st, y)[h] for y in Y) for h in range(st['n'])]
        r = dict(T0=Tp[0], R=int(sum(Tp)))
        if 'B' in methods:
            try:
                if strong and nf_set(st):
                    # the planted blocks are real enumerated blocks: join over that same list with the planted target
                    sB = b_solve_np(st, Tp, cap, blocks=BL)[0]
                else:
                    sB = b_solve_np(st, Tp, cap)[0] if strong else b_solve(st, Tp, cap)[0]
                r['B'] = canon(Y) in {canon(Z) for Z in sB}
            except (OverflowError, MemoryError): r['B'] = None
        if 'A' in methods:
            try:
                sA, _, _ = a_solve(key, Tp, cap); r['A'] = canon(Y) in {canon(Z) for Z in sA}
            except Exception: r['A'] = None
        if 'C' in methods:
            r['C'] = c_solve(st, Tp, csecs, strong=strong)[0]
        pl.append(r)
    rec['planted'] = pl
    rec['phase'] = 'final'
    print(json.dumps(rec), flush=True)
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
