#!/usr/bin/env python3
r"""[M] Independent re-implementation of the elliptic-character seed ansatz of report/openai_4p3_answer.md (Proposal 2),
for the review report/openai_4p3_review.md.

The ansatz, its reduced conditions (C), (K), (S) and the example parameters are by the OpenAI model.  This file is new
(Claude Code review) and is written from the text of the answer only: it imports no repository code, in particular
nothing from theory/openai_4p3_*.py.  The bordered construction that the cores feed is Kumar's form of EKLO Theorem 7.

Definition (answer, Proposal 2 section 1).  chi = quadratic character of F_p, chi(0) = 0; f_l(v) = v (v - 1)(v - l).
A seed is (l, eps, mask, u): l in F_p \ {0, 1}, eps = +-1, mask in 0..7, u != 0, and
    b(k) = eps chi(f_l(u k))           if u k is not a root,
    b(k) = -1 if bit j of mask is set, else +1,   if u k is root j of (0, 1, l).
A switched seed adds a root position r with b(r) = -1 and sets b'(r) = +1.  Row 0 has four translates of a seed b_0;
row i = 1, 2, 3 has b'_i in cell q = 0 and translates of b_i in cells 1, 2, 3 (cell shifts h_iq, h_i0 = 0).
Seed gate (K), t = 1..(p-1)/2:  P_0(t) + sum_i (3 P_i(t) + P'_i(t)) / 4 = -3.
Phase gate (S), d != 0, t in Z_p:  sum_{i,q} PCF(u_iq, u_{i,q+d}; t + h_iq - h_{i,q+d}) = 4.
Coordinates of a core: q * p + k, q in 0..3 (Q = Z2^2 with XOR, or Z4).

Subcommands
  gate     --primes 11,17,19   exhaustive seed gate: every parameter tuple, distinct vectors, exact 2 + 2 join; counts
                               all solutions (distinct-vector multisets and ordered parameter tuples)
  cores    --p 5 --q v4        exhaustive core search inside the ansatz: every gate solution, every pair of seed
                               sequences realising it (up to a cap), every shift (2 + 2 join, exact verification)
  examples                     rebuild the answer's five example cores from the parameters stored in its logs
  estimate --primes ...        lattice-Gaussian estimate of the number of gate solutions [E]
"""
import argparse, itertools, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TARGETS = [149, 179, 257, 271, 331, 359, 367, 389, 457, 461, 467, 509, 521, 523, 571, 599]


def is_prime(p):
    return p > 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def chi_table(p):
    c = -np.ones(p, np.int64)
    c[0] = 0
    c[(np.arange(1, p) ** 2) % p] = 1
    return c


def seed(p, lam, eps, mask, u, chi=None):
    chi = chi_table(p) if chi is None else chi
    k = np.arange(p, dtype=np.int64)
    v = (u * k) % p
    b = eps * chi[(v * (v - 1) % p) * ((v - lam) % p) % p]
    for j, r in enumerate((0, 1, lam)):
        b[v == r] = -1 if (mask >> j) & 1 else 1
    return b


def all_seeds(p):
    """every parameter tuple: (params, b, roots positions)"""
    chi = chi_table(p)
    out = []
    for lam in range(2, p):
        for eps in (1, -1):
            for mask in range(8):
                for u in range(1, p):
                    b = seed(p, lam, eps, mask, u, chi)
                    roots = [int(k) for k in range(p) if (u * k) % p in (0, 1, lam)]
                    out.append(((lam, eps, mask, u), b, roots))
    return out


def paf_half(B):
    """PAF of the rows of B (N x p) at t = 1..(p-1)/2, exact integers"""
    p = B.shape[1]
    return np.stack([(B * np.roll(B, -t, axis=1)).sum(1) for t in range(1, (p + 1) // 2)], 1)


def seed_vectors(p):
    """ordinary vectors P_0 (with parameter multiplicity) and switched vectors D = (3P + P')/4"""
    S = all_seeds(p)
    Bm = np.array([s[1] for s in S])
    P = paf_half(Bm)
    ordv = {}
    for (par, b, roots), pv in zip(S, P):
        ordv.setdefault(tuple(pv.tolist()), []).append(par)
    sw_rows, sw_par = [], []
    for (par, b, roots), pv in zip(S, P):
        for r in roots:
            if b[r] == -1:
                b2 = b.copy(); b2[r] = 1
                sw_rows.append(b2); sw_par.append((par, r, pv))
    swv = {}
    if sw_rows:
        P2 = paf_half(np.array(sw_rows))
        for (par, r, pv), pv2 in zip(sw_par, P2):
            num = 3 * pv + pv2
            assert not (num % 4).any()
            swv.setdefault(tuple((num // 4).tolist()), []).append((par, r))
    return ordv, swv, len(S), len(sw_rows)


def gate(p, want_list=True, cap=200000):
    t0 = time.time()
    ordv, swv, nseeds, nsw = seed_vectors(p)
    h = (p - 1) // 2
    T = np.full(h, -3, np.int64)
    A = np.array(list(ordv.keys()), np.int64).reshape(-1, h)
    Bv = np.array(list(swv.keys()), np.int64).reshape(-1, h)
    wa = np.array([len(v) for v in ordv.values()], np.int64)
    wb = np.array([len(v) for v in swv.values()], np.int64)
    rng = np.random.default_rng(12345)
    coef = rng.integers(1, 2 ** 62, size=h, dtype=np.int64)       # linear hash, exact under int64 wraparound
    ha, hb = A @ coef, Bv @ coef
    hT = int(T @ coef)
    # pairs (D_i, D_j), i <= j, keyed by hash; ordered weights for the parameter-tuple count
    nb = len(Bv)
    I, J = np.triu_indices(nb)
    hp = hb[I] + hb[J]
    o = np.argsort(hp, kind='stable')
    hp, I, J = hp[o], I[o], J[o]
    sols = set()
    ordered = 0
    hTa = np.full(nb, hT, np.int64)
    for a in range(len(A)):
        need = hTa - ha[a] - hb                                  # int64 arrays: exact wraparound
        lo, hi = np.searchsorted(hp, need, 'left'), np.searchsorted(hp, need, 'right')
        for c in np.flatnonzero(hi > lo):
            for x in range(lo[c], hi[c]):
                i, j = int(I[x]), int(J[x])
                if np.array_equal(A[a] + Bv[c] + Bv[i] + Bv[j], T):
                    sols.add((a, tuple(sorted((int(c), i, j)))))
    # ordered parameter tuples: for each multiset, number of orderings of rows 1..3 times weights
    for a, trip in sols:
        perms = len(set(itertools.permutations(trip)))
        ordered += int(wa[a]) * perms * int(np.prod([wb[x] for x in trip]))
    # each multiset was found once per choice of the singled-out row c; dedup handled by the set
    res = dict(p=p, half_lags=h, parameter_seeds=nseeds, switched_parameter_seeds=nsw,
               ordinary_vectors=len(A), switched_vectors=nb, switched_pairs=int(len(hp)),
               gate_solutions_distinct_vector_multisets=len(sols), gate_solutions_ordered_parameter_tuples=ordered,
               status='empty' if not sols else 'nonempty', secs=round(time.time() - t0, 1))
    if want_list:
        res['solutions'] = [dict(P0=A[a].tolist(), D=[Bv[x].tolist() for x in trip],
                                 P0_params=ordv[tuple(A[a].tolist())][:3],
                                 D_params=[swv[tuple(Bv[x].tolist())][:3] for x in trip]) for a, trip in sorted(sols)][:50]
    return res, (ordv, swv, A, Bv, sols)


# ------------------------------------------------------------------------------------------------ phases and cores
def partner(q):
    return [[(u ^ d) if q == 'v4' else (u + d) % 4 for u in range(4)] for d in range(4)]


def pcf_table(v, w):
    """C(s) = sum_k v(k) w(k + s), s in Z_p"""
    p = len(v)
    return np.array([int((v * np.roll(w, -s)).sum()) for s in range(p)], np.int64)


def row_vectors(cells, q):
    """for one row with base cells u_0..u_3 (cell q = translate of u_q by h_q, h_0 = 0): the (S)-vector for every
    (h_1, h_2, h_3), shape (p^3, 3p), indexed by d = 1..3 and t in Z_p"""
    p = len(cells[0])
    par = partner(q)
    C = {(a, b): pcf_table(cells[a], cells[b]) for a in range(4) for b in range(4) if a != b}
    H = np.array(list(itertools.product(range(p), repeat=3)), np.int64)
    H = np.concatenate([np.zeros((len(H), 1), np.int64), H], 1)          # h_0 = 0
    t = np.arange(p)
    out = np.zeros((len(H), 3, p), np.int64)
    for d in range(1, 4):
        for a in range(4):
            b = par[d][a]
            idx = (t[None, :] + H[:, a:a + 1] - H[:, b:b + 1]) % p
            out[:, d - 1, :] += C[(a, b)][idx]
    return H, out.reshape(len(H), 3 * p)


def build_core(p, rows_cells, shifts):
    xs = []
    for cells, h in zip(rows_cells, shifts):
        x = np.concatenate([np.roll(cells[qq], h[qq]) for qq in range(4)])     # cell q: u_q(k - h_q)
        xs.append(x.astype(np.int64))
    return xs


def core_check(xs, p, q):
    """direct integer PAF on G = Q x Z_p at every lag (no formula)"""
    par = partner(q)
    X = [x.reshape(4, p) for x in xs]
    for d in range(4):
        for t in range(p):
            v = sum(int((Xi[np.arange(4)] * np.roll(Xi[par[d]], -t, axis=1)).sum()) for Xi in X)
            want = (16 * p if t == 0 else -12) if d == 0 else 4
            if v != want:
                return False
    return True


class PhaseHits:
    """all hash-matches of the 2 + 2 shift join for fixed cells (a linear 64-bit hash is exact under wraparound, so no
    solution is missed); hits are verified exactly when fetched"""

    def __init__(self, p, q, rows_cells, coef):
        self.HV = [row_vectors(c, q) for c in rows_cells]
        self.T = np.full(3 * p, 4, np.int64)
        hT = np.int64(int(self.T @ coef))
        h = [hv[1] @ coef for hv in self.HV]
        self.n = len(h[0])
        s01 = (h[0][:, None] + h[1][None, :]).ravel()
        s23 = (h[2][:, None] + h[3][None, :]).ravel()
        self.o = np.argsort(s23, kind='stable')
        s23s = s23[self.o]
        need = hT - s01
        self.lo = np.searchsorted(s23s, need, 'left'); self.hi = np.searchsorted(s23s, need, 'right')
        self.cnt = self.hi - self.lo
        self.cum = np.cumsum(self.cnt)
        self.total = int(self.cum[-1]) if len(self.cum) else 0

    def fetch(self, k):
        """the k-th hit (0-based), verified: shift lists or None for a hash collision"""
        x = int(np.searchsorted(self.cum, k, 'right'))
        before = int(self.cum[x - 1]) if x else 0
        y = int(self.o[self.lo[x] + (k - before)])
        i0, i1 = divmod(x, self.n); i2, i3 = divmod(y, self.n)
        HV = self.HV
        V = HV[0][1][i0] + HV[1][1][i1] + HV[2][1][i2] + HV[3][1][i3]
        if not np.array_equal(V, self.T):
            return None
        return [HV[0][0][i0].tolist(), HV[1][0][i1].tolist(), HV[2][0][i2].tolist(), HV[3][0][i3].tolist()]


def realisations(p, ordv, swv, A, Bv, sol, chi):
    """distinct (row cells) realisations of a gate solution: row 0 = 4 copies of b_0; row i = (b'_i, b_i, b_i, b_i)"""
    a, trip = sol
    r0 = {}
    for par in ordv[tuple(A[a].tolist())]:
        b = seed(p, *par, chi=chi)
        r0.setdefault(tuple(b.tolist()), par)
    rs = []
    for x in trip:
        d = {}
        for par, r in swv[tuple(Bv[x].tolist())]:
            b = seed(p, *par, chi=chi)
            b2 = b.copy(); b2[r] = 1
            d.setdefault((tuple(b2.tolist()), tuple(b.tolist())), (par, r))
        rs.append(d)
    return r0, rs


def combo_cells(combo):
    (b0, par0), *rest = combo
    cells = [[np.array(b0)] * 4]
    meta = [dict(params=list(par0))]
    for (b2b, (par, r)) in rest:
        b2, b = b2b
        cells.append([np.array(b2), np.array(b), np.array(b), np.array(b)])
        meta.append(dict(params=list(par), root=int(r)))
    return cells, meta


def cores(p, q, mode='exist', seconds=540, pick_seed=2026, exclude=()):
    """mode 'exist': every gate solution and every realisation in a fixed order, stop at the first core (existence in
    the ansatz is decided when the loop completes or a core is found; a time cap means UNKNOWN).
    mode 'random': seeded uniform choice of a realisation (gate solution weighted by its number of realisations),
    then a uniform choice among all its shift solutions; repeated until a core not in `exclude` is found."""
    t0 = time.time()
    gres, (ordv, swv, A, Bv, sols) = gate(p, want_list=False)
    chi = chi_table(p)
    coef = np.random.default_rng(99).integers(1, 2 ** 62, size=3 * p, dtype=np.int64)
    sols = sorted(sols)
    reals = [realisations(p, ordv, swv, A, Bv, s, chi) for s in sols]
    sizes = [len(r0) * int(np.prod([len(d) for d in rs])) for r0, rs in reals]
    out = dict(p=p, q=q, mode=mode, gate={k: v for k, v in gres.items() if k != 'solutions'},
               gate_solutions=len(sols), realisations_per_solution=sizes, realisations_total=int(sum(sizes)),
               combos_tried=0, status=None)
    found = []
    if mode == 'exist':
        for s, (r0, rs) in zip(sols, reals):
            for combo in itertools.product(list(r0.items()), *[list(d.items()) for d in rs]):
                if time.time() - t0 > seconds:
                    out['status'] = 'capped (unknown)'
                    break
                out['combos_tried'] += 1
                cells, meta = combo_cells(combo)
                ph = PhaseHits(p, q, cells, coef)
                for k in range(ph.total):
                    sh = ph.fetch(k)
                    if sh is not None:
                        xs = build_core(p, cells, sh)
                        assert core_check(xs, p, q), 'phase join returned a non-core'
                        found.append(dict(rows=meta, shifts=sh, sequences=[x.tolist() for x in xs],
                                          shift_hits_in_combo=ph.total))
                        break
                if found:
                    out['status'] = 'core found'
                    break
            if out['status']:
                break
        if out['status'] is None:
            out['status'] = 'no core in the ansatz (exhaustive)'
    else:
        rng = np.random.default_rng(pick_seed)
        w = np.array(sizes, float) / sum(sizes)
        while time.time() - t0 < seconds and not found:
            si = int(rng.choice(len(sols), p=w))
            r0, rs = reals[si]
            lists = [list(r0.items())] + [list(d.items()) for d in rs]
            combo = [L[int(rng.integers(len(L)))] for L in lists]
            out['combos_tried'] += 1
            cells, meta = combo_cells(combo)
            ph = PhaseHits(p, q, cells, coef)
            if not ph.total:
                continue
            for _ in range(50):
                k = int(rng.integers(ph.total))
                sh = ph.fetch(k)
                if sh is None:
                    continue
                xs = build_core(p, cells, sh)
                key = tuple(np.concatenate(xs).tolist())
                if key in exclude:
                    continue
                assert core_check(xs, p, q), 'phase join returned a non-core'
                found.append(dict(rows=meta, shifts=sh, sequences=[x.tolist() for x in xs],
                                  shift_hits_in_combo=ph.total, gate_solution_index=si))
                break
        out['status'] = 'core picked' if found else 'capped (no core picked)'
    out['secs'] = round(time.time() - t0, 1)
    return out, found


# ------------------------------------------------------------------------------------------------ the answer's examples
def examples():
    """rebuild the answer's five example cores from the parameter objects stored in logs/openai_4p3/core_p*_*.json
    (data only), using this file's seed definition, and compare with the stored sequences"""
    res = []
    for p, q in ((3, 'v4'), (3, 'z4'), (5, 'v4'), (5, 'z4'), (7, 'v4')):
        path = os.path.join(ROOT, 'logs/openai_4p3/core_p%d_%s.json' % (p, q))
        d = json.load(open(path))
        chi = chi_table(p)
        cells = []
        shifts = []
        for i, par in enumerate(d['parameters']):
            b = seed(p, par['lam'], par['sign'], par['mask'], par['scale'], chi)
            if i == 0:
                cells.append([b] * 4)
            else:
                r = par['root']
                assert b[r] == -1 and (par['scale'] * r) % p in (0, 1, par['lam'])
                b2 = b.copy(); b2[r] = 1
                cells.append([b2, b, b, b])
            shifts.append(par['shifts'])
        xs = build_core(p, cells, shifts)
        same = all(np.array_equal(x, np.array(s)) for x, s in zip(xs, d['sequences']))
        ok = core_check(xs, p, q)
        bad = [x.copy() for x in xs]; bad[1][2] *= -1
        res.append(dict(p=p, q=q, order=16 * p + 12, rebuilt_equals_stored=same, core_ok=ok,
                        corrupted_rejected=not core_check(bad, p, q), name='answer example p=%d %s' % (p, q),
                        sequences=[x.tolist() for x in xs], parameters=d['parameters']))
    return res


# ------------------------------------------------------------------------------------------------ estimate [E]
def sample_K_vectors(p, N, rng):
    """N uniform parameter quadruples (row 0 ordinary, rows 1..3 switched with a uniform negative root; tuples with no
    negative root are resampled): the (K) left-hand side, t = 1..(p-1)/2; also the fraction of switched-seed tuples
    with >= 1 negative root (for the count of parameter tuples)"""
    chi = chi_table(p)
    h = (p - 1) // 2
    V = np.zeros((N, h), np.int64)

    def rand_seed():
        lam = int(rng.integers(2, p)); eps = int(rng.choice((1, -1))); mask = int(rng.integers(8))
        u = int(rng.integers(1, p))
        return (lam, eps, mask, u)
    B0 = np.array([seed(p, *rand_seed(), chi=chi) for _ in range(N)])
    V += paf_half(B0)
    tries = 0
    for i in range(3):
        Bs, B2s = [], []
        while len(Bs) < N:                                   # uniform over (seed, negative root) pairs
            tries += 1
            par = rand_seed()
            b = seed(p, *par, chi=chi)
            negs = [k for k in range(p) if (par[3] * k) % p in (0, 1, par[0]) and b[k] == -1]
            if rng.random() >= len(negs) / 3:
                continue
            r = negs[int(rng.integers(len(negs)))]
            b2 = b.copy(); b2[r] = 1
            Bs.append(b); B2s.append(b2)
        P = paf_half(np.array(Bs)); P2 = paf_half(np.array(B2s))
        V += (3 * P + P2) // 4
    return V


def n_param_tuples(p):
    """number of ordered parameter quadruples: row 0: (p-2) * 2 * 8 * (p-1); rows 1..3: sum over seeds of the number
    of negative roots = (p-2) * 2 * (p-1) * sum_mask popcount(mask) = (p-2)(p-1) * 2 * 12"""
    n0 = (p - 2) * 2 * 8 * (p - 1)
    n1 = (p - 2) * 2 * (p - 1) * 12
    return n0 * n1 ** 3


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('gate'); g.add_argument('--primes', required=True); g.add_argument('--out', required=True)
    c = sub.add_parser('cores'); c.add_argument('--p', type=int, required=True); c.add_argument('--q', required=True)
    c.add_argument('--out', required=True); c.add_argument('--mode', choices=('exist', 'random'), default='exist')
    c.add_argument('--seconds', type=float, default=540)
    c.add_argument('--pick-seed', type=int, default=2026)
    c.add_argument('--exclude', default='', help='JSON core file(s) to exclude from the random pick (comma list)')
    e = sub.add_parser('examples'); e.add_argument('--out', required=True)
    s = sub.add_parser('estimate'); s.add_argument('--primes', required=True); s.add_argument('--out', required=True)
    s.add_argument('--N', type=int, default=4000); s.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    if a.cmd == 'gate':
        res = dict(script='theory/4p3_review_elliptic.py', command=' '.join(sys.argv), seed=None, records=[])
        for p in map(int, a.primes.split(',')):
            assert is_prime(p)
            r, _ = gate(p)
            res['records'].append(r)
            print(json.dumps({k: v for k, v in r.items() if k != 'solutions'}), flush=True)
            with open(a.out, 'w') as f:
                json.dump(res, f, indent=1)
    elif a.cmd == 'cores':
        excl = set()
        for fn in [x for x in a.exclude.split(',') if x]:
            d = json.load(open(fn))
            excl.add(tuple(np.concatenate([np.array(s) for s in d['sequences']]).tolist()))
        out, found = cores(a.p, a.q, mode=a.mode, seconds=a.seconds, pick_seed=a.pick_seed, exclude=excl)
        out.update(script='theory/4p3_review_elliptic.py', command=' '.join(sys.argv))
        out['cores'] = [dict(name='elliptic p=%d %s (%s)' % (a.p, a.q, a.mode), p=a.p, q=a.q,
                             source='theory/4p3_review_elliptic.py cores --mode %s --pick-seed %d' % (a.mode, a.pick_seed),
                             parameters=dict(rows=x['rows'], shifts=x['shifts']), sequences=x['sequences'],
                             shift_hits_in_combo=x['shift_hits_in_combo'], equals_excluded=False) for x in found]
        with open(a.out, 'w') as f:
            json.dump(out, f, indent=1)
        print(json.dumps({k: v for k, v in out.items() if k not in ('cores', 'gate', 'realisations_per_solution')}))
    elif a.cmd == 'examples':
        res = examples()
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
        print(json.dumps([{k: v for k, v in r.items() if k not in ('sequences', 'parameters')} for r in res]))
    elif a.cmd == 'estimate':
        sys.path.insert(0, HERE)
        import q523_gauss as QG                               # the repository's lattice Gaussian (estimate only)
        res = dict(script='theory/4p3_review_elliptic.py estimate', command=' '.join(sys.argv), seed=a.seed,
                   N=a.N, records=[])
        for p in map(int, a.primes.split(',')):
            t0 = time.time()
            rng = np.random.default_rng(a.seed + p)
            V = sample_K_vectors(p, a.N, rng)
            T = np.full((p - 1) // 2, -3, np.int64)
            gb = QG.gauss_bits(V, T, nlat=min(400, a.N - 1))
            lp = math.log2(n_param_tuples(p))
            rec = dict(p=p, equations=(p - 1) // 2, log2_param_quadruples=round(lp, 2), gauss=gb,
                       log2_expected=round(lp - gb['bits'], 2) if np.isfinite(gb.get('bits', np.nan)) else None,
                       secs=round(time.time() - t0, 1))
            res['records'].append(rec)
            print(json.dumps(rec), flush=True)
            with open(a.out, 'w') as f:
                json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
