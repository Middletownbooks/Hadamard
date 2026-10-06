"""Independent SAT cross-check for Williamson quadruples (pre-registration addendum S).

Shares no code with wil_search.py. Variables: v[X][j] = (x_j == -1) for the half j = 0..(n-1)/2 of symmetric block X.
Williamson condition per lag s: total number of disagreements x_i != x_{i+s} over the 4 blocks and i in Z_n equals 2n.

Usage:
  wil_sat.py enum n                  # S0: enumerate all normalised quadruples (blocking clauses), print count
  wil_sat.py inc n p q k [--shard i/K] [--out F]   # SAT+CAS: both compressions as assumptions, one solver per case
  wil_sat.py case n k [--z13 [--q Q] [--shard i/K]]  # solve row-sum case k, optionally split by Z_Q compression reps
"""
import sys, json, time, math, itertools
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153


def rowsum_cases(n):
    out = []
    for a in range(1, int(math.isqrt(4 * n)) + 1, 2):
        for b in range(1, a + 1, 2):
            for c in range(1, b + 1, 2):
                r = 4 * n - a * a - b * b - c * c
                if r <= 0:
                    continue
                d = math.isqrt(r)
                if d * d == r and d % 2 == 1 and d <= c:
                    out.append((a, b, c, d))
    return out


def half(i, n):
    i %= n
    return min(i, n - i)


class Enc:
    def __init__(self, n, case, prod=False):
        self.n, self.case = n, case
        self.h = (n + 1) // 2
        self.pool = IDPool()
        self.cnf = CNF()
        self.v = [[self.pool.id(("v", X, j)) for j in range(self.h)] for X in range(4)]
        self._rows()
        self._lags()
        if prod:
            self._product()
            self._pairs()

    def _product(self):
        """Williamson (1944): a_j b_j c_j d_j is the same for every j != 0. Equivalently the number of -1 entries
        among the four blocks at position j has fixed parity; with W the total number of -1 entries, that parity
        is (W - n) mod 2. Encoded as 4-variable XOR constraints (8 clauses each)."""
        n = self.n
        W = sum((n - s) // 2 for s in self.case)
        par = (W - n) % 2
        for j in range(1, self.h):
            vs = [self.v[X][j] for X in range(4)]
            for signs in itertools.product([0, 1], repeat=4):
                if sum(signs) % 2 != par:     # forbid assignments whose parity differs
                    self.cnf.append([(-v if b else v) for v, b in zip(vs, signs)])
        # Williamson's theorem also fixes position 0: a_0 b_0 c_0 d_0 = -eps (opposite parity)
        vs = [self.v[X][0] for X in range(4)]
        for signs in itertools.product([0, 1], repeat=4):
            if sum(signs) % 2 == par:
                self.cnf.append([(-v if b else v) for v, b in zip(vs, signs)])

    def _pairs(self):
        """Redundant second form of the definition: PAF_X(s) = n - 4 w_X + 4 #{i : x_i = x_{i+s} = -1}, so
        sum_X #{i : both -1} = W - n for every s (W = total number of -1 entries). AND auxiliaries + cardinality."""
        n = self.n
        W = sum((n - s) // 2 for s in self.case)
        for s in range(1, self.h):
            lits = []
            for X in range(4):
                for i in range(n):
                    a, b = half(i, n), half(i + s, n)
                    if a == b:
                        lits.append(self.v[X][a])          # x_a = x_a = -1  <=>  v_a
                    else:
                        key = (min(a, b), max(a, b))
                        t = self.pool.id(("and", X, key))
                        p_, q_ = self.v[X][key[0]], self.v[X][key[1]]
                        if not hasattr(self, "_andset"): self._andset = set()
                        if (X, key) not in self._andset:
                            self.cnf.extend([[-t, p_], [-t, q_], [t, -p_, -q_]]); self._andset.add((X, key))
                        lits.append(t)
            # each unordered pair occurs twice and each self-pair once; duplicates go through equal copies
            cnt, L = {}, []
            for l in lits:
                if l in cnt:
                    y = self.pool.id(("pdup", s, l, cnt[l])); cnt[l] += 1
                    self.cnf.extend([[-y, l], [y, -l]]); L.append(y)
                else:
                    cnt[l] = 1; L.append(l)
            self.card_eq(L, W - n)

    def card_eq(self, lits, k):
        if k < 0 or k > len(lits):
            self.cnf.append([])          # impossible
            return
        enc = CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.totalizer)
        self.cnf.extend(enc.clauses)

    def _rows(self):
        n = self.n
        for X, s in enumerate(self.case):
            w = (n - s) // 2                      # number of -1 entries in the full block
            v0 = w % 2                            # x_0 is counted once, the others twice
            self.cnf.append([self.v[X][0] if v0 else -self.v[X][0]])
            self.card_eq(self.v[X][1:], (w - v0) // 2)

    def _lags(self):
        n = self.n
        for s in range(1, self.h):
            mult = {}
            for i in range(n):
                a, b = half(i, n), half(i + s, n)
                if a != b:
                    key = (min(a, b), max(a, b))
                    mult[key] = mult.get(key, 0) + 1
            assert set(mult.values()) == {2}, "unexpected multiplicity"
            lits = []
            for X in range(4):
                for (a, b) in mult:
                    t = self.pool.id(("x", X, a, b))
                    p, q = self.v[X][a], self.v[X][b]
                    self.cnf.extend([[-t, p, q], [-t, -p, -q], [t, -p, q], [t, p, -q]])
                    lits.append(t)
            self.card_eq(lits, n)                 # 2n disagreements, each half-pair counted twice

    def add_z13(self, comp):
        """comp[X][r] for r = 0..(q-1)/2: sum of x_i over i = r mod q (q = 13 by default)."""
        n, q = self.n, len(comp[0]) * 2 - 1
        for X in range(4):
            for r in range(len(comp[X])):
                lits = []
                for i in range(n):
                    if i % q == r:
                        j = half(i, n)
                        lits.append(self.v[X][j])
                # duplicate half-variables (r = 0) are passed through equal aux copies
                seen, L = {}, []
                for l in lits:
                    if l in seen:
                        y = self.pool.id(("dup", X, r, l, seen[l]))
                        self.cnf.extend([[-y, l], [y, -l]])
                        seen[l] += 1
                        L.append(y)
                    else:
                        seen[l] = 1
                        L.append(l)
                m = len(L)                          # number of -1 entries = (m - sum)/2
                self.card_eq(L, (m - comp[X][r]) // 2)

    def decode(self, model):
        ms = set(l for l in model if l > 0)
        n = self.n
        Q = []
        for X in range(4):
            hv = [(-1 if self.v[X][j] in ms else 1) for j in range(self.h)]
            Q.append([hv[half(i, n)] for i in range(n)])
        return Q

    def block(self, model):
        ms = set(l for l in model if l > 0)
        return [(-l if l in ms else l) for X in range(4) for l in self.v[X]]


def verify(Q):
    n = len(Q[0])
    ok = all(all(x in (1, -1) for x in X) for X in Q) and all(X[i] == X[(-i) % n] for X in Q for i in range(n))
    for s in range(1, n):
        if sum(sum(X[i] * X[(i + s) % n] for i in range(n)) for X in Q) != 0:
            return False
    return ok


def enum_all(n):
    total = 0
    for case in rowsum_cases(n):
        e = Enc(n, case)
        with Cadical153(bootstrap_with=e.cnf.clauses) as S:
            while S.solve():
                m = S.get_model()
                Q = e.decode(m)
                assert verify(Q), "SAT model fails verification"
                total += 1
                S.add_clause(e.block(m))
    return total


def comp13_reps(n, case, q=13):
    """Z_q compression quadruples via the independent stage-A implementation, canonical under (Z_q)*/+-1."""
    sys.path.insert(0, "scripts")
    from wil_stageA_alt import quads
    Qs = quads(n, q, case)
    units = [u for u in range(1, (q + 1) // 2) if math.gcd(u, q) == 1]
    def full(hv):
        return [hv[min(i, q - i)] for i in range(q)]
    def act(quad, u):
        out = []
        for hv in quad:
            f = full(hv)
            g = [f[(u * i) % q] for i in range(q)]
            out.append(tuple(g[: (q + 1) // 2]))
        return tuple(out)
    reps = set()
    for quad in Qs:
        reps.add(min(act(quad, u) for u in units))
    return sorted(reps)



def unary(cnf, pool, lits, tag):
    """Bidirectional totalizer: returns u[1..m] with u[k] <-> (at least k of lits are true)."""
    def build(L, key):
        if len(L) == 1:
            return [None, L[0]]
        mid = len(L) // 2
        a, b = build(L[:mid], key + "a"), build(L[mid:], key + "b")
        na, nb = len(a) - 1, len(b) - 1
        c = [None] + [pool.id((tag, key, k)) for k in range(1, na + nb + 1)]
        for i in range(na + 1):
            for j in range(nb + 1):
                k = i + j
                if k >= 1:          # a>=i and b>=j  ->  c>=i+j
                    cl = [c[k]]
                    if i: cl.append(-a[i])
                    if j: cl.append(-b[j])
                    cnf.append(cl)
                if k + 1 <= na + nb:  # a<=i and b<=j  ->  c<=i+j
                    cl = [-c[k + 1]]
                    if i + 1 <= na: cl.append(a[i + 1])
                    if j + 1 <= nb: cl.append(b[j + 1])
                    cnf.append(cl)
        return c
    return build(list(lits), "")


class IncCase(Enc):
    """Williamson constraints for one row-sum case, plus unary counters for every coset of Z_p and Z_q,
    so that any pair of compression quadruples can be imposed by assumptions on one incremental solver."""
    def __init__(self, n, case, p, q, prod=False):
        super().__init__(n, case, prod)
        self.p, self.q = p, q
        self.ctr = {}
        for m in (p, q):
            for X in range(4):
                for r in range((m + 1) // 2):
                    lits, seen = [], {}
                    for i in range(n):
                        if i % m == r:
                            l = self.v[X][half(i, n)]
                            if l in seen:
                                y = self.pool.id(("dup2", m, X, r, l, seen[l])); seen[l] += 1
                                self.cnf.extend([[-y, l], [y, -l]]); l = y
                            else:
                                seen[l] = 1
                            lits.append(l)
                    self.ctr[(m, X, r)] = (len(lits), unary(self.cnf, self.pool, lits, ("u", m, X, r)))

    def assumptions(self, comp_p, comp_q):
        A = []
        for m, comp in ((self.p, comp_p), (self.q, comp_q)):
            for X in range(4):
                for r, val in enumerate(comp[X]):
                    size, u = self.ctr[(m, X, r)]
                    k = (size - val) // 2          # number of -1 entries in this coset
                    if k >= 1: A.append(u[k])
                    if k + 1 <= size: A.append(-u[k + 1])
        return A


def comp_reps(n, case, m):
    return comp13_reps(n, case, m)


def run_inc(n, p, q, k, shard="0/1", log=None, prod=False):
    case = rowsum_cases(n)[k]
    t0 = time.time()
    Rp, Rq = comp_reps(n, case, p), comp_reps(n, case, q)
    e = IncCase(n, case, p, q, prod)
    si, sk = map(int, shard.split("/"))
    sols, done = [], 0
    with Cadical153(bootstrap_with=e.cnf.clauses) as S:
        for ip, cp in enumerate(Rp):
            if ip % sk != si:
                continue
            for cq in Rq:
                A = e.assumptions(cp, cq)
                while S.solve(assumptions=A):
                    m_ = S.get_model(); Q = e.decode(m_)
                    assert verify(Q), "SAT model fails verification"
                    sols.append([X[: (n + 1) // 2] for X in Q])
                    S.add_clause(e.block(m_))
                done += 1
            if log:
                print(json.dumps({"case": case, "ip": ip, "done": done, "sols": len(sols), "secs": round(time.time() - t0, 1)}), flush=True)
    return {"n": n, "case": case, "p": p, "q": q, "shard": shard, "Rp": len(Rp), "Rq": len(Rq),
            "product_condition": prod, "combos_this_shard": done, "sols": sols, "N_found": len(sols), "secs": round(time.time() - t0, 1)}


import numpy as _np
from pysat.engines import Propagator
from pysat.solvers import Cadical195


class PSDProp(Propagator):
    """Programmatic (CAS) check: a complete block must have PSD <= 4n at every frequency, and two complete blocks
    must have PSD sum <= 4n. Violations become blocking clauses over exactly the blocks involved. Decisions fill
    blocks in order A, B, C, D so the checks fire early."""
    def __init__(self, enc):
        super().__init__()
        self.e, self.n, self.h = enc, enc.n, enc.h
        self.var2 = {}
        for X in range(4):
            for j, v in enumerate(enc.v[X]):
                self.var2[v] = (X, j)
        self.val = {}                      # var -> bool (True: x = -1)
        self.trail, self.lim = [], []
        self.cnt = [0] * 4
        self.pending = []
        self.bound = 4 * n_ if (n_ := enc.n) else 0
        k = _np.arange(1, self.h)
        self.C = _np.cos(2 * _np.pi * _np.outer(k, _np.arange(1, self.h)) / self.n)   # [freq, j>=1]
        self.psd_cache = {}
        self.checks = self.blocks_cut = self.pairs_cut = 0

    def setup_observe(self, solver):
        for v in self.var2:
            solver.observe(v)

    def on_assignment(self, lit, fixed=False):
        v = abs(lit)
        if v in self.var2 and v not in self.val:
            self.val[v] = lit > 0
            self.trail.append(v)
            X = self.var2[v][0]
            self.cnt[X] += 1
            if self.cnt[X] == self.h:
                self._check(X)

    def on_new_level(self):
        self.lim.append(len(self.trail))

    def on_backtrack(self, to):
        while len(self.lim) > to:
            m = self.lim.pop()
            while len(self.trail) > m:
                v = self.trail.pop()
                del self.val[v]
                self.cnt[self.var2[v][0]] -= 1

    def _psd(self, X):
        x = _np.array([-1.0 if self.val[v] else 1.0 for v in self.e.v[X]])
        return (x[0] + 2 * self.C @ x[1:]) ** 2

    def _lits(self, X):
        return [(-v if self.val[v] else v) for v in self.e.v[X]]    # clause literals falsified now

    def _check(self, X):
        self.checks += 1
        P = self._psd(X)
        if _np.any(P > self.bound + 1e-6):
            self.blocks_cut += 1
            self.pending.append(self._lits(X)); return
        for Y in range(4):
            if Y != X and self.cnt[Y] == self.h:
                if _np.any(P + self._psd(Y) > self.bound + 1e-6):
                    self.pairs_cut += 1
                    self.pending.append(self._lits(X) + self._lits(Y)); return

    def check_model(self, model):
        return not self.pending

    def decide(self):
        for X in range(4):
            if self.cnt[X] < self.h:
                for v in self.e.v[X]:
                    if v not in self.val:
                        return -v
        return 0

    def propagate(self):
        return []

    def provide_reason(self, lit):
        return []

    def has_clause(self):
        return bool(self.pending)

    def add_clause(self):
        return self.pending.pop() if self.pending else []


def run_inc2(n, p, q, k, shard="0/1", log=None, limit_secs=None, sample=None, seed=1):
    """As run_inc, with every CNF-expressible Williamson relation plus the PSD propagator."""
    import random
    case = rowsum_cases(n)[k]
    t0 = time.time()
    Rp, Rq = comp_reps(n, case, p), comp_reps(n, case, q)
    e = IncCase(n, case, p, q, prod=True)
    P = PSDProp(e)
    si, sk = map(int, shard.split("/"))
    combos = [(ip, iq) for ip in range(len(Rp)) if ip % sk == si for iq in range(len(Rq))]
    if sample:
        random.seed(seed); combos = random.sample(combos, min(sample, len(combos)))
    sols, done, times = [], 0, []
    with Cadical195(bootstrap_with=e.cnf.clauses) as S:
        S.connect_propagator(P); P.setup_observe(S)
        for ip, iq in combos:
            A = e.assumptions(Rp[ip], Rq[iq])
            ts = time.time()
            while S.solve(assumptions=A):
                m_ = S.get_model(); Q = e.decode(m_)
                assert verify(Q), "SAT model fails verification"
                sols.append([X[: (n + 1) // 2] for X in Q])
                S.add_clause(e.block(m_))
            times.append(time.time() - ts); done += 1
            if log and done % 500 == 0:
                print(json.dumps({"case": case, "done": done, "of": len(combos), "sols": len(sols),
                                  "secs": round(time.time() - t0, 1)}), flush=True)
            if limit_secs and time.time() - t0 > limit_secs:
                break
    times.sort()
    return {"n": n, "case": case, "p": p, "q": q, "shard": shard, "Rp": len(Rp), "Rq": len(Rq),
            "combos_total": len(combos), "combos_done": done, "sols": sols, "N_found": len(sols),
            "solve_median": times[len(times) // 2] if times else None, "solve_mean": sum(times) / max(1, len(times)),
            "solve_max": times[-1] if times else None, "psd_checks": P.checks, "psd_block_cuts": P.blocks_cut,
            "psd_pair_cuts": P.pairs_cut, "secs": round(time.time() - t0, 1)}

if __name__ == "__main__":
    mode, n = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    if mode == "enum":
        print(json.dumps({"n": n, "N_sat": enum_all(n), "secs": round(time.time() - t0, 1)}), flush=True)
    elif mode == "inc":
        p, q, k = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        shard = sys.argv[sys.argv.index("--shard") + 1] if "--shard" in sys.argv else "0/1"
        out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
        res = run_inc(n, p, q, k, shard, log=True, prod="--prod" in sys.argv)
        print(json.dumps({kk: v for kk, v in res.items() if kk != "sols"}), flush=True)
        if out: json.dump(res, open(out, "w"))
    elif mode == "inc2":
        p, q, k = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        g = lambda f, d=None, t=str: t(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
        res = run_inc2(n, p, q, k, g("--shard", "0/1"), log=True, limit_secs=g("--limit", None, float),
                       sample=g("--sample", None, int))
        print(json.dumps({kk: v for kk, v in res.items() if kk != "sols"}), flush=True)
        if g("--out"): json.dump(res, open(g("--out"), "w"))
    elif mode == "case":
        k = int(sys.argv[3]); case = rowsum_cases(n)[k]
        if "--z13" not in sys.argv:
            e = Enc(n, case)
            with Cadical153(bootstrap_with=e.cnf.clauses) as S:
                r = S.solve()
                res = {"n": n, "case": case, "split": None, "sat": r, "secs": round(time.time() - t0, 1)}
                if r:
                    Q = e.decode(S.get_model()); res["verified"] = verify(Q); res["Q"] = Q
            print(json.dumps(res), flush=True)
        else:
            q = int(sys.argv[sys.argv.index("--q") + 1]) if "--q" in sys.argv else 13
            reps = comp13_reps(n, case, q)
            shard = sys.argv[sys.argv.index("--shard") + 1] if "--shard" in sys.argv else "0/1"
            si, sk = map(int, shard.split("/"))
            done = sat = 0
            for idx, comp in enumerate(reps):
                if idx % sk != si:
                    continue
                e = Enc(n, case); e.add_z13(comp)
                with Cadical153(bootstrap_with=e.cnf.clauses) as S:
                    r = S.solve()
                done += 1
                if r:
                    sat += 1
                    print(json.dumps({"SAT": True, "case": case, "comp": comp}), flush=True)
                if done % 200 == 0:
                    print(json.dumps({"progress": done, "of_shard": len(range(si, len(reps), sk)), "sat": sat,
                                      "secs": round(time.time() - t0, 1)}), flush=True)
            print(json.dumps({"n": n, "case": case, "split": "z13", "reps": len(reps), "shard": shard,
                              "done": done, "sat": sat, "secs": round(time.time() - t0, 1)}), flush=True)
