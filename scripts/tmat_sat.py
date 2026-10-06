"""SAT+CAS search for K-invariant T-matrices of prime order t (K <= (Z/t)* of odd order k).

GS form: four +-1 circulants A,B,C,D (A,B,C,D) = H4 (T1..T4), K-invariant, with
  (1) a_j b_j c_j d_j = +1 for all j        (product condition; XOR of the four '-1' bits = 0)
  (2) sum_X PAF_X(s) = 0 for s != 0          (exactly 2t disagreements x_i != x_{i+s} over the 4 blocks)
  (3) row sums (a,b,c,d) with a^2+b^2+c^2+d^2 = 4t   (normalised: |a|>=|b|>=|c|>=|d|, a,b,c > 0; negating two
      blocks and permuting blocks preserve (1),(2))
CAS (IPASIR-UP propagator): complete blocks must satisfy PSD_X(w) <= 4t; pairs and triples of complete blocks
their PSD sums <= 4t, at every frequency class.
Usage:
  tmat_sat.py solve t k [--all] [--limit SECS]   # first solution (or all, with blocking clauses) over all row-sum cases
"""
import sys, json, time, math, itertools
import numpy as np
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical195
from pysat.engines import Propagator
sys.path.insert(0, "scripts")
from tmat_search import subgroup


def orbit_data(t, k):
    K = subgroup(t, k)
    seen, orbits = set(), []
    for a in range(t):
        if a in seen: continue
        o = sorted({a * h % t for h in K}); seen.update(o); orbits.append(o)
    oid = [0] * t
    for i, o in enumerate(orbits):
        for x in o: oid[x] = i
    reps, seen = [], set()
    for o in orbits[1:]:
        if o[0] in seen: continue
        reps.append(o[0]); seen.update(o); seen.update((-x) % t for x in o)
    return K, orbits, oid, reps


def rowsum_cases(t):
    out = []
    r = int(math.isqrt(4 * t))
    for a in range(1, r + 1, 2):
        for b in range(1, a + 1, 2):
            for c in range(1, b + 1, 2):
                rem = 4 * t - a * a - b * b - c * c
                if rem <= 0: continue
                d = math.isqrt(rem)
                if d * d == rem and d % 2 == 1 and d <= c:
                    out.append((a, b, c, d)); out.append((a, b, c, -d))
    return out


class TEnc:
    def __init__(self, t, k, case):
        self.t, self.k, self.case = t, k, case
        self.K, self.orbits, self.oid, self.reps = orbit_data(t, k)
        self.m = len(self.orbits)
        self.pool, self.cnf = IDPool(), CNF()
        self.v = [[self.pool.id(("v", X, o)) for o in range(self.m)] for X in range(4)]
        self._product(); self._rows(); self._lags()

    def card_eq(self, lits, bound):
        if bound < 0 or bound > len(lits):
            self.cnf.append([]); return
        self.cnf.extend(CardEnc.equals(lits=lits, bound=bound, vpool=self.pool, encoding=EncType.totalizer).clauses)

    def _copies(self, lits, tag):
        cnt, out = {}, []
        for l in lits:
            if l in cnt:
                y = self.pool.id(("dup", tag, l, cnt[l])); cnt[l] += 1
                self.cnf.extend([[-y, l], [y, -l]]); out.append(y)
            else:
                cnt[l] = 1; out.append(l)
        return out

    def _product(self):
        for o in range(self.m):
            vs = [self.v[X][o] for X in range(4)]
            for bits in itertools.product((0, 1), repeat=4):
                if sum(bits) % 2 == 1:          # forbid odd number of -1 entries (assignment 'bits' -> clause)
                    self.cnf.append([(-v if b else v) for v, b in zip(vs, bits)])

    def _rows(self):
        for X, s in enumerate(self.case):
            wneg = (self.t - s) // 2           # number of -1 entries
            lits = []
            for o, orb in enumerate(self.orbits):
                lits += [self.v[X][o]] * len(orb)
            self.card_eq(self._copies(lits, ("row", X)), wneg)

    def _lags(self):
        t, oid = self.t, self.oid
        for s in self.reps:
            lits = []
            for X in range(4):
                for x in range(t):
                    a, b = oid[x], oid[(x + s) % t]
                    if a == b: continue
                    key = (min(a, b), max(a, b))
                    y = self.pool.id(("x", X, key))
                    if not hasattr(self, "_xs"): self._xs = set()
                    if (X, key) not in self._xs:
                        p, q = self.v[X][key[0]], self.v[X][key[1]]
                        self.cnf.extend([[-y, p, q], [-y, -p, -q], [y, -p, q], [y, p, -q]]); self._xs.add((X, key))
                    lits.append(y)
            self.card_eq(self._copies(lits, ("lag", s)), 2 * t)

    def decode(self, model):
        ms = set(l for l in model if l > 0)
        Q = np.ones((4, self.t), dtype=np.int64)
        for X in range(4):
            for o, orb in enumerate(self.orbits):
                if self.v[X][o] in ms: Q[X, orb] = -1
        return Q

    def block(self, model):
        ms = set(l for l in model if l > 0)
        return [(-l if l in ms else l) for X in range(4) for l in self.v[X]]


def verify_gs_product(Q):
    t = Q.shape[1]
    if not (np.prod(Q, axis=0) == 1).all(): return False
    F = np.fft.fft(Q, axis=1); tot = np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real).astype(np.int64)
    return bool(tot[0] == 4 * t and (tot[1:] == 0).all())


def to_T(Q):
    H = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])
    T = (H @ Q) // 4
    return T


class TProp(Propagator):
    def __init__(self, e):
        super().__init__()
        self.e, self.t, self.m = e, e.t, e.m
        self.var2 = {v: (X, o) for X in range(4) for o, v in enumerate(e.v[X])}
        self.val, self.trail, self.lim = {}, [], []
        self.cnt = [0] * 4
        self.pending = []
        self.bound = 4 * e.t
        # frequency class reps = lag reps (same orbit structure); Fourier of orbit indicator at each rep
        w = np.array(e.reps)
        S = np.zeros((len(w), e.m), dtype=complex)
        for o, orb in enumerate(e.orbits):
            S[:, o] = np.exp(2j * np.pi * np.outer(w, orb) / e.t).sum(1)
        self.S = S
        self.checks = self.cuts = 0

    def setup_observe(self, solver):
        for v in self.var2: solver.observe(v)

    def on_assignment(self, lit, fixed=False):
        v = abs(lit)
        if v in self.var2 and v not in self.val:
            self.val[v] = lit > 0; self.trail.append(v)
            X = self.var2[v][0]; self.cnt[X] += 1
            if self.cnt[X] == self.m: self._check(X)

    def on_new_level(self):
        self.lim.append(len(self.trail))

    def on_backtrack(self, to):
        while len(self.lim) > to:
            mk = self.lim.pop()
            while len(self.trail) > mk:
                v = self.trail.pop(); del self.val[v]; self.cnt[self.var2[v][0]] -= 1

    def _psd(self, X):
        x = np.array([-1.0 if self.val[v] else 1.0 for v in self.e.v[X]])
        return np.abs(self.S @ x) ** 2

    def _lits(self, X):
        return [(-v if self.val[v] else v) for v in self.e.v[X]]

    def _check(self, X):
        self.checks += 1
        P = self._psd(X)
        if np.any(P > self.bound + 1e-6):
            self.cuts += 1; self.pending.append(self._lits(X)); return
        done = [Y for Y in range(4) if Y != X and self.cnt[Y] == self.m]
        for r in (1, 2):
            for Ys in itertools.combinations(done, r):
                tot = P + sum(self._psd(Y) for Y in Ys)
                if np.any(tot > self.bound + 1e-6):
                    self.cuts += 1
                    self.pending.append(self._lits(X) + [l for Y in Ys for l in self._lits(Y)]); return

    def has_clause(self):
        return bool(self.pending)

    def add_clause(self):
        return self.pending.pop()

    def check_model(self, model):
        return not self.pending

    def decide(self):
        for X in range(4):
            if self.cnt[X] < self.m:
                for v in self.e.v[X]:
                    if v not in self.val: return -v
        return 0


def solve(t, k, all_sols=False, limit=None, log=True):
    t0 = time.time(); out = {"t": t, "k": k, "cases": [], "solutions": []}
    for case in rowsum_cases(t):
        e = TEnc(t, k, case); P = TProp(e); ns = 0; tc = time.time()
        with Cadical195(bootstrap_with=e.cnf.clauses) as S:
            S.connect_propagator(P); P.setup_observe(S)
            while S.solve():
                Q = e.decode(S.get_model())
                assert verify_gs_product(Q), "model fails verification"
                ns += 1; out["solutions"].append({"case": case, "T": to_T(Q).tolist()})
                if log: print(json.dumps({"FOUND": True, "t": t, "k": k, "case": case, "secs": round(time.time() - t0, 1)}), flush=True)
                if not all_sols: break
                S.add_clause(e.block(S.get_model()))
                if limit and time.time() - t0 > limit: break
        out["cases"].append({"case": case, "sols": ns, "secs": round(time.time() - tc, 1), "psd_checks": P.checks, "psd_cuts": P.cuts})
        if log: print(json.dumps(out["cases"][-1]), flush=True)
        if out["solutions"] and not all_sols: break
        if limit and time.time() - t0 > limit:
            out["timeout"] = True; break
    out["secs"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    mode, t, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    g = lambda f, d=None, ty=str: ty(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
    if mode == "solve":
        res = solve(t, k, "--all" in sys.argv, g("--limit", None, float))
        print(json.dumps({kk: v for kk, v in res.items() if kk != "solutions"} | {"n_solutions": len(res["solutions"])}))
        if g("--out"): json.dump(res, open(g("--out"), "w"))
