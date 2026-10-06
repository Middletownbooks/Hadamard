"""Twisted-symmetry screen for T-matrices of prime order t.

Ansatz: U = <u> <= (Z/t)* of order d, g = (pi, eps) in B4 (signed permutations of the 4 rows) with g^d = 1, and
    T_{pi(i)}(u x) = eps_i T_i(x)      i.e.  label(u x) = g(label(x)),  label = (row, sign) in 8 values.
Each nonzero U-orbit carries one free base label; position 0 must carry a g-fixed label. Sum of PAFs is invariant,
so the T-condition is one equation per class of <U, -1> on nonzero lags.
For each (t, d, conjugacy class of g): orbits, classes, exact row-sum feasibility (ILP), Gaussian log2 E
(Monte Carlo mean/cov of the class sums), and log2 E per symmetry class (/ ((t-1)/d * |C_B4(g)|)).
Usage: tmat_twist.py screen t [t ...]      |   tmat_twist.py validate
"""
import sys, json, math, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

# ---- B4 ----
def b4_elements():
    out = []
    for p in itertools.permutations(range(4)):
        for e in itertools.product((1, -1), repeat=4):
            out.append((p, e))
    return out

def g_label_map(g):
    p, e = g   # label l = r | (4 if sign<0); row r -> p[r], sign *= e[r]
    m = []
    for l in range(8):
        r, s = l & 3, (-1 if l & 4 else 1)
        r2, s2 = p[r], s * e[r]
        m.append(r2 | (4 if s2 < 0 else 0))
    return tuple(m)

def compose(a, b):  # label maps: (a o b)(l) = a[b[l]]
    return tuple(a[b[l]] for l in range(8))

def order_of(mp):
    k, cur, idm = 1, mp, tuple(range(8))
    while cur != idm:
        cur = compose(mp, cur); k += 1
    return k

def cycle_type(g):
    p, e = g; seen, ct = set(), []
    for r in range(4):
        if r in seen: continue
        L, sgn, x = 0, 1, r
        while x not in seen:
            seen.add(x); sgn *= e[x]; x = p[x]; L += 1
        ct.append((L, sgn))
    return tuple(sorted(ct))

B4 = b4_elements()
LMAP = {g: g_label_map(g) for g in B4}
def centralizer_size(g):
    mg = LMAP[g]
    return sum(1 for h in B4 if compose(LMAP[h], mg) == compose(mg, LMAP[h]))

def class_reps():
    reps = {}
    for g in B4:
        reps.setdefault(cycle_type(g), g)
    return reps

# ---- group / orbits ----
def prime_factors(n):
    f, p = set(), 2
    while p * p <= n:
        while n % p == 0: f.add(p); n //= p
        p += 1
    if n > 1: f.add(n)
    return f

def gen_of_order(t, d):
    for a in range(2, t):
        h = pow(a, (t - 1) // d, t)
        if pow(h, d, t) == 1 and all(pow(h, d // q, t) != 1 for q in prime_factors(d)):
            return h
    if d == 1: return 1
    raise ValueError

def setup(t, d, g):
    u = gen_of_order(t, d)
    mg = LMAP[g]
    pows = [tuple(range(8))]
    for j in range(1, d): pows.append(compose(mg, pows[-1]))
    seen, orbits = {0}, []
    for x in range(1, t):
        if x in seen: continue
        o = [x * pow(u, j, t) % t for j in range(d)]; seen.update(o); orbits.append(o)
    G = np.array(pows)                                    # d x 8 : label image under g^j
    fix0 = [l for l in range(8) if mg[l] == l]
    # lag classes: orbits of <u, -1> on nonzero
    seen, reps = set(), []
    for x in range(1, t):
        if x in seen: continue
        cl = {(sg * x * pow(u, j, t)) % t for j in range(d) for sg in (1, -1)}; seen.update(cl); reps.append(x)
    return u, orbits, G, fix0, np.array(reps)

ROWS = np.array([l & 3 for l in range(8)]); SGN = np.array([-1 if l & 4 else 1 for l in range(8)])

def orbit_rowvec(G, base):
    v = np.zeros(4, np.int64)
    for j in range(G.shape[0]):
        l = G[j, base]; v[ROWS[l]] += SGN[l]
    return v

def rowsum_feasible(t, orbits, G, fix0):
    V = sorted({tuple(orbit_rowvec(G, b)) for b in range(8)})
    Vm = np.array(V).T                                    # 4 x |V|
    nb = len(orbits); r = int(math.isqrt(t)); targets = []
    for q in itertools.product(range(-r, r + 1), repeat=4):
        if sum(x * x for x in q) == t: targets.append(np.array(q))
    for l0 in fix0:
        c0 = np.zeros(4, np.int64); c0[ROWS[l0]] = SGN[l0]
        A = np.vstack([Vm, np.ones((1, Vm.shape[1]))])
        for q in targets:
            b = np.concatenate([q - c0, [nb]]).astype(float)
            res = milp(c=np.zeros(Vm.shape[1]), constraints=LinearConstraint(A, b, b),
                       integrality=np.ones(Vm.shape[1]), bounds=Bounds(0, nb))
            if res.status == 0:
                return True
    return False

def sample_F(t, orbits, G, fix0, reps, N, rng, batch=20000):
    pos = np.array(orbits)                                # nb x d
    out = []
    for s in range(0, N, batch):
        n = min(batch, N - s)
        base = rng.integers(0, 8, (n, len(orbits)))
        lab = np.empty((n, t), np.int64)
        lab[:, 0] = rng.choice(fix0, n)
        for j in range(G.shape[0]):
            lab[:, pos[:, j]] = G[j][base]
        tot = np.zeros((n, t))
        for r in range(4):
            x = np.where(ROWS[lab] == r, SGN[lab], 0).astype(float)
            X = np.fft.rfft(x, axis=1); tot += np.fft.irfft(np.abs(X) ** 2, n=t, axis=1)
        out.append(np.rint(tot[:, reps]).astype(np.int64))
    return np.vstack(out)

def estimate(t, d, g, N=100000, seed=1):
    u, orbits, G, fix0, reps = setup(t, d, g)
    rec = {"t": t, "d": d, "g": cycle_type(g), "orbits": len(orbits) + 1, "classes": len(reps), "fix0": len(fix0)}
    if not fix0:
        rec.update(feasible=False, why="no g-fixed label for position 0"); return rec
    if not rowsum_feasible(t, orbits, G, fix0):
        rec.update(feasible=False, why="row sums"); return rec
    F = sample_F(t, orbits, G, fix0, reps, N, np.random.default_rng(seed))
    C = F.shape[1]
    mu = F.mean(0); S = np.cov(F.T).reshape(C, C) if C > 1 else np.array([[F[:, 0].var()]])
    d_ = F - F[0]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d_[:, c])))) for c in range(C)])
    on = all((-F[0, c]) % steps[c] == 0 for c in range(C))
    sign, logdet = np.linalg.slogdet(S)
    if sign <= 0 or not on:
        rec.update(feasible=bool(on), why="degenerate covariance" if on else "zero off lattice", log2E=None); return rec
    q = mu @ np.linalg.solve(S, mu)
    log2P = (np.log(steps).sum() - 0.5 * C * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
    nconf = len(orbits) * 3 + math.log2(len(fix0))
    log2E = nconf + log2P
    sym = math.log2((t - 1) / d * centralizer_size(g))
    rec.update(feasible=True, log2E=round(log2E, 1), n_bits=round(nconf - log2E, 1), per_class=round(log2E - sym, 1))
    return rec

def screen(t, Nsamp=100000):
    reps = class_reps()
    for d in sorted(x for x in range(2, t) if (t - 1) % x == 0):
        if (t - 1) // d < 3: continue
        for ct, g in sorted(reps.items()):
            if d % order_of(LMAP[g]) != 0: continue
            try:
                r = estimate(t, d, g, Nsamp)
            except Exception as ex:
                r = {"t": t, "d": d, "g": ct, "error": repr(ex)}
            print(json.dumps(r), flush=True)

# ---- brute-force validation (small t) ----
def brute(t, d, g):
    u, orbits, G, fix0, reps = setup(t, d, g)
    pos = np.array(orbits); cnt = 0
    for l0 in fix0:
        for base in itertools.product(range(8), repeat=len(orbits)):
            lab = np.empty(t, np.int64); lab[0] = l0
            for j in range(G.shape[0]): lab[pos[:, j]] = G[j][list(base)]
            T = np.zeros((4, t)); T[ROWS[lab], np.arange(t)] = SGN[lab]
            F = np.fft.fft(T, axis=1); tot = np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real)
            if (tot[1:] == 0).all(): cnt += 1
    return cnt

if __name__ == "__main__":
    if sys.argv[1] == "screen":
        for t in map(int, sys.argv[2:]): screen(t)
    elif sys.argv[1] == "validate":
        reps = class_reps()
        for t, d in [(13, 6), (19, 6), (31, 6), (37, 6), (43, 6), (29, 4), (37, 4), (19, 3), (31, 5)]:
            for ct, g in sorted(reps.items()):
                if d % order_of(LMAP[g]) != 0: continue
                u, orbits, G, fix0, rp = setup(t, d, g)
                if 8 ** len(orbits) * max(1, len(fix0)) > 3e6: continue
                b = brute(t, d, g); e = estimate(t, d, g, 50000)
                print(json.dumps({"t": t, "d": d, "g": ct, "orbits": len(orbits) + 1, "brute": b,
                                  "brute_log2": round(math.log2(b), 1) if b else None,
                                  "est_log2E": e.get("log2E"), "feasible": e.get("feasible"), "why": e.get("why")}), flush=True)
