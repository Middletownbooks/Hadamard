"""Twisted-multiplier screen for b=0 Goethals-Seidel cores of prime order t (four +-1 circulants, sum PAF = 4t*delta).

Ansatz: U = <u> <= (Z/t)* of order d, g = (pi, eps) a signed permutation of the blocks (B4) with g^d = 1, and
    X_{pi(i)}(u x) = eps_i X_i(x)     (labels: 16 = sign pattern of (A,B,C,D) at a position; label(ux) = g(label(x))).
g = id is the usual multiplier ansatz.  Position 0 must carry a g-fixed label.  One equation per class of <U,-1>.
Per (t, d, conjugacy class of g): orbits, classes, exact block-sum feasibility (ILP vs all (a,b,c,d) with
sum a^2 = 4t), Gaussian log2 E, n bits = log2(#configs) - log2 E, log2 E per symmetry class.
Usage: gs_twist.py screen t [d ...]  |  gs_twist.py validate
"""
import sys, json, math, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, "scripts")
from tmat_twist import B4, cycle_type, gen_of_order, class_reps

def lab_map(g):
    p, e = g; m = []
    for l in range(16):
        v = [(-1 if (l >> X) & 1 else 1) for X in range(4)]
        w = [0] * 4
        for X in range(4): w[p[X]] = e[X] * v[X]
        m.append(sum(1 << X for X in range(4) if w[X] < 0))
    return tuple(m)

def comp(a, b): return tuple(a[b[l]] for l in range(16))
def order16(mp):
    k, cur, idm = 1, mp, tuple(range(16))
    while cur != idm: cur = comp(mp, cur); k += 1
    return k
LM = {g: lab_map(g) for g in B4}
def centralizer(g):
    mg = LM[g]; return sum(1 for h in B4 if comp(LM[h], mg) == comp(mg, LM[h]))

SIGNS = np.array([[(-1 if (l >> X) & 1 else 1) for X in range(4)] for l in range(16)])   # 16 x 4

def setup(t, d, g):
    u = gen_of_order(t, d); mg = LM[g]
    pows = [tuple(range(16))]
    for j in range(1, d): pows.append(comp(mg, pows[-1]))
    seen, orbits = {0}, []
    for x in range(1, t):
        if x in seen: continue
        o = [x * pow(u, j, t) % t for j in range(d)]; seen.update(o); orbits.append(o)
    G = np.array(pows)
    fix0 = [l for l in range(16) if mg[l] == l]
    seen, reps = set(), []
    for x in range(1, t):
        if x in seen: continue
        seen.update({(sg * x * pow(u, j, t)) % t for j in range(d) for sg in (1, -1)}); reps.append(x)
    return orbits, G, fix0, np.array(reps)

def rowsum_feasible(t, orbits, G, fix0):
    V = sorted({tuple(SIGNS[G[:, b]].sum(0)) for b in range(16)})
    Vm = np.array(V).T; nb = len(orbits); r = int(math.isqrt(4 * t))
    r -= (r % 2 == 0)                       # block sums are odd
    targets = [np.array(q) for q in itertools.product(range(-r, r + 1, 2), repeat=4) if sum(x * x for x in q) == 4 * t]
    A = np.vstack([Vm, np.ones((1, Vm.shape[1]))])
    for l0 in fix0:
        c0 = SIGNS[l0]
        for q in targets:
            b = np.concatenate([q - c0, [nb]]).astype(float)
            res = milp(c=np.zeros(Vm.shape[1]), constraints=LinearConstraint(A, b, b),
                       integrality=np.ones(Vm.shape[1]), bounds=Bounds(0, nb))
            if res.status == 0: return True, [int(x) for x in q]
    return False, None

def sample_F(t, orbits, G, fix0, reps, N, rng, batch=20000):
    pos = np.array(orbits); out = []
    for s in range(0, N, batch):
        n = min(batch, N - s)
        base = rng.integers(0, 16, (n, len(orbits)))
        lab = np.empty((n, t), np.int64); lab[:, 0] = rng.choice(fix0, n)
        for j in range(G.shape[0]): lab[:, pos[:, j]] = G[j][base]
        tot = np.zeros((n, t))
        for X in range(4):
            x = SIGNS[lab, X].astype(float)
            F = np.fft.rfft(x, axis=1); tot += np.fft.irfft(np.abs(F) ** 2, n=t, axis=1)
        out.append(np.rint(tot[:, reps]).astype(np.int64))
    return np.vstack(out)

def estimate(t, d, g, N=100000, seed=1):
    orbits, G, fix0, reps = setup(t, d, g)
    rec = {"t": t, "d": d, "g": cycle_type(g), "orbits": len(orbits) + 1, "classes": len(reps), "fix0": len(fix0)}
    if not fix0: rec.update(feasible=False, why="no g-fixed label at 0"); return rec
    ok, q = rowsum_feasible(t, orbits, G, fix0)
    if not ok: rec.update(feasible=False, why="row sums"); return rec
    rec["example_rowsums"] = q
    F = sample_F(t, orbits, G, fix0, reps, N, np.random.default_rng(seed)); C = F.shape[1]
    mu = F.mean(0); S = np.atleast_2d(np.cov(F.T))
    d_ = F - F[0]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d_[:, c])))) for c in range(C)])
    on = all((-F[0, c]) % steps[c] == 0 for c in range(C))
    sign, logdet = np.linalg.slogdet(S)
    if sign <= 0 or not on:
        rec.update(feasible=bool(on), why="degenerate cov" if on else "zero off lattice"); return rec
    q_ = mu @ np.linalg.solve(S, mu)
    log2P = (np.log(steps).sum() - 0.5 * C * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q_) / math.log(2)
    nconf = 4 * len(orbits) + math.log2(len(fix0)); log2E = nconf + log2P
    sym = math.log2((t - 1) / d * centralizer(g))
    rec.update(feasible=True, log2E=round(log2E, 1), n_bits=round(nconf - log2E, 1), per_class=round(log2E - sym, 1))
    return rec

def brute(t, d, g, batch=1 << 16):
    orbits, G, fix0, reps = setup(t, d, g); pos = np.array(orbits); nb = len(orbits); cnt = 0
    total = 16 ** nb
    for l0 in fix0:
        for s in range(0, total, batch):
            idx = np.arange(s, min(total, s + batch))
            base = np.stack([(idx >> (4 * o)) & 15 for o in range(nb)], 1)
            lab = np.empty((len(idx), t), np.int64); lab[:, 0] = l0
            for j in range(G.shape[0]): lab[:, pos[:, j]] = G[j][base]
            tot = np.zeros((len(idx), t))
            for X in range(4):
                x = SIGNS[lab, X].astype(float); F = np.fft.rfft(x, axis=1)
                tot += np.fft.irfft(np.abs(F) ** 2, n=t, axis=1)
            cnt += int((np.abs(tot[:, 1:]).max(1) < 0.5).sum())
    return cnt

def twisted_reps():
    return {ct: g for ct, g in class_reps().items()}

if __name__ == "__main__":
    if sys.argv[1] == "validate":
        for t, d in [(37, 9), (61, 15), (109, 27), (73, 9), (31, 3), (43, 7), (37, 3), (61, 5)]:
            for ct, g in sorted(twisted_reps().items()):
                if d % order16(LM[g]) != 0: continue
                orbits, G, fix0, reps = setup(t, d, g)
                if 16 ** len(orbits) * max(1, len(fix0)) > 2 ** 22: continue
                b = brute(t, d, g); e = estimate(t, d, g, 50000)
                print(json.dumps({"t": t, "d": d, "g": ct, "orbits": len(orbits) + 1, "brute": b,
                                  "brute_log2": round(math.log2(b), 1) if b else None, "est_log2E": e.get("log2E"),
                                  "feasible": e.get("feasible"), "why": e.get("why")}), flush=True)
    elif sys.argv[1] == "screen":
        t = int(sys.argv[2]); ds = [int(x) for x in sys.argv[3:]] or [x for x in range(2, t) if (t - 1) % x == 0 and (t - 1) // x >= 3]
        for d in ds:
            for ct, g in sorted(twisted_reps().items()):
                if d % order16(LM[g]) != 0: continue
                print(json.dumps(estimate(t, d, g)), flush=True)
