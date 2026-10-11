r"""Stage 1 by simulated annealing, for the y-screen (report/paley_523_answer.md, section 5) [M for every y printed].

Find y_k in the letter set L on Z_F (k = 0..3) with sum_k PAF(y_k)(h) = -4q for h != 0 and a fixed lag-0 composition
(n_a entries with |y| = a), as in theory/q523_stage1.py.  The composition is kept by the moves (sign flips and swaps of
two entries anywhere), so lag 0 holds throughout; the energy is sum_{h=1}^{F/2} (sum_k PAF_k(h) + 4q)^2.  At energy 0
the row sums satisfy sum_k r_k^2 = 4, so they are (2, 0, 0, 0) up to order and sign; blocks are reordered and negated
to (2, 0, 0, 0).  Every y printed is re-verified exactly (all lags, rows) and canonicalised (per block: the least of its
translates, reversals; this is a symmetry of the slice, see section 5) so that duplicates can be dropped.
usage: python3 theory/q523_ysa.py F q s --comp a1:n1,... [--runs 200] [--iters 3000000] [--seed 0] [--out DIR]
       [--want 4] [--secs 600] [--T0 16q^2] [--T1 2]"""
import sys, os, json, time, math
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_stage1 import letters

@njit(cache=True)
def _paf(Y, F, H):
    P = np.zeros((Y.shape[0], H + 1), np.int64)
    for k in range(Y.shape[0]):
        for h in range(1, H + 1):
            s = 0
            for g in range(F): s += Y[k, g] * Y[k, (g + h) % F]
            P[k, h] = s
    return P

@njit(cache=True)
def _set(Y, S, F, H, k, g, new):
    # change Y[k, g] to new and update the lag sums S (h = 1..H): dPAF(h) = (new - old)(y_{g+h} + y_{g-h})
    d = new - Y[k, g]
    if d == 0: return
    for h in range(1, H + 1):
        S[h] += d * (Y[k, (g + h) % F] + Y[k, (g - h) % F])
    Y[k, g] = new

@njit(cache=True)
def _energy(S, tgt, H):
    e = 0
    for h in range(1, H + 1):
        r = S[h] - tgt
        e += r * r
    return e

@njit(cache=True)
def sa(Y0, F, tgt, iters, T0, T1, seed):
    """Metropolis on E with sign flips and (signed) swaps; returns the final state and energy."""
    np.random.seed(seed)
    Y = Y0.copy(); H = F // 2; nb = Y.shape[0]
    P = _paf(Y, F, H); S = np.zeros(H + 1, np.int64)
    for h in range(1, H + 1):
        for k in range(nb): S[h] += P[k, h]
    E = _energy(S, tgt, H)
    lr = math.log(T1 / T0)
    for it in range(iters):
        if E == 0: break
        temp = T0 * math.exp(lr * it / iters)
        if np.random.random() < 0.4:
            k = np.random.randint(nb); g = np.random.randint(F); old = Y[k, g]
            _set(Y, S, F, H, k, g, -old); En = _energy(S, tgt, H)
            if En <= E or np.random.random() < math.exp((E - En) / temp): E = En
            else: _set(Y, S, F, H, k, g, old)
        else:
            k1 = np.random.randint(nb); g1 = np.random.randint(F); k2 = np.random.randint(nb); g2 = np.random.randint(F)
            if k1 == k2 and g1 == g2: continue
            a = Y[k1, g1]; b = Y[k2, g2]
            sg = 1 if np.random.random() < 0.5 else -1
            na = sg * b; nbv = sg * a                                  # (a, b) -> (+-b, +-a): the multiset of |y| is kept
            if na == a and nbv == b: continue
            _set(Y, S, F, H, k1, g1, na); _set(Y, S, F, H, k2, g2, nbv)
            En = _energy(S, tgt, H)
            if En <= E or np.random.random() < math.exp((E - En) / temp): E = En
            else:
                _set(Y, S, F, H, k2, g2, b); _set(Y, S, F, H, k1, g1, a)
    return Y, E

def canon_block(r):
    F = len(r); best = None
    for rr in (list(r), list(r[::-1])):
        for t in range(F):
            c = tuple(rr[t:] + rr[:t])
            if best is None or c < best: best = c
    return best

def canonical(Y):
    """(2,0,0,0) row order and sign, then each block canonicalised; zero-row blocks sorted (and up to negation)."""
    Y = [list(r) for r in Y]
    rs = [sum(r) for r in Y]; i = next(j for j in range(4) if abs(rs[j]) == 2)
    b0 = Y[i] if rs[i] == 2 else [-x for x in Y[i]]
    rest = [Y[j] for j in range(4) if j != i]
    rest = [min(canon_block(r), canon_block([-x for x in r])) for r in rest]
    return [list(canon_block(b0))] + [list(r) for r in sorted(rest)]

def verify(Y, F, q):
    m = F * q; T0 = 4 * m - 4 * (q - 1)
    P = [sum(sum(Y[k][g] * Y[k][(g + h) % F] for g in range(F)) for k in range(4)) for h in range(F)]
    return P[0] == T0 and all(P[h] == -4 * q for h in range(1, F)) and [sum(r) for r in Y] == [2, 0, 0, 0]

def initial(F, q, s, comp, rng):
    L, f = letters(q, s)
    mags = []
    for a, n in comp.items(): mags += [a] * n
    assert len(mags) == 4 * F, (len(mags), 4 * F)
    mags = np.array(mags); rng.shuffle(mags)
    sg = rng.choice([-1, 1], 4 * F)
    return (mags * sg).reshape(4, F).astype(np.int64)

def main(F, q, s, comp, runs, iters, seed, out, want, secs, T0=None, T1=2.0):
    rng = np.random.default_rng(seed); tgt = -4 * q; found = {}; t0 = time.time()
    T0 = T0 or float(16 * q * q)
    for r in range(runs):
        if len(found) >= want or time.time() - t0 > secs: break
        Y0 = initial(F, q, s, comp, rng)
        Y, E = sa(Y0, F, tgt, iters, T0, T1, int(rng.integers(1 << 31)))
        rec = dict(run=r, E=int(E), secs=round(time.time() - t0, 1))
        if E == 0:
            Yc = canonical(Y.tolist()); ok = verify(Yc, F, q)
            key = json.dumps(Yc)
            rec.update(verified=ok, new=key not in found)
            if ok and key not in found:
                found[key] = Yc
                if out:
                    tag = "_".join(f"{a}-{n}" for a, n in sorted(comp.items()))
                    fn = os.path.join(out, f"y_sa_{tag}_s{seed}_{len(found)}.json")
                    json.dump(dict(F=F, q=q, s=s, comp={str(a): n for a, n in comp.items()}, y=Yc, verified=True,
                                   method="q523_ysa", seed=seed, run=r), open(fn, "w"))
                    rec["file"] = fn
        print(json.dumps(rec), flush=True)
    return found

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--runs": "200", "--iters": "3000000", "--seed": "0", "--comp": "", "--out": "",
                             "--want": "4", "--secs": "600", "--T0": "0", "--T1": "2"}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    comp = {int(t.split(':')[0]): int(t.split(':')[1]) for t in opt["--comp"].split(',')}
    if opt["--out"]: os.makedirs(opt["--out"], exist_ok=True)
    main(F, q, s, comp, int(opt["--runs"]), int(opt["--iters"]), int(opt["--seed"]), opt["--out"], int(opt["--want"]),
         float(opt["--secs"]), float(opt["--T0"]) or None, float(opt["--T1"]))

def repair(Y, F, q, s, free, time_limit=20.0, workers=1, seed=0):
    """Exact repair of a near-miss: CP-SAT over the blocks in `free` (their multiset of |y| kept), the others fixed."""
    from ortools.sat.python import cp_model
    L, f = letters(q, s); Y = [list(map(int, r)) for r in Y]; H = F // 2; tgt = -4 * q
    fixed = [k for k in range(4) if k not in free]
    R = [tgt - sum(sum(Y[k][g] * Y[k][(g + h) % F] for g in range(F)) for k in fixed) for h in range(H + 1)]
    M = cp_model.CpModel()
    b = {(k, g, i): M.NewBoolVar("") for k in free for g in range(F) for i in range(len(L))}
    y = {(k, g): M.NewIntVar(-max(L), max(L), "") for k in free for g in range(F)}
    for k in free:
        for g in range(F):
            M.AddExactlyOne([b[k, g, i] for i in range(len(L))]); M.Add(y[k, g] == sum(L[i] * b[k, g, i] for i in range(len(L))))
    for a in sorted({abs(x) for x in L}):
        n = sum(1 for k in free for g in range(F) if abs(Y[k][g]) == a)
        M.Add(sum(b[k, g, i] for k in free for g in range(F) for i in range(len(L)) if abs(L[i]) == a) == n)
    sq = cp_model.Domain.FromValues(sorted({x * z for x in L for z in L}))
    for h in range(1, H + 1):
        terms = []
        for k in free:
            for g in range(F):
                p = M.NewIntVarFromDomain(sq, ""); M.AddMultiplicationEquality(p, [y[k, g], y[k, (g + h) % F]]); terms.append(p)
        M.Add(sum(terms) == R[h])
    for k in free:
        for g in range(F): M.AddHint(y[k, g], Y[k][g])
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    st = S.StatusName(S.Solve(M))
    if st in ("FEASIBLE", "OPTIMAL"):
        Z = [r[:] for r in Y]
        for k in free:
            for g in range(F): Z[k][g] = S.Value(y[k, g])
        return st, Z
    return st, None
