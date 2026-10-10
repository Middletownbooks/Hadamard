r"""Positive control for the stage-1 CP-SAT models of theory/p571_p19cpsat.py, generalised to the Paley(p) ansatz on
Z_n (b = 1): the published v = 419 core (data/core_v419_b1.txt) is an exact Paley(11) core on Z38 (theory/
paley_anatomy419.py).  Its gadget class (W, eps) is extracted, and both formulations must be FEASIBLE on it:
  model A: sum_k [PAF(A_k)(h) + (p-1) X_k(h)] = -4p - (p-1)^2 s(h), row sums r_k - (p-1) t_k, Lemma 6 at psi
           (sum_k What_k (2 Ahat_k + (p-2) What_k) in {sum_k uhat_k^2}), XOR literals for A A';
  model y: sum_k PAF(y_k)(h) = -4p with y = A + (p-1) W (pure S1).
Also: the published A, plugged in, satisfies both models (all Booleans fixed).
usage: python3 theory/p571_p19control419.py"""
import sys, os, json, time
import numpy as np
from ortools.sat.python import cp_model
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))

def model(W, eps, n, p, R, kind='A', fixA=None):
    m = cp_model.CpModel()
    b = [[m.NewBoolVar('') if W[k][g] == 0 else None for g in range(n)] for k in range(4)]
    if fixA is not None:
        for k in range(4):
            for g in range(n):
                if W[k][g] == 0: m.Add(b[k][g] == (1 if fixA[k][g] > 0 else 0))
    A = lambda k, g: int(eps[k][g % n] * W[k][g % n]) if W[k][g % n] != 0 else 2 * b[k][g % n] - 1
    yv = lambda k, g: int((eps[k][g % n] + p - 1) * W[k][g % n]) if W[k][g % n] != 0 else 2 * b[k][g % n] - 1
    val = A if kind == 'A' else yv
    t = [int(sum(W[k])) for k in range(4)]
    for k in range(4):
        m.Add(sum(val(k, g) for g in range(n)) == (R[k] - (p - 1) * t[k] if kind == 'A' else R[k]))
    xc = {}
    def pr(k, g1, g2):
        g1 %= n; g2 %= n
        f1, f2 = W[k][g1] != 0, W[k][g2] != 0
        if f1 and f2: return val(k, g1) * val(k, g2)
        if f1: return val(k, g1) * (2 * b[k][g2] - 1)
        if f2: return val(k, g2) * (2 * b[k][g1] - 1)
        key = (k, min(g1, g2), max(g1, g2))
        if key not in xc:
            z = m.NewBoolVar(''); m.AddBoolXOr([b[k][g1], b[k][g2], z.Not()]); xc[key] = z
        return 1 - 2 * xc[key]
    s = [sum(int((np.asarray(W[k]) * np.roll(W[k], -h)).sum()) for k in range(4)) for h in range(n)]
    for h in range(1, n // 2 + 1):
        terms = [pr(k, g, g + h) for k in range(4) for g in range(n)]
        if kind == 'A':
            for k in range(4):
                for g in np.nonzero(W[k])[0]:
                    terms.append((p - 1) * int(W[k][g]) * A(k, g + h)); terms.append((p - 1) * int(W[k][g]) * A(k, g - h))
            m.Add(sum(terms) == -4 * p - (p - 1) ** 2 * s[h])
        else:
            m.Add(sum(terms) == -4 * p)
    if kind == 'A' and n % 2 == 0:
        sg = [1 - 2 * (g % 2) for g in range(n)]
        What = [int(sum(sg[g] * W[k][g] for g in range(n))) for k in range(4)]
        expr = sum(What[k] * (2 * sum(sg[g] * A(k, g) for g in range(n)) + (p - 2) * What[k]) for k in range(4))
        Ks = [int(np.count_nonzero(W[k])) for k in range(4)]
        sets = [{q * q for q in range(-(n - K), n - K + 1) if (q - (n - K)) % 2 == 0} for K in Ks]
        tot = {0}
        for st_ in sets: tot = {x + y_ for x in tot for y_ in st_}
        m.AddLinearExpressionInDomain(expr, cp_model.Domain.FromValues(sorted(tot)))
    return m

def main():
    from paley_anatomy419 import load
    from paley_algebra import components
    v, p, n = 419, 11, 38; m_ = n * p
    X = load(); shifts = [0, 208, 208, 208]; Y = []
    for x, sh in zip(X, shifts):
        xs = np.roll(x, -sh); y = np.zeros((n, p), np.int64)
        for i in range(m_): y[i % n, i % p] = xs[i]
        Y.append(y)
    comps = [components(y, p) for y in Y]
    A = np.array([c[0] for c in comps]); W = np.array([c[2] for c in comps])
    eps = np.where(W != 0, np.sign(A * W), 0)
    yk = A + (p - 1) * W; R = [int(r) for r in yk.sum(1)]
    out = dict(K=int((W != 0).sum()), Nminus=int((eps == -1).sum()), rows=R)
    for kind in ('A', 'y'):
        for label, fix in (('published A fixed', A), ('free', None)):
            mm = model(W, eps, n, p, R, kind=kind, fixA=fix)
            solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = 600; solver.parameters.num_search_workers = 4
            t0 = time.time(); st = solver.Solve(mm)
            out[f'model_{kind}_{label}'] = dict(status=solver.StatusName(st), secs=round(time.time() - t0, 1))
            print(json.dumps({f'model_{kind}_{label}': out[f'model_{kind}_{label}']}), flush=True)
    print(json.dumps(out), flush=True)

if __name__ == "__main__":
    main()
