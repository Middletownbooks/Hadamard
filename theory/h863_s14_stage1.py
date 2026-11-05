r"""[M if decided] Exact stage-1 (H-level, G/Z43) gate of the S14 sibling at v = 863 (census key fbf4f510b2566552:
the S7 host plus tau = (q1, q2, u, t) -> (q2, q1, -u, t); set-up of both hosts by OpenAI Codex, construction Kumar /
EKLO).  Independent of the Codex scripts (theory/h3452_compress.py, h3452_crosscheck.py) and of bt_engine.

H = Z2^2 x Z5 (cells q*5 + u, q = 2 q1 + q2); tau-orbits of cells: (00,0), (11,0), {(00,1),(00,4)}, {(00,2),(00,3)},
{(11,1),(11,4)}, {(11,2),(11,3)}, {(01,u),(10,-u)} for u = 0..4: 11 per block.  y on an orbit takes values in
A = {e + 7 S} = {+-1, +-13, +-15, +-27, +-29, +-41, +-43}.
Q-level (enumerated exactly here, up to block order, per-block sign and the translation by (1,1), which swaps a and c):
  per block (a, b, c) = (z(00), z(01) = z(10), z(11)); sum a^2 + 2b^2 + c^2 = 872, sum 2(ab + bc) = 860,
  sum 2ac + 2b^2 = 860.
H-level: sum_i PAF(y_i)(h) = 2936 at 0, -516 at (0, u != 0), 172 at q != 0 (all 20 lags), with the Q-compressions of y
fixed to the Q-level solution.  One CP-SAT model per Q-level solution (products of orbit values as exact product
variables).  FEASIBLE solutions are re-checked by a direct PAF evaluation on the explicit 20-cell sequences.
usage: python3 theory/h863_s14_stage1.py [secs] [workers] [first] [last] [out.jsonl]"""
import sys, os, json, time, itertools as it
from collections import defaultdict
from ortools.sat.python import cp_model

A = [s * v for v in (1, 13, 15, 27, 29, 41, 43) for s in (1, -1)]
SW = [0, 2, 1, 3]
TAU = [SW[k // 5] * 5 + (-(k % 5)) % 5 for k in range(20)]
ORB = sorted({tuple(sorted({k, TAU[k]})) for k in range(20)})
OIDX = {k: n for n, o in enumerate(ORB) for k in o}
assert len(ORB) == 11


def q_solutions():
    lst = []
    for a in range(-29, 30, 2):
        for b in range(-21, 22, 2):
            for c in range(-29, 30, 2):
                e0 = a * a + 2 * b * b + c * c
                if e0 <= 872: lst.append((a, b, c, e0, 2 * (a * b + b * c), 2 * a * c + 2 * b * b))
    pairs = defaultdict(list)
    for x in lst:
        for y in lst:
            s = (x[3] + y[3], x[4] + y[4], x[5] + y[5])
            if s[0] <= 872: pairs[s].append((x[:3], y[:3]))
    def canon(bl):
        out = []
        for (a, b, c) in bl:
            if a + 2 * b + c < 0: a, b, c = -a, -b, -c
            out.append((a, b, c))
        return tuple(sorted(out))
    sols = set()
    for s, L in pairs.items():
        t = (872 - s[0], 860 - s[1], 860 - s[2])
        for l in L:
            for r in pairs.get(t, ()):
                bl = [l[0], l[1], r[0], r[1]]
                sols.add(min(canon(bl), canon([(c, b, a) for (a, b, c) in bl])))
    return sorted(sols)


def hpaf(y, h0):
    q0, u0 = divmod(h0, 5)
    return sum(y[k] * y[(((k // 5) ^ q0) * 5) + (k % 5 + u0) % 5] for k in range(20))


def target(h0):
    q0, u0 = divmod(h0, 5)
    return 2936 if h0 == 0 else (-516 if q0 == 0 else 172)


def solve(z, secs, workers):
    M = cp_model.CpModel()
    dom = cp_model.Domain.FromValues(A)
    Y = [[M.NewIntVarFromDomain(dom, 'y%d_%d' % (i, o)) for o in range(11)] for i in range(4)]
    for i, (a, b, c) in enumerate(z):
        cell = lambda k: Y[i][OIDX[k]]
        M.Add(sum(cell(0 * 5 + u) for u in range(5)) == a)
        M.Add(sum(cell(1 * 5 + u) for u in range(5)) == b)
        M.Add(sum(cell(3 * 5 + u) for u in range(5)) == c)
    prods = {}
    def prod(i, o1, o2):
        key = (i, min(o1, o2), max(o1, o2))
        if key not in prods:
            v = M.NewIntVar(-1849, 1849, ''); M.AddMultiplicationEquality(v, [Y[i][key[1]], Y[i][key[2]]])
            prods[key] = v
        return prods[key]
    done = set()
    for h0 in range(20):
        q0, u0 = divmod(h0, 5)
        cls = (q0 if q0 in (0, 3) else 1, min(u0, (-u0) % 5))
        if cls in done: continue
        done.add(cls)
        terms = defaultdict(int)
        for i in range(4):
            for k in range(20):
                k2 = (((k // 5) ^ q0) * 5) + (k % 5 + u0) % 5
                terms[(i, min(OIDX[k], OIDX[k2]), max(OIDX[k], OIDX[k2]))] += 1
        M.Add(sum(m * prod(*key) for key, m in terms.items()) == target(h0))
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = secs; S.parameters.num_workers = workers
    t = time.time(); st = S.StatusName(S.Solve(M))
    out = dict(z=z, rows=[a + 2 * b + c for a, b, c in z], status=st, secs=round(time.time() - t, 1),
               classes=len(done))
    if st in ('OPTIMAL', 'FEASIBLE'):
        ys = [[S.Value(Y[i][OIDX[k]]) for k in range(20)] for i in range(4)]
        ok = all(sum(hpaf(y, h0) for y in ys) == target(h0) for h0 in range(20))
        ok &= all(y[k] == y[TAU[k]] for y in ys for k in range(20)) and all(v in A for y in ys for v in y)
        out.update(y=ys, recheck=ok)
    return out


if __name__ == '__main__':
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 60
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    Z = q_solutions()
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    last = int(sys.argv[4]) if len(sys.argv) > 4 else len(Z)
    outp = sys.argv[5] if len(sys.argv) > 5 else None
    print(json.dumps(dict(q_solutions=len(Z))), flush=True)
    for n in range(first, min(last, len(Z))):
        r = solve(Z[n], secs, workers); r['index'] = n
        s = json.dumps(r); print(s, flush=True)
        if outp:
            with open(outp, 'a') as fo: fo.write(s + '\n')
