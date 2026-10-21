r"""Independent check (shares no code with theory/bt_*) that a compressed level of a bordered-GS structure has no
solution [M]: CP-SAT on the compressed system.

Structure: G = Z_{m_1} x ... x Z_{m_k} (first coordinate slowest), K given by generators, Gamma per coset of K
(given as the value 4 sum_{chi in T} chi), one symmetry matrix list S.  For the subgroup N (generators), the
compression y_i(u) = sum_{g in u} x_i(g) on G/N, u running over cosets.  Constraints:
- each y_i(u) is a sum over the S-orbits of G inside u of +-(orbit points in u): domain from exact subset sums,
  and y_i is constant on the S-orbits of G/N;
- sum_i y_i(u) y_i(u + w) summed over u = F_N(w) for every coset w, F_N(w) = sum_{g in w} F(g);
- the row sums sum_u y_i(u) are given (a sorted pattern; blocks are interchangeable and negatable).
INFEASIBLE means: no core of the structure has these row sums.
usage: python3 theory/bt_xcheck_level.py SPEC.json  (keys: mods, Kgens, gamma (per coset, cosets numbered by the
       smallest element), gens, Ngens, rows (list of patterns))
"""
import sys, json, itertools, time
import numpy as np
from ortools.sat.python import cp_model


def main():
    spec = json.load(open(sys.argv[1]))
    mods = spec['mods']; k = len(mods)
    n = int(np.prod(mods))
    co = np.array(list(itertools.product(*[range(m) for m in mods])), dtype=np.int64)   # first slowest
    marr = np.array(mods)
    idx = lambda c: int(np.ravel_multi_index(tuple(np.mod(c, marr)), mods))

    def span(gens):
        H = {0}; fr = [0]
        while fr:
            x = fr.pop()
            for g in gens:
                y = idx(co[x] + np.array(g))
                if y not in H:
                    H.add(y); fr.append(y)
        return sorted(H)
    K = span(spec['Kgens']); N = span(spec['Ngens'])
    coset = -np.ones(n, np.int64); q = 0
    for g in range(n):
        if coset[g] < 0:
            for h in K:
                coset[idx(co[g] + co[h])] = q
            q += 1
    gamma = spec['gamma']
    F = np.array([-gamma[coset[g]] for g in range(n)], np.int64); F[0] = 4 * n
    # cosets of N
    cn = -np.ones(n, np.int64); nb = 0
    for g in range(n):
        if cn[g] < 0:
            for h in N:
                cn[idx(co[g] + co[h])] = nb
            nb += 1
    FN = np.zeros(nb, np.int64)
    for g in range(n):
        FN[cn[g]] += F[g]
    rep = [int(np.flatnonzero(cn == u)[0]) for u in range(nb)]
    addb = [[int(cn[idx(co[rep[u]] + co[rep[w]])]) for w in range(nb)] for u in range(nb)]
    # S-orbits on G
    perms = [np.array([idx(np.array(M) @ co[g]) for g in range(n)]) for M in spec['gens']]
    orb = -np.ones(n, np.int64); no = 0
    for g in range(n):
        if orb[g] < 0:
            stack = [g]; orb[g] = no
            while stack:
                x = stack.pop()
                for P in perms:
                    y = int(P[x])
                    if orb[y] < 0:
                        orb[y] = no; stack.append(y)
            no += 1
    # per coset u of N: multiset of (orbit, count in u); domain of y(u)
    dom = []
    for u in range(nb):
        cnt = {}
        for g in np.flatnonzero(cn == u):
            cnt[int(orb[g])] = cnt.get(int(orb[g]), 0) + 1
        sums = {0}
        for o, c in cnt.items():
            sums = {s + c for s in sums} | {s - c for s in sums}
        dom.append(sorted(sums))
    # S-orbits on G/N (y constant on them)
    borb = -np.ones(nb, np.int64); nbo = 0
    for u in range(nb):
        if borb[u] < 0:
            stack = [u]; borb[u] = nbo
            while stack:
                x = stack.pop()
                for P in perms:
                    y = int(cn[P[rep[x]]])
                    if borb[y] < 0:
                        borb[y] = nbo; stack.append(y)
            nbo += 1
    out = []
    for rows in spec['rows']:
        t0 = time.time()
        m = cp_model.CpModel()
        Y = [[m.NewIntVarFromDomain(cp_model.Domain.FromValues(dom[u]), f'y{i}_{u}') for u in range(nb)] for i in range(4)]
        for i in range(4):
            for u in range(nb):
                ru = int(np.flatnonzero(borb == borb[u])[0])
                if ru != u:
                    m.Add(Y[i][u] == Y[i][ru])
            m.Add(sum(Y[i]) == rows[i])
        lim = max(max(abs(v) for v in d) for d in dom)
        for w in range(nb):
            terms = []
            for i in range(4):
                for u in range(nb):
                    p = m.NewIntVar(-lim * lim, lim * lim, '')
                    m.AddMultiplicationEquality(p, [Y[i][u], Y[i][addb[u][w]]])
                    terms.append(p)
            m.Add(sum(terms) == int(FN[w]))
        sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = float(spec.get('secs', 600)); sol.parameters.num_workers = 1
        st = sol.Solve(m)
        res = dict(rows=rows, status=sol.StatusName(st), secs=round(time.time() - t0, 1), nb=nb, orbits=no)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            res['y'] = [[sol.Value(Y[i][u]) for u in range(nb)] for i in range(4)]
        print(json.dumps(res), flush=True)
        out.append(res)


if __name__ == '__main__':
    main()
