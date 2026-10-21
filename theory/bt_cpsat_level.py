r"""CP-SAT feasibility of compressed levels of census structures [M when INFEASIBLE; FEASIBLE gives a level solution].

Independent of theory/bt_engine and bt_chain (own group arithmetic, orbits, cosets, targets).  For a census record
(mods, K generators, Gamma per coset, symmetry generators) and an invariant subgroup N (or every subgroup generated
by a subset of the cyclic factors' natural subgroups), the compressed system on G/N:
  y_i(u) = sum over the S-orbits of G inside the coset u of +-(points of the orbit in u)   (exact domains),
  y_i constant on the S-orbits of G/N,
  sum_i sum_u y_i(u) y_i(u + w) = F_N(w) for every coset w   (F_N(w) = sum_{g in w} F(g)),
with symmetry breaking 0 <= r_1 <= r_2 <= r_3 <= r_4 (blocks share S, so they can be reordered and negated).
INFEASIBLE: no core of the structure exists [M].
usage: python3 theory/bt_cpsat_level.py RECORDS.jsonl OUT.jsonl [--keys k1,..] [--secs 600] [--min-logE x]
"""
import sys, os, json, time, itertools, argparse
import numpy as np
from ortools.sat.python import cp_model


class G0:
    def __init__(self, mods):
        self.mods = list(mods); self.n = int(np.prod(mods))
        self.co = np.array(list(itertools.product(*[range(m) for m in mods])), dtype=np.int64)
        self.marr = np.array(mods)

    def idx(self, c):
        return int(np.ravel_multi_index(tuple(np.mod(c, self.marr)), self.mods))

    def span(self, gens):
        H = {0}; fr = [0]
        while fr:
            x = fr.pop()
            for g in gens:
                y = self.idx(self.co[x] + np.array(g))
                if y not in H:
                    H.add(y); fr.append(y)
        return sorted(H)


def level_system(rec, Ngens):
    G = G0(rec['Pmods'] + rec['Omods']); n = G.n
    K = G.span(rec['Kgens']); N = G.span(Ngens)
    coset = -np.ones(n, np.int64); q = 0
    for g in range(n):
        if coset[g] < 0:
            for h in K:
                coset[G.idx(G.co[g] + G.co[h])] = q
            q += 1
    gamma = rec['gamma'] if 'gamma' in rec else ([12] + [-4] * 3)
    F = np.array([-gamma[coset[g]] for g in range(n)], np.int64); F[0] = 4 * n
    cn = -np.ones(n, np.int64); nb = 0
    for g in range(n):
        if cn[g] < 0:
            for h in N:
                cn[G.idx(G.co[g] + G.co[h])] = nb
            nb += 1
    FN = np.bincount(cn, weights=F, minlength=nb).astype(np.int64)
    rep = [int(np.flatnonzero(cn == u)[0]) for u in range(nb)]
    addb = [[int(cn[G.idx(G.co[rep[u]] + G.co[rep[w]])]) for w in range(nb)] for u in range(nb)]
    perms = [np.array([G.idx(np.array(M) @ G.co[g]) for g in range(n)]) for M in rec['gens']]
    orb = -np.ones(n, np.int64); no = 0
    for g in range(n):
        if orb[g] < 0:
            st = [g]; orb[g] = no
            while st:
                x = st.pop()
                for P in perms:
                    y = int(P[x])
                    if orb[y] < 0:
                        orb[y] = no; st.append(y)
            no += 1
    dom = []
    for u in range(nb):
        cnt = {}
        for g in np.flatnonzero(cn == u):
            cnt[int(orb[g])] = cnt.get(int(orb[g]), 0) + 1
        sums = {0}
        for c in cnt.values():
            sums = {s + c for s in sums} | {s - c for s in sums}
        dom.append(sorted(sums))
    borb = -np.ones(nb, np.int64); nbo = 0
    for u in range(nb):
        if borb[u] < 0:
            st = [u]; borb[u] = nbo
            while st:
                x = st.pop()
                for P in perms:
                    y = int(cn[P[rep[x]]])
                    if borb[y] < 0:
                        borb[y] = nbo; st.append(y)
            nbo += 1
    return dict(nb=nb, FN=FN, addb=addb, dom=dom, borb=borb, orbits=no)


def solve(sysd, secs=600, workers=1):
    nb, FN, addb, dom, borb = sysd['nb'], sysd['FN'], sysd['addb'], sysd['dom'], sysd['borb']
    m = cp_model.CpModel()
    Y = [[m.NewIntVarFromDomain(cp_model.Domain.FromValues(dom[u]), '') for u in range(nb)] for i in range(4)]
    for i in range(4):
        for u in range(nb):
            ru = int(np.flatnonzero(borb == borb[u])[0])
            if ru != u:
                m.Add(Y[i][u] == Y[i][ru])
    R = [sum(Y[i]) for i in range(4)]
    m.Add(R[0] >= 0)
    for i in range(3):
        m.Add(R[i] <= R[i + 1])
    lim = max(max(abs(v) for v in d) for d in dom)
    sq = {}
    for w in range(nb):
        terms = []
        for i in range(4):
            for u in range(nb):
                v = addb[u][w]
                key = (i, min(u, v), max(u, v))
                if key not in sq:
                    p = m.NewIntVar(-lim * lim, lim * lim, '')
                    m.AddMultiplicationEquality(p, [Y[i][u], Y[i][v]])
                    sq[key] = p
                terms.append(sq[key])
        m.Add(sum(terms) == int(FN[w]))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = float(secs); s.parameters.num_workers = workers
    st = s.Solve(m)
    out = dict(status=s.StatusName(st), secs=round(s.WallTime(), 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        out['y'] = [[s.Value(Y[i][u]) for u in range(nb)] for i in range(4)]
    return out


def natural_subgroups(rec):
    """subgroups N <= K generated by multiples of the coordinate axes (invariant under per-axis maps and, for equal
    2-power axes, usually under the coupled automorphisms; invariance is checked)"""
    mods = rec['Pmods'] + rec['Omods']; k = len(mods)
    G = G0(mods)
    K = set(G.span(rec['Kgens']))
    perms = [np.array([G.idx(np.array(M) @ G.co[g]) for g in range(G.n)]) for M in rec['gens']]
    divs = []
    for m_ in mods:
        divs.append([d for d in range(1, m_ + 1) if m_ % d == 0])
    out = []
    for choice in itertools.product(*divs):
        gens = []
        for i, d in enumerate(choice):
            v = [0] * k; v[i] = d % mods[i]; gens.append(v)
        N = G.span(gens)
        if len(N) in (1, G.n) or not set(N) <= K or len(N) == len(K):
            continue
        Ns = set(N)
        if not all(set(int(P[x]) for x in N) == Ns for P in perms):
            continue
        out.append((len(N), gens))
    uniq = {}
    for sz, gens in out:
        key = tuple(sorted(G.span(gens)))
        uniq[key] = (sz, gens)
    return sorted(uniq.values(), key=lambda t: -t[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('records'); ap.add_argument('out')
    ap.add_argument('--keys', default=None); ap.add_argument('--secs', type=float, default=600)
    ap.add_argument('--min-logE', type=float, default=-1e9); ap.add_argument('--max-nb', type=int, default=60)
    a = ap.parse_args()
    recs = []
    for f in a.records.split(','):
        for L in open(f):
            d = json.loads(L)
            if a.keys:
                if d.get('key') in a.keys.split(','):
                    recs.append(d)
            elif d.get('verdict', '').startswith('passes') and (d.get('log2E_Q') or -1e9) >= a.min_logE:
                recs.append(d)
    recs.sort(key=lambda d: -(d.get('log2E_Q') or -1e9))
    done = set()
    if os.path.exists(a.out):
        for L in open(a.out):
            done.add(json.loads(L)['key'])
    for rec in recs:
        if rec['key'] in done:
            continue
        if 'Kgens' not in rec:                            # k4 census record (bt_census_run): rebuild K gens
            PG = G0(rec['Pmods'])
            rec['Kgens'] = [list(map(int, PG.co[g])) + [0] * len(rec['Omods']) for g in rec['K2'] if g] + \
                [[0] * len(rec['Pmods']) + [int(i == j) for i in range(len(rec['Omods']))] for j in range(len(rec['Omods']))]
        res = dict(key=rec['key'], n=rec['n'], fam=rec.get('fam', 'k4'), order=rec['order'], r=rec['r'], e=rec['e'],
                   log2E_Q=rec.get('log2E_Q'), levels=[])
        verdict = 'undecided (no level proved infeasible)'
        for sz, gens in natural_subgroups(rec):
            sysd = level_system(rec, gens)
            if sysd['nb'] > a.max_nb:
                res['levels'].append(dict(N=sz, nb=sysd['nb'], status='skipped (too fine)'))
                continue
            r = solve(sysd, secs=a.secs)
            res['levels'].append(dict(N=sz, nb=sysd['nb'], coords=int(sysd['borb'].max()) + 1, status=r['status'], secs=r['secs']))
            if r['status'] == 'INFEASIBLE':
                verdict = f'no core [M] (CP-SAT: level G/N, |N| = {sz}, infeasible)'
                break
        res['verdict'] = verdict
        with open(a.out, 'a') as f:
            f.write(json.dumps(res) + '\n')
        print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
