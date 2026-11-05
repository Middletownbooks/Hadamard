"""Bounded, list-free feasibility experiments for the Z43 compression at v=863.

Uses bt_engine only for structure arithmetic; certificates are checked independently
by bt_cpsat_level.level_system (which shares no bt_engine code).
FEASIBLE certifies a compressed level, NOT a core or Hadamard matrix.
UNKNOWN is inconclusive. INFEASIBLE with fixed z excludes only that z fibre.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time
import numpy as np
from ortools.sat.python import cp_model
import bt_engine as E
from bt_cpsat_level import level_system


def structure(kind):
    mods = (4, 5, 43) if kind == 'cyclic' else (2, 2, 5, 43)
    M = np.eye(len(mods), dtype=int)
    M[-1, -1] = 4
    if kind == 's14':
        M[:2, :2] = [[0, 1], [1, 0]]
        M[-2, -2] = 4
    kg = np.eye(len(mods), dtype=int)[-2:].tolist()
    st = E.BStruct(mods, kg, 'k', [M.tolist()])
    rec = dict(Pmods=list(mods[:-2]), Omods=[5, 43], Kgens=kg,
               gens=[M.tolist()], gamma=[12, -4, -4, -4])
    return st, rec


def domains(qt):
    result = []
    for u in range(len(qt.bsize[0])):
        d = {0}
        for w in qt.w[0][qt.omap[0] == u]:
            d = {v + int(w) for v in d} | {v - int(w) for v in d}
        result.append(sorted(d))
    return result


def verify(y, rec, ngens):
    """All lags, exact integers, independent quotient construction and domains."""
    sd = level_system(rec, ngens)
    y = np.asarray(y, dtype=np.int64)
    assert y.shape == (4, sd['nb'])
    for row in y:
        for u, val in enumerate(row):
            assert int(val) in sd['dom'][u]
        for u in range(sd['nb']):
            assert np.all(row[sd['borb'] == sd['borb'][u]] == row[u])
    paf = np.array([sum(int(row @ row[np.array(sd['addb'])[:, w]]) for row in y)
                    for w in range(sd['nb'])])
    assert np.array_equal(paf, sd['FN']), (paf - sd['FN']).tolist()
    return True


def solve(st, qt, secs, seed, workers, z=None, max_abs=None, hint=None, fixed_blocks=()):
    d = domains(qt)
    if max_abs is not None:
        d = [[v for v in vals if abs(v) <= max_abs] for vals in d]
        assert all(d)
    model = cp_model.CpModel()
    Y = [[model.new_int_var_from_domain(cp_model.Domain.from_values(v), f'y{i}_{u}')
          for u, v in enumerate(d)] for i in range(4)]
    if hint is not None:
        for i in range(4):
            for u, rep in enumerate(qt.breps[i]):
                model.add_hint(Y[i][u], int(hint[i][rep]))
                if i in fixed_blocks:
                    model.add(Y[i][u] == int(hint[i][rep]))
    # Share products across every lag and use reduced orbit coordinates.
    products = {}
    lim = max(abs(v) for vals in d for v in vals)
    for lag in qt.cls_rep:
        terms = Counter()
        for u in range(qt.nb):
            a, b = sorted((int(qt.borb[0][u]), int(qt.borb[0][qt.add[u, lag]])))
            terms[a, b] += 1
        total = []
        for i in range(4):
            for (a, b), multiplicity in terms.items():
                key = i, a, b
                if key not in products:
                    vals = sorted({x * v for x in d[a] for v in d[b]}) if a != b else sorted({x*x for x in d[a]})
                    p = model.new_int_var_from_domain(cp_model.Domain.from_values(vals), f'p{i}_{a}_{b}')
                    if a == b:
                        model.add_allowed_assignments([Y[i][a], p], [(x, x*x) for x in d[a]])
                    else:
                        model.add_multiplication_equality(p, [Y[i][a], Y[i][b]])
                    products[key] = p
                total.append(multiplicity * products[key])
        model.add(sum(total) == int(qt.FN[lag]))
    if z is not None:
        cos = st.coset[qt.rep]
        for i in range(4):
            for c in range(st.q):
                model.add(sum(Y[i][int(qt.borb[i][u])] for u in range(qt.nb) if cos[u] == c) == int(z[i][c]))
    else:
        rows = [sum(int(qt.bsize[i][u])*v for u, v in enumerate(Y[i])) for i in range(4)]
        model.add(rows[0] >= 0)
        for i in range(3):
            model.add(rows[i] <= rows[i+1])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = secs
    solver.parameters.num_workers = workers
    solver.parameters.random_seed = seed
    status = solver.solve(model)
    out = dict(status=solver.status_name(status), seconds=solver.wall_time,
               branches=solver.num_branches, conflicts=solver.num_conflicts,
               compressed_positions=qt.nb, variables_per_block=len(d),
               domain=d, seed=seed, workers=workers, z=z, max_abs=max_abs,
               fixed_blocks=list(fixed_blocks))
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        out['y'] = [[solver.value(Y[i][int(qt.borb[i][u])]) for u in range(qt.nb)] for i in range(4)]
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kind', choices=['cyclic', 's7', 's14'], default='s14')
    p.add_argument('--secs', type=float, default=60)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--z-index', type=int)
    p.add_argument('--max-abs', type=int)
    p.add_argument('--hint')
    p.add_argument('--fixed-blocks', default='')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    st, rec = structure(a.kind)
    ngens = [[0]*(st.G.k-1)+[1]]
    qt = st.quotient(st.G.subgroup(ngens))
    z = None
    if a.z_index is not None:
        entries = [json.loads(l) for l in open('logs/bt/border_census.jsonl')]
        entry = next(e for e in entries if e['key'] == 'fbf4f510b2566552')
        z = entry['detail'][a.z_index]['z']
    hint = json.load(open(a.hint))['y'] if a.hint else None
    fixed = tuple(int(i) for i in a.fixed_blocks.split(',') if i)
    assert not fixed or hint is not None
    out = solve(st, qt, a.secs, a.seed, a.workers, z, a.max_abs, hint, fixed)
    out['hint_file'] = a.hint
    out.update(kind=a.kind, spec=rec, Ngens=ngens, stage='Z43 compression only', order=3452)
    if 'y' in out:
        out['independent_exact_check'] = verify(out['y'], rec, ngens)
        if z is not None:
            assert np.array_equal(np.array(out['y']).reshape(4,4,5).sum(2), z)
    path = Path(a.out); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('domain','y','spec')}), flush=True)


if __name__ == '__main__':
    main()
