#!/usr/bin/env python3
r"""[M] Review of Claim 1' (the small-repair obstruction) of report/openai_4p3_answer.md: independent recomputation of
the minimum-change table (Proposal 1 section 4) and of the per-index cyclotomic gate (section 3).

The claim, the cost function and the table are by the OpenAI model.  This script is new (Claude Code review) and imports
no repository code.  It uses the normal form of Claim 1 (proved, and checked by theory/4p3_review_nf.py): after
permuting and negating rows, the Q-compression is one constant row (a, a, a, a) and three rows (t+2, t, t, t), a and
t_j odd, m = a^2 + sum_j t_j (t_j + 1).

Cost model (the answer's).  A cell of odd length p with sum +-1 ("balanced") needs at least d(z) = (|z| - 1) / 2 sign
changes to reach the odd sum z.  Total cost of a form = 4 d(a) + sum_j (3 d(t_j) + d(t_j + 2)); cell cap = the largest
d over the 16 cells.  Both minimised over every normal form.

Cyclotomic gate.  Index e | p - 1, class size s = (p - 1)/e: the unmodified sums are b = eta + s (e - 2j), eta = +-1,
j = 0..e negative classes; with at most c sign changes, z = b + 2w with -min(c, (p + b)/2) <= w <= min(c, (p - b)/2).
The least c for which a normal form has all 16 entries attainable is computed per e.  The attainable-sum formula is
checked against explicit cells (primitive root, every class-sign vector, every zero sign, every flip set of size <= c)
for small p.

usage: python3 theory/4p3_review_repair.py --out logs/4p3_review/repair_table.json
"""
import argparse, itertools, json, math, sys, time

TARGETS = [149, 179, 257, 271, 331, 359, 367, 389, 457, 461, 467, 509, 521, 523, 571, 599]
SMALL = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
ANSWER_TABLE = {149: (29, 4), 179: (29, 4), 257: (37, 5), 271: (38, 5), 331: (41, 5), 359: (43, 5), 367: (39, 6),
                389: (41, 6), 457: (47, 6), 461: (47, 6), 467: (43, 6), 509: (44, 7), 521: (49, 6), 523: (51, 6),
                571: (53, 6), 599: (61, 7)}


def forms(m):
    """every normal form (a > 0, ordered (t1, t2, t3)), a and t odd"""
    L = math.isqrt(m) + 3
    ts = [t for t in range(-L - 1, L + 1) if t % 2 and t * (t + 1) <= m]
    out = []
    for a in range(1, math.isqrt(m) + 1, 2):
        for t1, t2, t3 in itertools.product(ts, repeat=3):
            if a * a + t1 * (t1 + 1) + t2 * (t2 + 1) + t3 * (t3 + 1) == m:
                out.append((a, t1, t2, t3))
    return out


def cells(f):
    a, *ts = f
    return [a] * 4 + [x for t in ts for x in (t + 2, t, t, t)]


def d(z):
    return (abs(z) - 1) // 2


def qcheck(f, m):
    """direct quotient correlations of the form on both Q (no normal-form formula)"""
    a, *ts = f
    rows = [[a] * 4] + [[t + 2, t, t, t] for t in ts]
    for q in ('v4', 'z4'):
        tot = [sum(r[u] * r[(u ^ dd) if q == 'v4' else (u + dd) % 4] for r in rows for u in range(4)) for dd in range(4)]
        if tot != [4 * m + 12, 4 * m, 4 * m, 4 * m]:
            return False
    return True


def attainable(p, e, c):
    s = (p - 1) // e
    S = set()
    for j in range(e + 1):
        for eta in (-1, 1):
            b = eta + s * (e - 2 * j)
            for w in range(-min(c, (p + b) // 2), min(c, (p - b) // 2) + 1):
                S.add(b + 2 * w)
    return S


def explicit_sums(p, e, c):
    """sums of explicit cells: class of k != 0 is log_g(k) mod e (g = least primitive root), zero sign eta at k = 0,
    every class-sign vector, then every flip set of size <= c"""
    g = next(g for g in range(2, p) if len({pow(g, k, p) for k in range(p - 1)}) == p - 1)
    lg = {pow(g, k, p): k for k in range(p - 1)}
    out = set()
    for signs in itertools.product((1, -1), repeat=e):
        for eta in (1, -1):
            cell = [eta] + [signs[lg[k] % e] for k in range(1, p)]
            for r in range(c + 1):
                for F in itertools.combinations(range(p), r):
                    out.add(sum(cell) - 2 * sum(cell[k] for k in F))
    return out


def table_row(p, indices=(2, 4, 6, 8, 12, 16), cmax=12):
    fs = forms(p)
    assert fs and all(qcheck(f, p) for f in fs[:50])
    tot = min(sum(d(z) for z in cells(f)) for f in fs)
    cap = min(max(d(z) for z in cells(f)) for f in fs)
    wt = min(fs, key=lambda f: sum(d(z) for z in cells(f)))
    wc = min(fs, key=lambda f: max(d(z) for z in cells(f)))
    per_e = {}
    for e in indices:
        if (p - 1) % e:
            continue
        for c in range(cmax + 1):
            S = attainable(p, e, c)
            w = next((f for f in fs if all(z in S for z in cells(f))), None)
            if w is not None:
                per_e[e] = dict(min_c=c, witness=list(w))
                break
        else:
            per_e[e] = dict(min_c=None, note=f'no form through c = {cmax}')
    return dict(p=p, N=16 * p + 12, forms=len(fs), min_total=tot, min_cap=cap, total_witness=list(wt),
                cap_witness=list(wc), all_forms_pass_direct_Q_check=all(qcheck(f, p) for f in fs), per_e=per_e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--host', default='logs/openai_4p3/compression_host.json',
                    help="the answer's log, read only to compare numbers")
    a = ap.parse_args()
    t0 = time.time()
    # 1. attainable-sum formula against explicit cells
    formula_checks = []
    for p in (5, 7, 13, 17, 19, 29, 31, 37):
        for e in (1, 2, 3, 4, 6):
            if (p - 1) % e:
                continue
            for c in (0, 1, 2):
                ok = attainable(p, e, c) == explicit_sums(p, e, c)
                formula_checks.append([p, e, c, ok])
    # 2. table
    rows = [table_row(p) for p in SMALL + TARGETS]
    try:
        host = {r['p']: r for r in json.load(open(a.host))['records']}
    except Exception:
        host = {}
    cmp = []
    for r in rows:
        h = host.get(r['p'])
        ans = ANSWER_TABLE.get(r['p'])
        cmp.append(dict(p=r['p'], ours=[r['min_total'], r['min_cap']],
                        answer_table=list(ans) if ans else None,
                        answer_log=[h['min_total'], h['min_per_cell']] if h else None,
                        per_e_equal_cap=all(v.get('min_c') == r['min_cap'] for v in r['per_e'].values())))
    agree_table = all(c['answer_table'] is None or c['answer_table'] == c['ours'] for c in cmp)
    agree_log = all(c['answer_log'] is None or c['answer_log'] == c['ours'] for c in cmp)
    res = dict(script='theory/4p3_review_repair.py', command=' '.join(sys.argv), seed=None,
               formula_vs_explicit_cells=formula_checks, formula_all_ok=all(x[3] for x in formula_checks),
               rows=rows, comparison=cmp, agree_with_answer_table=agree_table, agree_with_answer_log=agree_log,
               targets_per_e_equal_cap=all(c['per_e_equal_cap'] for c in cmp if c['p'] in TARGETS),
               secs=round(time.time() - t0, 1))
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    for c in cmp:
        if c['p'] in TARGETS:
            r = next(x for x in rows if x['p'] == c['p'])
            print(c['p'], c['ours'], c['answer_table'], {e: v['min_c'] for e, v in r['per_e'].items()})
    print(json.dumps(dict(formula_all_ok=res['formula_all_ok'], agree_with_answer_table=agree_table,
                          agree_with_answer_log=agree_log, targets_per_e_equal_cap=res['targets_per_e_equal_cap'],
                          secs=res['secs'])))


if __name__ == '__main__':
    main()
