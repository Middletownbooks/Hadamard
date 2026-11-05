#!/usr/bin/env python3
r"""[M] Review of Claim 1 (k4 compression normal form, by the OpenAI model, report/openai_4p3_answer.md): cross-check
against the repository's own Q-level solution lists.

The Q-level lists and joins are the census code of this repository (theory/bt_gates.Gates.q_level on structures built
by theory/bt_census_run.build; theory/bt_engine, theory/bt_chain; Claude Code).  Those lists come from an exact join
that knows nothing about the normal form, so every listed Q-level solution must be in normal form if Claim 1 holds.
The normal-form test and the Q-level correlation test below are new code (Claude Code review); nothing is imported
from theory/openai_4p3_*.py.

Per structure (census key, or the 1916 / 3964 controls): every Q-level solution (rows, quadruple) is expanded to full
Q-vectors on the cosets of K, checked by direct quotient correlations on the structure's own Q addition table
(sum_i PAF = 4m + 12 at 0, 4m off 0), and tested for the normal form: exactly one constant row and three rows with three
equal entries and a fourth differing by 2, with m = a^2 + sum t_j (t_j + 1).

usage: python3 theory/4p3_review_qxcheck.py --out logs/4p3_review/qxcheck.json [--seconds 900]
"""
import argparse, glob, json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_gates as BG, bt_census_run as CR, bt_chain as C


def nf_of(rows, m):
    """normal-form test of 4 full Q-vectors: returns (ok, description)"""
    const = [r for r in rows if len(set(r)) == 1]
    sing = []
    for r in rows:
        if len(set(r)) == 1:
            continue
        vals = sorted(r)
        if vals[0] == vals[2] and vals[3] == vals[2] + 2:
            sing.append(vals[0])                       # t, exceptional t + 2
        elif vals[1] == vals[3] and vals[0] == vals[1] - 2:
            sing.append(-vals[3])                      # negated row: three -u, one -u + 2
        else:
            return False, 'row not constant or singleton: %s' % (r,)
    if len(const) != 1 or len(sing) != 3:
        return False, 'shape count const=%d sing=%d' % (len(const), len(sing))
    a = abs(const[0][0])
    if a * a + sum(t * (t + 1) for t in sing) != m:
        return False, 'scalar identity fails'
    return True, dict(a=a, t=sorted(sing))


def qcorr_ok(rows, add, m):
    q = len(add)
    tot = [sum(r[u] * r[add[u][d]] for r in rows for u in range(q)) for d in range(q)]
    want = [4 * m + 12] + [4 * m] * (q - 1)
    return tot == want


def check_struct(name, st, g=None, patterns=None):
    t0 = time.time()
    g = g or BG.Gates(st)
    qs = g.q_level(patterns=patterns)
    qQ = g.qQ
    add = qQ.add.tolist()
    m = len(st.Kel)
    nsol = nf_ok = corr_ok = 0
    forms = {}
    bad = []
    for rows, sol, L in qs:
        for s in sol:
            z = [qQ.full_from_reduced(i, L[i][0][s[i]]).astype(int).tolist() for i in range(4)]
            nsol += 1
            c = qcorr_ok(z, add, m)
            corr_ok += c
            ok, desc = nf_of(z, m)
            nf_ok += ok
            if ok:
                k = (desc['a'], tuple(desc['t']))
                forms[k] = forms.get(k, 0) + 1
            elif len(bad) < 5:
                bad.append(dict(z=z, why=desc))
    return dict(structure=name, G=list(st.G.mods), K_order=m, Q_add=add, q_solutions=nsol, q_corr_ok=corr_ok,
                normal_form_ok=nf_ok, forms={'a=%d t=%s' % (k[0], list(k[1])): v for k, v in sorted(forms.items())},
                violations=bad, secs=round(time.time() - t0, 1))


def census_records(keys):
    out = {}
    for f in sorted(glob.glob('logs/bt/census/k4*.jsonl')):
        for L in open(f):
            try:
                r = json.loads(L)
            except Exception:
                continue
            if r.get('key') in keys and r['key'] not in out:
                out[r['key']] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seconds', type=float, default=900)
    ap.add_argument('--keys', default='2e93be9718324148')
    ap.add_argument('--auto', type=int, default=8, help='also take this many more passing k4 structures at 839/863 '
                    '(odd |K|), alternating Q type, smallest r first')
    a = ap.parse_args()
    t0 = time.time()
    keys = [k for k in a.keys.split(',') if k]
    if a.auto:
        cand = []
        for f in sorted(glob.glob('logs/bt/census/k4*.jsonl')):
            for L in open(f):
                try:
                    r = json.loads(L)
                except Exception:
                    continue
                if r.get('verdict') == 'passes Q-level' and r.get('n') in (839, 863) and r['key'] not in keys:
                    cand.append(r)
        cand.sort(key=lambda r: (r['r'], r['key']))
        byq = {}
        for r in cand:
            byq.setdefault((r['n'], r['Qtype']), []).append(r)
        order = []
        groups = sorted(byq)
        while len(order) < a.auto and any(byq[gk] for gk in groups):
            for gk in groups:
                if byq[gk] and len(order) < a.auto:
                    order.append(byq[gk].pop(0)['key'])
        keys += order
    recs = census_records(set(keys))
    res = dict(script='theory/4p3_review_qxcheck.py', command=' '.join(sys.argv), seed=None, records=[])
    # controls: 3964 (Kumar's structure) and the 1916 split structure (theory/bt_controls.py)
    st = E.control_3964()
    res['records'].append(dict(check_struct('3964 control (Kumar structure)', st), rows_pattern='all'))
    print(json.dumps({k: v for k, v in res['records'][-1].items() if k not in ('Q_add', 'violations')}), flush=True)
    A2 = np.array([[1, 1], [1, 0]])
    g12a = E.block_diag(A2, [[16]]); g12b = E.block_diag(np.eye(2, dtype=int), [[13]])
    split = E.BStruct((2, 2, 119), [(0, 0, 1)], 'k', [[g12a, g12b], [g12b], [g12b], [g12a, g12b]], name='1916-split')
    res['records'].append(check_struct('1916 split structure (theory/bt_controls.py)', split))
    print(json.dumps({k: v for k, v in res['records'][-1].items() if k not in ('Q_add', 'violations')}), flush=True)
    for k in keys:
        if time.time() - t0 > a.seconds:
            res['stopped_at_time_limit'] = True
            break
        r = recs.get(k)
        if r is None:
            res['records'].append(dict(structure=k, error='key not found'))
            continue
        st = CR.build(r)
        try:
            out = check_struct('census %s n=%d %s S%d r=%d' % (k, r['n'], r['Qtype'], r['order'], r['r']), st)
        except MemoryError as ex:
            out = dict(structure=k, error='Q-level truncated: %s' % ex)
        res['records'].append(out)
        print(json.dumps({x: v for x, v in out.items() if x not in ('Q_add', 'violations')}), flush=True)
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
    good = [r for r in res['records'] if 'q_solutions' in r]
    res['total_solutions'] = sum(r['q_solutions'] for r in good)
    res['all_in_normal_form'] = all(r['normal_form_ok'] == r['q_solutions'] == r['q_corr_ok'] for r in good)
    res['secs'] = round(time.time() - t0, 1)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(dict(structures=len(good), total_solutions=res['total_solutions'],
                          all_in_normal_form=res['all_in_normal_form'], secs=res['secs'])))


if __name__ == '__main__':
    main()
