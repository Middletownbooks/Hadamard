r"""Measured cost of the next exact gate for census survivors [M sizes; the cost is the work of the exact join].

For a structure and each invariant subgroup N <= K (N != K, coarse first), every Q-level solution (all row patterns)
is refined to G/N block by block; the refined list of block i has (at most, before the per-character filter) the
product over the Q-cells of the per-cell list sizes (theory/bt_chain.refine_list, counted without listing).  The exact
2 + 2 join at G/N costs L_a L_d + L_b L_c pair operations (sizes sorted), summed over the Q-level solutions.  The
cheapest level is the cost of continuing the exact chain one step; levels whose lists fit were already run by
theory/bt_ladder.py.
usage: python3 theory/bt_cost.py OUT.jsonl CENSUS.jsonl[,...] [--top 25] [--min-logE -10] [--maxq 64] [--keys k1,k2]
"""
import sys, os, json, math, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_gates as BG, bt_census_run as R, bt_census2 as C2, bt_ladder as LD
import hashlib


def nkey(N):
    """identifier of an invariant subgroup (several can share an order, e.g. the two axes of Z_11^2 at 971)"""
    return hashlib.sha1(np.sort(np.asarray(N, np.int64)).tobytes()).hexdigest()[:12]


def cell_sizes(qc, yc, qf, i, cvals):
    """per coarse cell: number of fine combinations meeting the cell value (refine_list, counted)"""
    cell, wt = C.nested_map(qc, qf, i)
    out = []
    for c in range(len(qc.bsize[i])):
        idx = np.flatnonzero(cell == c)
        tgt = int(yc[c])
        cnt = {0: 1}
        for u in idx:
            nxt = {}
            for s, m in cnt.items():
                for v in cvals[u]:
                    t = s + wt[u] * v
                    nxt[t] = nxt.get(t, 0) + m
            cnt = nxt
        out.append(cnt.get(tgt, 0))
    return out


def cost(st, maxq=64):
    g = BG.Gates(st)
    qs = g.q_level(patterns=LD.patterns_for(st))
    sols = [(rows, s, L) for rows, sol, L in qs for s in sol]
    res = dict(q_solutions=len(sols), levels=[])
    if not sols:
        return res
    rng = np.random.default_rng(1)
    pick = sols if len(sols) <= maxq else [sols[k] for k in rng.choice(len(sols), maxq, replace=False)]
    scale = len(sols) / len(pick)
    K = st.Kel
    for N in sorted([N for N in st.invariant_subgroups() if 1 < len(N) < len(K) and np.isin(N, K).all()],
                    key=lambda N: -len(N)):
        qA = st.quotient(N)
        cvA = [C.coord_values(qA, i) for i in range(4)]
        tot = 0.0; mx = 0.0
        for rows, s, L in pick:
            sz = []
            for i in range(4):
                cs = cell_sizes(g.qQ, L[i][0][s[i]], qA, i, cvA[i])
                sz.append(float(np.prod([float(x) for x in cs])))
            sz.sort()
            tot += sz[0] * sz[3] + sz[1] * sz[2]
            mx = max(mx, sz[3])
        tot *= scale
        res['levels'].append(dict(N=int(len(N)), Nkey=nkey(N), coords=[len(b) for b in qA.bsize],
                                  log2_pairs=round(math.log2(tot), 1) if tot > 0 else None,
                                  log2_maxlist=round(math.log2(mx), 1) if mx > 0 else None))
    fin = [l['log2_pairs'] for l in res['levels'] if l['log2_pairs'] is not None]
    res['log2_cost_next'] = min(fin) if fin else None
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('census')
    ap.add_argument('--top', type=int, default=25); ap.add_argument('--min-logE', type=float, default=-10)
    ap.add_argument('--maxq', type=int, default=64)
    ap.add_argument('--keys', default=None, help='only these structure keys (any estimate)')
    a = ap.parse_args()
    only = set(a.keys.split(',')) if a.keys else None
    done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            done.add(json.loads(x)['key'])
    alpha = {}
    if os.path.exists('logs/bt/census/thm_alpha.jsonl'):
        for x in open('logs/bt/census/thm_alpha.jsonl'):
            d = json.loads(x); alpha[d['key']] = d['dead']
    import glob as _g                                        # structures an exact gate already closed
    for f in _g.glob('logs/bt/gates_run*.jsonl') + _g.glob('logs/bt/ladder_run*.jsonl'):
        for x in open(f):
            d = json.loads(x)
            if d.get('verdict', '').startswith('no core'):
                alpha[d['key']] = True
    by = {}
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x)
            if only is not None:
                if d['key'] in only and d.get('verdict', '').startswith('passes'):
                    by.setdefault((d['n'], d.get('fam', 'k4')), []).append(d)
                continue
            if d.get('verdict', '').startswith('passes') and (d.get('log2E_Q') or -1e9) >= a.min_logE \
                    and not alpha.get(d['key']):
                by.setdefault((d['n'], d.get('fam', 'k4')), []).append(d)
    for k in sorted(by):
        for d in sorted(by[k], key=lambda d: -(d.get('log2E_Q') or -1e9))[:a.top]:
            if d['key'] in done:
                continue
            st = C2.build(d) if 'fam' in d else R.build(d)
            t = time.time()
            try:
                res = cost(st, maxq=a.maxq)
            except Exception as ex:
                res = dict(error=f'{type(ex).__name__}: {ex}')
            res.update(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), order=d['order'], r=d['r'], e=d['e'],
                       log2E_Q=d.get('log2E_Q'), secs=round(time.time() - t, 1))
            with open(a.out, 'a') as fo:
                fo.write(json.dumps(res) + '\n')
            print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'e', 'log2E_Q', 'q_solutions',
                                                      'log2_cost_next', 'levels', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
