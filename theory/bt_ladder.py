r"""Gate ladder: exact compression gates of a structure at every invariant N <= K, coarse to fine [M].

For every invariant subgroup N <= K (N != K) in order of increasing |G/N|, every Q-level solution (all row patterns)
is refined to G/N (product over the Q-cells, per-character bounds) and joined exactly.  The first level with no
solution decides the structure: no core [M].  Levels whose refined lists exceed the cap are recorded as such and
skipped.  Output: one JSON line per structure with the per-level counts.
usage: python3 theory/bt_ladder.py CENSUS.jsonl[,...] OUT.jsonl [--maxlen 4e6] [--part i/W] [--min-logE -1e9]
"""
import sys, os, json, math, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_gates as BG, bt_census_run as R, bt_census2 as C2


def patterns_for(st):
    triv_in_T = abs(sum(st.gamma) / st.q / 4 - 1) < 1e-9
    tot = 4 * st.t if triv_in_T else 4 * st.nH
    par = st.G.n % 2; lim = math.isqrt(tot); out = []
    for a in range(par, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = tot - a * a - b * b - c * c
                if d2 >= c * c and math.isqrt(d2) ** 2 == d2 and math.isqrt(d2) % 2 == par:
                    out.append((a, b, c, math.isqrt(d2)))
    return out


def ladder(st, maxlen=4 * 10 ** 6, tmax=600):
    t0 = time.time()
    g = BG.Gates(st, maxlen=maxlen, maxpairs=2e9)
    qs = g.q_level(patterns=patterns_for(st))
    res = dict(q_solutions=int(sum(len(q[1]) for q in qs)), levels=[])
    if not qs:
        res['verdict'] = 'no core [M] (Q-level empty)'
        return res
    K = st.Kel
    levels = [N for N in st.invariant_subgroups() if 1 < len(N) < len(K) and np.isin(N, K).all()]
    levels.sort(key=lambda N: -len(N))                     # coarse (large N) first
    for N in levels:
        if time.time() - t0 > tmax:
            res['levels'].append(dict(N=int(len(N)), note='time'))
            continue
        qA = st.quotient(N)
        chA = C.char_matrix(qA)
        cvA = [C.coord_values(qA, i) for i in range(4)]
        tot = 0; capped = False; checked = 0
        for rows, sol, L in qs:
            for s in sol:
                try:
                    LA, solA, noA, mult = g.refine_level(qA, (rows, s, L), chA=chA, cvA=cvA)
                except MemoryError:
                    capped = True; continue
                checked += 1
                tot += noA
                if tot:
                    break
            if tot:
                break
        # the level passes as soon as one Q-level solution refines; it is empty only if every one was checked
        rec = dict(N=int(len(N)), coords=len(qA.bsize[0]), solutions=int(tot) if (tot or not capped) else None,
                   capped=capped and not tot, checked=checked)
        res['levels'].append(rec)
        if not rec['capped'] and tot == 0:
            res['verdict'] = f'no core [M] (level G/N, |N| = {len(N)}, empty)'
            res['secs'] = round(time.time() - t0, 1)
            return res
    res['verdict'] = 'passes every listable level'
    res['secs'] = round(time.time() - t0, 1)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('census'); ap.add_argument('out')
    ap.add_argument('--maxlen', type=float, default=4e6); ap.add_argument('--part', default='0/1')
    ap.add_argument('--min-logE', type=float, default=-1e9)
    a = ap.parse_args()
    pi, pw = map(int, a.part.split('/'))
    done = set()
    if os.path.exists(a.out):
        for L in open(a.out):
            done.add(json.loads(L)['key'])
    recs = []
    for f in a.census.split(','):
        for L in open(f):
            d = json.loads(L)
            if d.get('verdict', '').startswith('passes') and (d.get('log2E_Q') or -1e9) >= a.min_logE:
                recs.append(d)
    recs.sort(key=lambda d: -(d.get('log2E_Q') or -1e9))
    for k, rec in enumerate(recs):
        if k % pw != pi or rec['key'] in done:
            continue
        st = C2.build(rec) if 'fam' in rec else R.build(rec)
        try:
            res = ladder(st, maxlen=int(a.maxlen))
        except Exception as e:
            res = dict(verdict=f'error: {type(e).__name__}: {e}')
        res.update(key=rec['key'], n=rec['n'], fam=rec.get('fam', 'k4'), T=rec.get('T'), Pmods=rec['Pmods'],
                   Omods=rec['Omods'], order=rec['order'], r=rec['r'], e=rec['e'], log2E_Q=rec.get('log2E_Q'))
        with open(a.out, 'a') as f:
            f.write(json.dumps(res) + '\n')
        print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'e', 'log2E_Q', 'verdict', 'levels', 'secs')}),
              flush=True)


if __name__ == '__main__':
    main()
