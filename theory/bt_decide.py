r"""Exact decision of small bordered-GS structures by the full compression chain [M].

For a census record: choose levels automatically, A = the finest invariant N_A < K whose refined lists fit, B = an
invariant N_B with |N_A cap N_B| minimal and the most new coordinates, C = the invariant N_C with the most coordinates
among the rest (final-join split); run theory/bt_search.Run over ALL row patterns (sorted, nonnegative: the blocks
share S, so block order and signs are free) with no randomisation cap.  A run that ends with event 'exhausted' and no
core decides the structure: no core [M].
usage: python3 theory/bt_decide.py CENSUS.jsonl[,...] OUTDIR [--keys k1,k2] [--rmax 45] [--hours h] [--part i/W]
"""
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_census_run as R, bt_search as BS


def choose_levels(st, cap=3 * 10 ** 6):
    subs = st.invariant_subgroups()
    K = st.Kel
    inK = [N for N in subs if 1 < len(N) < len(K) and np.isin(N, K).all()]
    inK.sort(key=len)
    if not inK:
        return None
    qQ = st.quotient(K); chQ = C.char_matrix(qQ)
    # pick A: finest whose refined lists fit for a probe Q-level solution (first pattern with solutions)
    A = inK[-1]
    for N in inK:
        A = N
        break
    others = [N for N in subs if 1 < len(N) < st.G.n and not np.array_equal(N, A)]

    def inter(X, Y):
        return len(np.intersect1d(X, Y))
    Bc = sorted(others, key=lambda N: (inter(N, A), len(N)))
    B = Bc[0] if Bc else None
    rest = [N for N in others if B is None or not np.array_equal(N, B)]
    Cc = sorted(rest, key=lambda N: (inter(N, np.intersect1d(A, B) if B is not None else A), len(N)))
    Cn = Cc[0] if Cc else None
    return dict(A=int(len(A)), B=int(len(B)) if B is not None else None, C=int(len(Cn)) if Cn is not None else None)


def spec_from_record(rec, levels):
    st = R.build(rec)
    Pmods, Omods = tuple(rec['Pmods']), tuple(rec['Omods'])
    PG = E.Grp(Pmods)
    Kgens = []
    for g in rec['K2']:
        if g:
            Kgens.append(list(map(int, PG.co[g])) + [0] * len(Omods))
    for j in range(len(Omods)):
        v = [0] * (len(Pmods) + len(Omods)); v[len(Pmods) + j] = 1; Kgens.append(v)
    return dict(name=f"{rec['n']} P{Pmods} O{Omods} S{rec['order']} r{rec['r']} key {rec['key']}",
                mods=list(Pmods + Omods), Kgens=Kgens, gens=rec['gens'], levels=levels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('census'); ap.add_argument('outdir')
    ap.add_argument('--keys', default=None); ap.add_argument('--rmax', type=int, default=45)
    ap.add_argument('--hours', type=float, default=0.5); ap.add_argument('--part', default='0/1')
    a = ap.parse_args()
    pi, pw = map(int, a.part.split('/'))
    recs = []
    for f in a.census.split(','):
        for L in open(f):
            d = json.loads(L)
            if a.keys:
                if d.get('key') in a.keys.split(','):
                    recs.append(d)
            elif d.get('verdict', '').startswith('passes') and d.get('r', 999) <= a.rmax:
                recs.append(d)
    recs.sort(key=lambda d: d['r'])
    os.makedirs(a.outdir, exist_ok=True)
    summ = os.path.join(a.outdir, f'decide_w{pi}.jsonl')
    done = set()
    if os.path.exists(summ):
        for L in open(summ):
            done.add(json.loads(L)['key'])
    for k, rec in enumerate(recs):
        if k % pw != pi or rec['key'] in done:
            continue
        st = R.build(rec)
        lv = choose_levels(st)
        out = dict(key=rec['key'], n=rec['n'], Pmods=rec['Pmods'], Omods=rec['Omods'], order=rec['order'],
                   r=rec['r'], e=rec['e'], levels=lv)
        t0 = time.time()
        if lv is None:
            out['verdict'] = 'no levels'
        else:
            spec = spec_from_record(rec, lv)
            d = os.path.join(a.outdir, rec['key'])
            try:
                run = BS.Run(spec, d, seed=1)
                run.run(hours=a.hours)
                ev = [json.loads(L) for L in open(run.logf)]
                last = [e for e in ev if e['event'] in ('exhausted', 'stop')][-1]
                out['levelA_solutions'] = sum(e.get('solutions', 0) for e in ev if e['event'] == 'levelA')
                out['levelB_combos'] = sum(e.get('combos', 0) for e in ev if e['event'] == 'levelB')
                out['finished_combos'] = sum(1 for e in ev if e['event'] == 'combo_done')
                out['cores'] = run.ncores
                out['verdict'] = ('no core [M] (chain exhausted)' if last['event'] == 'exhausted' and run.ncores == 0
                                  else ('CORES FOUND' if run.ncores else 'undecided (time)'))
            except MemoryError as e:
                out['verdict'] = f'undecided (cap: {e})'
        out['secs'] = round(time.time() - t0, 1)
        with open(summ, 'a') as f:
            f.write(json.dumps(out) + '\n')
        print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
