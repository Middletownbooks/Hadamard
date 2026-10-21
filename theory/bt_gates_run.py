r"""Run the affordable exact gates of census survivors, cheapest first [M].

Input: the next-gate costs (theory/bt_cost.py output: per structure and invariant N <= K, log2 of the 2 + 2 join pairs
summed over the Q-level solutions, and log2 of the largest refined list).  Every (structure, level) with
log2_pairs <= --cap that no ladder run has decided is run: each Q-level solution (all row patterns) is refined to G/N
and joined exactly (theory/bt_gates.Gates.refine_level).  A level with no solution over all Q-level solutions closes
the structure: no core [M].  Otherwise the number of level solutions (and their total multiplicity) is recorded.
One JSON line per (structure, level); restart-safe.
usage: python3 theory/bt_gates_run.py COSTS.jsonl OUT.jsonl CENSUS.jsonl[,...] [--cap 34] [--part i/W] [--keys k1,k2]
"""
import sys, os, json, math, time, argparse, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_gates as BG, bt_census_run as R, bt_census2 as C2, bt_ladder as LD
import bt_cost as CO


def run_level(st, NA, cap_pairs, maxlen):
    g = BG.Gates(st, maxlen=maxlen, maxpairs=cap_pairs)
    qs = g.q_level(patterns=LD.patterns_for(st))
    qA = st.quotient(NA)
    chA = C.char_matrix(qA); cvA = [C.coord_values(qA, i) for i in range(4)]
    tot = 0; mult = 0.0; per = []
    for rows, sol, L in qs:
        w = BG.pattern_weight(tuple(abs(x) for x in rows)) if not st.split else 1
        ps = 0
        for s in sol:
            LA, solA, noA, m = g.refine_level(qA, (rows, s, L), chA=chA, cvA=cvA)
            ps += noA; tot += w * noA; mult += w * m
        per.append(dict(rows=list(map(int, rows)), solutions=int(ps)))
    return dict(solutions_weighted=int(tot), log2_mult=round(math.log2(mult), 2) if mult > 0 else None, per_pattern=per)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('costs'); ap.add_argument('out'); ap.add_argument('census')
    ap.add_argument('--cap', type=float, default=34); ap.add_argument('--part', default='0/1')
    ap.add_argument('--keys', default=None, help='restrict to these structure keys (comma-separated, or @FILE)')
    a = ap.parse_args()
    if a.keys and a.keys.startswith('@'):                    # @FILE: keys separated by commas or newlines
        a.keys = open(a.keys[1:]).read().replace('\n', ',').strip(',')
    only = set(a.keys.split(',')) if a.keys else None
    pi, pw = map(int, a.part.split('/'))
    recs = {}
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x); recs[d['key']] = d
    decided = set()
    for f in glob.glob('logs/bt/ladder_run*.jsonl'):
        for x in open(f):
            d = json.loads(x)
            if d.get('verdict', '').startswith('no core'):
                decided.add(d['key'])
    done = set(); old_done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            d = json.loads(x)
            if d.get('verdict', '').startswith('capped'):
                continue                                        # retried (the caps may have changed)
            done.add((d['key'], d.get('Nkey') or d['N']))
            if not d.get('Nkey'):
                old_done.add((d['key'], d['N']))           # recorded before subgroup keys: valid if N is unique
            if d.get('verdict', '').startswith('no core'):
                decided.add(d['key'])
    jobs = []
    for x in open(a.costs):
        c = json.loads(x)
        for l in c.get('levels', []):
            if l.get('log2_pairs') is not None and l['log2_pairs'] <= a.cap:
                jobs.append((l['log2_pairs'], c['key'], l['N'], l['log2_maxlist'], l.get('Nkey')))
    jobs.sort()
    for k, (lp, key, Nsz, lml, nk) in enumerate(jobs):
        if k % pw != pi or key in decided or (key, nk or Nsz) in done or key not in recs or (only and key not in only):
            continue
        rec = recs[key]
        st = C2.build(rec) if 'fam' in rec else R.build(rec)
        cands = [N for N in st.invariant_subgroups() if len(N) == Nsz and np.isin(N, st.Kel).all()]
        if (key, Nsz) in old_done and len(cands) == 1:
            continue
        NA = next((N for N in cands if CO.nkey(N) == nk), cands[0]) if nk else cands[0]
        t = time.time()
        try:
            # per-cell partial lists can exceed the final list bound, so the list cap has a floor of 4e6
            res = run_level(st, NA, cap_pairs=2.0 ** (a.cap + 1), maxlen=int(min(max(2 ** (lml + 1), 4e6), 6e7)))
            res['verdict'] = 'no core [M] (level empty)' if res['solutions_weighted'] == 0 else 'passes'
        except MemoryError as e:
            res = dict(verdict=f'capped ({e})')
        res.update(key=key, n=rec['n'], fam=rec.get('fam', 'k4'), order=rec['order'], r=rec['r'], e=rec['e'],
                   log2E_Q=rec.get('log2E_Q'), N=int(Nsz), Nkey=CO.nkey(NA), log2_pairs=lp, secs=round(time.time() - t, 1))
        if res['verdict'].startswith('no core'):
            decided.add(key)
        with open(a.out, 'a') as fo:
            fo.write(json.dumps(res) + '\n')
        print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'log2E_Q', 'N', 'log2_pairs', 'verdict',
                                                  'solutions_weighted', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
