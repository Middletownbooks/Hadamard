r"""First-compression-level gate and stratified count for a census structure [M gate, E count].

For a structure (census record or spec) and each invariant subgroup N_A <= K whose refined lists fit, every Q-level
solution (all row patterns) is refined to G/N_A and joined exactly; the number of level-A solutions and their
multiplicity are exact [M].  The stratified count sums, over (a sample of) the level-A solutions, the multiplicity times
the lattice Gaussian of the remaining class equations from blocks sampled uniformly inside their level-A fibres,
read at the true target [E]; ordered signed row vectors are weighted by the pattern weight.
usage: python3 theory/bt_level1.py CENSUS.jsonl KEY [--NA size] [--N 2000] [--maxA 40] [--out FILE]
"""
import sys, os, json, math, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C, bt_gates as BG, bt_census_run as R


def level_candidates(st, maxlen=3 * 10 ** 6):
    """invariant subgroups N <= K (N != K, N != 0 allowed) ordered by decreasing |G/N|"""
    out = []
    for N in st.invariant_subgroups():
        if len(N) < len(st.Kel) and np.isin(N, st.Kel).all() and len(N) > 1:
            out.append(N)
    return sorted(out, key=len)


def run_level(st, NA, N=2000, maxA=40, seed=1, maxlen=3 * 10 ** 6, patterns=None, log=None):
    g = BG.Gates(st, maxlen=maxlen)
    qA = st.quotient(NA)
    chA = C.char_matrix(qA)
    cvA = [C.coord_values(qA, i) for i in range(4)]
    qs = g.q_level(patterns=patterns)
    rng = np.random.default_rng(seed)
    res = dict(NA=int(len(NA)), coordsA=[len(b) for b in qA.bsize], q_solutions=int(sum(len(q[1]) for q in qs)))
    tot_sol = 0; tot_mult = 0.0; est = 0.0; samples = []
    per_pattern = []
    for rows, sol, L in qs:
        w = BG.pattern_weight(tuple(abs(x) for x in rows)) if not st.split else 1
        pat_sol = 0; pat_est = 0.0
        for s in sol:
            try:
                LA, solA, noA, mult = g.refine_level(qA, (rows, s, L), chA=chA, cvA=cvA)
            except MemoryError as e:
                res['cap'] = str(e)
                return res
            pat_sol += noA; tot_sol += w * noA; tot_mult += w * mult
            if noA == 0:
                continue
            k = min(noA, max(1, maxA // max(1, len(qs))))
            idx = rng.choice(noA, k, replace=False)
            vals = []
            for a in solA[idx]:
                ys = [LA[i][0][a[i]] for i in range(4)]
                lm = sum(math.log2(LA[i][1][a[i]]) for i in range(4))
                r = BG.stratified_bits(st, qA, ys, N=N, seed=int(rng.integers(1 << 30)))
                v = lm - r['bits'] if np.isfinite(r['bits']) else -math.inf
                vals.append(v)
                samples.append(dict(rows=list(map(int, rows)), log2E=round(v, 2) if np.isfinite(v) else None,
                                    mahal=r.get('mahal'), rank=r.get('rank')))
            vals = np.array(vals)
            m = np.mean(np.where(np.isfinite(vals), 2.0 ** vals, 0.0))
            pat_est += noA * m
        est += w * pat_est
        per_pattern.append(dict(rows=list(map(int, rows)), levelA_solutions=int(pat_sol),
                                log2E=round(math.log2(pat_est), 2) if pat_est > 0 else None))
    res.update(levelA_solutions_weighted=int(tot_sol), log2_mult=round(math.log2(tot_mult), 2) if tot_mult else None,
               log2E_A=round(math.log2(est), 2) if est > 0 else None, per_pattern=per_pattern,
               samples=sorted(samples, key=lambda d: -(d['log2E'] if d['log2E'] is not None else -1e9))[:8])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('census'); ap.add_argument('key')
    ap.add_argument('--NA', type=int, default=None)
    ap.add_argument('--N', type=int, default=2000); ap.add_argument('--maxA', type=int, default=40)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    rec = None
    for f in a.census.split(','):
        for L in open(f):
            d = json.loads(L)
            if d.get('key') == a.key:
                rec = d
    assert rec is not None
    st = R.build(rec)
    cands = level_candidates(st)
    if a.NA is not None:
        cands = [N for N in cands if len(N) == a.NA]
    for NA in cands[::-1]:
        t = time.time()
        res = run_level(st, NA, N=a.N, maxA=a.maxA)
        res.update(key=a.key, n=rec['n'], r=st.r[0], e=st.e, order=rec['order'], secs=round(time.time() - t, 1))
        print(json.dumps(res), flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(res) + '\n')


if __name__ == '__main__':
    main()
