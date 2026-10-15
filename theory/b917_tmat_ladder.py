r"""Same-shape control ladder for the T-matrix route of H(3668): cyclotomic T-matrices of order t (K = the order-5
subgroup of F_t^*) times Williamson matrices of order 7 -> H(28 t) [M, P, E].

For each t = 1 (mod 10) prime:
  rowsum   the exact row-sum condition (three row sums = 0, one = +-1 (mod 5), sum of squares t) [P, enumerated];
  estimate scripts/tmat_estimate.py (log2 E of labelled solutions) [E];
  exact    with --all: the complete count over every support S (2 <= |S| <= m - 2) by exact 2+2 joins
           (theory/b917_tmat_census.count_S); by the T1 <-> T2 / sign symmetry each labelled solution with both
           sides nonempty is counted once per 16, so labelled = 16 * count [P];
  search   time to the first solution with scripts/tmat_kt.py (random supports);
  build    Cooper-Wallis with W(7) (theory/b917_cw), then the independent scripts/verify_h3668.py on the result.
Log: logs/b917/tmat_ladder.jsonl.
usage: python3 theory/b917_tmat_ladder.py --t 41,61,101 [--all 41,61] [--build]"""
import sys, os, json, time, math, itertools, subprocess
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
from tmat_kt import orbit_setup
from b917_tmat_census import count_S
import b917_cw as CW


def rowsum_solutions(t, k):
    """all |s| patterns: three s = 0 (mod k), one s = +-1 (mod k), sum s^2 = t"""
    out = set()
    lim = int(math.isqrt(t))
    for s4 in range(0, lim + 1):
        if s4 % k not in (1, k - 1): continue
        rest = t - s4 * s4
        if rest % (k * k): continue
        r = rest // (k * k)
        for a in range(0, int(math.isqrt(r)) + 1):
            for b in range(a, int(math.isqrt(r)) + 1):
                c2 = r - a * a - b * b
                if c2 < b * b: continue
                c = int(math.isqrt(c2))
                if c * c == c2: out.add((s4, k * a, k * b, k * c))
    return sorted(out)


def exact_total(t, k):
    K, orbits, reps = orbit_setup(t, k); m = len(orbits)
    tot = 0; t0 = time.time(); per_size = {}
    for size in range(2, m - 1):
        for S in itertools.combinations(range(m), size):
            n = len(count_S(t, orbits, reps, list(S), 0, 7, 10 ** 9))
            tot += n; per_size[size] = per_size.get(size, 0) + n
    return dict(orbits=m, census_total=tot, labelled=16 * tot, log2_labelled=round(math.log2(16 * tot), 2) if tot else None,
                per_size=per_size, secs=round(time.time() - t0, 1))


def first_solution(t, k, seed=1, nS=400):
    t0 = time.time()
    p = subprocess.run([sys.executable, 'scripts/tmat_kt.py', str(t), str(k), '--nS', str(nS), '--seed', str(seed), '--buckets', '200'],
                       capture_output=True, text=True, cwd=ROOT)
    recs = [json.loads(l) for l in p.stdout.splitlines() if l.startswith('{"S_index"')]
    hit = [r for r in recs if r.get('FOUND')]
    return (np.array(hit[0]['T'], np.int64) if hit else None), dict(supports_tried=len(recs), secs=round(time.time() - t0, 1))


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    ts = [int(x) for x in opt('--t', '41,61,101').split(',')]
    alls = {int(x) for x in opt('--all', '').split(',') if x}
    out = os.path.join(ROOT, 'logs/b917/tmat_ladder.jsonl')
    for t in ts:
        k = 5
        rec = dict(t=t, k=k, rowsum_patterns=rowsum_solutions(t, k))
        est = subprocess.run([sys.executable, 'tmat_estimate.py', str(t), str(k), '200000'], capture_output=True, text=True,
                             cwd=os.path.join(ROOT, 'scripts')).stdout.strip().splitlines()
        e = json.loads(est[-1]) if est else {}
        rec.update(estimate_log2E=e.get('log2E'), estimate_rowsum_feasible=e.get('rowsum_feasible'), orbits=e.get('orbits'))
        if t in alls:
            rec['exact'] = exact_total(t, k)
        if rec['rowsum_patterns']:
            T, info = first_solution(t, k)
            rec['search'] = info
            if T is not None and '--build' in a:
                CW.check_T(T)
                X = CW.cooper_wallis(T, CW.williamson(7))
                Xc = CW.to_cyclic(X)
                v = 7 * t
                seq = os.path.join(ROOT, f'logs/b917/seq_ladder_H{4 * v}_Z{v}.txt')
                with open(seq, 'w') as fh:
                    for lab, x in zip('ABCD', Xc): fh.write(f"{lab}=" + ''.join('+' if z > 0 else '-' for z in x) + "\n")
                p = subprocess.run([sys.executable, 'scripts/verify_h3668.py', seq, f'/tmp/ladder_H{4 * v}.npz', '--slice', str(t), str((t - 1) // 5), '--umult', '6'],
                                   capture_output=True, text=True, cwd=ROOT)
                last = [json.loads(l) for l in p.stdout.splitlines() if l.startswith('{')]
                rec['build'] = dict(order=4 * v, hadamard=last[-1].get('hadamard') if last else None,
                                    slice_invariant=last[-1].get('slice', {}).get('invariant') if last else None)
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
