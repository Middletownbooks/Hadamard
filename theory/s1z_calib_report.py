r"""Summary of the z-conditioned stage-1 calibration (logs/s1z/calib/*.jsonl, logs/s1z/calib_347.jsonl) [M vs E].

Per rung: the symmetry order (theory/s1z_lib.log2_symmetry), the estimated and exact totals of ordered stage-1
solutions, both also per symmetry class (minus log2 of the order), the per-multiset agreement, and the error of a total
estimated from K multisets drawn with probability proportional to ordered * N0 (importance sampling, 2000 bootstrap
draws), against the all-multiset estimate and against the exact total.
usage: python3 theory/s1z_calib_report.py [OUT.json]"""
import sys, os, json, math, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import log2_symmetry

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(f):
    L = [json.loads(l) for l in open(f)]
    head = L[0]; ms = [x for x in L if x.get('stage') == 'ms']; summ = [x for x in L if x.get('stage') == 'summary']
    return head, ms, (summ[0] if summ else None)


def main():
    files = [os.path.join(ROOT, 'logs/s1z/calib_347.jsonl')] + sorted(glob.glob(os.path.join(ROOT, 'logs/s1z/calib/calib_*.jsonl')))
    rows = []
    rng = np.random.default_rng(5)
    for f in files:
        head, ms, summ = load(f)
        if not summ or not ms: 
            rows.append(dict(rung=os.path.basename(f)[6:-6], v=head['cfg']['v'], multisets=head['multisets'], note='no compression quadruple (stage 1 empty [M])'))
            continue
        cfg = SD.cfg_from_json(json.dumps(head['cfg']))
        sym = log2_symmetry(cfg)['log2_sym']
        have_exact = all('exact' in m for m in ms)
        W = np.array([m['ordered'] * 2.0 ** m['log2N0'] for m in ms])
        P = np.array([2.0 ** (m['log2E'] - m['log2N0']) if m.get('log2E') is not None and np.isfinite(m['log2E']) else 0.0 for m in ms])
        est_all = float((W * P).sum())
        r = dict(rung=os.path.basename(f)[6:-6], v=head['cfg']['v'], L1=head['L1'], rest=head['rest_lags'], multisets=len(ms),
                 log2_sym=sym, log2_est=round(math.log2(est_all), 2) if est_all > 0 else None)
        if have_exact:
            X = np.array([m['ordered'] * m['exact'] for m in ms], float); ex = float(X.sum())
            r['log2_exact'] = round(math.log2(ex), 2) if ex > 0 else None
            r['exact_zero_multisets'] = int(sum(1 for m in ms if m['exact'] == 0))
            r['est_positive_but_exact_zero'] = int(sum(1 for m in ms if m['exact'] == 0 and m.get('log2E') is not None and np.isfinite(m['log2E']) and m['log2E'] > 0))
            if ex > 0 and est_all > 0: r['err_bits'] = round(math.log2(est_all) - math.log2(ex), 2)
            r['est_per_sym'] = round(math.log2(est_all) - sym, 2) if est_all > 0 else None
            r['exact_per_sym'] = round(math.log2(ex) - sym, 2) if ex > 0 else None
            nsol = sum(m['ordered'] * m['exact'] for m in ms); nsym = sum(m['ordered'] * m.get('sym_all', 0) for m in ms)
            r['sym_share'] = round(nsym / nsol, 4) if nsol else None
            # per multiset log-ratio where both positive
            lr = [m['log2E'] - m['log2_exact'] for m in ms if m['exact'] > 0 and m.get('log2E') is not None and np.isfinite(m['log2E'])]
            if lr: r['per_ms_err_bits'] = [round(float(np.min(lr)), 2), round(float(np.median(lr)), 2), round(float(np.max(lr)), 2)]
        else:
            r['est_per_sym'] = round(math.log2(est_all) - sym, 2) if est_all > 0 else None
        # K-subsample bootstrap of the importance estimator
        if est_all > 0 and len(ms) > 1:
            pw = W / W.sum(); res = {}
            for K in (5, 10, 20, 40):
                vals = []
                for _ in range(2000):
                    idx = rng.choice(len(ms), size=K, replace=True, p=pw)
                    e = W.sum() * P[idx].mean()
                    vals.append(math.log2(e) if e > 0 else -99)
                vals = np.array(vals) - math.log2(est_all)
                res[K] = [round(float(np.percentile(vals, 5)), 2), round(float(np.median(vals)), 2), round(float(np.percentile(vals, 95)), 2)]
            r['K_bootstrap_vs_all'] = res
        rows.append(r)
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'logs/s1z/calib/calib_summary.json')
    json.dump(rows, open(out, 'w'), indent=1)
    for r in rows: print(json.dumps(r))


if __name__ == '__main__':
    main()
