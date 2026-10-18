r"""Tally of a theory/m2_prod.py run on the 347 control [M counts, E model].

Unique S-solutions D-checked (the saved S-value index arrays of every middle value, deduplicated across workers by
their bytes), unique cores (by the sha of the core), S-solutions per core against the model 2^D_def (D_def = 16.14
at the true target, theory/m2_levelcal.py), a 90% Poisson interval for the rate, the S-join excess per middle value
(log2 found / predicted), and the wall-clock per core.
usage: python3 theory/m2_tally347.py TAG [--out FILE]"""
import sys, os, json, glob, math
import numpy as np
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D_DEF = 16.14


def poisson_ci(k, conf=0.90):
    """exact two-sided interval for a Poisson mean given k events (chi-square quantiles via scipy if present)"""
    try:
        from scipy.stats import chi2
        a = 1 - conf
        lo = 0.0 if k == 0 else chi2.ppf(a / 2, 2 * k) / 2
        hi = chi2.ppf(1 - a / 2, 2 * k + 2) / 2
        return lo, hi
    except Exception:
        return max(0.0, k - 1.645 * math.sqrt(k)), k + 1.645 * math.sqrt(k + 1)


def main(tag, out=None):
    d = os.path.join(ROOT, 'logs/m2/prod347')
    recs = []
    for fn in sorted(glob.glob(os.path.join(d, f'{tag}_w*.jsonl'))):
        for l in open(fn):
            try: recs.append(json.loads(l))
            except Exception: pass
    done = {r['mid']: r for r in recs if r.get('stage') == 'middle_done' and 'S_solutions' in r}
    prog = defaultdict(int)
    for r in recs:
        if r.get('stage') == 'dprogress': prog[r['mid']] = max(prog[r['mid']], r['next'])
    seen = set(); total = 0; checked_unique = 0; dup = 0
    for fn in sorted(glob.glob(os.path.join(d, f'{tag}_w*_*.npy'))):
        mid = os.path.basename(fn)[:-4].split('_')[-1]
        S = np.load(fn)
        n_ck = len(S) if mid in done else prog.get(mid, 0)
        for i in range(n_ck):
            h = S[i].tobytes(); total += 1
            if h in seen: dup += 1; continue
            seen.add(h); checked_unique += 1
    cores = {}
    for r in recs:
        if 'CORE' in r: cores.setdefault(r['sha'], r)
    sj = [r for r in recs if r.get('stage') == 'sjoin' and r.get('excess_bits') is not None]
    ex = np.array([r['excess_bits'] for r in sj]) if sj else np.zeros(0)
    jsecs = sum(r.get('join_secs', 0) for r in recs if r.get('stage') == 'sjoin')
    dsecs = sum(r.get('dcheck_secs', 0) for r in done.values())
    k = len(cores)
    lo, hi = poisson_ci(k)
    res = dict(tag=tag, S_solutions_checked=total, unique=checked_unique, duplicates_across_workers=dup,
               middles_done=len(done), cores_unique=k, hadamard_1388=sum(1 for c in cores.values() if c.get('hadamard_1388')),
               identical_to_published=sum(1 for c in cores.values() if c.get('identical_to_published')),
               log2_S_per_core=round(math.log2(checked_unique / k), 2) if k else None,
               log2_S_per_core_90=[round(math.log2(checked_unique / hi), 2), round(math.log2(checked_unique / lo), 2) if lo > 0 else None],
               model_log2_S_per_core=D_DEF, expected_cores_model=round(checked_unique * 2.0 ** -D_DEF, 2),
               sjoin_excess_bits=dict(n=int(len(ex)), mean=round(float(ex.mean()), 3) if len(ex) else None,
                                      min=round(float(ex.min()), 3) if len(ex) else None,
                                      max=round(float(ex.max()), 3) if len(ex) else None),
               cpu_secs_join=round(jsecs, 1), cpu_secs_dcheck=round(dsecs, 1),
               ms_per_S_solution=round(1000 * (jsecs + dsecs) / max(1, total), 3),
               cpu_minutes_per_core=round((jsecs + dsecs) / 60 / k, 2) if k else None)
    print(json.dumps(res, indent=1))
    if out: json.dump(res, open(out, 'w'), indent=1)
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[a.index('--out') + 1] if '--out' in a else None)
