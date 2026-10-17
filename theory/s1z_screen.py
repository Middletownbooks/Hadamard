r"""z-conditioned stage-1 screen over the b >= 3 configurations the exact batches skipped (more than 26 G'-orbits)
[counts exact; remaining lags E].  Replaces the never-run theory/ab_s1_zscreen.py, whose rows came from screen A's list
(truncated at 40 quadruples, empty when screen A gives up) and whose totals extrapolated from K first-found or random
compression quadruples without their multiplicities.

Per configuration (theory/s1z_lib.ZOrb):
  * rows: every achievable row sum in a quadruple with sum r^2 = 4v; compressions with the Parseval and per-character
    bounds; classes under the per-block symmetries; class multisets with sum PAF_Q = T_Q; the exact lag-0 count N0 of
    each multiset and the number of ordered compression quadruples it stands for.  No multiset with N0 > 0 => stage 1
    is empty [M].
  * stage-1 total [E]: sum over multisets of ordered * N0 * 2^-bits (lattice Gaussian of the remaining lag classes,
    theory/s1z_lib.ZOrb.estimate, finite-sample log-det correction).  All multisets when there are at most --all of them,
    otherwise --K multisets drawn with probability proportional to ordered * N0 (importance sample; the total is
    N0_total times the mean Gaussian probability).  A multiset whose four block lists are small (fewer than --exact-max
    y's per block) is counted exactly instead (the Gaussian is unreliable there: theory/s1z_calib.py, the v = 131 rung).
Restart-safe JSONL (configurations already in OUT are skipped); time cap per configuration.
usage: python3 theory/s1z_screen.py OUT.jsonl [--part i/W] [--N 3000] [--K 40] [--all 60] [--exact-max 30000]
                                    [--secs 3600] [--cfgs FILE (JSON list of cfgs; default: the 112 skipped)]"""
import sys, os, json, math, time, signal
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import ZOrb
from s1z_budget_batch import skipped_configs, cfgkey


class TO(Exception): pass


def screen(cfg, N=3000, K=40, nall=60, exact_max=30000, seed=1, log=None):
    t0 = time.time()
    c2 = SD.cfg_from_json(json.dumps(cfg))
    Z = ZOrb(c2)
    rows = Z.rows(); C = Z.compressions(rows); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
    cache = {}
    for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
    qd = [q for q in qd if q['n0'] > 0]
    rec = dict(L1=Z.L1, rest_lags=len(Z.rest), letters=sorted({x for L in Z.letters for x in L}), rows=len(rows),
               compressions=len(C), classes=len(cls), multisets=len(qd), ordered_zquads=sum(q['ordered'] for q in qd))
    if not qd:
        rec.update(stage1_empty=True, reason='no compression quadruple meets the Q-level and lag-0 equations [M]')
        return rec
    W = np.array([float(q['ordered']) * float(q['n0']) for q in qd]); N0tot = float(W.sum())
    rec['log2_N0_total'] = round(math.log2(N0tot), 3)
    rng = np.random.default_rng(seed)
    if len(qd) <= nall: pick = list(range(len(qd))); mode = 'all'
    else: pick = rng.choice(len(qd), size=K, replace=True, p=W / N0tot).tolist(); mode = 'importance'
    vals = []; nexact = 0; ndeg = 0; nzero_exact = 0; ests = []
    for t in pick:
        q = qd[t]; zq = [cls[i]['rep'] for i in q['m']]
        sizes = [sum(cls[i]['E'].values()) for i in q['m']]
        if max(sizes) <= exact_max:
            ndeg += 1
            try:
                ex = Z.exact_count(zq)[0]; nexact += 1; nzero_exact += ex == 0
                p = ex / q['n0']; ests.append(dict(i=t, exact=ex))
            except (RuntimeError, MemoryError):
                r = Z.estimate(zq, N=N, seed=seed + t); p = 2.0 ** -r['bits'] if r.get('bits') is not None and np.isfinite(r['bits']) else 0.0
                ests.append(dict(i=t, log2E=r.get('log2E'), degenerate=True))
        else:
            r = Z.estimate(zq, N=N, seed=seed + t)
            p = 2.0 ** -r['bits'] if r.get('bits') is not None and np.isfinite(r['bits']) else 0.0
            ests.append(dict(i=t, log2E=r.get('log2E'), rank=r.get('rank'), mahal=r.get('mahal')))
        vals.append((t, p))
        if log: log(dict(stage='ms', i=t, m=list(q['m']), **ests[-1], secs=round(time.time() - t0, 1)))
    if mode == 'all':
        tot = sum(qd[t]['ordered'] * qd[t]['n0'] * p for t, p in vals)
    else:
        tot = N0tot * float(np.mean([p for _, p in vals]))
    rec.update(mode=mode, evaluated=len(pick), exact_multisets=nexact, exact_zero=nzero_exact,
               log2_stage1_total=round(math.log2(tot), 3) if tot > 0 else None)
    if mode == 'importance':
        lp = np.array([p for _, p in vals]); se = lp.std(ddof=1) / math.sqrt(len(lp)) if len(lp) > 1 else float('nan')
        rec['rel_se'] = round(float(se / lp.mean()), 3) if lp.mean() > 0 else None
    rec['secs'] = round(time.time() - t0, 1)
    return rec


def main():
    a = sys.argv[1:]; outp = a[0]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    N = int(opt('--N', '3000')); K = int(opt('--K', '40')); nall = int(opt('--all', '60'))
    exact_max = int(opt('--exact-max', '30000')); secs = int(opt('--secs', '3600'))
    import resource
    memgb = float(opt('--memgb', '3.5'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    cfgs = json.load(open(opt('--cfgs', ''))) if '--cfgs' in a else skipped_configs()
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfgkey']))
            except Exception: pass
    cfgs = sorted(cfgs, key=lambda c: (c.get('L1', 0), c['v']))
    cfgs = [c for i, c in enumerate(cfgs) if i % nW == part]
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for cfg in cfgs:
        key = cfgkey(cfg)
        if json.dumps(key) in done: continue
        rec = dict(cfgkey=key); t0 = time.time(); signal.alarm(secs)
        try:
            rec.update(screen(cfg, N=N, K=K, nall=nall, exact_max=exact_max))
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
        signal.alarm(0); rec['secs'] = round(time.time() - t0, 1)
        with open(outp, 'a') as f: f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == '__main__':
    main()
