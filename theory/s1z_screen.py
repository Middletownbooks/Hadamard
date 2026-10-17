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
--steal: take the configurations left in OUT from the end of the whole list instead of a fixed part, re-reading OUT
before each one and claiming it with a file in OUT.claims/ (process id and start time; claims of dead processes are
taken over), so finished workers relieve the busy ones.
usage: python3 theory/s1z_screen.py OUT.jsonl [--part i/W | --steal] [--N 3000] [--K 40] [--all 60] [--exact-max 30000]
                                    [--secs 3600] [--cfgs FILE (JSON list of cfgs; default: the 112 skipped)]"""
import sys, os, json, math, time, signal, hashlib
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
    # log2 weights (N0 can pass 2^1024: 643 with L1 = 128)
    lw = np.array([math.log2(q['ordered']) + math.log2(q['n0']) for q in qd]); wmax = float(lw.max())
    W = 2.0 ** (lw - wmax); lN0tot = wmax + math.log2(float(W.sum()))
    rec['log2_N0_total'] = round(lN0tot, 3)
    rng = np.random.default_rng(seed)
    if len(qd) <= nall: pick = list(range(len(qd))); mode = 'all'
    else: pick = rng.choice(len(qd), size=K, replace=True, p=W / W.sum()).tolist(); mode = 'importance'
    vals = []; nexact = 0; ndeg = 0; nzero_exact = 0; ests = []        # vals: (t, log2 p), p = Pr(remaining lags | lag 0)
    lg = lambda bits: -bits if bits is not None and np.isfinite(bits) else -math.inf
    for t in pick:
        q = qd[t]; zq = [cls[i]['rep'] for i in q['m']]
        sizes = [sum(cls[i]['E'].values()) for i in q['m']]
        if max(sizes) <= exact_max:
            ndeg += 1
            try:
                ex = Z.exact_count(zq)[0]; nexact += 1; nzero_exact += ex == 0
                lp = math.log2(ex) - math.log2(q['n0']) if ex else -math.inf; ests.append(dict(i=t, exact=ex))
            except (RuntimeError, MemoryError):
                r = Z.estimate(zq, N=N, seed=seed + t); lp = lg(r.get('bits'))
                ests.append(dict(i=t, log2E=r.get('log2E'), degenerate=True))
        else:
            r = Z.estimate(zq, N=N, seed=seed + t); lp = lg(r.get('bits'))
            ests.append(dict(i=t, log2E=r.get('log2E'), rank=r.get('rank'), mahal=r.get('mahal')))
        vals.append((t, lp))
        if log: log(dict(stage='ms', m=list(q['m']), **ests[-1], secs=round(time.time() - t0, 1)))
    def lsum(xs):
        xs = [x for x in xs if x > -math.inf]
        if not xs: return -math.inf
        m = max(xs); return m + math.log2(sum(2.0 ** (x - m) for x in xs))
    if mode == 'all':
        ltot = lsum([float(lw[t]) + lp for t, lp in vals])
    else:
        ltot = lN0tot + lsum([lp for _, lp in vals]) - math.log2(len(vals))
    rec.update(mode=mode, evaluated=len(pick), exact_multisets=nexact, exact_zero=nzero_exact,
               log2_stage1_total=round(ltot, 3) if ltot > -math.inf else None)
    if mode == 'importance':
        lps = np.array([lp for _, lp in vals]); m = lps.max()
        lp = 2.0 ** (lps - m) if m > -math.inf else np.zeros(len(lps))
        se = lp.std(ddof=1) / math.sqrt(len(lp)) if len(lp) > 1 else float('nan')
        rec['rel_se'] = round(float(se / lp.mean()), 3) if lp.mean() > 0 else None
    rec['secs'] = round(time.time() - t0, 1)
    return rec


def read_done(outp):
    """configurations with a record in OUT; records with an error are redone (once per process)"""
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try:
                r = json.loads(l)
                if 'error' not in r: done.add(json.dumps(r['cfgkey']))
            except Exception: pass
    return done


def proc_start(pid):
    try: return open(f'/proc/{pid}/stat').read().rsplit(')', 1)[1].split()[19]
    except OSError: return None


def claim(cdir, key):
    """claim a configuration for this process; a claim whose process is gone is taken over"""
    path = os.path.join(cdir, hashlib.md5(key.encode()).hexdigest()[:16] + '.claim')
    for _ in range(2):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, json.dumps(dict(pid=os.getpid(), start=proc_start(os.getpid()), cfgkey=key)).encode()); os.close(fd)
            return path
        except FileExistsError:
            try: c = json.load(open(path))
            except Exception: c = {}
            if c.get('pid') is not None and proc_start(c['pid']) == c.get('start'): return None
            try: os.remove(path)
            except OSError: pass
    return None


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
    done = read_done(outp)
    cfgs = sorted(cfgs, key=lambda c: (c.get('L1', 0), c['v']))
    steal = '--steal' in a; cdir = outp + '.claims'
    if steal: os.makedirs(cdir, exist_ok=True); cfgs = cfgs[::-1]
    else: cfgs = [c for i, c in enumerate(cfgs) if i % nW == part]
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    tried = set()
    while True:
        cl = None
        if steal:
            done = read_done(outp) | tried
            for c in cfgs:
                k = json.dumps(cfgkey(c))
                if k not in done and (cl := claim(cdir, k)): cfg = c; break
            if cl is None: break
        else:
            todo = [c for c in cfgs if json.dumps(cfgkey(c)) not in done]
            if not todo: break
            cfg = todo[0]
        key = cfgkey(cfg); done.add(json.dumps(key)); tried.add(json.dumps(key))
        rec = dict(cfgkey=key); t0 = time.time(); signal.alarm(secs)
        try:
            rec.update(screen(cfg, N=N, K=K, nall=nall, exact_max=exact_max))
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:200]
        signal.alarm(0); rec['secs'] = round(time.time() - t0, 1)
        with open(outp, 'a') as f: f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
        if cl:
            try: os.remove(cl)
            except OSError: pass


if __name__ == '__main__':
    main()
