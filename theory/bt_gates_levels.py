r"""Exact gates at chosen levels, with a time limit instead of a cost cap [M].

The pair-count bound of theory/bt_cost.py can overestimate a gate by many orders of magnitude (at 669 a level bounded
by 2^37 pairs ran in 8 s).  This runner takes the levels by the order of N (every S-invariant N <= K of that order;
'auto': every invariant N <= K, coarsest level first) and the structures by key or by family, and runs
theory/bt_gates_run.run_level under an alarm.  A level with no solution closes the structure: no core [M].  One JSON
line per (structure, level), restart-safe; output in the format of bt_gates_run (read by theory/bt_summary.py).
A level capped under a shorter time limit counts as untried under a longer one, so a pass at --tmax 30 followed by
one at --tmax 300 retries only the levels the first pass could not finish.  Likewise a level whose allocation was
refused under an address-space limit (ulimit -v; recorded as 'capped (memory limit N MB)', unlike the list caps of
bt_chain) counts as untried under a larger limit, and a level over a list cap ('refined list x > N') counts as
untried under a larger --maxlen when its recorded size fits.  --part i/W runs every W-th key from the i-th, for W parallel
runners sharing one output file; --retry-only skips the keys whose records are all final under the current limits.
usage: python3 theory/bt_gates_levels.py OUT.jsonl CENSUS.jsonl[,...] --levels 111,74|auto
       (--keys k1,k2|@FILE | --fam n:fam) [--tmax 900] [--maxlen 2e7] [--seen-all] [--part i/W] [--retry-only]
"""
import sys, os, re, json, time, math, argparse, signal, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census_run as R, bt_census2 as C2, bt_cost as CO, bt_gates_run as GR
import bt_summary as SU


class Timeout(Exception):
    pass


FIRED = [False]


def _alarm(signum, frame):
    FIRED[0] = True
    raise Timeout()


def _mem_limit_mb():
    lim = resource.getrlimit(resource.RLIMIT_AS)[0]
    return None if lim == resource.RLIM_INFINITY else lim // 2 ** 20


def _tried(d, tmax, mem_mb, maxlen=2e7):
    """A recorded level counts as tried unless it was capped under a shorter time limit than tmax, under a smaller
    address-space limit than mem_mb (None: unlimited), or under a smaller list cap than maxlen.  An allocation failure
    recorded without its limit (older runners) counts as untried only under a limit of 8 GiB or more."""
    v = d.get('verdict', '')
    m = re.match(r'capped \(time limit ([\d.]+) s\)', v)
    if m:
        return float(m.group(1)) >= tmax
    m = re.match(r'capped \(memory limit (\d+) MB\)', v)
    if m:
        return mem_mb is not None and int(m.group(1)) >= mem_mb
    m = re.match(r'capped \(refined list (\S+) > (\d+)\)', v)
    if m:                                       # rerun only a list that fits under the larger cap
        return not (int(m.group(2)) < maxlen and float(m.group(1)) <= maxlen)
    if re.match(r'capped \((?:refined cell list|level list) > \d+\)', v):
        return True                             # size not recorded: on 739 k4 such reruns took 7-200 s and all re-capped
    if v.startswith(('capped (Unable to allocate', 'capped (Allocation failed')):
        return not (mem_mb is not None and mem_mb >= 8192)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('census')
    ap.add_argument('--levels', required=True); ap.add_argument('--keys', default=None)
    ap.add_argument('--fam', default=None, help='n:fam, every alive structure of that order and family')
    ap.add_argument('--tmax', type=float, default=900); ap.add_argument('--maxlen', type=float, default=2e7)
    ap.add_argument('--seen-all', action='store_true',
                    help='also skip (structure, level) pairs recorded in any logs/bt/gates_run*.jsonl')
    ap.add_argument('--part', default=None, help='i/W: every W-th key from the i-th (W parallel runners)')
    ap.add_argument('--retry-only', action='store_true',
                    help='skip keys that have records, all of them final under the current limits')
    a = ap.parse_args()
    mem_mb = _mem_limit_mb()
    levels = None if a.levels == 'auto' else [int(x) for x in a.levels.split(',')]
    recs = {}
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x)
            if 'key' in d:
                recs[d['key']] = d
    if a.keys:
        ks = open(a.keys[1:]).read().replace('\n', ',') if a.keys.startswith('@') else a.keys
        keys = [k for k in ks.split(',') if k]
    else:
        n, fam = a.fam.split(':')
        D = SU.Data()
        keys = [d['key'] for d in sorted(D.groups().get((int(n), fam), []), key=lambda d: -SU.e_of(d))
                if D.status(d) == 'alive']
    if a.part:
        i, W = (int(x) for x in a.part.split('/'))
        keys = keys[i::W]
    done = set(); dead = set(); seen = set(); untried = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            d = json.loads(x)
            seen.add(d['key'])
            (done if _tried(d, a.tmax, mem_mb, a.maxlen) else untried).add((d['key'], d.get('Nkey')))
            if d.get('verdict', '').startswith('no core'):
                dead.add(d['key'])
    if a.seen_all:                                             # levels another runner already tried or decided
        import glob
        for f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(a.out)), 'gates_run*.jsonl')):
            for x in open(f):
                d = json.loads(x)
                seen.add(d['key'])
                if d.get('Nkey'):
                    (done if _tried(d, a.tmax, mem_mb, a.maxlen) else untried).add((d['key'], d['Nkey']))
                if d.get('verdict', '').startswith('no core'):
                    dead.add(d['key'])
    if a.retry_only:                                           # keys never run, or with a level to rerun
        retry = {k for k, nk in untried if (k, nk) not in done}
        keys = [k for k in keys if k not in seen or k in retry]
    signal.signal(signal.SIGALRM, _alarm)
    for key in keys:
        if key in dead or key not in recs:
            continue
        rec = recs[key]
        try:
            st = C2.build(rec) if 'fam' in rec else R.build(rec)
        except MemoryError:                # memory still held from an earlier level: resume in a fresh process
            sys.exit(EXIT_FRESH)
        lv = levels
        if lv is None:                                           # every invariant N inside K, coarsest level first
            lv = sorted({len(N) for N in st.invariant_subgroups() if 1 < len(N) < len(st.Kel)
                         and np.isin(N, st.Kel).all()}, reverse=True)
        for Nsz in lv:
            if key in dead:
                break
            for NA in [N for N in st.invariant_subgroups() if len(N) == Nsz and np.isin(N, st.Kel).all()]:
                nk = CO.nkey(NA)
                if (key, nk) in done:
                    continue
                t = time.time()
                FIRED[0] = False
                real = False
                signal.alarm(int(a.tmax))
                try:
                    try:
                        # the alarm cannot interrupt a compiled join, so joins that could not finish in time are
                        # refused up front: pair cap = measured rate (about 3e6 pairs/s) x time limit
                        res = GR.run_level(st, NA, cap_pairs=max(2.0 ** 30, 3e6 * a.tmax), maxlen=int(a.maxlen))
                        res['verdict'] = 'no core [M] (level empty)' if res['solutions_weighted'] == 0 else 'passes'
                    except MemoryError as e:                     # a list cap of bt_chain, or the address-space limit
                        msg = str(e)
                        real = not msg or msg.startswith(('Unable to allocate', 'Allocation failed'))
                        res = dict(verdict=f'capped (memory limit {mem_mb} MB)' if real and mem_mb else f'capped ({e})')
                    except Exception:
                        if not FIRED[0]:                         # inside numba the Timeout arrives as a SystemError
                            raise
                        res = dict(verdict=f'capped (time limit {a.tmax:.0f} s)')
                    finally:
                        signal.alarm(0)
                except Timeout:      # the alarm fired as the level ended, in an except clause or the finally (Oct 28)
                    signal.alarm(0)
                    real = False
                    res = dict(verdict=f'capped (time limit {a.tmax:.0f} s)')
                res.update(key=key, n=rec['n'], fam=rec.get('fam', 'k4'), order=rec['order'], r=rec['r'], e=rec['e'],
                           log2E_Q=rec.get('log2E_Q'), N=int(Nsz), Nkey=nk, log2_pairs=None,
                           secs=round(time.time() - t, 1))
                with open(a.out, 'a') as fo:
                    fo.write(json.dumps(res) + '\n')
                done.add((key, nk))
                print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'log2E_Q', 'N', 'verdict',
                                                          'solutions_weighted', 'secs')}), flush=True)
                if res['verdict'].startswith('no core'):
                    dead.add(key)
                    break
                if real:                    # after a refused allocation the memory can stay held (every later
                    sys.exit(EXIT_FRESH)    # allocation then fails at once): resume in a fresh process instead
                if res['verdict'].startswith('capped (time limit') and res['secs'] > a.tmax + 10:
                    sys.exit(EXIT_FRESH)    # a compiled join overran the alarm; its memory can stay held too


EXIT_FRESH = 75                             # closure.sh restarts a runner that exits with this code


if __name__ == '__main__':
    main()
