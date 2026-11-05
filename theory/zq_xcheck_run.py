r"""CP-SAT gate on the compressed levels of the alive Z_q structures [M] (theory/bt_xcheck_level.py, which shares no
code with theory/bt_*, run as a subprocess on a JSON spec).  For every S-invariant N <= K (1 < |N| < |K|, coarsest
first) and every Q-level row pattern of the structure, the level-G/N system is posed to CP-SAT with a time limit;
a level whose every row pattern is INFEASIBLE closes the structure.  One JSON line per (structure, level),
restart-safe.  The engine (bt_engine) is used only to list the invariant subgroups.
usage: python3 theory/zq_xcheck_run.py TAG [--secs 60] [--keys k1,k2]    (TAG = fam_n, e.g. z7all_1169)
"""
import sys, os, json, time, argparse, subprocess, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census2 as C2, zq_summary as SU

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def gens_of(G, N):
    """a small generating set of the subgroup N (element indices) of G"""
    N = set(int(x) for x in N); gens = []; span = {0}
    for x in sorted(N, key=lambda x: -G.order(x)):
        if x in span:
            continue
        gens.append(x)
        new = set(span)
        while True:
            add = {int(G.idx(G.co[a] + G.co[g])) for a in new for g in gens} | new
            if add == new:
                break
            new = add
        span = new
        if span == N:
            break
    return [G.co[g].tolist() for g in gens]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tag'); ap.add_argument('--secs', type=float, default=60)
    ap.add_argument('--keys', default=None)
    a = ap.parse_args()
    C, Gt = SU.load(a.tag)
    out = f'{ROOT}/logs/zq/xcheck/{a.tag}.jsonl'
    done = set(); dead = set()
    if os.path.exists(out):
        for L in open(out):
            d = json.loads(L)
            if d['status'] == 'dead' or d['secs_limit'] >= a.secs:
                done.add((d['key'], d['N']))
            if d['status'] == 'dead':
                dead.add(d['key'])
    keys = a.keys.split(',') if a.keys else [k for k, d in C.items() if SU.status(d, Gt.get(k, [])) == 'alive']
    for k in keys:
        if k in dead:
            continue
        rec = C[k]; st = C2.build(rec); G = st.G
        Ns = [N for N in st.invariant_subgroups() if 1 < len(N) < len(st.Kel) and np.isin(N, st.Kel).all()]
        Ns.sort(key=lambda N: -len(N))
        for N in Ns:
            gens = gens_of(G, N)
            Nid = f'{len(N)}:' + ';'.join(','.join(map(str, g)) for g in gens)
            if (k, Nid) in done:
                continue
            spec = dict(mods=rec['Pmods'] + rec['Omods'], Kgens=rec['Kgens'], gamma=rec['gamma'], gens=rec['gens'],
                        Ngens=gens, rows=rec['q_patterns'], secs=a.secs)
            with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
                json.dump(spec, f); fn = f.name
            t0 = time.time()
            try:
                p = subprocess.run([sys.executable, os.path.join(ROOT, 'theory/bt_xcheck_level.py'), fn],
                                   capture_output=True, text=True, timeout=a.secs * (len(rec['q_patterns']) + 1) + 600)
                res = [json.loads(L) for L in p.stdout.splitlines() if L.startswith('{')]
            except subprocess.TimeoutExpired:
                res = []
            os.unlink(fn)
            sts = [r['status'] for r in res]
            status = 'dead' if res and len(res) == len(rec['q_patterns']) and all(s == 'INFEASIBLE' for s in sts) \
                else ('feasible' if any(s in ('OPTIMAL', 'FEASIBLE') for s in sts) else 'undecided')
            d = dict(key=k, n=rec['n'], fam=rec['fam'], order=rec['order'], r=rec['r'], log2E_Q=rec.get('log2E_Q'),
                     N=Nid, Nsize=len(N), status=status, per_pattern=res, secs_limit=a.secs,
                     secs=round(time.time() - t0, 1))
            with open(out, 'a') as f:
                f.write(json.dumps(d) + '\n')
            print(json.dumps({x: d[x] for x in ('n', 'fam', 'order', 'r', 'log2E_Q', 'Nsize', 'status', 'secs')}),
                  flush=True)
            if status == 'dead':
                break


if __name__ == '__main__':
    main()
