r"""Census of the Z_q all-character bordered family at the composite open orders [P screens, M gates, E counts]
(Kumar's construction, search/NOTES.md section 18, https://github.com/vatskishan03/hadamard, Apache-2.0; this driver is
new and reuses theory/bt_census2.py, bt_thm_rigid.py and bt_engine.py).

Family z{q}all (q odd, Q = Z_q; z9x3all: Q = Z_3^2): G abelian of order n - q, K of index q, T = all characters of Q,
target sum_i PAF_i = -4q on K \ 0 and 0 off K, border of size 4q.  Per structure, cheapest first: coset-Parseval for T
[P], Lemma R [P] (theory/bt_thm_rigid.rigid), the Q-level gate [M] (every row pattern with sum r^2 = 4q), the Q-level
stratified count [E].  No prefilter (every structure in the r range is decided by these gates or passes them).
One JSON line per structure, restart-safe.
usage: python3 theory/zq_census.py OUT.jsonl --fam z5all --ns 1115 [--rmin 12] [--rmax 120] [--part i/W]
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census2 as C2, bt_thm_rigid as TR


def evaluate(rec, N=1200):
    t0 = time.time()
    st = C2.build(rec)
    cp = st.coset_parseval()
    if any(d['dead'] for d in cp):
        out = dict(rec); out['key'] = C2.key_of(rec); out['r'] = st.r[0]; out['e'] = st.e
        out['cp'] = [(d['p'], round(d['margin'], 3)) for d in cp]
        out['verdict'] = 'no core [P] (coset-Parseval)'; out['secs'] = round(time.time() - t0, 1)
        return out
    Nr = TR.rigid(st)
    if Nr is not None:
        out = dict(rec); out['key'] = C2.key_of(rec); out['r'] = st.r[0]; out['e'] = st.e
        out['verdict'] = f'no core [P] (Lemma R, |N| = {Nr})'; out['secs'] = round(time.time() - t0, 1)
        return out
    return C2.evaluate(rec, N=N, prefilter=None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--fam', required=True); ap.add_argument('--ns', required=True)
    ap.add_argument('--rmin', type=int, default=12); ap.add_argument('--rmax', type=int, default=120)
    ap.add_argument('--N', type=int, default=1200); ap.add_argument('--part', default='0/1')
    ap.add_argument('--list-only', action='store_true')
    ap.add_argument('--cap', type=float, default=2 ** 20)
    a = ap.parse_args()
    pi, pw = map(int, a.part.split('/'))
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            try: done.add(json.loads(line)['key'])
            except Exception: pass
    k = 0
    for n in map(int, a.ns.split(',')):
        t0 = time.time(); recs = []
        for ip in range(C2.n_two_options(n, a.fam, cap=int(a.cap))):
            cache = f'{a.out}.structs_{n}_{a.rmin}_{a.rmax}_{int(a.cap)}_p{ip}.json'
            if os.path.exists(cache):
                part = json.load(open(cache))
            else:
                part = list(C2.structures(n, a.fam, a.rmin, a.rmax, cap=int(a.cap), only=ip))
                with open(cache + '.tmp', 'w') as fo:
                    json.dump(part, fo)
                os.replace(cache + '.tmp', cache)
            recs += part
        print(json.dumps(dict(n=n, fam=a.fam, structures=len(recs), secs=round(time.time() - t0, 1))), flush=True)
        if a.list_only:
            continue
        for rec in recs:
            k += 1
            if (k - 1) % pw != pi or C2.key_of(rec) in done:
                continue
            try:
                res = evaluate(rec, N=a.N)
            except Exception as e:
                res = dict(rec, key=C2.key_of(rec), verdict=f'error: {type(e).__name__}: {e}')
            with open(a.out, 'a') as f:
                f.write(json.dumps(res) + '\n')
            print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'Pmods', 'Omods', 'order', 'r', 'e', 'verdict',
                                                     'q_solutions', 'log2E_Q', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
