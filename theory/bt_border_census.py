r"""Border-problem status of census survivors [M]: for every Q-level solution (all row patterns) of a structure,
is there a border (X, Rh, Ch) for its Q-compressions?  (theory/bt_border.border_for_compressions; frames: Sylvester
frames with column permutations for q = 4, t = 3; Paley H12 frames for q = 3, t = 3; for q = 8, t = 7 the balanced
frames H_4 (x) H_8[1:] with the admissible set enumerated through block-balanced sign vectors,
theory/bt_border.border_balanced.)

A border depends only on the Q-compressions, so this decides the border problem for every core the structure could
have: "all" means every Q-level class has a border, "none" means no core of the structure can be completed (with these
frames; a failure is a failure of the frame search, not a proof).
usage: python3 theory/bt_border_census.py OUT.jsonl CENSUS.jsonl[,...] [--top 30] [--min-logE -10] [--max-frames 3000]
"""
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_gates as BG, bt_border as BB, bt_census_run as R, bt_census2 as C2, bt_ladder as LD


def border_status(st, max_frames=3000, tmax=2.0, maxcls=40):
    g = BG.Gates(st)
    qs = g.q_level(patterns=LD.patterns_for(st), allow_truncated=True)
    qQ = g.qQ
    assert (qQ.pi == st.coset).all(), "coset numbering of the Q quotient differs from QData"
    seen = {}; out = []
    for rows, sol, L in qs:
        for s in sol:
            z = np.array([qQ.full_from_reduced(i, L[i][0][s[i]]) for i in range(4)], np.int64)
            key = tuple(sorted(tuple(sorted(map(int, abs(z[i])))) for i in range(4)))   # coarse class key
            if key in seen:
                continue
            if len(seen) >= maxcls:
                break
            t = time.time()
            if st.q == 8 and st.t == 7:
                res = BB.border_balanced(st, z, max_frames=min(max_frames, 40), tmax=30.0)
            else:
                res = BB.border_for_compressions(st, z, max_frames=max_frames, tmax=tmax)
            seen[key] = res is not None
            out.append(dict(rows=list(map(int, rows)), z=z.tolist(), border=res is not None,
                            frame=None if res is None else int(res[3]), secs=round(time.time() - t, 1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('census')
    ap.add_argument('--top', type=int, default=30); ap.add_argument('--min-logE', type=float, default=-10)
    ap.add_argument('--max-frames', type=int, default=3000)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        for x in open(a.out):
            done.add(json.loads(x)['key'])
    alpha = {}
    if os.path.exists('logs/bt/census/thm_alpha.jsonl'):
        for x in open('logs/bt/census/thm_alpha.jsonl'):
            d = json.loads(x); alpha[d['key']] = d['dead']
    recs = []
    for f in a.census.split(','):
        for x in open(f):
            d = json.loads(x)
            if d.get('verdict', '').startswith('passes') and (d.get('log2E_Q') or -1e9) >= a.min_logE \
                    and not alpha.get(d['key']):
                recs.append(d)
    by_n = {}
    for d in sorted(recs, key=lambda d: -(d.get('log2E_Q') or -1e9)):
        by_n.setdefault((d['n'], d.get('fam', 'k4')), []).append(d)
    for key, lst in sorted(by_n.items()):
        for d in lst[:a.top]:
            if d['key'] in done:
                continue
            st = C2.build(d) if 'fam' in d else R.build(d)
            if not ((st.q == 4 and st.t == 3) or (st.q == 3 and st.t == 3) or (st.q == 8 and st.t == 7)):
                res = dict(status='no frames implemented for (q, t) = (%d, %d)' % (st.q, st.t))
            else:
                t = time.time()
                cl = border_status(st, max_frames=a.max_frames)
                nb = sum(c['border'] for c in cl)
                res = dict(status='all' if nb == len(cl) else ('none' if nb == 0 else 'some'), classes=len(cl),
                           with_border=nb, detail=cl[:6], secs=round(time.time() - t, 1))
            res.update(key=d['key'], n=d['n'], fam=d.get('fam', 'k4'), order=d['order'], r=d['r'], e=d['e'],
                       log2E_Q=d.get('log2E_Q'))
            with open(a.out, 'a') as fo:
                fo.write(json.dumps(res) + '\n')
            print(json.dumps({x: res.get(x) for x in ('n', 'fam', 'order', 'r', 'log2E_Q', 'status', 'classes',
                                                      'with_border', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
