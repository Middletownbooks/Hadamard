r"""The k = 6 border problem (q = 6, t = 5, T = all nontrivial characters of Q): a control on Alpoge's 347 template [M].

The 347 template (data/template347.npz: compressions, frames K and Kc, corner C, and the GS layout of the published
H(1388)) has non-monomial compressions on Q = Z_2 x Z_3.
1. In its own layout, with Q in (Z_2, Z_3) coordinates, the corner is reproduced (theory/ab_assemble.forced_corner),
   and either template frame used as the frame Rh gives |S_R| = 20 = 4t: the admissible set is exactly the corner.
2. The frame family of Lemma 1 (theory/bt_border.frames_lemma1: H_2 (x) Paley H(12), every choice of the row s, then
   random column permutations within blocks) gives no border for the same compressions.  Its largest admissible set
   is far below 20, both in the 347 layout and in Kumar's layout (theory/bt_border.gs_array).
So, as at k = 8 (theory/bt_border_k8.py), the simple frame families do not decide the k = 6 border problem: borders
exist (347), but they need frames of template type.  The k = 6 survivors' Q-level classes are not monomial, so the
Singer-type search of theory/ab_singer.py does not apply to them either; their border status is open.
S_R by a meet in the middle over the four blocks (theory/bt_border.admissible_balanced), about 5 ms per frame.
usage: python3 theory/bt_border_k6.py [--frames 3000] [--survivors 3] [--out logs/bt/border_k6.jsonl]
"""
import sys, os, json, time, itertools, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, ab_assemble as AS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
q, t = 6, 5


def search(Mc, frames, tmax=5.0):
    best = 0; found = None; k = -1
    for k, F in enumerate(frames):
        S = BB.admissible_balanced(Mc, F, t, q)
        best = max(best, len(S))
        if len(S) >= 4 * t and BB.clique(S, 4 * t, tmax=tmax) is not None:
            found = k
            break
    return best, found, k + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--frames', type=int, default=3000)
    ap.add_argument('--survivors', type=int, default=0, help='also test the top N alive k6 structures per order')
    ap.add_argument('--out', default=os.path.join(ROOT, 'logs/bt/border_k6.jsonl'))
    a = ap.parse_args()
    T = np.load(os.path.join(ROOT, 'data/template347.npz'), allow_pickle=True)
    z = T['comp'].astype(np.int64); L = json.loads(str(T['layout']))
    K, Kc = T['K'].astype(np.int64), T['Kc'].astype(np.int64)
    out = []
    # 1. the template in its own layout
    Mc = AS.compressed_core(z, (2, 3), L)
    Cf = AS.forced_corner(z, (2, 3), L, K, Kc)
    rec = dict(case='347 layout, Q = Z2 x Z3', corner_reproduced=bool(Cf is not None and (Cf == T['C']).all()),
               template_frames_admissible={nm: len(BB.admissible_balanced(Mc, F, t, q)) for nm, F in (('K', K), ('Kc', Kc))})
    t0 = time.time()
    best, found, tried = search(Mc, BB.frames_lemma1(q, count=a.frames, seed=5))
    rec.update(lemma1_frames=tried, lemma1_max_admissible=best, lemma1_border_frame=found, secs=round(time.time() - t0, 1))
    print(json.dumps(rec), flush=True); out.append(rec)
    # 2. Kumar's layout (the one used by theory/bt_border.assemble), same compressions in Z6 coordinates (CRT)
    I = np.arange(q)
    add = (I[:, None] + I[None, :]) % q; neg = (-I) % q
    crt = np.array([3 * (x % 2) + (x % 3) for x in range(q)])     # Z6 element x -> (x mod 2, x mod 3) index
    z6 = z[:, crt]
    Mk = BB.gs_array(list(z6), add, neg)
    rec = dict(case="Kumar's layout, Q = Z6")
    rec['template_frames_admissible'] = {}
    for nm, F0 in (('K', K), ('Kc', Kc)):
        F6 = F0.reshape(4 * t, 4, q)[:, :, crt].reshape(4 * t, 4 * q)
        rec['template_frames_admissible'][nm] = max(len(BB.admissible_balanced(Mk, F6 * np.repeat(np.array(sg), q)[None, :], t, q))
                                                    for sg in itertools.product((1, -1), repeat=4))
    t0 = time.time()
    best, found, tried = search(Mk, BB.frames_lemma1(q, count=a.frames, seed=5))
    rec.update(lemma1_frames=tried, lemma1_max_admissible=best, lemma1_border_frame=found, secs=round(time.time() - t0, 1))
    print(json.dumps(rec), flush=True); out.append(rec)
    # 3. the survivors: template frames (block signs, random column permutations within blocks) on their classes
    if a.survivors:
        out += survivors(K, Kc, crt, a.survivors, a.frames)
    with open(a.out, 'w') as fo:
        for r in out:
            fo.write(json.dumps(r) + '\n')


def template_variants(K, Kc, crt, count, seed=1):
    rng = np.random.default_rng(seed)
    base = [F0.reshape(4 * t, 4, q)[:, :, crt].reshape(4 * t, 4 * q) for F0 in (K, Kc)]
    k = 0
    while k < count:
        F = base[k % 2] * np.repeat(rng.choice([-1, 1], 4), q)[None, :]
        cols = np.concatenate([J * q + (np.arange(q) if k < 32 else rng.permutation(q)) for J in range(4)])
        yield F[:, cols]; k += 1


def survivors(K, Kc, crt, top, frames):
    import bt_census2 as C2, bt_gates as BG, bt_ladder as LD, bt_summary as SU
    D = SU.Data(); out = []
    for (n, fam), Ls in sorted(D.groups().items()):
        if fam != 'k6':
            continue
        for d in sorted([d for d in Ls if D.status(d) == 'alive'], key=lambda d: -SU.e_of(d))[:top]:
            st = C2.build(d); g = BG.Gates(st); qQ = g.qQ
            Q = BB.QData(st)
            seen = set(); res = []
            for rows, sol, Lq in g.q_level(patterns=LD.patterns_for(st)):
                for s in sol:
                    z = np.array([qQ.full_from_reduced(i, Lq[i][0][s[i]]) for i in range(4)], np.int64)
                    key = tuple(sorted(tuple(sorted(map(int, abs(z[i])))) for i in range(4)))
                    if key in seen or len(seen) >= 6:
                        continue
                    seen.add(key)
                    Mc = BB.gs_array(list(z), Q.add, Q.neg)
                    t0 = time.time()
                    best, found, tried = search(Mc, template_variants(K, Kc, crt, frames))
                    res.append(dict(z=z.tolist(), max_admissible=best, border_frame=found, frames=tried,
                                    secs=round(time.time() - t0, 1)))
            rec = dict(case='k6 survivor', key=d['key'], n=n, order=d['order'], r=d['r'], log2E_Q=d.get('log2E_Q'),
                       classes=len(res), with_border=sum(x['border_frame'] is not None for x in res),
                       max_admissible=max([x['max_admissible'] for x in res], default=0), detail=res)
            print(json.dumps({x: rec[x] for x in ('case', 'n', 'order', 'r', 'log2E_Q', 'classes', 'with_border',
                                                   'max_admissible')}), flush=True)
            out.append(rec)
    return out


if __name__ == '__main__':
    main()
