r"""Census of k = 4 (and k4one) structures at n = 971 on G = P x Z_11^2 (|P| = 8) with abelian symmetry groups whose
Z_11^2-part lies in a maximal torus of GL(2, 11): the split torus D = {diag(l, m)} (C_10 x C_10), the Singer torus
C_120 (the multiplicative group of F_121 acting on F_121 = F_11^2), and Z x U = {l [[1, t], [0, 1]]} (C_10 x C_11)
[M orbit counts and gates, E counts; evaluation by theory/bt_census2.evaluate].
usage: python3 theory/bt_census971.py OUT.jsonl [--fam k4] [--rmin 12] [--rmax 120] [--part i/W]
"""
import sys, os, json, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_census as BC, bt_census2 as C2

P11 = 11


def singer_matrix():
    """companion matrix of a primitive quadratic over F_11 (order 120)"""
    for a in range(P11):
        for b in range(1, P11):
            M = np.array([[0, (-b) % P11], [1, (-a) % P11]], np.int64)          # x^2 + a x + b
            X = np.eye(2, dtype=np.int64); o = 0
            for k in range(1, 121):
                X = (X @ M) % P11
                if (X == np.eye(2, dtype=np.int64)).all():
                    o = k; break
            if o == 120:
                return M
    raise RuntimeError


def torus(kind):
    I = np.eye(2, dtype=np.int64)
    if kind == 'split':
        return [np.array([[l, 0], [0, m]], np.int64) for l in range(1, 11) for m in range(1, 11)]
    if kind == 'singer':
        s = singer_matrix(); out = []; X = I.copy()
        for _ in range(120):
            out.append(X.copy()); X = (X @ s) % P11
        return out
    if kind == 'zu':
        return [(l * np.array([[1, t], [0, 1]], np.int64)) % P11 for l in range(1, 11) for t in range(11)]
    raise ValueError(kind)


def structures(fam, rmin, rmax):
    F = C2.FAMS[fam]
    t = F['t']
    m = 971 - t
    assert m == 968
    for popt in C2.two_options_q(3, 4):
        P = popt['P']
        mods = P.mods + (11, 11)
        G = E.Grp(mods)
        Kgens = [tuple(int(c) for c in P.G.co[g]) + (0, 0) for g in popt['KP'] if g] + \
                [(0,) * len(P.mods) + (1, 0), (0,) * len(P.mods) + (0, 1)]
        Kel = G.subgroup(Kgens)
        Kmask = np.zeros(G.n, bool); Kmask[Kel] = True
        cos = -np.ones(G.n, np.int64); qn = 0
        for g in range(G.n):
            if cos[g] < 0:
                cos[G.idx(G.co[g] + G.co[Kel])] = qn; qn += 1
        O = E.Grp((11, 11))
        for label, gam in C2.gamma_options(fam, G, Kel, cos, qn):
            F_tgt = -gam[cos]
            pm = popt['mats']
            Pp = [G.perm_of(E.block_diag(M, np.eye(2, dtype=np.int64))) for M in pm]
            okP = [i for i, p in enumerate(Pp) if (F_tgt[p] == F_tgt).all()]
            TP = BC.AutTable([popt['perms'][i] for i in okP])
            for kind in ('split', 'singer', 'zu'):
                Ts = torus(kind)
                TO = BC.AutTable([O.perm_of(M) for M in Ts])
                res = BC.enumerate_S(TP, TO, rmin, rmax, mode='all' if TP.n * TO.n <= 3000 else 'cyc+odd')
                seen = set()
                for S, gens in res:
                    key = BC.conj_canonical(TP, TO, S)
                    if key in seen:
                        continue
                    seen.add(key)
                    mats = [E.block_diag(pm[okP[ga]], Ts[gb]).tolist() for ga, gb in gens]
                    yield dict(n=971, fam=fam, T=label + f':{kind}', Pmods=list(P.mods), Omods=[11, 11],
                               Kgens=[list(k) for k in Kgens], gamma=gam.tolist(), gens=mats,
                               order=len(S), r=BC.orbit_count(TP, TO, S))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--fam', default='k4')
    ap.add_argument('--rmin', type=int, default=12); ap.add_argument('--rmax', type=int, default=120)
    ap.add_argument('--part', default='0/1'); ap.add_argument('--list-only', action='store_true')
    a = ap.parse_args()
    pi, pw = map(int, a.part.split('/'))
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            try: done.add(json.loads(line)['key'])
            except Exception: pass
    recs = list(structures(a.fam, a.rmin, a.rmax))
    print(f'971 fam {a.fam}: {len(recs)} structures', flush=True)
    if a.list_only:
        return
    for k, rec in enumerate(recs):
        if k % pw != pi or C2.key_of(rec) in done:
            continue
        try:
            res = C2.evaluate(rec)
        except Exception as e:
            res = dict(rec, key=C2.key_of(rec), verdict=f'error: {type(e).__name__}: {e}')
        with open(a.out, 'a') as f:
            f.write(json.dumps(res) + '\n')
        print(json.dumps({x: res.get(x) for x in ('T', 'Pmods', 'order', 'r', 'e', 'rule', 'verdict', 'q_solutions',
                                                 'log2E_Q', 'log2E_whole', 'secs')}), flush=True)


if __name__ == '__main__':
    main()
