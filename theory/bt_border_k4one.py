r"""Border status of the k4one family (Q of order 4, T = {1} and two nontrivial characters, t = 3) [M].

V_T is the space of block vectors with Fourier support in T.  On one block it is spanned by the +-1 vectors c with
<c, chi> = 0 for every character chi of Q outside T: three of them, c_0 = 1 and c_1, c_2 (for Q = Z_2^2 the two
characters in T, for Q = Z_4 the vectors (1, 1, -1, -1) and (1, -1, -1, 1), which span the order-4 characters).
A frame row restricted to a block is then +-c_k for some k, and the frames searched are
    Rh[(i, k), (J, .)] = H_4[i, J] s_J c_{sigma_J(k)}      (sigma_J a permutation of {0, 1, 2}, sigma_0 = 1; s_J = +-1)
whenever these 12 rows are orthogonal (6^3 * 2^4 candidates).  For each frame S_R (theory/bt_border.admissible, all
2^11 sign vectors) is searched for 12 orthogonal vectors, and a border is accepted only after the exact checks
(Ch in {+-1}, Ch Ch^T = I_4 (x) Gamma_Q).  With sigma_J = 1 for every block (the frames H_4 (x) C) no class of the
863 hosts has a border; with the block permutations every class tested has one.
usage: python3 theory/bt_border_k4one.py [--min-logE 20] [--maxcls 8] [--out logs/bt/border_k4one.jsonl]
"""
import sys, os, json, glob, itertools, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_census2 as C2, bt_gates as BG, bt_ladder as LD, bt_border as BB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H4 = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])


def t_basis(st, Q):
    """three orthogonal +-1 vectors spanning V_T on one block (cosets in QData numbering)"""
    q = st.q
    I = np.arange(q)
    # characters of Q as complex vectors: from the addition table, via a generator search over q-th roots of unity
    chars = []
    for vals in itertools.product(range(q), repeat=q):
        c = np.exp(2j * np.pi * np.array(vals) / q)
        if all(abs(c[Q.add[a, b]] - c[a] * c[b]) < 1e-9 for a in I for b in I):
            chars.append(c)
    out_T = [c for c in chars if abs(np.dot(st.gamma, np.conj(c)) / (4 * q)) < 1e-9]   # characters outside T
    cand = [np.array(v) for v in itertools.product((1, -1), repeat=q)
            if all(abs(np.dot(np.array(v), np.conj(c))) < 1e-9 for c in out_T)]
    for B in itertools.combinations(cand, 3):
        M = np.array(B)
        if (M @ M.T == q * np.eye(3)).all():
            return M
    return None


def frames(Cm):
    perms = list(itertools.permutations(range(3)))
    for sg in itertools.product(perms, repeat=3):
        for sc in itertools.product((1, -1), repeat=4):
            Hs = H4 * np.array(sc)[None, :]
            Rh = np.zeros((12, 16), np.int64)
            for i in range(4):
                for k in range(3):
                    for J in range(4):
                        kk = k if J == 0 else sg[J - 1][k]
                        Rh[3 * i + k, 4 * J:4 * J + 4] = Hs[i, J] * Cm[kk]
            if (Rh @ Rh.T == 16 * np.eye(12, dtype=np.int64)).all():
                yield Rh


def border(Mc, Gm, Cm, t):
    for Rh in frames(Cm):
        S = BB.admissible(Mc, Rh, t)
        if len(S) < 4 * t:
            continue
        Xc = BB.clique(S, 4 * t, tmax=2)
        if Xc is None:
            continue
        X = Xc.T.copy(); num = -(Mc @ Rh.T @ X)
        if (num % (4 * t)).any():
            continue
        Ch = num // (4 * t)
        if (np.abs(Ch) == 1).all() and (Ch @ Ch.T == np.kron(np.eye(4, dtype=np.int64), Gm)).all() \
                and (X @ X.T == 4 * t * np.eye(4 * t, dtype=np.int64)).all():
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-logE', type=float, default=20); ap.add_argument('--maxcls', type=int, default=8)
    ap.add_argument('--out', default=os.path.join(ROOT, 'logs/bt/border_k4one.jsonl'))
    a = ap.parse_args()
    recs = []
    for f in glob.glob(os.path.join(ROOT, 'logs/bt/census2/k4one_w*.jsonl')):
        for x in open(f):
            d = json.loads(x)
            if (d.get('log2E_Q') or -1e9) >= a.min_logE:
                recs.append(d)
    out = []
    for d in sorted(recs, key=lambda d: -d['log2E_Q']):
        t0 = time.time()
        st = C2.build(d); Q = BB.QData(st)
        Cm = t_basis(st, Q)
        g = BG.Gates(st); qs = g.q_level(patterns=LD.patterns_for(st), allow_truncated=True); qQ = g.qQ
        seen = set(); res = []
        for rows, sol, L in qs:
            for s in sol:
                z = np.array([qQ.full_from_reduced(i, L[i][0][s[i]]) for i in range(4)], np.int64)
                key = tuple(sorted(tuple(sorted(map(int, abs(z[i])))) for i in range(4)))
                if key in seen or len(seen) >= a.maxcls:
                    continue
                seen.add(key)
                res.append(border(BB.gs_array(list(z), Q.add, Q.neg), Q.Gmat, Cm, st.t))
        rec = dict(key=d['key'], n=d['n'], fam='k4one', Pmods=d['Pmods'], T=d.get('T'), order=d['order'], r=d['r'],
                   log2E_Q=d['log2E_Q'], classes=len(res), with_border=int(sum(res)),
                   status='all' if res and all(res) else ('none' if not any(res) else 'some'),
                   secs=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True)
        out.append(rec)
    with open(a.out, 'w') as fo:
        for r in out:
            fo.write(json.dumps(r) + '\n')


if __name__ == '__main__':
    main()
