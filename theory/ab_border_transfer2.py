r"""Border transfer from the published templates (347: c = 6, 479: c = 4, 359: c = 8) to every open order [P, M].

[P] The forced corner C = -K Cb^T Kc^T / (4c) is linear in the compression; for each template the kernel of
z -> K Cb(z)^T Kc^T, block by block, is one-dimensional, spanned by w (w = 1 for 347, a real character of Q for 479
and 359: the templates' F_k are the natural cores twisted by that character).  Moving the compression along w keeps
the corner: z_k -> z_k + s_k w (s_k even keeps the coset-sum parity, odd flips it, matching h' = h + 1 mod 2), and the
natural row sum a_k = <z_k, w> moves to a_k + c s_k.  So the template's border fits order v' (same c) whenever
sum_k (a_k + c s_k)^2 = 4v' has a solution with every s_k = h' - h (mod 2).  Each solution is checked exactly here:
the forced corner of the moved compression is +-1 and C C^T = 4b I.
usage: python3 theory/ab_border_transfer2.py"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS
from ab_lib import AG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN = [537, 571, 599, 643, 659, 669, 719, 739, 751, 789, 823, 839, 859, 863, 883, 907, 917, 919, 933, 947, 955, 971, 991]

def load(N):
    T = np.load(os.path.join(ROOT, f'data/template{N}_ab.npz'), allow_pickle=True)
    return dict(K=T['K'], Kc=T['Kc'], C=T['C'], comp=T['comp'], layout=json.loads(str(T['layout'])),
                Qmods=tuple(T['Qmods'].tolist()), qg=T['qg'], mods=tuple(T['mods'].tolist()))

def kernel_dir(T):
    import sympy
    c = AG(T['Qmods']).n; out = []
    for k in range(4):
        cols = []
        for qq in range(c):
            D = np.zeros((4, c), np.int64); D[k, qq] = 1
            Cb = AS.compressed_core(D, T['Qmods'], T['layout'], T['qg'])
            cols.append((T['K'] @ Cb.T @ T['Kc'].T).ravel())
        ns = sympy.Matrix(np.array(cols).T.tolist()).nullspace()
        assert len(ns) == 1
        w = np.array([int(x) for x in ns[0]], np.int64); w = w // np.gcd.reduce(w)
        out.append(w)
    return out

if __name__ == "__main__":
    res = []
    for N, c, hmod in ((1916, 4, 119), (1388, 6, 57), (1436, 8, 44)):
        T = load(N); W = kernel_dir(T); b = c - 1
        a = [int(T['comp'][k] @ W[k]) for k in range(4)]
        for v in OPEN:
            if (v + 1) % c: continue
            h = (v + 1) // c - 1; par = (h - hmod) % 2
            sols = []
            rng = range(-12, 13)
            for s in itertools.product(rng, repeat=4):
                if any((x - par) % 2 for x in s): continue
                r = [ak + c * sk for ak, sk in zip(a, s)]
                if sum(x * x for x in r) == 4 * v: sols.append((s, r))
            ok = 0
            for (s, r) in sols[:50]:
                Z = T['comp'] + np.array([s[k] * W[k] for k in range(4)])
                Cf = AS.forced_corner(Z, T['Qmods'], T['layout'], T['K'], T['Kc'], T['qg'])
                ok += int(Cf is not None and (Cf @ Cf.T == 4 * b * np.eye(4 * b, dtype=np.int64)).all())
            rec = dict(template=N, c=c, v=v, h=h, row_solutions=len(sols), checked=min(50, len(sols)), corner_ok=ok,
                       example_rows=sols[0][1] if sols else None)
            print(json.dumps(rec), flush=True); res.append(rec)
    json.dump(res, open(os.path.join(ROOT, 'logs/ab/border_transfer.json'), 'w'), indent=1)
