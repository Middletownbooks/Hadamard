r"""Map a Singer-coordinate template (theory/ab_singer.py) to the GS layout of the cores and verify it there [M].

In Singer coordinates Cb_S = Q_S (x) I_c.  The layout's compressed core for z_k = a_k delta_0 (347 layout, Q = Z_c or
Z_2^m) is Cb_L[(I, g), (J, h)] = Q_L[I, J] X_IJ[g, h] with X = I (sigma = h - g) or rho (g + h, -g - h), rho = negation.
We find signed block permutations with Q_L[pi1(I), pi2(J)] = s1(I) s2(J) Q_S[I, J] (any arrangement a of the profile
whose Q_L is two-sided equivalent to Q_S), point labellings lambda^r_I, lambda^c_J : {inf, Gamma} -> Q with
X_IJ[lambda^r_I(x), lambda^c_J(x')] = [x = x'] (propagated along the nonzero blocks), and map
    K_L[:, (pi2 J, lambda^c(x))] = s2(J) K_S[:, (J, x)],   Kc_L[:, (pi1 J, lambda^r(x))] = s1(J) Kc_S[:, (J, x)].
The result is checked with ab_assemble.forced_corner (C must equal the Singer corner) and the frame conditions.
usage: python3 theory/ab_singer_layout.py TEMPLATE.npz --Q 10 | 2,2,2 [--out FILE]"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT = AS.load_template347()['layout']


def layout_Q(a):
    Q = np.zeros((4, 4), np.int64); X = {}
    for key, L in LAYOUT.items():
        I, J = int(key[0]), int(key[1])
        Q[I, J] = L['rho_sign'] * L['gam_sign'] * a[L['k']]
        X[(I, J)] = 'id' if L['sigma'] == 'h-g' else 'rho'
    return Q, X


def to_layout(N, Qs, Qmods, b):
    c = b + 1; nb = 4 * b
    Kc_S = N[nb:, :nb].T; K_S = N[:nb, nb:]; C = N[:nb, :nb]
    prof = sorted(np.abs(Qs[0]).tolist())
    G = AS.AG(Qmods); assert G.n == c
    neg = G.idx(-G.co)
    ops = [(p, sg) for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
    for perm in sorted(set(itertools.permutations(prof))):
        for sgn in itertools.product((1, -1), repeat=4):
            a = [s * x for s, x in zip(sgn, perm)]
            if any(x == 0 and s < 0 for s, x in zip(sgn, perm)): continue
            QL, X = layout_Q(a)
            for (p1, s1) in ops:
                # rows: QL[p1[I], :] = s1[I] * (Qs[I, :] permuted by p2 with signs s2)
                R = np.array([s1[I] * QL[p1[I]] for I in range(4)])      # R[I, J'] = s1 QL[p1 I, J']
                for (p2, s2) in ops:
                    if all(R[I, p2[J]] == s2[J] * Qs[I, J] for I in range(4) for J in range(4)):
                        return a, (p1, s1, p2, s2), finish(a, QL, X, p1, s1, p2, s2, K_S, Kc_S, C, Qmods, b, neg)
    return None, None, None


def finish(a, QL, X, p1, s1, p2, s2, K_S, Kc_S, C, Qmods, b, neg):
    c = b + 1; nb = 4 * b
    # point labellings: lambda_r[I][x], lambda_c[J][x] for x in 0..c-1 (0 = inf, 1 + t = t in Gamma)
    lam_r = {0: np.arange(c)}; lam_c = {}
    changed = True
    while changed:
        changed = False
        for (I, J), t in X.items():
            if QL[I, J] == 0: continue
            mp = (lambda v: v) if t == 'id' else (lambda v: neg[v])
            if I in lam_r and J not in lam_c: lam_c[J] = mp(lam_r[I]); changed = True
            if J in lam_c and I not in lam_r:
                inv = np.argsort(neg) if t == 'rho' else np.arange(c)
                lam_r[I] = lam_c[J] if t == 'id' else neg[lam_c[J]]; changed = True
    assert len(lam_r) == 4 and len(lam_c) == 4
    K_L = np.zeros_like(K_S); Kc_L = np.zeros_like(Kc_S)
    for J in range(4):
        Jl = p2[J]
        K_L[:, Jl * c + lam_c[Jl]] = s2[J] * K_S[:, J * c:(J + 1) * c]
        Il = p1[J]
        Kc_L[:, Il * c + lam_r[Il]] = s1[J] * Kc_S[:, J * c:(J + 1) * c]
    z = np.zeros((4, c), np.int64)
    for k in range(4): z[k, 0] = a[k]
    Cf = AS.forced_corner(z, Qmods, LAYOUT, K_L, Kc_L)
    ok = Cf is not None and bool((Cf == C).all())
    ok &= bool((K_L @ K_L.T == 4 * c * np.eye(nb, dtype=np.int64)).all() and (Kc_L @ Kc_L.T == 4 * c * np.eye(nb, dtype=np.int64)).all())
    ok &= all((K_L[:, J * c:(J + 1) * c].sum(1) == 0).all() and (Kc_L[:, J * c:(J + 1) * c].sum(1) == 0).all() for J in range(4))
    return dict(ok=ok, K=K_L, Kc=Kc_L, C=C, z=z)


def main():
    a_ = sys.argv[1:]
    f = a_[0]; Qm = tuple(int(x) for x in a_[a_.index('--Q') + 1].split(','))
    d = np.load(f)
    N = d['N'].astype(np.int64); Qs = d['Q'].astype(np.int64)
    b = (N.shape[0] - 4) // 8
    a, ops, res = to_layout(N, Qs, Qm, b)
    out = dict(file=f, Qmods=list(Qm), a=a, found_map=res is not None, verified=bool(res and res['ok']))
    if res and res['ok'] and '--out' in a_:
        o = a_[a_.index('--out') + 1]
        np.savez_compressed(o, K=res['K'].astype(np.int8), Kc=res['Kc'].astype(np.int8), C=res['C'].astype(np.int8),
                            z=res['z'], a=np.array(a), Qmods=np.array(Qm), layout=json.dumps(LAYOUT))
        out['saved'] = o
    print(json.dumps(out))


if __name__ == "__main__":
    main()
