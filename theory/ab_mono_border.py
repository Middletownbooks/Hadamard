r"""Borders for monomial Q-compressions from the Singer templates [M].

A monomial compression quadruple has z_k = s_k + a_k*delta_{q_k} (a_k = 0 for a constant block) with
sum a_k^2 = 4b; for c = 6 the profile is {4, 2, 0, 0} (report section 3.2).  In the core's GS layout (347 layout) the
compressed core is Cb[(I, g), (J, h)] = Q_L[I, J] * [h = X_IJ(g)] (plus block constants, which do not matter on balanced
frames), with Q_L[I, J] = sign(I, J) * a_{k(I, J)} and X_IJ the spike map of block (I, J):
    sigma = h - g : h = g + q,   g + h : h = q - g,   -g - h : h = -q - g,   g - h : h = g - q   (q = q_{k(I, J)}).
A Singer template (theory/ab_singer.py, Singer coordinates, Cb_S = Q_S (x) I_c) is carried over by signed block
permutations with Q_L[p1 I, p2 J] = s1(I) s2(J) Q_S[I, J] and point labellings lambda^r_I, lambda^c_J with
lambda^c_J = X_IJ o lambda^r_I on every nonzero block (propagated per component of the nonzero-block graph; the
compositions around cycles are consistent for every choice of spike positions, report section 3.2).  The resulting
K, Kc are checked with ab_assemble.forced_corner on the actual z (C must be the Singer corner) and the frame conditions.
usage (library): border_for(z, Qmods) -> dict(K, Kc, C, template) or None"""
import sys, os, glob, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS
from ab_lib import AG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT = AS.load_template347()['layout']
_TPL = {}


def singer_templates(c):
    if c not in _TPL:
        out = []
        fs = sorted(glob.glob(os.path.join(ROOT, f'logs/ab/singer/template_c{c}_cyc_*.npz')))
        fs += sorted(glob.glob(os.path.join(ROOT, f'logs/ab/singer/c{c}_control/template_c{c}_cyc_*.npz')))
        for f in fs:
            if f.endswith('_layout.npz'): continue
            d = np.load(f)
            out.append((os.path.relpath(f, ROOT), d['N'].astype(np.int64), d['Q'].astype(np.int64)))
        _TPL[c] = out
    return _TPL[c]


def spikes(z):
    """per block (base s_k, spike a_k, position q_k or None); None if z is not monomial"""
    out = []
    for zk in np.asarray(z):
        vals, cnt = np.unique(zk, return_counts=True)
        if len(vals) == 1: out.append((int(vals[0]), 0, None)); continue
        if len(vals) != 2 or cnt.min() != 1: return None
        base = int(vals[np.argmax(cnt)]); q = int(np.nonzero(zk != base)[0][0])
        out.append((base, int(zk[q]) - base, q))
    return out


def layout_Q(a):
    Q = np.zeros((4, 4), np.int64); typ = {}
    for key, L in LAYOUT.items():
        I, J = int(key[0]), int(key[1])
        Q[I, J] = L['rho_sign'] * L['gam_sign'] * a[L['k']]
        typ[(I, J)] = (L['sigma'], L['k'])
    return Q, typ


def x_map(G, sigma, q):
    co = G.co; cq = G.co[q]
    if sigma == 'h-g': return G.idx(co + cq)
    if sigma == 'g+h': return G.idx(cq - co)
    if sigma == '-g-h': return G.idx(-cq - co)
    if sigma == 'g-h': return G.idx(co - cq)
    raise ValueError(sigma)


def _signed_perm_match(QL, QS):
    ops = [(p, sg) for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
    for (p1, s1) in ops:
        R = np.array([s1[I] * QL[p1[I]] for I in range(4)])
        for (p2, s2) in ops:
            if all(R[I, p2[J]] == s2[J] * QS[I, J] for I in range(4) for J in range(4)):
                return p1, s1, p2, s2
    return None


def border_for(z, Qmods, first_only=True):
    """K, Kc, C bordering the monomial compression z (4 x c, Q index = ravel of Qmods), or None"""
    G = AG(tuple(Qmods)); c = G.n; b = c - 1; nb = 4 * b
    sp = spikes(z)
    if sp is None: return None
    a = [s[1] for s in sp]; qpos = [s[2] if s[2] is not None else 0 for s in sp]
    if sum(x * x for x in a) != 4 * b: return None
    QL, typ = layout_Q(a)
    X = {}
    for (I, J), (sigma, k) in typ.items():
        if QL[I, J] != 0: X[(I, J)] = x_map(G, sigma, qpos[k])
    for (fname, N, QS) in singer_templates(c):
        m = _signed_perm_match(QL, QS)
        if m is None: continue
        p1, s1, p2, s2 = m
        lam_r = {}; lam_c = {}; bad = False
        while (len(lam_r) < 4 or len(lam_c) < 4) and not bad:
            if len(lam_r) < 4: lam_r[min(set(range(4)) - set(lam_r))] = np.arange(c)
            else: lam_c[min(set(range(4)) - set(lam_c))] = np.arange(c)
            changed = True
            while changed and not bad:
                changed = False
                for (I, J), Xm in X.items():
                    if I in lam_r:
                        want = Xm[lam_r[I]]
                        if J not in lam_c: lam_c[J] = want; changed = True
                        elif not (lam_c[J] == want).all(): bad = True; break
                    elif J in lam_c:
                        lam_r[I] = np.argsort(Xm)[lam_c[J]]; changed = True
        if bad: continue
        K_S = N[:nb, nb:]; Kc_S = N[nb:, :nb].T; C = N[:nb, :nb]
        K_L = np.zeros_like(K_S); Kc_L = np.zeros_like(Kc_S)
        for J in range(4):
            K_L[:, p2[J] * c + lam_c[p2[J]]] = s2[J] * K_S[:, J * c:(J + 1) * c]
            Kc_L[:, p1[J] * c + lam_r[p1[J]]] = s1[J] * Kc_S[:, J * c:(J + 1) * c]
        Cf = AS.forced_corner(np.asarray(z, np.int64), tuple(Qmods), LAYOUT, K_L, Kc_L)
        ok = Cf is not None and bool((Cf == C).all())
        ok &= bool((K_L @ K_L.T == 4 * c * np.eye(nb, dtype=np.int64)).all() and (Kc_L @ Kc_L.T == 4 * c * np.eye(nb, dtype=np.int64)).all())
        ok &= all((K_L[:, J * c:(J + 1) * c].sum(1) == 0).all() and (Kc_L[:, J * c:(J + 1) * c].sum(1) == 0).all() for J in range(4))
        if ok: return dict(K=K_L, Kc=Kc_L, C=C, template=fname, a=a, q=qpos)
    return None


def coverage(c, Qmods):
    """which signed arrangements of the monomial profile have a border (spikes at random positions)"""
    G = AG(tuple(Qmods)); prof = {6: (4, 2, 0, 0), 4: (2, 2, 2, 0)}[c]
    rng = np.random.default_rng(0); res = {}
    for perm in sorted(set(itertools.permutations(prof))):
        for sg in itertools.product((1, -1), repeat=4):
            if any(x == 0 and s < 0 for s, x in zip(sg, perm)): continue
            a = [s * x for s, x in zip(sg, perm)]
            z = np.zeros((4, c), np.int64)
            for k in range(4):
                z[k] += 1
                if a[k]: z[k, rng.integers(c)] += a[k]
            res[tuple(a)] = border_for(z, Qmods) is not None
    return res


if __name__ == "__main__":
    for Qm in ((2, 3), (3, 2), (6,)):
        r = coverage(6, Qm)
        print('c = 6, Q', Qm, ':', sum(r.values()), 'of', len(r), 'signed arrangements bordered;',
              'missing:', [k for k, v in r.items() if not v][:8])
