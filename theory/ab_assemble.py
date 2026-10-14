r"""Assemble a full b >= 3 Hadamard matrix from a core and an inflated border [M].

Layout (report/bge3_construction.md, data/template347.npz): G = Q x H (split, Q = the first nQ axes), m = |G|,
c = |Q|, off = 4b.  The 4v x 4v matrix is
    [[ C               , K inflated by h ],
     [ Kc^T inflated   , core GS array   ]]
with  H[p, off + J m + g] = K[p, J c + q(g)],  H[off + J m + g, p] = Kc[p, J c + q(g)],  H[:off, :off] = C,
and core sub-block (I, J)[g, h] = sign_IJ * F_k(sigma_IJ(g, h)), sigma in {h-g, g+h, g-h, -g-h} (layout of the
published 347 matrix; the same layout fits 479 and 359 by bge3_construction.md).  The corner is forced:
C = -K Cb^T Kc^T / (4c), Cb = the compressed GS array of the Q-compressions (scripts/bge3_border_check.py);
it depends only on the deviations of the compressions from their block means.
usage: python3 theory/ab_assemble.py --check-template      (rebuild H(1388) from data/template347.npz)"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import AG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_template347():
    T = np.load(os.path.join(ROOT, 'data/template347.npz'), allow_pickle=True)
    return dict(C=T['C'].astype(np.int64), K=T['K'].astype(np.int64), Kc=T['Kc'].astype(np.int64),
                comp=T['comp'].astype(np.int64), layout=json.loads(str(T['layout'])), F=T['F'].astype(np.int64))


def sigma_index(G, name):
    co = G.co
    if name == 'h-g': return G.idx(co[None, :, :] - co[:, None, :])
    if name == 'g+h': return G.idx(co[:, None, :] + co[None, :, :])
    if name == 'g-h': return G.idx(co[:, None, :] - co[None, :, :])
    if name == '-g-h': return G.idx(-co[:, None, :] - co[None, :, :])
    raise ValueError(name)


def compressed_core(xbar, Qmods, layout, qg=None):
    """4c x 4c GS array of the Q-compressions (same layout as the core), as in scripts/bge3_border_check.py;
    coset-constant character twists rho, gam (layout 'arr') are included (qg: coset index of each element of G)"""
    Q = AG(Qmods); c = Q.n
    M = np.zeros((4 * c, 4 * c), np.int64)
    for key, L in layout.items():
        I_, J_ = int(key[0]), int(key[1]); sg = L['rho_sign'] * L['gam_sign']
        blk = sg * np.asarray(xbar[L['k']])[sigma_index(Q, L['sigma'])]
        if L.get('rho') == 'arr':
            rho = np.array(L['rho_arr']); gam = np.array(L['gam_arr'])
            rq = np.array([rho[qg == q][0] for q in range(c)]); gq = np.array([gam[qg == q][0] for q in range(c)])
            blk = blk * np.outer(rq, gq)
        M[I_ * c:(I_ + 1) * c, J_ * c:(J_ + 1) * c] = blk
    return M


def forced_corner(xbar, Qmods, layout, K, Kc, qg=None):
    c = AG(Qmods).n
    Cb = compressed_core(xbar, Qmods, layout, qg)
    num = -(K @ Cb.T @ Kc.T)
    if (num % (4 * c)).any(): return None
    C = num // (4 * c)
    return C if (np.abs(C) == 1).all() else None


def assemble(F, mods, nQ, layout, K, Kc, C, qg=None, Qmods=None):
    """qg/Qmods: coset index of every element of G and the quotient Q (for non-split G); default Q = first nQ axes"""
    G = AG(mods); m = G.n
    if qg is None:
        Qm = AG(mods[:nQ]); qg = Qm.idx(G.co[:, :nQ])
    else:
        Qm = AG(Qmods)
    c = Qm.n
    off = K.shape[0]; N = off + 4 * m
    H = np.zeros((N, N), np.int8)
    H[:off, :off] = C
    for J in range(4):
        H[:off, off + J * m: off + (J + 1) * m] = K[:, J * c + qg]
        H[off + J * m: off + (J + 1) * m, :off] = Kc[:, J * c + qg].T
    for key, L in layout.items():
        I_, J_ = int(key[0]), int(key[1]); sg = L['rho_sign'] * L['gam_sign']
        blk = sg * F[L['k']][sigma_index(G, L['sigma'])]
        if L.get('rho') == 'arr': blk = blk * np.outer(np.array(L['rho_arr']), np.array(L['gam_arr']))
        H[off + I_ * m: off + (I_ + 1) * m, off + J_ * m: off + (J_ + 1) * m] = blk
    return H


def is_hadamard(H):
    Hi = H.astype(np.int64)
    N = Hi.shape[0]
    return bool(np.array_equal(Hi @ Hi.T, N * np.eye(N, dtype=np.int64)))


def compressions(F, mods, nQ, qg=None, Qmods=None):
    G = AG(mods)
    if qg is None:
        Qm = AG(mods[:nQ]); qg = Qm.idx(G.co[:, :nQ])
    else:
        Qm = AG(Qmods)
    Z = np.zeros((4, Qm.n), np.int64)
    for k in range(4): np.add.at(Z[k], qg, F[k])
    return Z


if __name__ == "__main__":
    if '--check-template' in sys.argv:
        T = load_template347()
        H = assemble(T['F'], (2, 3, 57), 2, T['layout'], T['K'], T['Kc'], T['C'])
        Z = compressions(T['F'], (2, 3, 57), 2)
        Cf = forced_corner(Z, (2, 3), T['layout'], T['K'], T['Kc'])
        print(json.dumps(dict(order=int(H.shape[0]), hadamard=is_hadamard(H), comp_matches=bool((Z == T['comp']).all()),
                              forced_corner_equals_C=bool(Cf is not None and (Cf == T['C']).all()))))
