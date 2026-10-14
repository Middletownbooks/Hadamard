r"""Border for a c = 4 (Q = Z_2^2) core, then the full H(4v) [M].

Given a core X (4 x m on G = Q x H, Q = the first two axes), the Q-compressions z and their GS array Cb (layout of the
479 template, untwisted), search frames: for a candidate Kc (479's own, and the Sylvester frames H4 (x) M with M the
three balanced rows of H4, under the 4!^4 block-internal column orders reduced to the 8 automorphisms of Z_2^2 per
block), enumerate the admissible rows {k block-balanced : -k Cb^T Kc^T/16 in {+-1}^12} (6^4 candidates) and look
for 12 mutually orthogonal ones (K).  Then C = -K Cb^T Kc^T/16, assemble with theory/ab_assemble.assemble and check
H H^T = 4v I.
usage (library): border_c4(X, mods) -> (H, K, Kc, C) or None"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS
from ab_border10 import admissible, clique

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def template1916():
    T = np.load(os.path.join(ROOT, 'data/template1916_ab.npz'), allow_pickle=True)
    return dict(K=T['K'].astype(np.int64), Kc=T['Kc'].astype(np.int64), layout=json.loads(str(T['layout'])))


def _autos():
    Qel = [(0, 0), (0, 1), (1, 0), (1, 1)]
    autos = []
    for A in itertools.product((0, 1), repeat=4):
        a, b, c, d = A
        if (a * d - b * c) % 2 == 0: continue
        for t in Qel:
            autos.append([Qel.index(((a * x + b * y + t[0]) % 2, (c * x + d * y + t[1]) % 2)) for x, y in Qel])
    return autos                                   # AGL(2, 2): 24 relabelings of Q = Z_2^2


def sylvester_frames():
    """H4 (x) M (M = the three balanced rows of H4) with Q relabelled inside blocks 1-3 by AGL(2,2) (24^3 frames)"""
    H4 = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
    base = np.kron(H4, H4[1:])
    autos = _autos()
    for p1, p2, p3 in itertools.product(range(len(autos)), repeat=3):
        cols = list(range(4))
        for J, pi in ((1, p1), (2, p2), (3, p3)): cols += [J * 4 + autos[pi][i] for i in range(4)]
        yield base[:, cols]


def border_c4(X, mods, nQ=2):
    layout = AS.load_template347()['layout']        # untwisted GS layout (valid on any abelian Q)
    Z = AS.compressions(X, mods, nQ)
    Cb = AS.compressed_core(Z, tuple(mods[:nQ]), layout)
    for Kc in sylvester_frames():
        S = admissible(Cb, Kc, 4)
        if len(S) < 24: continue
        K, nS, deg = clique(S, 12, tmax=20)
        if K is None: continue
        num = -(K @ Cb.T @ Kc.T)
        if (num % 16).any(): continue
        C = num // 16
        if not (np.abs(C) == 1).all(): continue
        H = AS.assemble(X, tuple(mods), nQ, layout, K, Kc, C)
        if AS.is_hadamard(H): return H, K, Kc, C
    return None


if __name__ == "__main__":
    # self-test on the 479 template's core, untwisted by the character (-1)^{q_2} (the template stores a twisted form)
    from ab_lib import AG
    Tz = np.load(os.path.join(ROOT, 'data/template1916_ab.npz'), allow_pickle=True)
    F = Tz['F'].astype(np.int64); mods = tuple(int(x) for x in Tz['mods'])
    G = AG(mods); F = F * ((-1) ** G.co[:, 1])[None, :]
    r = border_c4(F, mods)
    print(json.dumps(dict(selftest_479=r is not None, order=None if r is None else int(r[0].shape[0]))))
