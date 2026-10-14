r"""Does a published frame (K, Kc) border a given Q-compression?  [M]
For a stage-1 solution y (on G' = Q x H'), its Q-compression z_k (k = 0..3) is tested against a template frame:
for every block order and per-block sign, the forced corner C = -K Cb(z)^T Kc^T / (4c) must be +-1 (then C C^T = 4b I
automatically and the full matrix is Hadamard once the core is completed; theory/ab_assemble.py).
Templates: data/template{1388,1916,1436}_ab.npz (347: c = 6, Q = Z_2 x Z_3; 479: c = 4, Z_2^2; 359: c = 8, Z_2^3).
usage: (library) fits(z, N) -> list of (perm, signs)"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_T = {}

def template(N):
    if N not in _T:
        T = np.load(os.path.join(ROOT, f'data/template{N}_ab.npz'), allow_pickle=True)
        _T[N] = dict(K=T['K'], Kc=T['Kc'], C=T['C'], comp=T['comp'], layout=json.loads(str(T['layout'])),
                     Qmods=tuple(T['Qmods'].tolist()), qg=T['qg'])
    return _T[N]

def fits(z, N, first_only=True):
    """z: 4 x c compressions in the template's Q indexing; returns the (perm, signs) giving a +-1 corner"""
    T = template(N); out = []
    for perm in itertools.permutations(range(4)):
        for sg in itertools.product((1, -1), repeat=4):
            Z = np.array([sg[i] * np.asarray(z[perm[i]]) for i in range(4)])
            Cf = AS.forced_corner(Z, T['Qmods'], T['layout'], T['K'], T['Kc'], T['qg'])
            if Cf is not None:
                out.append((perm, sg))
                if first_only: return out
    return out
