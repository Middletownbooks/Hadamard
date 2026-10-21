r"""Control: Kumar's 2148 core (Q = Z_3, T = all characters; hadamard/2148/certificate.json in
https://github.com/vatskishan03/hadamard, Apache-2.0) through our engine and border code [M]: the core condition on
G = Z_3 x Z_178, its Q-compressions, a border from our Paley-H12 frames (theory/bt_border.border_for_compressions),
the assembled matrix, H H^T = 2148 I.
usage: python3 theory/bt_control2148.py /path/to/vatskishan03/hadamard/hadamard/2148/certificate.json
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_border as BB

cert = json.load(open(sys.argv[1]))
mods = tuple(cert['group'])
st = E.BStruct(mods, [tuple(g) for g in cert['subgroup_generators']], 'all',
               [np.array(M) for M in cert['symmetries']], name='2148')
xs = [np.array([1 if ch == '+' else -1 for ch in s], np.int64) for s in cert['sequences']]
out = dict(control='2148', mods=list(mods), q=st.q, t=st.t, r=st.r, e=st.e)
out['invariant'] = all(bool((x[st.bperms[0][0]] == x).all()) for x in xs)
out['core'] = st.check_core(xs)
out['rows'] = [int(x.sum()) for x in xs]
z = BB.compressions(st, xs)
out['compressions'] = z.tolist()
res = BB.border_for_compressions(st, z, max_frames=20000, tmax=5.0)
out['border_found'] = res is not None
if res is not None:
    X, Rh, Ch, k = res
    H = BB.assemble(st, xs, X, Rh, Ch)
    out['frame_index'] = int(k)
    out['order'] = int(H.shape[0])
    out['hadamard'] = BB.is_hadamard(H)
print(json.dumps(out))
