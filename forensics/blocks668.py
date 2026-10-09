"""Identify the 16 core blocks of the announced H(668): circulant f(j-i), back-circulant f(i+j), with optional
row/column sign vectors; recover sequences and compare with the repo core (data/core_v167_b1.txt) and the Lean variant."""
import numpy as np, sys
sys.path.insert(0, '/home/user/Hadamard/scripts')
H = np.load('H668_announced.npy')
n = 166
M = H[4:, 4:]
def blk(p, q): return M[p*n:(p+1)*n, q*n:(q+1)*n]
I = np.arange(n)
def classify(B):
    out = []
    for kind in ('circ', 'back'):
        idx = (I[None, :] - I[:, None]) % n if kind == 'circ' else (I[None, :] + I[:, None]) % n
        f = np.zeros(n, dtype=int); cnt = np.zeros(n, dtype=int)
        # direct test (no sign twists): constant along the diagonals
        vals = {}
        consistent = True
        for i in range(n):
            for j in range(n):
                k = idx[i, j]
                if k in vals and vals[k] != B[i, j]:
                    consistent = False; break
                vals[k] = B[i, j]
            if not consistent: break
        if consistent:
            out.append((kind, np.array([vals[k] for k in range(n)])))
    return out
seqs = {}
for p in range(4):
    row = []
    for q in range(4):
        c = classify(blk(p, q))
        row.append(c[0][0] if c else 'other')
        if c: seqs[(p, q)] = c[0][1]
    print('block row', p, row)
np.save('seqs668.npy', np.array([seqs[(0, q)] for q in range(4)]))
# compare first block-row sequences with repo core (up to shift, reversal, negation)
from hadlib import Core
core = Core('/home/user/Hadamard/data/core_v167_b1.txt')
def equiv(x, y):
    for s in range(n):
        for rev in (False, True):
            z = np.roll(y[::-1] if rev else y, s)
            if np.array_equal(x, z): return f'shift {s}{" reversed" if rev else ""}'
            if np.array_equal(x, -z): return f'-shift {s}{" reversed" if rev else ""}'
    return None
for q in range(4):
    x = seqs[(0, q)]
    m = [(k, equiv(x, core.X[k])) for k in range(4)]
    print(f'announced block (0,{q}) sum={x.sum():+d}: matches repo core blocks:', [(k, e) for k, e in m if e])
