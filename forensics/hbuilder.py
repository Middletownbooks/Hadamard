"""The b=1 'h' builder (bordered Goethals-Seidel), inferred from the authors' H(668) and checked on H(716), H(1772);
and the exact rederivation of an original from its normalised copy.

Layout: 4 border rows/columns, then 4x4 core blocks of size m = v-1.  Block (p,q) = sign * T(f_k(e*y + c)) with T
circulant f(j-i) or back-circulant f(i+j) and f_k one of the four row-0 sequences of the first block row (the records).
The border is record-independent: a fixed 4x4 corner, and border rows/columns constant (+-1) on each core block.
Since normalisation is hoa = diag(H[:,0]) H diag(H[0,0] H[0,:]) and H's first row and column are part of the fixed
border, the original is recovered exactly from a normalised copy:  H = diag(col0) hoa diag(H[0,0] row0)."""
import os, sys, gzip
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gbuilder import infer_layout as infer_core_layout, build as build_core
from load_originals import load_any
def border_pattern(H, m):
    """corner (4x4) and per-block signs of the 4 border rows (over core column-blocks) and columns (over core row-blocks)"""
    corner = H[:4, :4].copy()
    rows = np.array([[H[i, 4 + q * m] for q in range(4)] for i in range(4)])
    cols = np.array([[H[4 + p * m, i] for p in range(4)] for i in range(4)])
    for i in range(4):
        for q in range(4):
            assert (H[i, 4 + q * m: 4 + (q + 1) * m] == rows[i, q]).all() and (H[4 + q * m: 4 + (q + 1) * m, i] == cols[i, q]).all()
    return corner, rows, cols
def build(seqs, layout, border, m):
    corner, rows, cols = border
    N = 4 + 4 * m; H = np.zeros((N, N), dtype=np.int64)
    H[:4, :4] = corner
    for i in range(4):
        for q in range(4):
            H[i, 4 + q * m: 4 + (q + 1) * m] = rows[i, q]; H[4 + q * m: 4 + (q + 1) * m, i] = cols[i, q]
    H[4:, 4:] = build_core(seqs, layout, m)
    return H
def records(H, m):
    return [H[4, 4 + q * m: 4 + (q + 1) * m].copy() for q in range(4)]
def scaled_layout(layout, m_from, m_to):
    return {k: (a, t, e, c if c < m_from // 2 else m_to - (m_from - c), s) for k, (a, t, e, c, s) in layout.items()}
def denormalise(hoa, border, m):
    """undo hoa = diag(H[:,0]) H diag(H[0,0] H[0,:]) using the builder's fixed first row and column"""
    corner, rows, cols = border
    row0 = np.concatenate([corner[0], np.repeat(rows[0], m)]); col0 = np.concatenate([corner[:, 0], np.repeat(cols[0], m)])
    return col0[:, None] * hoa * (corner[0, 0] * row0)[None, :]
def load_hoa(N):
    return np.loadtxt(gzip.open(os.path.join(HERE, '..', 'data', f'hadamard_{N}.csv.gz')), delimiter=',', dtype=np.int64)
H668 = load_any('matrix_H668_v167_h_b1'); M668 = 166
LAYOUT = {k: v for k, v in infer_core_layout(H668[4:, 4:], M668).items()}
BORDER = border_pattern(H668, M668)
def rederive(N):
    v = N // 4; m = v - 1
    hoa = load_hoa(N)
    H = denormalise(hoa, BORDER, m)
    lay = scaled_layout(LAYOUT, M668, m)
    rec = records(H, m)
    exact = np.array_equal(build(rec, lay, BORDER, m), H)
    return H, rec, exact
if __name__ == '__main__':
    print('h layout from H(668): block (p,q) = sign * type(f_k(e*y + c)) -> (k, type, e, c, sign)')
    for p in range(4): print('  ', [LAYOUT[(p, q)] for q in range(4)])
    for name, v in [('matrix_H668_v167_h_b1', 167), ('matrix_H716_v179_h_b1', 179), ('matrix_H1772_v443_h_b1', 443)]:
        H = load_any(name); m = v - 1
        same = np.array_equal(build(records(H, m), scaled_layout(LAYOUT, M668, m), BORDER, m), H)
        print(f'{name}: rebuilt exactly from its own records with the h layout + fixed border: {same}')
    print('rederivation from the normalised hoa64 copy (validation where the original is known):')
    for N, name in [(668, 'matrix_H668_v167_h_b1'), (716, 'matrix_H716_v179_h_b1'), (1772, 'matrix_H1772_v443_h_b1')]:
        H, rec, exact = rederive(N)
        print(f'   N={N}: result is an exact h-builder matrix: {exact};  equals the authors\' original: {np.array_equal(H, load_any(name))}')
    H, rec, exact = rederive(1676)
    ok = np.array_equal(H @ H.T, 1676 * np.eye(1676, dtype=np.int64))
    print(f'   N=1676: result is an exact h-builder matrix: {exact};  HH^T = NI: {ok}')
    if exact:
        with open(os.path.join(HERE, 'originals', 'seq_v419_h_rederived.txt'), 'w') as fh:
            fh.write('v=419 N=1676 func=h b=1 m=418 (rederived from data/hadamard_1676.csv.gz by forensics/hbuilder.py)\n')
            for k, x in zip('ABCD', rec): fh.write(k + '=' + ''.join('+' if t > 0 else '-' for t in x) + '\n')
        print('   wrote originals/seq_v419_h_rederived.txt')
