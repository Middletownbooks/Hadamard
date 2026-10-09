"""Infer the b=0 'g' builder layout from an original matrix and rebuild matrices from record files (seq_v*_g.txt).
Each block (p,q) is s * T(f_k(e*y + c)) with T circulant (f(j-i)) or back-circulant (f(i+j)), f_k one of the four
row-0 sequences of the first block row (A circulant; B,C,D back-circulant as stored), e = +-1, c a shift."""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def circ(f):
    m = len(f); I = np.arange(m); return f[(I[None, :] - I[:, None]) % m]
def back(f):
    m = len(f); I = np.arange(m); return f[(I[None, :] + I[:, None]) % m]
def infer_layout(H, m):
    I = np.arange(m)
    seqs = [H[0, q * m:(q + 1) * m] for q in range(4)]
    lay = {}
    for p in range(4):
        for q in range(4):
            B = H[p * m:(p + 1) * m, q * m:(q + 1) * m]; r0 = B[0]
            typ = 'circ' if np.array_equal(B, circ(r0)) else ('back' if np.array_equal(B, back(r0)) else None)
            hit = None
            for k in range(4):
                for e in (1, -1):
                    for c in range(m):
                        g = seqs[k][(e * I + c) % m]
                        if np.array_equal(r0, g): hit = (k, typ, e, c, 1)
                        elif np.array_equal(r0, -g): hit = (k, typ, e, c, -1)
                        if hit: break
                    if hit: break
                if hit: break
            lay[(p, q)] = hit
    return lay
def build(seqs, lay, m):
    I = np.arange(m); H = np.zeros((4 * m, 4 * m), dtype=np.int64)
    for (p, q), (k, typ, e, c, s) in lay.items():
        g = seqs[k][(e * I + c) % m]
        H[p * m:(p + 1) * m, q * m:(q + 1) * m] = s * (circ(g) if typ == 'circ' else back(g))
    return H
def read_seq(path):
    lines = [l.strip() for l in open(path) if l.strip()]
    return lines[0], [np.array([1 if c == '+' else -1 for c in l.split('=', 1)[1]]) for l in lines[1:5]]
if __name__ == '__main__':
    from load_originals import load_any
    H = load_any('matrix_H1244_v311_g_b0'); m = 311
    lay = infer_layout(H, m)
    print('g-builder layout from H(1244): block (p,q) = sign * type(f_k(e*y + c)) -> (k, type, e, c, sign)')
    for p in range(4): print('  ', [lay[(p, q)] for q in range(4)])
    assert np.array_equal(build([H[0, q*m:(q+1)*m] for q in range(4)], lay, m), H)
    H2 = load_any('matrix_H1964_v491_g_b0'); m2 = 491
    lay2 = infer_layout(H2, m2)
    norm = lambda L, mm: {k: (a, b, e, (c if c < mm // 2 else c - mm), s) for k, (a, b, e, c, s) in L.items()}
    print('H(1964) has the same layout (shifts compared as signed residues):', norm(lay2, m2) == norm(lay, m))
    for f in ['seq_v223_g.txt', 'seq_v283_g.txt']:
        hdr, seqs = read_seq(os.path.join(HERE, 'originals', f)); mm = len(seqs[0])
        lay_m = {k: (a, b, e, c if c < m // 2 else mm - (m - c), s) for k, (a, b, e, c, s) in lay.items()}
        Hb = build(seqs, lay_m, mm); N = 4 * mm
        ok = np.array_equal(Hb @ Hb.T, N * np.eye(N, dtype=np.int64))
        np.save(os.path.join(HERE, 'originals', f.replace('.txt', '_built.npy')), Hb)
        print(f'{f} ("{hdr}"): rebuilt with the g layout, HH^T = NI: {ok}')
