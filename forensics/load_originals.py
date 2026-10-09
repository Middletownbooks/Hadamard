"""Load the user-supplied original matrices, verify them, and compare with the normalised hoa64 copies in data/."""
import re, gzip, glob, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def load(path):
    if not os.path.exists(path) and os.path.exists(path + '.gz'):
        lines = gzip.open(path + '.gz', 'rt').read().split('\n')
    else:
        lines = open(path).read().split('\n')
    hdr = lines[0]
    rows = [l.strip() for l in lines[1:] if l.strip()]
    H = np.array([[1 if ch == '+' else -1 for ch in r] for r in rows], dtype=np.int64)
    meta = dict(re.findall(r'(\w+)=(\w+)', hdr))
    return H, meta, hdr
def header(name):
    """header line of an original by base name, from .txt or .txt.gz"""
    f = os.path.join(HERE, 'originals', name + '.txt')
    return (open(f) if os.path.exists(f) else gzip.open(f + '.gz', 'rt')).readline().strip()
def originals():
    """base names of all matrix originals present (.txt or .txt.gz), sorted by order"""
    names = {os.path.basename(p).split('.txt')[0] for p in glob.glob(os.path.join(HERE, 'originals', 'matrix_H*.txt*'))}
    return sorted(names, key=lambda s: int(re.search(r'_H(\d+)_', s).group(1)))
def load_any(name):
    """matrix by base name (e.g. 'matrix_H668_v167_h_b1'): cached .npy if present, else parse .txt / .txt.gz"""
    npy = os.path.join(HERE, 'originals', name + '.npy')
    if os.path.exists(npy): return np.load(npy)
    H, _, _ = load(os.path.join(HERE, 'originals', name + '.txt')); np.save(npy, H); return H
if __name__ == '__main__':
    for f in sorted(sorted(set(p[:-3] if p.endswith('.gz') else p for p in glob.glob(os.path.join(HERE, 'originals', 'matrix_H*.txt*')))), key=lambda s: int(re.search(r'_H(\d+)_', s).group(1))):
        H, meta, hdr = load(f)
        N = H.shape[0]
        ok = H.shape == (N, N) and np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64))
        normalised = bool((H[0] == 1).all() and (H[:, 0] == 1).all())
        np.save(f.replace('.txt', '.npy'), H)
        hoa = np.loadtxt(gzip.open(os.path.join(HERE, '..', 'data', f'hadamard_{N}.csv.gz')), delimiter=',', dtype=np.int64)
        Hn = H * H[:, :1]; Hn = Hn * Hn[:1, :]
        print(f'{os.path.basename(f)}: "{hdr}"  HH^T=NI: {ok}  normalised: {normalised}  '
              f'normalising it gives the hoa64 copy: {np.array_equal(Hn, hoa)}')
    ann = os.path.join(HERE, 'H668_announced.npy')
    if os.path.exists(ann):
        print('668 original identical to the formal-conjectures announced matrix:',
              np.array_equal(load_any('matrix_H668_v167_h_b1'), np.load(ann)))
    else:
        print('(H668_announced.npy not present: run decode668.py on the formal-conjectures Lean file to compare)')
