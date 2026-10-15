r"""Second construction of H(3668) [M]: the Kronecker-form Cooper-Wallis array of scripts/cw_build.build (the code that
built H(2884)), fed with the T-matrices of order 131 (logs/b917/T131_k5.npy) and the order-7 Williamson quadruple of
theory/b917_cw.williamson(7).  The blocks are block-circulant over Z_131 x Z_7 (Kronecker products of circulants), the
array uses the reversal g -> -g on Z_131 x Z_7; H H^T = 3668 I is checked exactly in int64 row blocks.  This matrix is
built from the T-matrices directly, without the CRT map to Z_917 and without the Goethals-Seidel array of
scripts/verify_h3668.py, so it is an independent check of the construction.
usage: python3 theory/b917_cw_kron.py [--T logs/b917/T131_k5.npy] [--save data/H3668_cwkron_packed.npz]"""
import sys, os, json, hashlib
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cw_build import build, check_T
from b917_cw import williamson


def exact_hadamard(H):
    N = H.shape[0]; H = H.astype(np.int64)
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True


def main():
    a = sys.argv[1:]
    opt = lambda n, d=None: a[a.index(n) + 1] if n in a else d
    T = np.load(os.path.join(ROOT, opt('--T', 'logs/b917/T131_k5.npy'))).astype(np.int64)
    check_T(T)                                            # cw_build's own check (integer PAF)
    W = williamson(7)
    H = build(T, W).astype(np.int8)
    assert set(np.unique(H).tolist()) == {-1, 1}
    ok = exact_hadamard(H)
    rec = dict(order=int(H.shape[0]), construction='Cooper-Wallis Kronecker array (scripts/cw_build.build)', hadamard=ok,
               sha256_of_packed_rows=hashlib.sha256(np.packbits(H > 0, axis=1).tobytes()).hexdigest())
    if ok and opt('--save'):
        np.savez_compressed(os.path.join(ROOT, opt('--save')), bits=np.packbits(H > 0, axis=1), N=H.shape[0])
        rec['saved'] = opt('--save')
    print(json.dumps(rec), flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
