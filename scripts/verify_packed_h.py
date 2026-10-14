"""Exact check of a packed +-1 matrix (np.packbits of H > 0 along rows, key 'H' or 'bits'): H H^T = N I in int64.
Shares no code with theory/.   usage: python3 scripts/verify_packed_h.py FILE.npz [N]"""
import hashlib, json, sys
import numpy as np

def main(path, N=None):
    Z = np.load(path)
    bits = Z['bits'] if 'bits' in Z.files else Z['H']
    N = int(N) if N else bits.shape[0]
    H = np.where(np.unpackbits(bits, axis=1)[:, :N] > 0, 1, -1).astype(np.int64)
    assert H.shape == (N, N)
    ok = True
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        ok &= bool(np.array_equal(G, E))
    print(json.dumps(dict(file=path, N=N, hadamard=ok, sha256_of_packed_rows=hashlib.sha256(np.packbits(H > 0, axis=1).tobytes()).hexdigest())))
    return ok

if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None) else 1)
