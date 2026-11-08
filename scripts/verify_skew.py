"""Independent verification of a skew Goethals-Seidel certificate (pattern of scripts/verify_h2092.py; shares no code
with theory/ or with Kumar's search code).

Certificate (JSON): {"v": v, "mods": [n1, ..., nk], "gens": [M1, M2, ...], "A": "+-...", "B": ..., "C": ..., "D": ...}
- G = Z_n1 x ... x Z_nk, element index = mixed radix of its coordinates (first coordinate most significant);
- each M is a k x k integer matrix acting on coordinate column vectors mod (n1, ..., nk) (an automorphism of G);
- A, B, C, D are strings of v signs, position = element index.
Checks, all in exact integer arithmetic:
 1. sum_X PAF_X(s) = 4v [s = 0] at every s in G, by the definition (sum over g of X(g) X(g + s));
 2. A(0) = +1 and A(-g) = -A(g) for g != 0;
 3. every block is invariant under every generator: X(M g) = X(g);
 4. the Goethals-Seidel array H with group-developed blocks X_{g,h} = X(h - g) and R_{g,h} = [h = -g]
      [[A, BR, CR, DR], [-BR, A, D^T R, -C^T R], [-CR, -D^T R, A, B^T R], [-DR, C^T R, -B^T R, A]]
    has entries +-1, H H^T = 4v I (float64 BLAS on +-1 entries: every entry of the product is an integer of size
    <= 4v < 2^53, so the product is exact) and H + H^T = 2I.
Writes the SHA-256 of the packed rows of H and optionally H as packed bits.
usage: python3 scripts/verify_skew.py CERT.json [OUT.npz]"""
import hashlib, json, sys
import numpy as np

def main(path, out=None):
    c = json.load(open(path))
    v, mods = int(c['v']), [int(n) for n in c['mods']]
    assert int(np.prod(mods)) == v and v % 2 == 1
    k = len(mods)
    coords = np.array(np.unravel_index(np.arange(v), mods)).T            # v x k, mixed radix (C order)
    modv = np.array(mods)
    idx = lambda C: np.ravel_multi_index(tuple((np.asarray(C) % modv).T), mods)
    X = {}
    for name in 'ABCD':
        s = c[name].strip(); assert len(s) == v and set(s) <= set('+-')
        X[name] = np.array([1 if ch == '+' else -1 for ch in s], dtype=np.int64)
    # 1. PAF by definition
    add = np.array([idx(coords + coords[s]) for s in range(v)])            # add[s, g] = index of g + s
    tot = np.zeros(v, dtype=np.int64)
    for name in 'ABCD':
        x = X[name]; tot += (x[None, :] * x[add]).sum(1)
    paf_ok = bool(tot[0] == 4 * v and np.all(tot[1:] == 0))
    # 2. skew first block
    neg = idx(-coords)
    A = X['A']
    skew_ok = bool(A[0] == 1 and all(A[neg[g]] == -A[g] for g in range(1, v)))
    # 3. invariance
    inv_ok = True
    for M in c['gens']:
        M = np.array(M, dtype=np.int64).reshape(k, k)
        img = idx(coords @ M.T)
        assert len(set(img.tolist())) == v, 'generator is not a bijection'
        for name in 'ABCD':
            inv_ok &= bool(np.all(X[name][img] == X[name]))
    rows = {n: int(X[n].sum()) for n in 'ABCD'}
    print(json.dumps(dict(cert=path, v=v, mods=mods, paf_ok=paf_ok, skew_first_block=skew_ok, invariant=inv_ok,
                          rowsums=rows)), flush=True)
    if not (paf_ok and skew_ok and inv_ok): return False
    # 4. the matrix
    sub = np.array([idx(coords - coords[g]) for g in range(v)])            # sub[g, h] = index of h - g
    dev = {n: X[n][sub] for n in 'ABCD'}                                    # X_{g,h} = X(h - g)
    R = np.zeros((v, v), dtype=np.int64); R[np.arange(v), neg] = 1
    Am, Bm, Cm, Dm = (dev[n] for n in 'ABCD')
    H = np.block([[Am, Bm @ R, Cm @ R, Dm @ R],
                  [-Bm @ R, Am, Dm.T @ R, -Cm.T @ R],
                  [-Cm @ R, -Dm.T @ R, Am, Bm.T @ R],
                  [-Dm @ R, Cm.T @ R, -Bm.T @ R, Am]]).astype(np.int8)
    N = 4 * v
    assert set(np.unique(H).tolist()) == {-1, 1}
    Hf = H.astype(np.float64)
    hh = bool(np.array_equal(Hf @ Hf.T, N * np.eye(N)))
    sk = bool(np.array_equal(H.astype(np.int64) + H.T.astype(np.int64), 2 * np.eye(N, dtype=np.int64)))
    bits = np.packbits(H > 0, axis=1)
    res = dict(N=N, HHt_eq_NI=hh, H_plus_Ht_eq_2I=sk, sha256_of_packed_rows=hashlib.sha256(bits.tobytes()).hexdigest())
    if out: np.savez_compressed(out, bits=bits, N=N); res['saved'] = out
    print(json.dumps(res), flush=True)
    return hh and sk

if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None) else 1)
