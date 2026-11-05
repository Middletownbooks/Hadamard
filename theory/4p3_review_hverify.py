#!/usr/bin/env python3
r"""[M] Standalone Hadamard verifier for the review of report/openai_4p3_answer.md (Claude Code).  Imports nothing from
the repository: no construction, border, assembly or OpenAI code.

Input: a .npy file holding a square matrix, or a JSON file with key "H".  Checks:
  1. square, every entry +1 or -1;
  2. H H^T = N I by exact integer arithmetic: rows bit-packed (bit = 1 for -1), <r_i, r_j> = N - 2 popcount(r_i xor r_j),
     every pair i <= j;
  3. the same Gram matrix by float64 BLAS (exact here: every partial sum is an integer of absolute value <= N < 2^53);
  4. a negative control: one entry flipped at a seeded random position must be rejected by check 2.
usage: python3 theory/4p3_review_hverify.py H.npy [--out result.json] [--seed 1]
"""
import argparse, hashlib, json, sys
import numpy as np


def load(path):
    if path.endswith('.npy'):
        return np.load(path)
    with open(path) as f:
        return np.array(json.load(f)['H'])


def gram_popcount_ok(H):
    N = H.shape[0]
    bits = np.packbits((H < 0).astype(np.uint8), axis=1)
    pad = (-bits.shape[1]) % 8
    if pad:
        bits = np.concatenate([bits, np.zeros((N, pad), np.uint8)], 1)
    W = bits.view(np.uint64)
    for i in range(N):
        x = np.bitwise_xor(W[i][None, :], W[i:])
        pc = np.bitwise_count(x).sum(1, dtype=np.int64)
        dots = N - 2 * pc
        want = np.zeros(N - i, np.int64)
        want[0] = N
        if not np.array_equal(dots, want):
            j = int(np.flatnonzero(dots != want)[0]) + i
            return False, dict(i=i, j=j, dot=int(dots[j - i]))
    return True, None


def verify(H):
    H = np.asarray(H)
    out = dict(shape=list(H.shape))
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        return dict(out, ok=False, why='not square')
    if not np.isin(H, (-1, 1)).all():
        return dict(out, ok=False, why='entries not +-1')
    H = H.astype(np.int8)
    N = H.shape[0]
    ok1, bad = gram_popcount_ok(H)
    Hf = H.astype(np.float64)
    G = Hf @ Hf.T
    ok2 = bool(np.array_equal(G, N * np.eye(N)))
    out.update(N=N, popcount_gram_ok=ok1, first_bad=bad, blas_gram_ok=ok2, ok=bool(ok1 and ok2),
               sha256_int8=hashlib.sha256(np.ascontiguousarray(H).tobytes()).hexdigest())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('matrix')
    ap.add_argument('--out')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    H = load(a.matrix)
    res = verify(H)
    res['input'] = a.matrix
    if res.get('ok'):
        rng = np.random.default_rng(a.seed)
        i, j = (int(v) for v in rng.integers(0, H.shape[0], 2))
        Hb = np.array(H, np.int8)
        Hb[i, j] *= -1
        okb, bad = gram_popcount_ok(Hb)
        res['negative_control'] = dict(flipped=[i, j], rejected=not okb)
        res['ok'] = bool(res['ok'] and not okb)
    res['command'] = ' '.join(sys.argv)
    text = json.dumps(res, indent=1)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(text + '\n')
    print(json.dumps({k: res[k] for k in ('N', 'ok', 'popcount_gram_ok', 'blas_gram_ok', 'negative_control', 'sha256_int8')
                      if k in res}))
    sys.exit(0 if res['ok'] else 1)


if __name__ == '__main__':
    main()
