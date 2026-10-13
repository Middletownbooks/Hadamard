r"""Exact alpha-congruence screen for generalised-Paley slices (b = 1) [P, M].

Slice: Z_F fixed pointwise, index-s classes on Z_q, f = (q - 1)/s.  On the coset g of block k, x = alpha_k(g) at z = 0
and c_j on the class C_j, so the stage-1 letter is y = alpha + f*sigma, sigma = sum_j c_j = s (mod 2).
[P] Write y = beta + 2f j with beta = alpha + f [s odd] and j an integer (sigma - [s odd] is even).  Then
    PAF(y)(h) = PAF(beta)(h) + 2f sum_g (beta_g j_{g+h} + j_g beta_{g+h}) + 4f^2 sum_g j_g j_{g+h}.
The middle sum is even: if beta is odd it is = sum_g (j_{g+h} + j_g) = 2 sum j (mod 2); if beta is even (s, f odd) each
term is even.  Hence for every lag h != 0
    sum_k PAF(beta_k)(h) = -4q          (mod 4f)
and rows sum_g beta_k(g) = r_k (mod 2f) (F even), rows (2, 0, 0, 0) after reordering and negating blocks.
The alpha_k are arbitrary +-1 sequences on Z_F, so this is a finite check: all 2^F sequences, keys (row mod 2f,
PAF(beta) lags 1..F/2 mod 4f), a key-level 2 + 2 join.  No solution => the slice is empty.
usage: python3 theory/gp_alpha_screen.py F q s [F q s ...]"""
import sys, json, time
import numpy as np

def keys(F, shift, mod, rmod):
    n = 1 << F
    A = (((np.arange(n, dtype=np.int64)[:, None] >> np.arange(F)[None, :]) & 1) * 2 - 1).astype(np.int32) + shift
    P = np.stack([(A * np.roll(A, -h, axis=1)).sum(1) for h in range(1, F // 2 + 1)], 1) % mod
    R = A.sum(1) % rmod
    return R, P

def screen(F, q, s):
    f = (q - 1) // s; mod = 4 * f; rmod = 2 * f; t0 = time.time(); shift = f if s % 2 else 0
    tgt = (-4 * q) % mod
    R, P = keys(F, shift, mod, rmod)
    H = P.shape[1]; base = np.array([mod ** i for i in range(H)], dtype=object)
    def enc(rows): return {tuple(r) for r in rows}
    K2 = enc(P[R == 2 % rmod]); K0 = enc(P[R == 0])
    T = np.full(H, tgt)
    S = {tuple((np.array(a) + np.array(b)) % mod) for a in K0 for b in K0} if len(K0) <= 6000 else None
    found = None
    if S is not None:
        for a in K2:
            for d in K0:
                need = tuple((T - np.array(a) - np.array(d)) % mod)
                if need in S: found = (a, d, need); break
            if found: break
    return dict(F=F, q=q, s=s, f=f, mod=mod, target=int(tgt), keys_row2=len(K2), keys_row0=len(K0),
                decided=S is not None, alpha_ok=(found is not None) if S is not None else None,
                secs=round(time.time() - t0, 1))

if __name__ == "__main__":
    a = list(map(int, sys.argv[1:]))
    for i in range(0, len(a), 3): print(json.dumps(screen(*a[i:i + 3])), flush=True)
