r"""The published v = 419 (b = 1) core is an exact Paley-coset core with p = 11 on G' = Z38, inside the order-3 multiplier
slice (report section 7.5) [M].

Checks, all exact:
  1. data/core_v419_b1.txt is a b = 1 core (sum PAF = 4m at 0, -4 elsewhere; row sums (2, 0, 0, 0)).
  2. After translating blocks B, C, D by 208 (= -210, i.e. by -1 in Z11 and by 18 in Z38; block A untranslated, which
     keeps every PAF), each block is alpha*delta + u*chi + w*eps on every Z11-coset of G = Z38 x Z11 (CRT order).
  3. The labels (alpha, u, w) on Z38 are invariant under g -> 7g (order 3), so the core lies in the "Paley(11) + mult:7"
     slice = the orbit slice of the order-15 multiplier group (squares mod 11) x <7 mod 38>, as the forensics found.
  4. Theorem 1 (S1, S2 real, X = 0), the norm identity, Lemma 4 data and Lemma 6 (at both real characters of Z38).
usage: python3 theory/paley_anatomy419.py"""
import sys, os, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_algebra import components
from paley_lemma6 import check as lemma6_check

def pafc(x): return np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64)
def ccc(f, g): return np.rint(np.fft.ifft(np.conj(np.fft.fft(f)) * np.fft.fft(g)).real).astype(np.int64)

def load():
    L = open(os.path.join(ROOT, 'data/core_v419_b1.txt')).read().split('\n')
    return [np.array([1 if c == '+' else -1 for c in l.split('=')[1].strip()], np.int64) for l in L[1:] if '=' in l and l[0] in 'ABCD']

def main():
    v, p, n = 419, 11, 38; m = n * p
    X = load()
    P = sum(pafc(x) for x in X)
    out = dict(core_ok=bool(P[0] == 4 * m and (P[1:] == -4).all()), rows=[int(x.sum()) for x in X])
    shifts = [0, 208, 208, 208]
    Y = []
    for x, s in zip(X, shifts):
        xs = np.roll(x, -s); y = np.zeros((n, p), np.int64)
        for i in range(m): y[i % n, i % p] = xs[i]
        Y.append(y)
    comps = [components(y, p) for y in Y]
    out['paley_coset_after_shifts'] = all(c[3] for c in comps)
    A = np.array([c[0] for c in comps]); U = np.array([c[1] for c in comps]); W = np.array([c[2] for c in comps])
    out['mult7_invariant'] = bool(all((M[:, (7 * np.arange(n)) % n] == M).all() for M in (A, U, W)))
    K_k = (W != 0).sum(1); Nm = int(((W != 0) & (W == -A)).sum()); K = int(K_k.sum())
    out.update(K=K, Nminus=Nm, Nplus=K - Nm, K_k=K_k.tolist(), norm_identity=(p + 1) * K - 4 * Nm == 4 * n - 4,
               t=W.sum(1).tolist(), rho=U.sum(1).tolist())
    # Theorem 1 on G' = Z38
    beta = A - W; y = A + (p - 1) * W
    TJ = np.full(n, -4 * p); TJ[0] = 4 * m - 4 * (p - 1)
    S1 = sum(pafc(y[k]) for k in range(4)); S2 = sum(pafc(beta[k]) + p * pafc(U[k]) for k in range(4))
    Xs = sum(ccc(U[k], beta[k]) - ccc(beta[k], U[k]) for k in range(4))
    out.update(S1=bool((S1 == TJ).all()), S2=bool(S2[0] == 4 * v and (S2[1:] == 0).all()), X0=bool((Xs == 0).all()))
    s = sum(pafc(W[k]) for k in range(4))
    out['lemma4_s_off0'] = sorted(set(s[1:].tolist())); out['sum_t2_eq_K_plus_sum_s'] = int((W.sum(1) ** 2).sum()) == K + int(s[1:].sum())
    # orbits of the multiplier 7 and which ones carry the gadget
    orbs, seen = [], set()
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * 7) % n
        seen.update(o); orbs.append(sorted(o))
    out['orbits'] = len(orbs)
    out['gadget_orbits'] = [[(o[0], len(o), int(W[k, o[0]]), 'Jprime' if W[k, o[0]] == -A[k, o[0]] else 'J') for o in orbs if W[k, o[0]] != 0] for k in range(4)]
    L6 = lemma6_check(v, 1, [y_ for y_ in Y], p)
    out['lemma6_all_real_chars'] = L6['all_identities']
    print(json.dumps(out))
    return out

if __name__ == "__main__":
    main()
