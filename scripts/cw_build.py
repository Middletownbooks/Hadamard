"""Cooper-Wallis: T-matrices of order t + Williamson matrices of order w -> Hadamard matrix of order 4tw.

Own implementation (numpy, exact integers). T-matrices data for t = 103 transcribed from
schneiderlo/hadamard-2060 (orbit seed under the multiplier subgroup {1,46,56} of (Z/103)*).
"""
import sys, json, hashlib, itertools
import numpy as np

T103_MULT = (1, 46, 56)
T103_STATE = "04037071214226636752235627311036350"

def t103():
    P, seen, orbits = 103, set(), []
    for a in range(P):
        if a in seen: continue
        o = sorted({a * h % P for h in T103_MULT}); seen.update(o); orbits.append(o)
    assert len(orbits) == len(T103_STATE)
    rows = np.zeros((4, P), dtype=np.int64)
    for s, o in zip(T103_STATE, orbits):
        v = int(s); rows[v & 3, o] = -1 if v & 4 else 1
    return rows

def paf(x):
    n = len(x); return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(n)])

def check_T(T):
    t = T.shape[1]
    assert set(np.unique(T)) <= {-1, 0, 1}
    assert (np.abs(T).sum(0) == 1).all(), "supports must partition Z_t"
    tot = sum(paf(r) for r in T)
    assert tot[0] == t and (tot[1:] == 0).all(), "sum of PAFs must be t*delta"

def williamson_brute(w):
    """All symmetric +-1 sequences of length w; first quadruple with sum PAF = 4w*delta."""
    h = (w - 1) // 2
    syms = []
    for bits in itertools.product((1, -1), repeat=h + 1):
        x = [bits[0]] + list(bits[1:]) + list(bits[1:][::-1])
        syms.append(np.array(x))
    pafs = [paf(x) for x in syms]
    for i, j, k, l in itertools.combinations_with_replacement(range(len(syms)), 4):
        tot = pafs[i] + pafs[j] + pafs[k] + pafs[l]
        if (tot[1:] == 0).all():
            return np.array([syms[i], syms[j], syms[k], syms[l]])
    raise ValueError("no Williamson quadruple")

def circ(x):
    n = len(x); return np.array([np.roll(x, i) for i in range(n)])  # row i = x shifted right by i

def build(T, W):
    t, w = T.shape[1], W.shape[1]
    Tc = [circ(r) for r in T]; A, B, C, D = [circ(r) for r in W]
    k = np.kron
    X = k(Tc[0], A) + k(Tc[1], B) + k(Tc[2], C) + k(Tc[3], D)
    Y = -k(Tc[0], B) + k(Tc[1], A) + k(Tc[2], D) - k(Tc[3], C)
    Z = -k(Tc[0], C) - k(Tc[1], D) + k(Tc[2], A) + k(Tc[3], B)
    V = -k(Tc[0], D) + k(Tc[1], C) - k(Tc[2], B) + k(Tc[3], A)
    m = t * w
    # R: reversal on Z_t x Z_w (g -> -g), index i*w + j
    perm = np.array([((-i) % t) * w + ((-j) % w) for i in range(t) for j in range(w)])
    R = np.zeros((m, m), dtype=np.int64); R[np.arange(m), perm] = 1
    H = np.block([[X, -Y @ R, -Z @ R, -V @ R],
                  [Y @ R, X, -V.T @ R, Z.T @ R],
                  [Z @ R, V.T @ R, X, -Y.T @ R],
                  [V @ R, -Z.T @ R, Y.T @ R, X]])
    return H

if __name__ == "__main__":
    w = int(sys.argv[1]); out = sys.argv[2] if len(sys.argv) > 2 else None
    T = t103(); check_T(T)
    W = williamson_brute(w)
    H = build(T, W).astype(np.int64)
    n = H.shape[0]
    assert set(np.unique(H)) == {-1, 1}
    G = H @ H.T
    ok = bool((G == n * np.eye(n, dtype=np.int64)).all())
    rec = {"order": n, "t": 103, "w": w, "williamson": W.tolist(), "HHt_eq_nI": ok,
           "sha256_rows": hashlib.sha256(((H < 0).astype(np.uint8)).tobytes()).hexdigest()}
    print(json.dumps(rec))
    if out and ok:
        np.save(out, (H < 0).astype(np.uint8))
