"""Is a known b=0 GS core (four +-1 circulants, sum PAF = 0 off zero) equivalent to T-matrix form,
i.e. do per-block shifts / reversals / negations make a_j b_j c_j d_j constant?
Equivalent to: canon(A o shift(B', s)) == canon(+-(C' o shift(D', e))) for some s, e and reversal choices,
where canon = lexicographically smallest rotation (rotations absorb the common shift of C', D')."""
import sys, json, itertools
import numpy as np

def load(path):
    blocks = {}
    for line in open(path):
        if "=" in line and line[0] in "ABCD" and line[1] == "=":
            blocks[line[0]] = np.array([1 if ch == "+" else -1 for ch in line[2:].strip()], dtype=np.int8)
    return [blocks[k] for k in "ABCD"]

def canon(x):
    n = len(x); xx = np.concatenate([x, x])
    rots = np.lib.stride_tricks.sliding_window_view(xx, n)[:n]
    best = min(range(n), key=lambda i: rots[i].tobytes())
    return rots[best].tobytes()

def keys(X, Y):
    out = set()
    for s in range(len(X)):
        p = X * np.roll(Y, s)
        out.add(canon(p)); out.add(canon(-p))
    return out

def check(path):
    A, B, C, D = load(path); n = len(A)
    paf = lambda x: np.array([int(np.dot(x, np.roll(x, -s))) for s in range(n)])
    tot = paf(A) + paf(B) + paf(C) + paf(D)
    assert (tot[1:] == 0).all(), "not a b=0 GS core"
    rev = lambda x: x[(-np.arange(n)) % n]
    hits = 0
    pairings = [((A, B), (C, D)), ((A, C), (B, D)), ((A, D), (B, C))]   # pairing is irrelevant but cheap
    for (X, Y), (Z, W) in pairings[:1]:
        for rY, rZ, rW in itertools.product((0, 1), repeat=3):
            Y2 = rev(Y) if rY else Y; Z2 = rev(Z) if rZ else Z; W2 = rev(W) if rW else W
            hits += len(keys(X, Y2) & keys(Z2, W2))
    return {"core": path, "n": n, "product_form_equivalent": hits > 0, "hits": hits}

if __name__ == "__main__":
    for p in sys.argv[1:]:
        print(json.dumps(check(p)), flush=True)
