"""Independent verification of b = 0 Goethals-Seidel cores found by the flat-pool k-tree, and an inequivalence test
against reference cores.  Shares no code with the search (pattern of scripts/verify_h3532.py: integer PAFs by dot
products, no FFT; the GS array built from circulants; H H^T = 4v I checked exactly in int64 row blocks).

Input: a solutions file of JSON lines with "blocks": {"A": "+-...", ...} (scripts/kt_join_slots*.py output), or a
text file with lines A=..., B=..., C=..., D=....  For every core:
  1. sum_k PAF_k = 4v at 0 and 0 elsewhere, exactly; row sums;
  2. the slice: invariant under x -> u x for some u of order d (d from --d);
  3. the GS array [[A, BR, CR, DR], [-BR, A, D^T R, -C^T R], [-CR, -D^T R, A, B^T R], [-DR, C^T R, -B^T R, A]] has
     H H^T = 4v I (with --hadamard; O(v^3), seconds at v ~ 300-900);
  4. equivalence class key under the standard core operations (cyclic shifts, negation and reversal of any block,
     permutations of the blocks, a common multiplier x(i) -> x(m i)): the multiset of the four PAF vectors is
     invariant under shifts, negation, reversal and permutation and is mapped by i -> m i under a multiplier, so the
     key = the least, over m in Z_v^*, of the sorted tuple of (PAF_k(m i))_i.  Equal keys are necessary for
     equivalence; different keys prove inequivalence under these operations [P].
usage: python3 scripts/verify_b0_core.py FILE --d D [--hadamard] [--ref REF1 [REF2 ...]]   (REF: .txt A=.. or .npy (4, v))"""
import json, sys
import numpy as np


def read_cores(path):
    if path.endswith('.npy'):
        return [np.load(path).astype(np.int64)]
    txt = open(path).read()
    if txt.lstrip().startswith('{'):
        out = []
        for l in txt.splitlines():
            if l.strip():
                b = json.loads(l)['blocks']
                out.append(np.array([[1 if c == '+' else -1 for c in b[k]] for k in 'ABCD'], dtype=np.int64))
        return out
    rows = [l.split('=', 1)[1].strip() for l in txt.splitlines() if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    return [np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)]


def paf(x):
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(len(x))], dtype=np.int64)


def circ(x):
    return np.array([np.roll(x, j) for j in range(len(x))], dtype=np.int64)


def is_hadamard(X):
    v = X.shape[1]; R = np.fliplr(np.eye(v, dtype=np.int64)); A, B, C, D = (circ(x) for x in X)
    H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    N = H.shape[0]
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T; E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True


def invariant_under(X, d):
    v = X.shape[1]; xs = np.arange(v)
    for u in range(2, v):
        if next(k for k in range(1, v) if pow(u, k, v) == 1) == d:
            return u, [bool(np.array_equal(x[(u * xs) % v], x)) for x in X]   # the order-d subgroup is unique (v prime)
    return None, []


def eq_key(X):
    v = X.shape[1]; P = [paf(x) for x in X]; xs = np.arange(v); best = None
    for m in range(1, v):
        k = tuple(sorted(tuple(p[(m * xs) % v].tolist()) for p in P))
        if best is None or k < best: best = k
    return best


def main():
    a = sys.argv[1:]
    d = int(a[a.index('--d') + 1]); had = '--hadamard' in a
    refs = a[a.index('--ref') + 1:] if '--ref' in a else []
    refkeys = {r: eq_key(c) for r in refs for c in read_cores(r)}
    keys = {}; ok_all = True
    for i, X in enumerate(read_cores(a[0])):
        v = X.shape[1]; P = sum(paf(x) for x in X)
        ok_paf = bool(P[0] == 4 * v and (P[1:] == 0).all())
        u, inv = invariant_under(X, d); k = eq_key(X)
        rec = dict(core=i, v=v, row_sums=X.sum(1).tolist(), paf_ok=ok_paf, u=u, u_invariant=inv,
                   class_id=keys.setdefault(k, len(keys)), equals_ref=[r for r, rk in refkeys.items() if rk == k])
        if had: rec['hadamard'] = is_hadamard(X)
        rec['ok'] = bool(ok_paf and all(inv) and rec.get('hadamard', True))
        ok_all &= rec['ok']; print(json.dumps(rec), flush=True)
    print(json.dumps(dict(cores=i + 1, classes=len(keys), all_ok=ok_all)))
    sys.exit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
