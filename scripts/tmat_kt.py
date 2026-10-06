"""T-matrices via support split + k-tree 4-list join.

Fix S = supp(T1) u supp(T2) (a union of K-orbits); T3, T4 live on S^c. On S write T1 = (P+R)/2, T2 = (P-R)/2 with
P, R arbitrary +-1 on S (K-invariant), likewise P', R' on S^c. Then
    sum_i PAF(T_i) = t*delta   <=>   PAF(P) + PAF(R) + PAF(P') + PAF(R') = 2t*delta,
four INDEPENDENT lists (P and R range over the same set). Solved per S with kt_generic.kt_join4.
Usage: tmat_kt.py t k [--known] [--nS N] [--seed s] [--buckets B] [--maxorb M]
"""
import sys, json, time
import numpy as np
sys.path.insert(0, "scripts")
from tmat_search import subgroup
from kt_generic import kt_join4


def orbit_setup(t, k):
    K = subgroup(t, k)
    seen, orbits = set(), []
    for a in range(t):
        if a in seen: continue
        o = sorted({a * h % t for h in K}); seen.update(o); orbits.append(o)
    reps, seen = [], set()
    for o in orbits[1:]:
        if o[0] in seen: continue
        reps.append(o[0]); seen.update(o); seen.update((-x) % t for x in o)
    return K, orbits, np.array(reps)


def side_list(t, orbits, Oset, reps):
    """All +-1 K-invariant vectors on the orbits in Oset (first orbit's sign fixed to +1), with PAFs at reps."""
    Oset = list(Oset); r = len(Oset)
    idx = np.arange(2 ** (r - 1), dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(r - 1)) & 1).astype(np.int8)
    signs = np.concatenate([np.ones((len(idx), 1), np.int8), 1 - 2 * bits], axis=1)   # n x r
    X = np.zeros((len(idx), t), np.float32)
    for j, o in enumerate(Oset): X[:, orbits[o]] = signs[:, [j]]
    P = np.empty((len(idx), len(reps)), np.int32)
    for s in range(0, len(idx), 1 << 16):
        F = np.fft.rfft(X[s:s + (1 << 16)], axis=1)
        P[s:s + (1 << 16)] = np.rint(np.fft.irfft(np.abs(F) ** 2, n=t, axis=1)[:, reps]).astype(np.int32)
    return signs, P


def tm_check(T):
    t = T.shape[1]
    if not (np.abs(T).sum(0) == 1).all(): return False
    F = np.fft.fft(T, axis=1); tot = np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real).astype(np.int64)
    return bool(tot[0] == t and (tot[1:] == 0).all())


def run_S(t, orbits, reps, Sset, nbuckets, seed, log):
    m = len(orbits); Sc = [o for o in range(m) if o not in Sset]
    if len(Sset) < 2 or len(Sc) < 2: return None
    sgS, PS = side_list(t, orbits, Sset, reps)
    sgC, PC = side_list(t, orbits, Sc, reps)
    found = []
    def verify(idx):
        T = np.zeros((4, t), np.int64)
        for (sg, Os, a, b) in ((sgS, Sset, 0, 1), (sgC, Sc, 2, 3)):
            pass
        P1, R1, P2, R2 = sgS[idx[0]], sgS[idx[1]], sgC[idx[2]], sgC[idx[3]]
        for j, o in enumerate(Sset):
            T[0, orbits[o]] = (int(P1[j]) + int(R1[j])) // 2; T[1, orbits[o]] = (int(P1[j]) - int(R1[j])) // 2
        for j, o in enumerate(Sc):
            T[2, orbits[o]] = (int(P2[j]) + int(R2[j])) // 2; T[3, orbits[o]] = (int(P2[j]) - int(R2[j])) // 2
        if tm_check(T): found.append(T); return True
        return False
    kt_join4([PS, PS, PC, PC], verify, nbuckets=nbuckets, max_solutions=1, log=log, seed=seed)
    return found[0] if found else None


if __name__ == "__main__":
    t, k = int(sys.argv[1]), int(sys.argv[2])
    g = lambda f, d=None, ty=str: ty(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
    nS, seed, nb = g("--nS", 1, int), g("--seed", 1, int), g("--buckets", 2000, int)
    K, orbits, reps = orbit_setup(t, k); m = len(orbits)
    rng = np.random.default_rng(seed); t0 = time.time()
    quiet = lambda s: print(s, flush=True) if ("kt_done" in s or "kt_setup" in s) else None
    print(json.dumps({"t": t, "k": k, "orbits": m, "classes": len(reps)}), flush=True)
    for i in range(nS):
        if "--known" in sys.argv:
            from cw_build import t103
            T = t103(); oid = {x: j for j, o in enumerate(orbits) for x in o}
            cands = [sorted({oid[x] for x in range(t) if T[a, x] != 0 or T[b, x] != 0}) for a, b in ((0, 1), (0, 2), (0, 3))]
            Sset = min(cands, key=lambda c: abs(2 * len(c) - m))          # most balanced pairing
        else:
            Sset = sorted(rng.choice(m, m // 2 + int(rng.integers(0, 2)) * (m % 2), replace=False).tolist())
        T = run_S(t, orbits, reps, Sset, nb, seed * 1000 + i, quiet)
        rec = {"S_index": i, "S_orbits": len(Sset), "found": T is not None, "secs": round(time.time() - t0, 1)}
        if T is not None:
            rec["T"] = T.tolist(); rec["FOUND"] = True
        print(json.dumps(rec), flush=True)
        if T is not None: break
