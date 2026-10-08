"""Enumerate generalized bordered-GS options (Q, Psi, h) for an order 4v, with a cheap border pre-check.

General architecture (report/bge3_answer.md section 6; original = Psi all nontrivial characters, b = c - 1):
  Q abelian of order c, Psi a conjugation-closed set of b nontrivial characters (b odd), v = c*h + b,
  core target  sum_k PAF_k(q, z) = 4v*delta - 4 * sum_{psi in Psi} psi(q),
  border rows/columns: every block-part is a +-1 vector on Q whose Fourier support lies in Psi.
Pre-check (necessary): +-1 vectors on Q with Fourier support in Psi exist AND span the b-dimensional Psi-space
(the border frame needs sum_i K_ij K_ij^T = 4c P_Psi, which has rank b).
Options are reduced modulo Aut(Q) (orbit representatives; using a generating sample of Aut can only leave duplicate
representatives, never lose an option).
usage: bge3_options.py v [v ...] [--cmax 24] [--bmax 23] [--hmin 3] [--skip K  (resume after K written options)]
output: one JSON line per option: v, c, Q, b, h, psi (list of character vectors), n_psi_vectors (negative = lower bound,\n        scan stopped early at full rank), rank, precheck"""
import sys, json, itertools, random
import numpy as np
SKIP = 0

def args():
    a = sys.argv[1:]; opt = {"--cmax": 24, "--bmax": 23, "--hmin": 3, "--skip": 0}
    vs = []
    i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = int(a[i + 1]); i += 2
        else: vs.append(int(a[i])); i += 1
    return vs, opt

def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield []; return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k): yield [k] + rest

def abelian_groups(c):
    per = [[[p ** e for e in lam] for lam in partitions(k)] for p, k in factor(c).items()]
    for combo in itertools.product(*per):
        yield sorted([x for part in combo for x in part], reverse=True)

def setup_group(Q):
    Q = list(Q); c = int(np.prod(Q)); Qm = np.array(Q)
    co = np.array(np.unravel_index(np.arange(c), Q)).T                    # elements (also used for characters)
    V = np.exp(2j * np.pi * ((co[:, None, :] * co[None, :, :]) / Qm).sum(-1))  # V[k, q] = psi_k(q)
    # automorphisms: images of the generators
    cand = [np.where((((co * Qm[i]) % Qm) == 0).all(1))[0] for i in range(len(Q))]
    autos = []
    for imgs in itertools.product(*cand):
        I = co[list(imgs)]; img = (co @ I) % Qm
        keys = {tuple(x) for x in img}
        if len(keys) == c:
            autos.append(np.ravel_multi_index(tuple(img.T), Q))            # permutation of Q elements
    return c, co, V, autos

def char_perm(V, perm):
    """character index permutation induced by q -> perm[q]: psi_k o alpha = psi_{k'}"""
    key = {tuple(np.round(V[k], 6)): k for k in range(len(V))}
    return np.array([key[tuple(np.round(V[k][perm], 6))] for k in range(len(V))])

def psi_vectors(V, psi, c):
    """count +-1 vectors y on Q with y^(chi) = 0 for chi not in Psi (incl. trivial), and the rank of their span"""
    out = [k for k in range(c) if k not in set(psi)]
    W = V[out]                                                              # rows: characters that must vanish
    cnt, basis = 0, np.zeros((0, c))
    for s in range(0, 2 ** c, 1 << 18):
        if basis.shape[0] == len(psi): return -cnt, len(psi)          # early exit: full rank reached (count is a lower bound)
        idx = np.arange(s, min(2 ** c, s + (1 << 18)), dtype=np.int64)
        Y = 1 - 2 * ((idx[:, None] >> np.arange(c)) & 1)
        ok = np.abs(Y @ W.conj().T).max(1) < 1e-6
        if ok.any():
            G = Y[ok].astype(float); cnt += len(G)
            if basis.shape[0] < len(psi):
                basis = np.vstack([basis, G[: 4 * c]])
                u, sv, vt = np.linalg.svd(basis, full_matrices=False); basis = vt[sv > 1e-8]
    return cnt, int(basis.shape[0])

def options(v, cmax, bmax, hmin):
    for c in range(2, cmax + 1, 2):
        for b in range(1, min(c - 1, bmax) + 1, 2):
            if (v - b) % c or (v - b) // c < hmin: continue
            h = (v - b) // c
            for Q in abelian_groups(c):
                c_, co, V, autos = setup_group(Q)
                nontriv = list(range(1, c))
                conj = {k: int(np.ravel_multi_index(tuple((-co[k]) % np.array(Q)), Q)) for k in nontriv}
                real = [k for k in nontriv if conj[k] == k]
                pairs = sorted({tuple(sorted((k, conj[k]))) for k in nontriv if conj[k] != k})
                gens = autos if len(autos) <= 64 else random.Random(0).sample(autos, 12)
                cps = [char_perm(V, p) for p in gens]
                seen = set()
                for nr in range(1, min(len(real), b) + 1, 2):
                    if (b - nr) % 2 or (b - nr) // 2 > len(pairs): continue
                    for R in itertools.combinations(real, nr):
                        for P in itertools.combinations(pairs, (b - nr) // 2):
                            psi = tuple(sorted(list(R) + [x for pr in P for x in pr]))
                            if psi in seen: continue
                            # orbit under the sampled automorphisms (closure)
                            orb, stack = {psi}, [psi]
                            while stack:
                                s_ = stack.pop()
                                for cp in cps:
                                    t_ = tuple(sorted(int(cp[k]) for k in s_))
                                    if t_ not in orb: orb.add(t_); stack.append(t_)
                            seen |= orb
                            options.count = getattr(options, "count", 0) + 1
                            if options.count <= SKIP: continue           # already written by a previous run
                            n, rk = psi_vectors(V, list(psi), c)
                            yield dict(v=v, c=c, Q=Q, b=b, h=h, psi=[co[k].tolist() for k in psi],
                                       orbit_size=len(orb), n_psi_vectors=n, rank=rk, precheck=bool(n != 0 and rk == b))

if __name__ == "__main__":
    vs, opt = args()
    SKIP = opt["--skip"]
    for v in vs:
        for rec in options(v, opt["--cmax"], opt["--bmax"], opt["--hmin"]):
            print(json.dumps(rec), flush=True)
