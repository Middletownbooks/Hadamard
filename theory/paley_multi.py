r"""Iterated ("multi-") Paley-coset ansatz: G = G'' x Z_q1 x ... x Z_qr, every block invariant under the multiplier group
QR_q1 x ... x QR_qr (each QR_qi acting on its own factor only).  r = 1 is the repo's Paley-coset ansatz; r = 2 is the
"doubly twisted" (9-term coset types) of report/prompt_paley_design.md Q1.

Orbits.  On Z_q the QR-orbits are {0}, QR, NQR (classes 0, 1, 2).  On prod Z_qi the orbits are products of classes, 3^r of
them, with sizes prod (1 or (q-1)/2).  A block is a +-1 value per (g'', orbit): 3^r bits per element of G''.

Reduced equations (exact).  For a lag (h'', s) the periodic autocorrelation depends only on h'' and the class vector of s:
    R(h'', S) = sum_{g''} sum_{O1, O2} x(g'', O1) x(g''+h'', O2) N_S[O1, O2],
    N_S[(a_i), (b_i)] = prod_i N_qi[a_i, b_i, S_i],   N_q[a, b, s] = #{z in O_a : z + s in O_b}  (cyclotomic numbers).
The GS condition sum_k PAF_k = T is equivalent to sum_k R_k = T on the |G''| * 3^r classes.  `validate` checks R against
FFT of the expanded sequences (every element of every class).

Field decomposition (proved in report/paley_design_answer.md section 2): span{delta, chi_q, eps_q} (x) Q = Q + Q(sqrt(q*)),
q* = (-1)^((q-1)/2) q, so the condition splits into 2^r independent "sum of norms" conditions over G'', one per subset of
the primes (the field Q(sqrt(q_i*), i in subset)).

usage: python3 theory/paley_multi.py validate
       python3 theory/paley_multi.py mitm v b g2mods... -- q1 q2 ...     (exhaustive 4-list join, small slices)
       python3 theory/paley_multi.py estimate v b g2mods... -- q1 q2 ... (Gaussian, row-sum conditioned)"""
import sys, os, json, itertools as it, math, time
from collections import Counter, defaultdict
import numpy as np

def legendre(p):
    chi = np.zeros(p, dtype=np.int64)
    for z in range(1, p): chi[z] = 1 if pow(z, (p - 1) // 2, p) == 1 else -1
    return chi

def gf(q):
    """GF(q), q = r^k (k <= 3): additive group Z_r^k (ravel index, last coordinate fastest) and the class of each element
    under the squares of GF(q)^*: 0 -> {0}, 1 -> nonzero squares, 2 -> non-squares."""
    r = next(t for t in range(2, q + 1) if q % t == 0); k = round(math.log(q, r)); assert r ** k == q
    if k == 1:
        chi = legendre(q); return (q,), np.where(np.arange(q) == 0, 0, np.where(chi == 1, 1, 2))
    # find a monic irreducible polynomial of degree k over F_r (no roots suffices for k <= 3)
    for c in it.product(range(r), repeat=k):
        if all((pow(x, k, r) + sum(c[i] * pow(x, i, r) for i in range(k))) % r for x in range(r)): break
    def mul(a, b):
        prod = [0] * (2 * k - 1)
        for i in range(k):
            for j in range(k): prod[i + j] = (prod[i + j] + a[i] * b[j]) % r
        for d in range(2 * k - 2, k - 1, -1):                 # reduce x^d = -sum c_i x^(d-k+i)
            t = prod[d]; prod[d] = 0
            for i in range(k): prod[d - k + i] = (prod[d - k + i] - t * c[i]) % r
        return tuple(prod[:k])
    elems = list(it.product(range(r), repeat=k))                 # tuples (a_0..a_{k-1}) = sum a_i x^i
    nz = [e for e in elems if any(e)]
    sq = {mul(e, e) for e in nz}
    # ravel index with last coordinate fastest; element tuple is (a_0, ..., a_{k-1})
    cl = np.zeros(q, np.int64)
    for e in elems:
        i = int(np.ravel_multi_index(e, (r,) * k))
        cl[i] = 0 if not any(e) else (1 if e in sq else 2)
    return (r,) * k, cl

def classes(q):
    return gf(q)[1]                                              # class index of each element of GF(q)'s additive group

def cyclo_N(q):
    """N[a, b, s] = #{z in O_a : z + rep_s in O_b}, reps of the classes 0, square, non-square (additive group of GF(q))"""
    mods, cl = gf(q); A = AG(mods)
    reps = [0, int(np.nonzero(cl == 1)[0][0]), int(np.nonzero(cl == 2)[0][0])]
    N = np.zeros((3, 3, 3), np.int64)
    for s, r in enumerate(reps):
        for z in range(q): N[cl[z], cl[A.add[z, r]], s] += 1
    # sanity: the count does not depend on the representative chosen
    for s in range(3):
        for r in np.nonzero(cl == s)[0][:3]:
            M = np.zeros((3, 3), np.int64)
            for z in range(q): M[cl[z], cl[A.add[z, r]]] += 1
            assert (M == N[:, :, s]).all(), ("classes not multiplier-invariant", q)
    return N

class AG:
    """finite abelian group Z_{n1} x ... (ravel index, last factor fastest)"""
    def __init__(s, mods):
        s.mods = tuple(int(t) for t in mods) if len(mods) else (1,)
        s.n = int(np.prod(s.mods))
        s.co = np.array(np.unravel_index(np.arange(s.n), s.mods)).T
        s.add = np.ravel_multi_index(tuple(np.moveaxis((s.co[:, None, :] + s.co[None, :, :]) % np.array(s.mods), -1, 0)), s.mods)
        s.neg = np.ravel_multi_index(tuple(np.moveaxis((-s.co) % np.array(s.mods), -1, 0)), s.mods)

class MP:
    def __init__(s, g2mods, primes):
        s.G2 = AG(g2mods); s.n2 = s.G2.n; s.primes = list(primes); s.r = len(primes)
        s.no = 3 ** s.r                                            # orbits per G''-element
        s.ocls = np.array(list(it.product(range(3), repeat=s.r)), np.int64).reshape(s.no, s.r)   # orbit -> class vector
        s.osize = np.array([math.prod(1 if c == 0 else (q - 1) // 2 for c, q in zip(oc, s.primes)) for oc in s.ocls], np.int64)
        Ns = [cyclo_N(q) for q in s.primes]
        # N_S[O1, O2] for each lag class S (same indexing as orbits)
        s.NS = np.zeros((s.no, s.no, s.no), np.int64)
        for S, sc in enumerate(s.ocls):
            for O1, c1 in enumerate(s.ocls):
                for O2, c2 in enumerate(s.ocls):
                    s.NS[O1, O2, S] = math.prod(int(N[a, b, t]) for N, a, b, t in zip(Ns, c1, c2, sc))
        s.m = s.n2 * math.prod(s.primes)

    def reduced_paf(s, x):
        """x: (N, n2, no) +-1 -> R: (N, n2, no) (lag (h'', class S))"""
        N = x.shape[0]; R = np.zeros((N, s.n2, s.no), np.int64)
        xf = x.astype(np.int64)
        for h in range(s.n2):
            B = xf[:, s.G2.add[:, h], :]                             # B[:, g] = x[:, g + h]
            M = np.einsum('ngi,ngj->nij', xf, B)                     # sum_g x(g, O1) x(g+h, O2)
            R[:, h, :] = np.einsum('nij,ijs->ns', M, s.NS)
        return R

    def rowsum(s, x):
        return (x * s.osize[None, None, :]).sum((1, 2))

    def fshape(s):
        return tuple(t for q in s.primes for t in gf(q)[0])

    def expand(s, xk):
        """one block (n2, no) -> full +-1 array over G'' x F_q1 x ... (additive groups; shape g2mods + factor mods)"""
        cls = [classes(q) for q in s.primes]
        full = np.zeros((s.n2,) + tuple(s.primes), np.int64)
        for z in it.product(*[range(q) for q in s.primes]):
            O = 0
            for c, zz in zip(cls, z): O = 3 * O + int(c[zz])
            full[(slice(None),) + z] = xk[:, O]
        return full.reshape(s.G2.mods + s.fshape())

    def full_paf(s, xk):
        f = s.expand(xk); F = np.fft.fftn(f)
        return np.rint(np.fft.ifftn(np.abs(F) ** 2).real).astype(np.int64)

    def target(s, b, Hmask2=None):
        """reduced target over (h'', S): b = 0: 0 off 0; b = 1: -4 off 0; b >= 3: -4b on H'' x (all classes), +4 else"""
        T = np.full((s.n2, s.no), -4 * b, np.int64)
        if b >= 3:
            T[:] = 4; T[Hmask2, :] = -4 * b
        T[0, 0] = 4 * s.m
        return T

def validate():
    rng = np.random.default_rng(0); ok = True; n = 0
    for g2, pr in [((2,), (3, 5)), ((2,), (5, 7)), ((3,), (5, 13)), ((1,), (3, 7, 5)), ((2, 2), (7,)), ((4,), (11, 3)),
                   ((2,), (9, 7)), ((1,), (9, 5)), ((2,), (27,)), ((1,), (25, 3))]:
        P = MP(g2, pr)
        x = rng.choice(np.array([-1, 1]), (6, P.n2, P.no))
        R = P.reduced_paf(x)
        cls = [classes(q) for q in pr]
        for k in range(6):
            F = P.full_paf(x[k]).reshape((P.n2,) + tuple(pr))
            for h in range(P.n2):
                for z in it.product(*[range(q) for q in pr]):
                    S = 0
                    for c, zz in zip(cls, z): S = 3 * S + int(c[zz])
                    ok &= bool(F[(h,) + z] == R[k, h, S]); n += 1
    print(json.dumps({"validate_multi_reduced_paf_vs_fft": bool(ok), "checks": n}))
    return ok

# ---------------- exhaustive 4-list join ----------------
def enumerate_blocks(P):
    nb = P.n2 * P.no
    assert nb <= 22, "too many bits per block for exhaustive enumeration"
    idx = np.arange(2 ** nb, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(nb)) & 1).astype(np.int8)
    x = (1 - 2 * bits).reshape(-1, P.n2, P.no)
    return x

def profiles(v, b, rows_avail):
    if b == 1: return [(2, 0, 0, 0)]
    pos = sorted(a for a in rows_avail if a >= 0)
    return [q for q in it.combinations_with_replacement(pos, 4) if sum(t * t for t in q) == 4 * v]

def mitm(v, b, g2mods, primes, Hmask2=None, max_sols=20, verbose=True, Tplanted=None):
    """Exact 4-list join.  Blocks are grouped by their compression y(g'') (sum over the Z_q1 x ... part); the compressed
    condition sum_k PAF_{G''}(y_k) = T_JJ (exact consequence) selects compatible 4-tuples of groups, and only those are
    joined on the full reduced vector.  Counts ordered 4-tuples with block A's row sum >= 0 (block negation symmetry)."""
    P = MP(g2mods, primes); assert P.m == v - b, (P.m, v, b)
    t0 = time.time()
    x = enumerate_blocks(P)
    Y = (x * P.osize[None, None, :]).sum(2)                         # (N, n2) compression onto G''
    T = P.target(b, Hmask2) if Tplanted is None else Tplanted
    TJJ = (T * P.osize[None, :]).sum(1)                            # compressed target over G''
    Tv = T.reshape(-1)
    groups = defaultdict(list)
    for i, yy in enumerate(map(tuple, Y.tolist())): groups[yy].append(i)
    keys = list(groups); Ky = np.array(keys, np.int64)
    def pafG(yv):
        return np.array([(yv * yv[P.G2.add[:, h]]).sum() for h in range(P.n2)], np.int64)
    Py = np.array([pafG(k) for k in Ky])                          # (G, n2)
    rowsum = Ky.sum(1)
    # compatible 4-tuples of compression groups: pairs (A,B) with A rowsum >= 0
    pair = defaultdict(list)
    for i in range(len(keys)):
        S = Py[i][None, :] + Py
        for j in np.nonzero((S <= TJJ[None, :] + 10 ** 9).all(1))[0]:
            pair[S[j].tobytes()].append((i, int(j)))
    quads = []
    for kb, lst in pair.items():
        need = (TJJ - np.frombuffer(kb, np.int64)).tobytes()
        if need in pair:
            for (i, j) in lst:
                if rowsum[i] < 0: continue
                for (k, l) in pair[need]: quads.append((i, j, k, l))
    out = {"v": v, "b": b, "G''": list(g2mods), "primes": list(primes), "bits_per_block": int(P.n2 * P.no),
           "compression_groups": len(keys), "compatible_group_quads": len(quads), "solutions": 0, "examples": []}
    if verbose: print(json.dumps({"enumerated": len(x), "groups": len(keys), "quads": len(quads),
                                  "secs": round(time.time() - t0, 1)}), flush=True)
    Rc = {}
    def Rg(gi):
        if gi not in Rc:
            ids = np.array(groups[keys[gi]], np.int64); Rc[gi] = (ids, P.reduced_paf(x[ids]).reshape(len(ids), -1))
        return Rc[gi]
    total = 0
    # linear 64-bit hash: key(R) = R . wts (mod 2^64), additive, so pair sums hash to sums of keys; hits verified exactly
    wts = np.random.default_rng(12345).integers(1, 2 ** 62, Tv.shape[0], dtype=np.int64).astype(np.uint64)
    def hkey(R): return (R.astype(np.uint64) * wts[None, :]).sum(1)          # uint64 wraparound
    tkey = (Tv.astype(np.uint64) * wts).sum()
    Kc = {}
    def Kg(gi):
        if gi not in Kc: Kc[gi] = hkey(Rg(gi)[1])
        return Kc[gi]
    byCD = defaultdict(list)
    for (i, j, k, l) in quads: byCD[(k, l)].append((i, j))
    for (k, l), ABs in byCD.items():
        (iC, RC), (iD, RD) = Rg(k), Rg(l)
        KCD = (Kg(k)[:, None] + Kg(l)[None, :]).reshape(-1)
        order = np.argsort(KCD, kind='stable'); Ks = KCD[order]
        for (i, j) in ABs:
            (iA, RA), (iB, RB) = Rg(i), Rg(j)
            need = (tkey - Kg(i)[:, None] - Kg(j)[None, :]).reshape(-1)
            lo = np.searchsorted(Ks, need, 'left'); hi = np.searchsorted(Ks, need, 'right')
            for jj in np.nonzero(hi > lo)[0]:
                a, bb = divmod(int(jj), len(RB))
                for h in order[lo[jj]:hi[jj]]:
                    c, d = divmod(int(h), len(RD))
                    if not (RA[a] + RB[bb] + RC[c] + RD[d] == Tv).all(): continue      # hash collision
                    total += 1
                    if len(out["examples"]) < max_sols:
                        out["examples"].append([x[iA[a]].tolist(), x[iB[bb]].tolist(), x[iC[c]].tolist(), x[iD[d]].tolist()])
    out["solutions"] = total; out["secs"] = round(time.time() - t0, 1)
    ok = True
    for ex in out["examples"]:
        tot = sum(P.full_paf(np.array(bk)) for bk in ex)
        ok &= bool((tot == full_target(P, b, Hmask2)).all())
    out["examples_verified_fft"] = ok if out["examples"] else None
    return out

def full_target(P, b, Hmask2=None):
    shp = P.G2.mods + P.fshape()
    T = np.full(shp, -4 * b, np.int64)
    if b >= 3:
        T[:] = 4
        Hm = np.asarray(Hmask2).reshape(P.G2.mods)
        T[Hm] = -4 * b
    T[(0,) * len(shp)] = 4 * P.m
    return T

# ---------------- Gaussian estimate (row-sum conditioned), for slices too large to enumerate ----------------
def estimate(v, b, g2mods, primes, N=60000, seed=1, Hmask2=None, maxprof=8):
    P = MP(g2mods, primes); rng = np.random.default_rng(seed)
    T = P.target(b, Hmask2).reshape(-1)
    nb = P.n2 * P.no
    sizes = np.tile(P.osize, P.n2)
    # exact row-sum counts by DP over the nb independent +-size terms
    cnt = Counter({0: 1})
    for sz in sizes.tolist():
        E = Counter()
        for a, c in cnt.items(): E[a + sz] += c; E[a - sz] += c
        cnt = E
    profs = profiles(v, b, cnt.keys())
    profs = sorted(profs, key=lambda q: -sum(math.log2(cnt[t]) for t in q))[:maxprof]
    stats = {}
    def samp(a):
        if a not in stats:
            res = []
            for _ in range(2000):
                xx = rng.choice(np.array([-1, 1], np.int8), (50000, P.n2, P.no))
                xx = xx[P.rowsum(xx) == a]; res.append(xx)
                if sum(len(t) for t in res) >= N: break
            xx = np.concatenate(res)[:N]
            stats[a] = None if len(xx) < 2000 else P.reduced_paf(xx).reshape(len(xx), -1)
        return stats[a]
    tot = 0.0; best = None
    for q in profs:
        Rs = [samp(a) for a in q]
        if any(r is None for r in Rs): continue
        mu = sum(r.mean(0) for r in Rs); Cv = sum(np.cov(r.T) for r in Rs); Rc = sum(r[:2000] for r in Rs)
        cols, seen = [], {}
        for j in range(1, len(T)):
            kk = Rc[:, j].tobytes()
            if kk in seen:
                if T[j] != T[seen[kk]]: cols = None; break
                continue
            seen[kk] = j; cols.append(j)
        if cols is None: continue
        w, V = np.linalg.eigh(Cv[np.ix_(cols, cols)]); kw = w > 1e-6 * w.max(); pr = V.T @ (T[cols] - mu[cols])
        if (~kw).any() and np.abs(pr[~kw]).max() > 1e-6 * max(1, np.abs(pr).max()): continue
        d = Rc[:, cols] - Rc[0, cols]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j])))) for j in range(len(cols))])
        if not all((T[cols][j] - Rc[0, cols][j]) % steps[j] == 0 for j in range(len(cols))): continue
        lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum() - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
        perms = len(set(it.permutations(q))); sgn = 2 ** sum(1 for t in q if t != 0)
        lE = sum(math.log2(cnt[t]) for t in q) + lp + math.log2(perms * sgn)
        tot += 2 ** lE
        if best is None or lE > best[1]: best = (list(q), round(lE, 1), round(-lp, 1), int(kw.sum()))
    rec = {"v": v, "b": b, "G''": list(g2mods), "primes": list(primes), "free_bits": 4 * nb}
    if best is None: rec.update(feasible=False)
    else: rec.update(log2E=round(math.log2(tot), 1), best_profile=best[0], constraint_bits=best[2], eqs=best[3])
    return rec

def parse(argv):
    i = argv.index('--'); a = list(map(int, argv[:i])); pr = list(map(int, argv[i + 1:]))
    return a[0], a[1], tuple(a[2:]), pr

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "validate": validate()
    elif cmd == "mitm":
        v, b, g2, pr = parse(sys.argv[2:])
        res = mitm(v, b, g2, pr)
        print(json.dumps({k: res[k] for k in res if k != "examples"}), flush=True)
        if res["solutions"]:
            fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs",
                              f"paley_multi_v{v}_b{b}_{'x'.join(map(str, g2))}_{'x'.join(map(str, pr))}.json")
            json.dump(res, open(fn, "w"))
    elif cmd == "planted":                                           # validation: target from 4 random blocks
        g2 = tuple(map(int, sys.argv[2].split(','))); pr = list(map(int, sys.argv[3].split(',')))
        P = MP(g2, pr); rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
        xs = rng.choice(np.array([-1, 1]), (4, P.n2, P.no)); xs[0] *= 1 if P.rowsum(xs[:1])[0] >= 0 else -1
        Tp = P.reduced_paf(xs).sum(0)
        res = mitm(P.m, 0, g2, pr, Tplanted=Tp, verbose=False)
        found = any(all((np.array(ex[k]) == xs[k]).all() for k in range(4)) for ex in res["examples"]) if res["examples"] else False
        print(json.dumps({"planted_G''": list(g2), "primes": pr, "solutions": res["solutions"], "compatible_group_quads": res["compatible_group_quads"],
                          "planted_found_among_first_examples": bool(found)}))
    elif cmd == "estimate":
        v, b, g2, pr = parse(sys.argv[2:])
        print(json.dumps(estimate(v, b, g2, pr)), flush=True)
