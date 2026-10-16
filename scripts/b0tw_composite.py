"""b = 0 (unbordered Goethals-Seidel) multiplier and twisted slices at the composite open orders v = p*q, on
G = Z_p x Z_q (= Z_v by CRT), with u = (a, mu) acting on both factors (scripts/b1tw_ab.py bit layout).

b = 0 target: sum_k PAF_k = 4v delta; block sums s_k odd with sum s_k^2 = 4v; 'BC' needs s_B = s_C, 'BCD' needs
s_B = s_C = s_D (signs are free per block: PAF is invariant under negation, so only |s_k| matters up to labelled
sign choices, which the DP counts).  Per slice: exact row-sum patterns (DP over orbit weights), exact quotient screens
onto Z_p and Z_q (compressed blocks: sum_k PAF_{y_k}(j) = 4v [j = 0] on the quotient, all lags), and the Gaussian
summed over the admissible row-sum patterns (one class dropped: its value is implied by the row sums).
usage: b0tw_composite.py OUT v p q [--rmax R] [--rmin R]"""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "theory"))
from abgroup import AbGroup
from b1tw_ab import AbSlice, cosets, letters
from b1tw_lib import count_signed, SignedSampler


def units(n): return [a for a in range(1, n) if math.gcd(a, n) == 1]


def order(a, n):
    k, x = 1, a % n
    while x != 1: x = x * a % n; k += 1
    return k


def sum_dist(w):
    """{sum: count} over all eps in {+-1}^n of sum w_i eps_i"""
    D = {0: 1}
    for x in w:
        E = {}
        for t, c in D.items():
            E[t + x] = E.get(t + x, 0) + c; E[t - x] = E.get(t - x, 0) + c
        D = E
    return D


def patterns(sl, v):
    """admissible labelled row-sum patterns (dict group -> sum of its representative block) with their counts"""
    gs = list(sl.groups)
    mult = {g: ({"A": 1, "D": 1, "C": 1, "B": {"id": 1, "BC": 2, "BCD": 3}[sl.tw]}[g]) for g in gs}
    dists = {g: {t: c for t, c in sum_dist(sl.gw[g][1]).items() if t % 2 and t * t * mult[g] <= 4 * v} for g in gs}
    out = []
    def rec(i, rem, cur, cnt):
        if i == len(gs):
            if rem == 0: out.append((dict(cur), cnt))
            return
        g = gs[i]
        for t, c in dists[g].items():
            e = mult[g] * t * t
            if e <= rem:
                cur[g] = t; rec(i + 1, rem - e, cur, cnt * c)
        cur.pop(g, None)
    rec(0, 4 * v, {}, 1)
    return out


def quotient_ok(sl, gens, v, pats, maxvec=300_000):
    G = sl.G; lab, creps, K = cosets(G, gens); Fn = len(creps); n = sl.n
    ucos = np.array([lab[sl.u[creps[b]]] for b in range(Fn)])
    addc = np.array([[lab[G.idx(G.coords[creps[a]] + G.coords[creps[b]])] for b in range(Fn)] for a in range(Fn)])
    tgt = np.zeros(Fn, np.int64); tgt[lab[0]] = 4 * v
    vpow = {"id": 1, "BC": 2, "BCD": 3}[sl.tw]
    contrib = {}
    for g in sl.groups:
        rep = {"A": 0, "B": 1, "C": 2, "D": 3}[g]; vp = vpow if g == "B" else 1
        vc = np.arange(Fn)
        for _ in range(vp): vc = ucos[vc]
        cl = -np.ones(Fn, np.int64); cells = []
        for b in range(Fn):
            if cl[b] >= 0: continue
            o, y = [b], int(vc[b])
            while y != b: o.append(y); y = int(vc[y])
            cl[o] = len(cells); cells.append(o)
        per = {}
        for i in sl.groups[g]:
            pos = np.flatnonzero(sl.idx[rep] == i); c = int(cl[lab[pos[0]]]); per.setdefault(c, []).append(len(pos) // len(cells[c]))
        cids = sorted(per); Ls = [letters(per[c]) for c in cids]
        tot = 1
        for x in Ls: tot *= len(x)
        if tot > maxvec: return None
        Y = np.zeros((tot, Fn), np.int64)
        grid = np.array(list(itertools.product(*Ls)), np.int64)
        for j, c in enumerate(cids): Y[:, cells[c]] = grid[:, [j]]
        blocks = [Y]
        if g == "B" and sl.tw in ("BC", "BCD"):
            uinv = np.argsort(ucos); blocks.append(Y[:, uinv])
            if sl.tw == "BCD": blocks.append(Y[:, uinv[uinv]])
        P = np.zeros((tot, Fn), np.int64)
        for Z in blocks:
            for j in range(Fn): P[:, j] += (Z * Z[:, addc[:, j]]).sum(1)
        contrib[g] = (Y.sum(1), P)
    W = np.random.default_rng(3).integers(1, 2 ** 62, Fn).astype(np.int64); gs = list(contrib)
    for t, _ in pats:
        sel = {g: contrib[g][1][contrib[g][0] == t[g]] for g in gs}
        if any(len(x) == 0 for x in sel.values()): continue
        acc = np.unique(sel[gs[0]] @ W)
        for g in gs[1:-1]:
            b = np.unique(sel[g] @ W)
            if len(acc) * len(b) > 50_000_000: return None
            acc = np.unique((acc[:, None] + b[None, :]).ravel())
        if np.isin(np.unique((tgt - sel[gs[-1]]) @ W), acc).any(): return True
    return False


def gaussian(sl, v, pats, N=4000, seed=1, top=40):
    """Gaussian over the `top` most numerous row-sum patterns (coverage = their share of all labelled patterns)"""
    rng = np.random.default_rng(seed); tot = 0.0; per = []
    drop = int(np.argmin(sl.csize)); keep = [c for c in range(sl.L) if c != drop]
    ps = sorted(pats, key=lambda t: -t[1]); allc = sum(c for _, c in ps); ps = ps[:top]
    cov = sum(c for _, c in ps) / allc
    for t, cnt in ps:
        beta = np.empty((N, sl.nbits), np.int8)
        for g, bits in sl.groups.items(): beta[:, bits] = SignedSampler(sl.gw[g][1], t[g]).sample(N, rng)
        F = sl.paf_classes(beta)[:, keep]
        steps = np.array([max(1, int(np.gcd.reduce(np.abs(F[:, c] - F[0, c])))) for c in range(len(keep))])
        if not all(((0 - F[0, c]) % steps[c]) == 0 for c in range(len(keep))): per.append((t, None)); continue
        mu = F.mean(0); S = np.atleast_2d(np.cov(F.T)); sign, logdet = np.linalg.slogdet(S)
        if sign <= 0: per.append((t, None)); continue
        q = mu @ np.linalg.solve(S, mu)
        lp = (np.log(steps).sum() - 0.5 * len(keep) * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
        tot += cnt * 2.0 ** lp; per.append((t, round(math.log2(cnt) + lp, 1)))
    return (math.log2(tot) if tot > 0 else None), per, cov


def main():
    a = sys.argv[1:]; OUT, v, p, q = a[0], int(a[1]), int(a[2]), int(a[3])
    opt = lambda f, d: int(a[a.index(f) + 1]) if f in a else d
    rmin, rmax = opt("--rmin", 6), opt("--rmax", 120)
    G = AbGroup([("C", p), ("C", q)]); assert G.n == v
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT): z = json.loads(l); done.add((z["v"], z["a"], z["mu"], z["tw"]))
    out = open(OUT, "a"); seen = set()
    for x in units(p):
        for mu in units(q):
            up = G.aut_perm([x, mu])
            # subgroup key
            S = set(); y = np.arange(v)
            d = 1; cur = up.copy()
            while not np.array_equal(cur, np.arange(v)): cur = up[cur]; d += 1
            key = frozenset(tuple(((pow(x, j, p)), pow(mu, j, q))) for j in range(d))
            if key in seen or d == 1: continue
            seen.add(key)
            for tw in ("id", "BC", "BCD"):
                if (tw == "BC" and d % 2) or (tw == "BCD" and d % 3) or (v, x, mu, tw) in done: continue
                sl = AbSlice(G, up, tw)
                if not (rmin <= sl.r <= rmax): continue
                t0 = time.time()
                rec = dict(v=v, a=x, mu=mu, d=d, tw=tw, orbits=sl.r, classes=sl.L, minus1_in_U=bool(sl.minus1_in_U),
                           group_bits={g: len(b) for g, b in sl.groups.items()})
                pats = patterns(sl, v)
                rec["patterns"] = len(pats)
                if not pats:
                    rec.update(alive=False, killed="row sums")
                else:
                    killed = None
                    for name, gens in (("Z_q", [int(G.idx([1, 0]))]), ("Z_p", [int(G.idx([0, 1]))])):
                        # quotient by the component generated: gens [1,0] generates Z_p, so G/H = Z_q, and vice versa
                        ok = quotient_ok(sl, gens, v, pats)
                        if ok is False: killed = name; break
                    rec.update(alive=killed is None, killed=killed)
                    if killed is None:
                        print(json.dumps(dict(stage="alive-before-gaussian", v=v, a=x, mu=mu, d=d, tw=tw, r=sl.r, L=sl.L,
                                              patterns=len(pats))), flush=True)
                        lE, per, cov = gaussian(sl, v, pats)
                        rec["pattern_coverage"] = round(cov, 3)
                        rec.update(log2E=None if lE is None else round(lE, 1),
                                   log2_configs=round(math.log2(sum(c for _, c in pats)), 1),
                                   n=None if lE is None else round(math.log2(sum(c for _, c in pats)) - lE, 1))
                rec["secs"] = round(time.time() - t0, 1)
                out.write(json.dumps(rec) + "\n"); out.flush(); print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
