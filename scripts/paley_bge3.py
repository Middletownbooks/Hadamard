"""Paley-coset ansatz for the bordered b >= 3 architecture (the form the published v = 359 core actually has).

Architecture (report/bge3_construction.md): v = c(h + 1) - 1, b = c - 1, m = v - b = c*h, G of order m, H <= G of
index c.  Core: GS quadruple over G with sum_k PAF_k = 4m at 0, -4b on H \ {0}, +4 off H.
Ansatz: G = G' x Z_p with p | h, p not dividing c (so Z_p <= H and H = H' x Z_p, H' <= G' of index c), G' abelian
(any invariant factors).  Per coset g' of Z_p: x(g', 0) = alpha, x(g', z) = beta*chi(z) (twisted) or gamma (constant).
Lag (h', s) lies in H iff h' in H', so the target is -4b on (H' x {0,QR,NQR}) \ (0,0) and +4 elsewhere.

Screens (all exact necessary conditions):
  h0      within-coset equations + row sums (paley_ansatz.exact_h0_screen_fast; sum r^2 = 4v for b >= 3)
  qnorm   compression onto Q = G'/H': per block the additive invariants (r, N = sum_q y_q^2, s, d), with
          y_q = coset-q sum.  Needs sum_k r_k^2 = 4v, sum_k N_k = 4(h + b), sum s = -4b, sum d = 0.
Estimate: row-sum-conditioned Gaussian on the reduced (h', class) equations.

usage: paley_bge3.py control                      (359 core: exactly in the ansatz, reduced equations = target)
       paley_bge3.py options v [v ...]            (list (c, p, G', Q) options)
       paley_bge3.py screen v [v ...]             (h0 + qnorm per option, JSON lines)
       paley_bge3.py screen1 v c p                (one option)
       paley_bge3.py estimate v c p mods...       (Gaussian estimate for one option, all index-c subgroups classes)"""
import sys, os, json, itertools as it, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_ansatz import legendre, corr_table, LABELS, LC, exact_h0_screen_fast

class AG:
    """finite abelian group Z_{n1} x ... x Z_{nk}; elements are ravel indices (last factor fastest)"""
    def __init__(s, mods):
        s.mods = np.array(mods, np.int64); s.n = int(np.prod(mods))
        s.co = np.array(np.unravel_index(np.arange(s.n), tuple(mods))).T
        s.add = s.idx(s.co[:, None, :] + s.co[None, :, :])                  # add[g, h] = g + h
    def idx(s, c):
        c = np.asarray(c) % s.mods
        return np.ravel_multi_index(tuple(np.moveaxis(c, -1, 0)), tuple(s.mods))

def factor(n):
    f, q = {}, 2
    while q * q <= n:
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield []; return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k): yield [k] + rest

def abelian_groups(n):
    """invariant-factor-free form: list of prime-power cyclic factors"""
    per = [[[q ** e for e in lam] for lam in partitions(k)] for q, k in factor(n).items()]
    for combo in it.product(*per):
        yield sorted([x for part in combo for x in part], reverse=True)

def index_c_subgroups(G, c):
    """kernels of surjections G -> Q, |Q| = c, Q abelian; returns {kernel (frozenset): Q factors}"""
    out = {}
    for Q in abelian_groups(c):
        Qg = AG(Q)
        # generator images: for factor i of G (order n_i) an element of Q of order dividing n_i
        cand = [[e for e in range(Qg.n) if ((Qg.co[e] * n_i) % Qg.mods == 0).all()] for n_i in G.mods]
        for imgs in it.product(*cand):
            I = Qg.co[list(imgs)]                                            # (k, dimQ)
            img = Qg.idx(G.co @ I)
            if len(set(img.tolist())) != c: continue
            ker = frozenset(np.nonzero(img == 0)[0].tolist())
            out.setdefault(ker, (Q, img))
    return out

def ker_class(G, ker, Q):
    """cheap isomorphism-class key for (G', H') (element-order multiset of H', Q type)"""
    ords = []
    for g in ker:
        o, x = 1, g
        while x != 0: x = G.add[x, g]; o += 1
        ords.append(o)
    return (tuple(sorted(ords)), tuple(Q))

def reduced_paf_g(lab, G, C):
    A = LC[lab]
    R = np.zeros((lab.shape[0], G.n, 3), np.int64)
    for h in range(G.n):
        B = A[:, G.add[:, h]]
        R[:, h, :] = np.einsum('ngi,ngj,ijc->nc', A, B, C)
    return R

def target(G, Hmask, b, m):
    T = np.where(Hmask[:, None], -4 * b, 4) * np.ones((1, 3), np.int64); T[0, 0] = 4 * m
    return T

# ---------------- positive control: the published 359 core ----------------
def control():
    L = open('data/core_v359_b7.txt').read().split('\n')
    X = [np.array([1 if ch == '+' else -1 for ch in l.split('=')[1].strip()]) for l in L if '=' in l and l[0] in 'ABCD']
    mods, p, b = (2, 2, 8, 11), 11, 7; m = 352; v = m + b
    G = AG(mods[:-1]); chi = legendre(p)
    labs, ok_form = [], True
    for x in X:
        Y = x.reshape(G.n, p); lab = []
        for g in range(G.n):
            a, rest = Y[g, 0], Y[g, 1:]
            if (rest == rest[0]).all(): lab.append(LABELS.index((int(a), 0, int(rest[0]))))
            elif (rest == rest[0] * chi[1:] * chi[1]).all(): lab.append(LABELS.index((int(a), 1, int(rest[0] * chi[1]))))
            else: ok_form = False; lab.append(-1)
        labs.append(lab)
    lab = np.array(labs)
    C, nqr = corr_table(p)
    R = reduced_paf_g(lab, G, C).sum(0)
    P = sum(np.rint(np.fft.ifftn(np.abs(np.fft.fftn(x.reshape(mods))) ** 2).real).astype(int) for x in X)
    Hmask = np.array([P.reshape(G.n, p)[g, 0] == -4 * b for g in range(G.n)]); Hmask[0] = True
    T = target(G, Hmask, b, m)
    ntw = int(sum(LABELS[l][1] for l in lab.ravel()))
    print(json.dumps({"control": 359, "in_ansatz": bool(ok_form), "twisted_cosets": ntw, "constant_cosets": 128 - ntw,
                      "H'_size": int(Hmask.sum()), "reduced_equals_target": bool((R == T).all()),
                      "h0_screen": exact_h0_screen_fast(v, b, G.n, p)["h0_screen"], "qnorm": qnorm_screen(v, b, G.n, p, 8)}))

# ---------------- options ----------------
def options(v):
    for c in range(4, v, 2):
        if (v + 1) % c: continue
        h = (v + 1) // c - 1; b = c - 1; m = c * h
        if h < 3: continue
        for p in sorted(factor(h)):
            if p < 3 or factor(m)[p] != 1 or c % p == 0: continue
            n = m // p
            if n % c: continue                                               # H' of index c in G' needs c | n'
            yield dict(v=v, c=c, b=b, h=h, m=m, p=p, n=n)

# ---------------- qnorm screen (exact, compression onto Q) ----------------
def coset_states(k, p):
    """reachable (y, s, d) over k cosets of Z_p (one H'-coset of G'): per coset as in exact_h0_screen_fast"""
    p1 = (p % 4 == 1)
    per = set()
    for a in (1, -1):
        for bt in (1, -1): per.add((a, -1, a * bt if p1 else 0))                       # twisted
        for g in (1, -1): per.add((a + (p - 1) * g, p - 2 + 2 * a * g, 0))             # constant
    S = {(0, 0, 0)}
    for _ in range(k): S = {(x[0] + y[0], x[1] + y[1], x[2] + y[2]) for x in S for y in per}
    return S

def qnorm_screen(v, b, n, p, c, cap=4_000_000):
    k = n // c; h = (v - b) // c; Nmax = 4 * (h + b)
    cs = [t for t in coset_states(k, p) if t[0] ** 2 <= Nmax]
    # block DP over c cosets: (r, N, s, d)
    S = {(0, 0, 0, 0)}
    for _ in range(c):
        S = {(x[0] + y[0], x[1] + y[0] ** 2, x[2] + y[1], x[3] + y[2]) for x in S for y in cs if x[1] + y[0] ** 2 <= Nmax}
        if len(S) > cap: return None                                          # undecided (too many states)
    # combine 4 blocks: sum r^2 = 4v, sum N = Nmax, sum s = -4b, sum d = 0.  (s, d) sets as bitsets:
    # bit (s + n_off)*W + (d + n_off), W = 8n + 1, so sums of up to 4 blocks never carry between fields.
    from collections import defaultdict
    W = 8 * n + 1 if p % 4 == 1 else 1                                      # d = 0 identically when p = 3 mod 4
    K = defaultdict(int)
    for r, N, s_, d in S: K[(r * r, N)] |= 1 << ((s_ + n) * W + (d + n))
    def bits(x):
        out = []
        while x:
            low = x & -x; out.append(low.bit_length() - 1); x ^= low
        return out
    keys = sorted(K); P = defaultdict(int)
    for i, k1 in enumerate(keys):
        for k2 in keys[i:]:
            if k1[0] + k2[0] > 4 * v or k1[1] + k2[1] > Nmax: continue
            a, bb = K[k1], K[k2]
            if bin(a).count("1") < bin(bb).count("1"): a, bb = bb, a
            acc = 0
            for j in bits(bb): acc |= a << j
            P[(k1[0] + k2[0], k1[1] + k2[1])] |= acc
    T = (-4 * b + 4 * n) * W + 4 * n
    for (r2, N), x in P.items():
        y = P.get((4 * v - r2, Nmax - N))
        if y and any(T - j >= 0 and x >> (T - j) & 1 for j in bits(y)): return True
    return False

def screen(v):
    for o in options(v):
        h0 = exact_h0_screen_fast(o["v"], o["b"], o["n"], o["p"])["h0_screen"]
        qn = qnorm_screen(o["v"], o["b"], o["n"], o["p"], o["c"]) if h0 else None
        print(json.dumps(dict(o, h0=h0, qnorm=qn)), flush=True)

# ---------------- Gaussian estimate (row-sum conditioned) ----------------
def gauss_estimate(v, c, p, mods, N=40000, seed=1, maxprof=6):
    b = c - 1; m = v - b; G = AG(mods); n = G.n; assert n * p == m
    C, nqr = corr_table(p); rng = np.random.default_rng(seed)
    rvals = np.array([LABELS[i][0] + (p - 1) * (0 if LABELS[i][1] else LABELS[i][2]) for i in range(8)])
    from collections import Counter
    per = Counter(rvals.tolist()); cnt = Counter({0: 1})
    for _ in range(n):
        E = Counter()
        for x, cx in cnt.items():
            for y, cy in per.items(): E[x + y] += cx * cy
        cnt = E
    A_ = sorted(a for a in cnt if a >= 0)
    profs = [q for q in it.combinations_with_replacement(A_, 4) if sum(x * x for x in q) == 4 * v]
    profs = sorted(profs, key=lambda q: -sum(math.log2(cnt[x]) for x in q))[:maxprof]
    subs = index_c_subgroups(G, c); classes = {}
    for ker, (Q, img) in subs.items(): classes.setdefault(ker_class(G, ker, Q), ker)
    out = []
    stats = {}
    def samp(a):
        if a not in stats:
            res = []
            for _ in range(400):
                lab = rng.integers(0, 8, (40000, n)); lab = lab[rvals[lab].sum(1) == a]; res.append(lab)
                if sum(len(x) for x in res) >= N: break
            lab = np.concatenate(res)[:N]
            stats[a] = None if len(lab) < 1000 else reduced_paf_g(lab, G, C).reshape(len(lab), -1)
        return stats[a]
    for key, ker in classes.items():
        Hmask = np.zeros(n, bool); Hmask[list(ker)] = True
        T = target(G, Hmask, b, m).reshape(-1)
        best, tot = None, 0.0
        for q in profs:
            Rs = [samp(a) for a in q]
            if any(r is None for r in Rs): continue
            mu = sum(r.mean(0) for r in Rs); Cv = sum(np.cov(r.T) for r in Rs); Rc = sum(r[:1500] for r in Rs)
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
            lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum()
                  - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
            perms = len(set(it.permutations(q))); signs = 2 ** sum(1 for x in q if x != 0)
            lE = sum(math.log2(cnt[x]) for x in q) + lp + math.log2(perms * signs)
            tot += 2 ** lE
            if best is None or lE > best[1]: best = (list(q), round(lE, 1), round(-lp, 1), int(kw.sum()))
        rec = dict(v=v, c=c, b=b, p=p, Gp=list(mods), Q=list(key[1]), H_orders=sorted(set(key[0])), free_bits=12 * n)
        if best is None: rec.update(feasible=False)
        else: rec.update(log2E=round(math.log2(tot), 1), profile=best[0], constraint_bits=best[2], eqs=best[3])
        out.append(rec); print(json.dumps(rec), flush=True)
    return out

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "control": control()
    elif cmd == "options":
        for v in map(int, sys.argv[2:]):
            for o in options(v): print(json.dumps(o))
    elif cmd == "screen":
        for v in map(int, sys.argv[2:]): screen(v)
    elif cmd == "screen1":                                                   # one option: v c p
        v, c, p = map(int, sys.argv[2:5])
        for o in options(v):
            if o["c"] == c and o["p"] == p:
                h0 = exact_h0_screen_fast(o["v"], o["b"], o["n"], o["p"])["h0_screen"]
                qn = qnorm_screen(o["v"], o["b"], o["n"], o["p"], o["c"]) if h0 else None
                print(json.dumps(dict(o, h0=h0, qnorm=qn)), flush=True)
    elif cmd == "estimate":
        v, c, p = map(int, sys.argv[2:5]); gauss_estimate(v, c, p, tuple(map(int, sys.argv[5:])))
