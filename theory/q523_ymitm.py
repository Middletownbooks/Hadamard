r"""Stage 1 by a Wagner 4-list join, for the y-screen (report/paley_523_answer.md, section 5) [M for every y printed].

Stage 1: y_k in L^F (k = 0..3) with sum_k PAF(y_k)(h) = -4q (h = 1..F/2), a fixed lag-0 composition and rows
(2, 0, 0, 0).  Each block gets a list of N random sequences with its own composition and row sum (so lag 0 and the rows
hold by construction).  Join: pick a middle value c on a set M of lags; left pairs (a, b) with PAF_a + PAF_b = c on M,
right pairs (c', d) with PAF_c' + PAF_d = T - c on M (hash look-ups); then the two pair lists are matched on the
remaining lags by a 64-bit linear hash, and every hit is re-verified exactly.  With N = 2^22 and |M| = 3 lags each
middle value costs about 2^24 operations and yields about one solution at 523 [E; measured counts are logged].
usage: python3 theory/q523_ymitm.py F q s --comp a1:n1,... [--N 22] [--mid 3] [--cents 50] [--seed 0] [--out DIR]
       [--secs 600] [--maxper 20] [--mode lags|psd]"""
import sys, os, json, time, math
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_stage1 import letters
from q523_ysa import canonical, verify

def sign_table(mags, rowsum):
    """all sign patterns of the |y| > 1 entries with the number of +-1 completions to the row sum (exact weights)"""
    big = [int(x) for x in mags if x > 1]; n1 = len(mags) - len(big); pats = []; wts = []; kps = []
    for m in range(1 << len(big)):
        sg = [big[i] if (m >> i) & 1 else -big[i] for i in range(len(big))]
        need = rowsum - sum(sg)
        if (need + n1) % 2 or abs(need) > n1: continue
        kp = (need + n1) // 2; pats.append(sg); wts.append(math.comb(n1, kp)); kps.append(kp)
    if not pats: return None
    w = np.array(wts, float)
    return (np.array(pats, np.int64).reshape(len(pats), len(big)), np.cumsum(w) / w.sum(), np.array(kps, np.int64),
            math.log2(w.sum()))

@njit(cache=True)
def gen_block(pats, cum, kps, F, N, H, seed):
    """N uniform random sequences with the given |y| multiset and row sum (sign pattern of the large entries drawn
    with its exact completion count, a random subset of the units positive, positions shuffled); (seqs, PAF[:, 1..H])"""
    np.random.seed(seed)
    X = np.zeros((N, F), np.int16); P = np.zeros((N, H), np.int32)
    nb = pats.shape[1]; n1 = F - nb; y = np.zeros(F, np.int64)
    for n in range(N):
        t = np.searchsorted(cum, np.random.random()); kp = kps[t]
        for i in range(nb): y[i] = pats[t, i]
        for i in range(n1): y[nb + i] = 1 if i < kp else -1
        for i in range(F - 1, 0, -1):
            j = np.random.randint(i + 1); u = y[i]; y[i] = y[j]; y[j] = u
        for i in range(F): X[n, i] = y[i]
        for h in range(1, H + 1):
            acc = 0
            for g in range(F): acc += y[g] * y[(g + h) % F]
            P[n, h - 1] = acc
    return X, P

@njit(cache=True)
def key3(P, idx, c, sgn, cls, xc):
    """packed key of (sgn * P[:, idx] + c) on the middle lags (18 bits each) and the residue class (cls ^ xc, 6 bits)"""
    n = P.shape[0]; K = np.empty(n, np.int64)
    for i in range(n):
        k = 0
        for t in range(idx.shape[0]):
            k = (k << 18) | ((sgn * P[i, idx[t]] + c[t] + (1 << 17)) & ((1 << 18) - 1))
        K[i] = (k << 6) | ((cls[i] ^ xc) & 63)
    return K

def rest_class(P, rest):
    """residue class mod 8 of the rest lags: bit t = 1 iff P[:, rest[t]] = 6 mod 8 (every block PAF is 2 mod 4)"""
    c = np.zeros(P.shape[0], np.int64)
    for t, h in enumerate(rest): c |= ((P[:, h] % 8) == 6).astype(np.int64) << t
    return c

def xor_conv(h1, h2):
    out = np.zeros(64)
    for u in range(64):
        for v in range(64): out[u ^ v] += h1[u] * h2[v]
    return out

@njit(cache=True)
def pair_join(KA, ordA, KB, PA, PB, rest, w, tgt_rest, sgn_out, maxp):
    """pairs (a, b) with KA[a] == KB[b]; returns (a, b, hash of sgn_out*(tgt_rest - PA - PB) on rest lags)"""
    nA = KA.shape[0]; outa = np.empty(maxp, np.int64); outb = np.empty(maxp, np.int64); outh = np.empty(maxp, np.uint64)
    sk = KA[ordA]; m = 0
    for b in range(KB.shape[0]):
        k = KB[b]
        lo = np.searchsorted(sk, k);
        while lo < nA and sk[lo] == k:
            a = ordA[lo]
            if m < maxp:
                hsh = np.uint64(0)
                for t in range(rest.shape[0]):
                    v = PA[a, rest[t]] + PB[b, rest[t]]
                    if sgn_out < 0: v = tgt_rest[t] - v
                    hsh += np.uint64(v + (1 << 20)) * w[t]
                outa[m] = a; outb[m] = b; outh[m] = hsh
            m += 1; lo += 1
    mm = min(m, maxp)
    return outa[:mm], outb[:mm], outh[:mm], m

def choose_c(PA, PB, PC, PD, cl, x, mid, T, rng, K=20000):
    """middle value realised on both sides with the residue classes fixed (left pairs of rest-class x, right pairs of
    rest-class x ^ 63): an exact collision on the middle lags between K random left and K random right pairs of those
    classes, nearest the centre of the two samples"""
    def side(P1, P2, c1, c2, xx, sgn):
        bucket = {}
        for i, cc in enumerate(c2.tolist()): bucket.setdefault(cc, []).append(i)
        bucket = {k: np.array(v) for k, v in bucket.items()}
        ia = rng.integers(len(P1), size=K); out = []
        for i in ia:
            bk = bucket.get(int(c1[i]) ^ xx)
            if bk is None: continue
            j = bk[rng.integers(len(bk))]; out.append(P1[i, mid].astype(np.int64) + P2[j, mid])
        V = np.array(out).reshape(-1, len(mid))
        return V if sgn > 0 else T[mid] - V
    Lv = side(PA, PB, cl[0], cl[1], x, 1); Rv = side(PC, PD, cl[2], cl[3], x ^ 63, -1)
    if len(Lv) == 0 or len(Rv) == 0: return None
    Rs = set(map(tuple, Rv.tolist())); ok = np.array([tuple(r) in Rs for r in Lv.tolist()])
    if not ok.any(): return None
    centre = (Lv.mean(0) + Rv.mean(0)) / 2; sd = np.sqrt((Lv.var(0) + Rv.var(0)) / 2) + 1
    d = (((Lv - centre) / sd) ** 2).sum(1) + rng.exponential(1.0, len(Lv)); d[~ok] = np.inf
    return Lv[int(np.argmin(d))].astype(np.int64)

def psd_feats(X):
    """integer PSD values of each sequence at the frequencies j = 9, 6, 3 of Z_18 (orders 2, 3, 6): A^2 with A the
    alternating sum, and |sum_r s_r w^r|^2 = (X^2 + 3 Y^2) / 4 for the Z_3 and Z_6 class sums (exact integers)"""
    X = X.astype(np.int64); F = X.shape[1]; g = np.arange(F)
    A = (X * np.where(g % 2 == 0, 1, -1)).sum(1)
    s3 = np.stack([X[:, g % 3 == r].sum(1) for r in range(3)], 1)
    x3 = 2 * s3[:, 0] - s3[:, 1] - s3[:, 2]; y3 = s3[:, 1] - s3[:, 2]
    s6 = np.stack([X[:, g % 6 == r].sum(1) for r in range(6)], 1)
    x6 = s6 @ np.array([2, 1, -1, -2, -1, 1]); y6 = s6 @ np.array([0, 1, 1, 0, -1, -1])
    return np.stack([A * A, (x3 * x3 + 3 * y3 * y3) // 4, (x6 * x6 + 3 * y6 * y6) // 4], 1)

def all_class(P):
    """residue class mod 8 of all lags 1..H (bit t = 1 iff PAF(t + 1) = 6 mod 8)"""
    c = np.zeros(P.shape[0], np.int64)
    for t in range(P.shape[1]): c |= ((P[:, t] % 8) == 6).astype(np.int64) << t
    return c

@njit(cache=True)
def lin_key(M, c, sgn, cls, xc, w):
    """64-bit linear hash of (sgn * M + c, cls ^ xc) per row (exact equality is re-checked downstream)"""
    n = M.shape[0]; K = np.empty(n, np.int64)
    for i in range(n):
        h = np.uint64(0)
        for t in range(M.shape[1]): h += np.uint64(sgn * M[i, t] + c[t] + (1 << 40)) * w[t]
        h += np.uint64((cls[i] ^ xc) + 1) * w[M.shape[1]]
        K[i] = np.int64(h >> np.uint64(1))
    return K

def xor_conv_n(h1, h2, nb):
    """XOR convolution of two class histograms over nb bits (via Walsh-Hadamard transforms)"""
    def wht(a):
        a = a.copy(); hlen = 1
        while hlen < len(a):
            for i in range(0, len(a), 2 * hlen):
                x = a[i:i + hlen].copy(); y = a[i + hlen:i + 2 * hlen].copy(); a[i:i + hlen] = x + y; a[i + hlen:i + 2 * hlen] = x - y
            hlen *= 2
        return a
    return wht(wht(h1) * wht(h2)) / len(h1)

def run_psd(F, q, s, comp, N_log2=20, cents=50, seed=0, out="", secs=600.0, log=print, maxper=20):
    """Variant: middle value = the integer PSD sums at j = 9, 6, 3 together with the mod-8 class of all lags (the left
    pair class x, the right x ^ all-ones); second level: all lags by a 64-bit hash; every hit verified exactly."""
    rng = np.random.default_rng(seed); H = F // 2; tgt = -4 * q; T = np.full(H, tgt, np.int64)
    PT = (4 * F * q - 4 * (q - 1)) + 4 * q                                    # sum_k PSD_k(j) for j != 0: T(0) - T(h)
    pool = []
    for a, n in comp.items(): pool += [a] * n
    assert len(pool) == 4 * F
    found = {}; t0 = time.time(); N = 1 << N_log2; rows = (2, 0, 0, 0); ones = (1 << H) - 1
    allr = np.arange(H, dtype=np.int64)
    while time.time() - t0 < secs:
        while True:
            pool_ = np.array(pool); rng.shuffle(pool_); mags = pool_.reshape(4, F)
            tabs = [sign_table(mags[k], rows[k]) for k in range(4)]
            if all(t is not None for t in tabs): break
        bc = [{int(a): int((np.abs(mags[k]) == a).sum()) for a in sorted(comp) if a > 1} for k in range(4)]
        lists = [gen_block(tabs[k][0], tabs[k][1], tabs[k][2], F, N, H, int(rng.integers(1 << 31))) for k in range(4)]
        P = [l[1] for l in lists]; M = [psd_feats(l[0]) for l in lists]; cl = [all_class(p) for p in P]
        hist = [np.bincount(c, minlength=1 << H).astype(float) for c in cl]
        Lh = xor_conv_n(hist[0], hist[1], H); Rh = xor_conv_n(hist[2], hist[3], H)
        score = np.array([Lh[x] * Rh[x ^ ones] for x in range(1 << H)])
        if score.max() <= 0.5:
            log(json.dumps(dict(block_comps=bc, note="residue classes incompatible mod 8"))); continue
        xs = np.argsort(-score)[:16]; xs = xs[score[xs] > 0.5]
        w = rng.integers(1, 2 ** 63, 4, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        wr = rng.integers(1, 2 ** 63, H, dtype=np.int64).astype(np.uint64) | np.uint64(1)
        for ci in range(cents):
            if time.time() - t0 > secs: break
            x = int(xs[ci % len(xs)])
            # middle value: exact collision of left PSD sums (class x pairs) with right needs PT - PSD_C - PSD_D
            def side(M1, M2, c1, c2, xx, n=20000):
                bucket = {}
                for i, cc in enumerate(c2.tolist()): bucket.setdefault(cc, []).append(i)
                ia = rng.integers(len(M1), size=n); V = []
                for i in ia:
                    bk = bucket.get(int(c1[i]) ^ xx)
                    if bk: V.append(M1[i] + M2[bk[rng.integers(len(bk))]])
                return np.array(V).reshape(-1, 3)
            Lv = side(M[0], M[1], cl[0], cl[1], x); Rv = PT - side(M[2], M[3], cl[2], cl[3], x ^ ones)
            if len(Lv) == 0 or len(Rv) == 0:
                log(json.dumps(dict(centre=ci, cls=x, note="empty class side"))); continue
            Rs = set(map(tuple, Rv.tolist())); ok = np.array([tuple(r) in Rs for r in Lv.tolist()])
            if not ok.any():
                log(json.dumps(dict(centre=ci, block_comps=bc, cls=x, note="no PSD collision"))); continue
            cands = Lv[ok]; c = cands[rng.integers(len(cands))].astype(np.int64)
            z = np.zeros(3, np.int64)
            KA = lin_key(M[0], z, 1, cl[0], 0, w); ordA = np.argsort(KA, kind='stable')
            KB = lin_key(M[1], c, -1, cl[1], x, w)
            a_, b_, hl, nl = pair_join(KA, ordA, KB, P[0], P[1], allr, wr, T, 1, 1 << 23)
            KC = lin_key(M[2], z, 1, cl[2], 0, w); ordC = np.argsort(KC, kind='stable')
            KD = lin_key(M[3], PT - c, -1, cl[3], x ^ ones, w)
            c_, d_, hr, nr = pair_join(KC, ordC, KD, P[2], P[3], allr, wr, T, -1, 1 << 23)
            ol = np.argsort(hl); orr = np.argsort(hr); hs = hl[ol]; hrs = hr[orr]
            i = np.searchsorted(hrs, hs, side='left'); j = np.searchsorted(hrs, hs, side='right')
            nsol = 0; new = 0; hit = np.nonzero(j > i)[0]; nhit = int((j - i)[hit].sum()); rng.shuffle(hit)
            for t in hit:
                if new >= maxper: break
                for u in range(i[t], j[t]):
                    a, b = a_[ol[t]], b_[ol[t]]; cc, dd = c_[orr[u]], d_[orr[u]]
                    Y = [lists[0][0][a].tolist(), lists[1][0][b].tolist(), lists[2][0][cc].tolist(), lists[3][0][dd].tolist()]
                    if not verify(Y, F, q): continue
                    nsol += 1; Yc = canonical(Y); key = json.dumps(Yc)
                    if key in found: continue
                    found[key] = Yc; new += 1
                    if out:
                        tag = "_".join(f"{a}-{n}" for a, n in sorted(comp.items()))
                        fn = os.path.join(out, f"y_mitm_{tag}_s{seed}_{len(found)}.json")
                        json.dump(dict(F=F, q=q, s=s, comp={str(a): n for a, n in comp.items()}, y=Yc, verified=True,
                                       method="q523_ymitm psd", seed=seed), open(fn, "w"))
            log(json.dumps(dict(centre=ci, block_comps=bc, cls=x, c=c.tolist(), left_pairs=int(nl), right_pairs=int(nr),
                                hash_hits=nhit, verified=nsol, new=new, total=len(found), secs=round(time.time() - t0, 1))))
    return found

def run(F, q, s, comp, N_log2=22, nmid=3, cents=50, seed=0, out="", secs=600.0, log=print, maxper=20):
    rng = np.random.default_rng(seed); H = F // 2; tgt = -4 * q; T = np.full(H, tgt, np.int64)
    L, f = letters(q, s)
    # split the global composition over the 4 blocks at random (each block F entries)
    pool = []
    for a, n in comp.items(): pool += [a] * n
    assert len(pool) == 4 * F
    found = {}; t0 = time.time(); N = 1 << N_log2
    rows = (2, 0, 0, 0)
    while time.time() - t0 < secs:
        # per-block split of the composition, accepted with probability ~ (number of sequences with the row sums)
        while True:
            pool_ = np.array(pool); rng.shuffle(pool_); mags = pool_.reshape(4, F)
            tabs = [sign_table(mags[k], rows[k]) for k in range(4)]
            if all(t is not None for t in tabs): break
        bc = [{int(a): int((np.abs(mags[k]) == a).sum()) for a in sorted(comp) if a > 1} for k in range(4)]
        lists = [gen_block(tabs[k][0], tabs[k][1], tabs[k][2], F, N, H, int(rng.integers(1 << 31))) for k in range(4)]
        tg = time.time() - t0
        # middle lags: nmid lags chosen at random among 1..H-1 (lag H = F/2 is doubled, keep it in the rest)
        mid = np.sort(rng.choice(H - 1, nmid, replace=False)).astype(np.int64)
        rest = np.array([h for h in range(H) if h not in set(mid.tolist())], np.int64)
        w = rng.integers(1, 2 ** 63, len(rest), dtype=np.int64).astype(np.uint64) | np.uint64(1)
        PA, PB, PC, PD = [l[1] for l in lists]
        # residue classes mod 8 of the rest lags: the four block classes must XOR to all-ones (sum = 4 mod 8 per lag);
        # the left pair class x is chosen to maximise (#left pairs of class x) * (#right pairs of class x ^ 63)
        cl = [rest_class(P, rest) for P in (PA, PB, PC, PD)]
        hist = [np.bincount(c, minlength=64).astype(float) for c in cl]
        Lh = xor_conv(hist[0], hist[1]); Rh = xor_conv(hist[2], hist[3])
        score = np.array([Lh[x] * Rh[x ^ 63] for x in range(64)])
        if score.max() == 0:
            log(json.dumps(dict(block_comps=bc, note="residue classes incompatible mod 8"))); continue
        xs = np.argsort(-score)[:8]; xs = xs[score[xs] > 0]
        for ci in range(cents):
            if time.time() - t0 > secs: break
            x = int(xs[ci % len(xs)])
            c = choose_c(PA, PB, PC, PD, cl, x, mid, T, rng)
            if c is None:
                log(json.dumps(dict(centre=ci, block_comps=bc, cls=x, note="no compatible middle value"))); continue
            zero = np.zeros(nmid, np.int64)
            KA = key3(PA, mid, zero, 1, cl[0], 0); ordA = np.argsort(KA, kind='stable')
            KB = key3(PB, mid, c, -1, cl[1], x)                            # PA = c - PB, clsA = clsB ^ x
            a_, b_, hl, nl = pair_join(KA, ordA, KB, PA, PB, rest, w, T[rest], 1, 1 << 23)
            KC = key3(PC, mid, zero, 1, cl[2], 0); ordC = np.argsort(KC, kind='stable')
            KD = key3(PD, mid, T[mid] - c, -1, cl[3], x ^ 63)               # PC = (T - c) - PD, clsC = clsD ^ x ^ 63
            c_, d_, hr, nr = pair_join(KC, ordC, KD, PC, PD, rest, w, T[rest], -1, 1 << 23)
            ol = np.argsort(hl); orr = np.argsort(hr)
            hs = hl[ol]; hrs = hr[orr]
            i = np.searchsorted(hrs, hs, side='left'); j = np.searchsorted(hrs, hs, side='right')
            nsol = 0; new = 0; hit = np.nonzero(j > i)[0]; nhit = int((j - i)[hit].sum()); rng.shuffle(hit)
            for t in hit:
                if new >= maxper: break
                for u in range(i[t], j[t]):
                    a, b = a_[ol[t]], b_[ol[t]]; cc, dd = c_[orr[u]], d_[orr[u]]
                    Y = [lists[0][0][a].tolist(), lists[1][0][b].tolist(), lists[2][0][cc].tolist(), lists[3][0][dd].tolist()]
                    if not verify(Y, F, q): continue
                    nsol += 1; Yc = canonical(Y); key = json.dumps(Yc)
                    if key in found: continue
                    found[key] = Yc; new += 1
                    if out:
                        tag = "_".join(f"{a}-{n}" for a, n in sorted(comp.items()))
                        fn = os.path.join(out, f"y_mitm_{tag}_s{seed}_{len(found)}.json")
                        json.dump(dict(F=F, q=q, s=s, comp={str(a): n for a, n in comp.items()}, y=Yc, verified=True,
                                       method="q523_ymitm", seed=seed), open(fn, "w"))
            log(json.dumps(dict(centre=ci, block_comps=bc, mid=(mid + 1).tolist(), cls=x, c=c.tolist(), left_pairs=int(nl),
                                right_pairs=int(nr), hash_hits=nhit, verified=nsol, new=new, total=len(found), gen_secs=round(tg, 1),
                                secs=round(time.time() - t0, 1))))
    return found

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "22", "--mid": "3", "--cents": "50", "--seed": "0", "--comp": "", "--out": "",
                             "--secs": "600", "--maxper": "20", "--mode": "lags"}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    comp = {int(t.split(':')[0]): int(t.split(':')[1]) for t in opt["--comp"].split(',')}
    if opt["--out"]: os.makedirs(opt["--out"], exist_ok=True)
    if opt["--mode"] == "psd":
        run_psd(F, q, s, comp, int(opt["--N"]), int(opt["--cents"]), int(opt["--seed"]), opt["--out"], float(opt["--secs"]),
                log=lambda x: print(x, flush=True), maxper=int(opt["--maxper"]))
        sys.exit()
    run(F, q, s, comp, int(opt["--N"]), int(opt["--mid"]), int(opt["--cents"]), int(opt["--seed"]), opt["--out"],
        float(opt["--secs"]), log=lambda x: print(x, flush=True), maxper=int(opt["--maxper"]))
