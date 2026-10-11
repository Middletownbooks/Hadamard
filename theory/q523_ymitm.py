r"""Stage 1 by a Wagner 4-list join, for the y-screen (report/paley_523_answer.md, section 5) [M for every y printed].

Stage 1: y_k in L^F (k = 0..3) with sum_k PAF(y_k)(h) = -4q (h = 1..F/2), a fixed lag-0 composition and rows
(2, 0, 0, 0).  Each block gets a list of N random sequences with its own composition and row sum (so lag 0 and the rows
hold by construction).  Join: pick a middle value c on a set M of lags; left pairs (a, b) with PAF_a + PAF_b = c on M,
right pairs (c', d) with PAF_c' + PAF_d = T - c on M (hash look-ups); then the two pair lists are matched on the
remaining lags by a 64-bit linear hash, and every hit is re-verified exactly.  With N = 2^22 and |M| = 3 lags each
middle value costs about 2^24 operations and yields about one solution at 523 [E; measured counts are logged].
usage: python3 theory/q523_ymitm.py F q s --comp a1:n1,... [--N 22] [--mid 3] [--cents 50] [--seed 0] [--out DIR]
       [--secs 600] [--maxper 20]"""
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
def key3(P, idx, c, sgn):
    """packed key of sgn*(P[:, idx]) + c (three lags, 20 bits each)"""
    n = P.shape[0]; K = np.empty(n, np.int64)
    for i in range(n):
        k = 0
        for t in range(idx.shape[0]):
            k = (k << 20) | ((sgn * P[i, idx[t]] + c[t] + (1 << 19)) & ((1 << 20) - 1))
        K[i] = k
    return K

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

def choose_c(PA, PB, PC, PD, mid, T, rng, K=20000, mod=8):
    """a middle value realised on both sides AND residue-compatible on the other lags: an exact collision on the middle
    lags between K random left sums (PA + PB) and K random right values T - (PC + PD) whose remaining lags agree
    mod `mod` (block compositions tie the PAF residues of different lags together, so left and right pair lists can be
    incompatible mod 8 even when the middle lags match); the collision nearest the centre of the two samples"""
    ia, ib, ic, id_ = (rng.integers(len(X), size=K) for X in (PA, PB, PC, PD))
    L9 = PA[ia].astype(np.int64) + PB[ib]; R9 = T - (PC[ic].astype(np.int64) + PD[id_])
    rest = np.array([h for h in range(PA.shape[1]) if h not in set(mid.tolist())])
    keyR = set((tuple(r[mid]), tuple(r[rest] % mod)) for r in R9.tolist() for r in [np.array(r)])
    ok = np.array([(tuple(l[mid]), tuple(l[rest] % mod)) in keyR for l in L9])
    if not ok.any(): return None
    Lv = L9[:, mid]; Rv = R9[:, mid]
    centre = (Lv.mean(0) + Rv.mean(0)) / 2; sd = np.sqrt((Lv.var(0) + Rv.var(0)) / 2) + 1
    d = (((Lv - centre) / sd) ** 2).sum(1) + rng.exponential(1.0, K); d[~ok] = np.inf
    return Lv[int(np.argmin(d))].astype(np.int64)

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
        for ci in range(cents):
            if time.time() - t0 > secs: break
            # middle value: a realised left sum (PA + PB on the middle lags) whose residue class mod 16 is also realised
            # by T - (PC + PD), closest to the midpoint of the two means (block compositions can fix residues mod 8)
            c = choose_c(PA, PB, PC, PD, mid, T, rng)
            if c is None:
                log(json.dumps(dict(centre=ci, block_comps=bc, note="no residue-compatible middle value"))); break
            zero = np.zeros(nmid, np.int64)
            KA = key3(PA, mid, zero, 1); ordA = np.argsort(KA, kind='stable')
            KB = key3(PB, mid, c, -1)                                      # PA = c - PB
            a_, b_, hl, nl = pair_join(KA, ordA, KB, PA, PB, rest, w, T[rest], 1, 1 << 23)
            KC = key3(PC, mid, zero, 1); ordC = np.argsort(KC, kind='stable')
            KD = key3(PD, mid, T[mid] - c, -1)                              # PC = (T - c) - PD
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
            log(json.dumps(dict(centre=ci, block_comps=bc, mid=(mid + 1).tolist(), c=c.tolist(), left_pairs=int(nl),
                                right_pairs=int(nr), hash_hits=nhit, verified=nsol, new=new, total=len(found), gen_secs=round(tg, 1),
                                secs=round(time.time() - t0, 1))))
    return found

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "22", "--mid": "3", "--cents": "50", "--seed": "0", "--comp": "", "--out": "",
                             "--secs": "600", "--maxper": "20"}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    comp = {int(t.split(':')[0]): int(t.split(':')[1]) for t in opt["--comp"].split(',')}
    if opt["--out"]: os.makedirs(opt["--out"], exist_ok=True)
    run(F, q, s, comp, int(opt["--N"]), int(opt["--mid"]), int(opt["--cents"]), int(opt["--seed"]), opt["--out"],
        float(opt["--secs"]), log=lambda x: print(x, flush=True), maxper=int(opt["--maxper"]))
