r"""Stage 2 of a staged c = 4 slice with large fibres (the s = 4 Paley-type configurations), summed over all
stage-1 solutions [E], with an exact sub-level count on the same y's for calibration [M].

Input: the cfg and the raw int8 stage-1 file of theory/ab_s1_chars.py (4 x L1 orbit values per y, canonical).
  1. H2(y) = log2 #(lifts of y), exactly, for every y (letter -> fibre-size table).  H2 depends only on how many
     letters of each fibre size y has, so the y's fall into H2 strata; every stratum is counted exactly and
     sampled uniformly (smallest random keys).
  2. Per sampled y, from one set of N uniform lifts of y (uniform pattern per G'-orbit = uniform over the fibre):
     - the full stage-2 lattice Gaussian of sum_k PAF_G at the lag classes (theory/q523_gauss.gauss_bits, as in
       theory/ab_gauss.stage2_sampled): C2 = its bits, log2 E(cores | y) = H2 - C2;
     - the same Gaussian at the Z_q level (compression along G', lag classes of Z_q under K and -1).
  3. The exact Z_q-level count for the same y [M]: per block a DP over the G'-orbits of the compressed vector
     (exact integer counts), its distribution of Z_q-PAF signatures, then the 4-fold convolution
     (sum of products of exact counts, accumulated in float64), against H2 - (Z_q-level Gaussian bits).
  4. E(cores) = sum over strata of (stratum size) x (sample mean of 2^(H2 - C2)).
usage: python3 theory/ab_s2_s4.py CFG.json YBIN OUT.jsonl [--per 3] [--seed 1] [--N 12000] [--exact-zq]"""
import sys, os, json, time, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs
from q523_gauss import gauss_bits


def zq_setup(sl):
    q = sl.q
    Mo = np.zeros((sl.r, q), np.int64); np.add.at(Mo, (sl.orb, sl.zco), 1)         # orbit o -> points at each t
    K = sorted({int(u[sl.qaxis]) % q for u in sl.S})
    kor = -np.ones(q, np.int64); kreps = []
    for t in range(q):
        if kor[t] < 0:
            for k in K: kor[(t * k) % q] = len(kreps)
            kreps.append(t)
    lcl = -np.ones(q, np.int64); lreps = []
    for t in range(q):
        if lcl[t] < 0:
            for k in K:
                for e in (1, -1): lcl[(e * t * k) % q] = len(lreps)
            lreps.append(t)
    Tq = np.zeros(q, np.int64); np.add.at(Tq, sl.zco, sl.T)
    return dict(q=q, Mo=Mo, kor=kor, kreps=kreps, lreps=lreps, lsize=np.bincount(lcl), Tq=Tq[lreps])


def zq_paf(Xb, lreps):
    q = Xb.shape[1]
    return np.stack([(Xb * np.roll(Xb, -t, axis=1)).sum(1) for t in lreps], 1)


def block_zq_dist(sl, Z, pats_k):
    """exact distribution of the Z_q-compressed vector of one block over its fibre: (values at the K-orbit reps,
    counts as int64)"""
    kr = Z['kreps']; nk = len(kr); R = 2 * sl.Gp.n + 1
    st = np.zeros((1, nk), np.int64); cnt = np.ones(1, np.int64)
    for j, P in enumerate(pats_k):
        os_ = [o for o, _ in sl.parts[j]]
        contrib = P.astype(np.int64) @ Z['Mo'][os_][:, kr]                     # (patterns, nk)
        nst = (st[:, None, :] + contrib[None, :, :]).reshape(-1, nk)
        ncnt = np.repeat(cnt, len(P))
        key = ((nst + sl.Gp.n) * (R ** np.arange(nk))).sum(1)
        u, inv = np.unique(key, return_inverse=True)
        c2 = np.zeros(len(u), np.int64); np.add.at(c2, inv.ravel(), ncnt)
        first = np.zeros(len(u), np.int64); first[inv.ravel()[::-1]] = np.arange(len(inv.ravel()))[::-1]
        st, cnt = nst[first], c2
    return st, cnt


def zq_exact(sl, Z, pats):
    """exact number of lifts (all 4 blocks) satisfying sum_k PAF_{Z_q}(compressed) = T compressed; float64 sum"""
    q = Z['q']; lreps = Z['lreps']; nl = len(lreps)
    dists = []
    for k in range(4):
        st, cnt = block_zq_dist(sl, Z, pats[k])
        Xb = st[:, Z['kor']]                                                    # full Z_q vectors
        P = zq_paf(Xb, lreps)[:, :nl - 1]                                       # the last class is implied by the row sum
        u, inv = np.unique(P, axis=0, return_inverse=True)
        c = np.zeros(len(u)); np.add.at(c, inv.ravel(), cnt.astype(np.float64))
        dists.append((u, c))
    lo = np.min([d[0].min(0) for d in dists], 0) * 4; hi = np.max([d[0].max(0) for d in dists], 0) * 4
    T = Z['Tq'][:nl - 1]
    lo = np.minimum(lo, T); hi = np.maximum(hi, T); Rr = (hi - lo + 1)
    assert float(np.prod(Rr.astype(np.float64))) < 2 ** 62
    mul = np.cumprod(np.concatenate([[1], Rr[:-1]])).astype(np.int64)
    enc = lambda S, sh: ((S - sh) * mul).sum(-1)
    off = lo // 4                                                               # per-block shift (sum of 4 shifts = lo)

    def conv(a, b):
        ka = enc(a[0], off); kb = enc(b[0], off)
        keys = (ka[:, None] + kb[None, :]).ravel(); vals = (a[1][:, None] * b[1][None, :]).ravel()
        u, inv = np.unique(keys, return_inverse=True); s = np.zeros(len(u)); np.add.at(s, inv.ravel(), vals)
        return u, s
    A = conv(dists[0], dists[1]); B = conv(dists[2], dists[3])
    need = enc(T[None, :], 4 * off)[0] - B[0]
    pos = np.searchsorted(A[0], need); pos = np.minimum(pos, len(A[0]) - 1); hit = A[0][pos] == need
    return float((A[1][pos[hit]] * B[1][hit]).sum())


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); YB = a[1]; OUT = a[2]
    opt = lambda f, d: type(d)(a[a.index(f) + 1]) if f in a else d
    per, seed, N = opt('--per', 3), opt('--seed', 1), opt('--N', 12000)
    sl = W.build(cfg); L1 = len(sl.breps); t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    pats_of = [{int(y): AGs.fibre_patterns(sl, j, int(y))[0] for y in sl.letters[j]} for j in range(L1)]
    lmin = min(min(sl.letters[j]) for j in range(L1)); lmax = max(max(sl.letters[j]) for j in range(L1))
    fsz = np.zeros((L1, lmax - lmin + 1), np.int64)
    for j in range(L1):
        for y, P in pats_of[j].items(): fsz[j, y - lmin] = len(P)
    sizes = sorted({int(x) for x in fsz.ravel() if x > 0})
    # 1. exact H2 strata over all y's, with a uniform sample per stratum (smallest random keys)
    RB = 4 * L1; n_all = os.path.getsize(YB) // RB
    rng = np.random.default_rng(seed); strata = defaultdict(int); keep = {}
    jj = np.arange(L1)
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1).astype(np.int64)
        F = fsz[jj[None, None, :], Y - lmin]                                    # fibre size per letter
        assert (F > 0).all()
        prof = np.stack([(F == z).sum((1, 2)) for z in sizes], 1)              # letters per fibre size
        pk = (prof * (100 ** np.arange(len(sizes)))).sum(1)
        u, inv, c = np.unique(pk, return_inverse=True, return_counts=True)
        r = rng.random(len(Y))
        for i, key in enumerate(u.tolist()):
            strata[key] += int(c[i]); sel = np.nonzero(inv.ravel() == i)[0]
            cand = [(float(r[x]), s0 + int(x)) for x in sel[np.argsort(r[sel])[:per]]]
            keep[key] = sorted(keep.get(key, []) + cand)[:per]
    def H2_of(key):
        prof = [(key // 100 ** i) % 100 for i in range(len(sizes))]
        return sum(n * math.log2(z) for n, z in zip(prof, sizes)), prof
    tot_lifts = math.log2(sum(cnt * 2 ** H2_of(k)[0] for k, cnt in strata.items()))
    log(stage='strata', y_total=n_all, strata=len(strata), fibre_sizes=sizes, log2_total_lifts=round(tot_lifts, 2),
        H2_range=[round(min(H2_of(k)[0] for k in strata), 2), round(max(H2_of(k)[0] for k in strata), 2)])
    # 2.-3. per sampled y
    Z = zq_setup(sl); reps, _ = AGs.lag_reps(sl.G, sl.perms); Tfull = np.asarray(sl.T)[reps]
    done = {}
    if os.path.exists(OUT):
        for l in open(OUT):
            z = json.loads(l)
            if 'idx' in z: done[z['idx']] = z
    fo = open(OUT, 'a')
    order = sorted(strata, key=lambda k: -H2_of(k)[0])
    for key in order:
        H2, prof = H2_of(key)
        for _, idx in keep[key]:
            if idx in done: continue
            y = np.fromfile(YB, np.int8, count=RB, offset=idx * RB).reshape(4, L1).astype(np.int64)
            Yf = y[:, np.asarray(sl.borb)]
            assert (sl.Gp.paf(Yf).sum(0) == np.asarray(sl.T1)).all()
            pats = [[pats_of[j][int(y[k, j])] for j in range(L1)] for k in range(4)]
            rs = np.random.default_rng(1000 + idx); V = np.zeros((N, len(reps)), np.int64); Vq = np.zeros((N, len(Z['lreps'])), np.int64)
            for c0 in range(0, N, 1000):
                c = min(1000, N - c0); acc = np.zeros((c, sl.G.n), np.int64); accq = np.zeros((c, len(Z['lreps'])), np.int64)
                for k in range(4):
                    lab = np.zeros((c, sl.r), np.int8)
                    for j in range(L1):
                        P = pats[k][j]; pk_ = P[rs.integers(len(P), size=c)]
                        for i, (o, _) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
                    acc += sl.G.paf(lab[:, sl.orb])
                    accq += zq_paf(lab.astype(np.int64) @ Z['Mo'], Z['lreps'])
                V[c0:c0 + c] = acc[:, reps]; Vq[c0:c0 + c] = accq
            gb = gauss_bits(V, Tfull); gq = gauss_bits(Vq, Z['Tq'])
            rec = dict(idx=idx, stratum=prof, stratum_size=strata[key], H2=round(H2, 3),
                       C2=gb.get('bits'), rank=gb.get('rank'), mahal=gb.get('mahal'), in_lattice=gb.get('in_lattice'),
                       log2E=round(H2 - gb['bits'], 2) if gb.get('in_lattice') and gb.get('bits') is not None and np.isfinite(gb['bits']) else None,
                       zq_bits=gq.get('bits'), zq_rank=gq.get('rank'), zq_log2E=round(H2 - gq['bits'], 3) if gq.get('in_lattice') else None)
            if '--exact-zq' in a:
                ex = zq_exact(sl, Z, pats); rec['zq_exact_log2'] = round(math.log2(ex), 3) if ex > 0 else None
            fo.write(json.dumps(rec) + "\n"); fo.flush(); done[idx] = rec
            log(stage='y', **{k_: rec[k_] for k_ in ('idx', 'H2', 'C2', 'log2E', 'zq_log2E') if k_ in rec}, zq_exact_log2=rec.get('zq_exact_log2'))
    # 4. the sum over strata
    E = 0.0; contrib = []
    for key in order:
        recs = [done[i] for _, i in keep[key] if i in done and done[i].get('log2E') is not None]
        if not recs: continue
        m = float(np.mean([2.0 ** r['log2E'] for r in recs])); E += strata[key] * m
        contrib.append((math.log2(strata[key] * m) if m > 0 else -math.inf, H2_of(key)[1], strata[key]))
    contrib.sort(reverse=True)
    zr = [r['zq_exact_log2'] - r['zq_log2E'] for r in done.values() if r.get('zq_exact_log2') is not None and r.get('zq_log2E') is not None]
    log(stage='done', y_total=n_all, log2_expected_cores=round(math.log2(E), 2) if E > 0 else None,
        top_strata=[dict(log2_contrib=round(c, 2), profile=p, size=s) for c, p, s in contrib[:5]],
        zq_exact_minus_gauss_bits=dict(n=len(zr), mean=round(float(np.mean(zr)), 3), min=round(float(np.min(zr)), 3), max=round(float(np.max(zr)), 3)) if zr else None)


if __name__ == "__main__":
    main()
