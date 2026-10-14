r"""Exact stage 2 for small b >= 3 (or b = 0) GP slices, given stage-1 solutions [M].

For a stage-1 y (4 x |G'|), block k's fibre is the product over G'-orbits j of the sign patterns on the parts of j
with signed size-sum y_k(rep_j) (theory/ab_gauss.fibre_patterns).  When every block fibre is small (<= ~2^16), the
whole stage 2 is one exact 2+2 hashed join of the full PAF over G (every lag class of G under the slice group and -1)
against the target T (4m at 0, -4b on H \ 0, +4 off H).  Every match is re-checked exactly (sum_k PAF_G = T).
Cores are saved to logs/ab/small/core_<tag>_<n>.npy (4 x |G|, group coordinates of the cfg).
usage: python3 theory/ab_small_exact.py '<cfg json>' Y.json [--tag T] [--max-fibre 65536]
       (Y.json: a list of stage-1 y's, as written by theory/ab_s1_one.py --out (key "examples"))"""
import sys, os, json, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def block_fibre(sl, yk, cap):
    """all x on G (rows) in the slice with compression y_k; None if larger than cap"""
    pats = []; tot = 1
    for j, bb in enumerate(sl.breps):
        P, sizes = AGs.fibre_patterns(sl, j, int(yk[bb]))
        if len(P) == 0: return np.zeros((0, sl.G.n), np.int8)
        pats.append(P); tot *= len(P)
        if tot > cap: return None
    idx = np.array(list(itertools.product(*[range(len(P)) for P in pats])), np.int64).reshape(-1, len(pats))
    lab = np.zeros((len(idx), sl.r), np.int8)
    for j, P in enumerate(pats):
        sel = P[idx[:, j]]
        for i, (o, s_) in enumerate(sl.parts[j]): lab[:, o] = sel[:, i]
    return lab[:, sl.orb]


def solve_y(sl, Y, reps, gw, cap=1 << 16):
    Xs = []
    for k in range(4):
        X = block_fibre(sl, Y[k], cap)
        if X is None: return None, [None] * 4
        if len(X) == 0: return [], [0] * 4
        Xs.append(X)
    m = [len(x) for x in Xs]
    h = []
    for X in Xs:
        H_ = np.zeros(len(X), np.int64)
        for c0 in range(0, len(X), 4096):
            H_[c0:c0 + 4096] = sl.G.paf(X[c0:c0 + 4096].astype(np.int64))[:, reps] @ gw
        h.append(H_)
    tw = np.int64(sl.T[reps] @ gw)
    HL = (h[0][:, None] + h[1][None, :]).ravel(); o = np.argsort(HL); HLs = HL[o]
    HR = (tw - (h[2][:, None] + h[3][None, :])).ravel()
    lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
    hits = []
    for r in np.nonzero(hi > lo)[0]:
        c_, d_ = divmod(int(r), m[3])
        for p in o[lo[r]:hi[r]]:
            a_, b_ = divmod(int(p), m[1])
            Xc = np.array([Xs[0][a_], Xs[1][b_], Xs[2][c_], Xs[3][d_]], np.int64)
            if (sl.G.paf(Xc).sum(0) == sl.T).all(): hits.append(Xc)
    return hits, m


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(a[0]); Ys = json.load(open(a[1]))
    if isinstance(Ys, dict): Ys = Ys.get('examples') or Ys.get('solutions') or []
    tag = a[a.index('--tag') + 1] if '--tag' in a else 'small'
    cap = int(a[a.index('--max-fibre') + 1]) if '--max-fibre' in a else 1 << 16
    sl = W.build(cfg)
    reps, _ = AGs.lag_reps(sl.G, sl.perms)
    gw = np.random.default_rng(5).integers(1, 2 ** 62, size=len(reps), dtype=np.int64)
    outdir = os.path.join(ROOT, 'logs/ab/small'); os.makedirs(outdir, exist_ok=True)
    ncore = 0; t0 = time.time()
    for i, Y in enumerate(Ys):
        Y = np.array(Y, np.int64)
        hits, m = solve_y(sl, Y, reps, gw, cap)
        rec = dict(i=i, fibres=m, cores=None if hits is None else len(hits))
        if hits:
            for Xc in hits:
                ncore += 1; np.save(os.path.join(outdir, f'core_{tag}_{ncore}.npy'), Xc.astype(np.int8))
        print(json.dumps(rec), flush=True)
    print(json.dumps(dict(stage="end", ys=len(Ys), cores=ncore, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    main()
