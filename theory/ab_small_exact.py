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


def block_labels(sl, yk, cap):
    """all orbit-label vectors (n x r, int8) of the slice with compression y_k; None if more than cap"""
    pats = []; tot = 1
    for j, bb in enumerate(sl.breps):
        P, sizes = AGs.fibre_patterns(sl, j, int(yk[bb]))
        if len(P) == 0: return np.zeros((0, sl.r), np.int8)
        pats.append(P); tot *= len(P)
        if tot > cap: return None
    idx = np.array(list(itertools.product(*[range(len(P)) for P in pats])), np.int64).reshape(-1, len(pats))
    lab = np.zeros((len(idx), sl.r), np.int8)
    for j, P in enumerate(pats):
        sel = P[idx[:, j]]
        for i, (o, s_) in enumerate(sl.parts[j]): lab[:, o] = sel[:, i]
    return lab


def block_fibre(sl, yk, cap):
    lab = block_labels(sl, yk, cap)
    return None if lab is None else lab[:, sl.orb]


_CACHE = {}


def block_data(sl, yk, reps, gw, cap):
    """cached per block vector y_k: (orbit labels, PAF-hash of every fibre element)"""
    key = np.asarray(yk, np.int64).tobytes()
    if key in _CACHE: return _CACHE[key]
    lab = block_labels(sl, yk, cap)
    if lab is None: res = (None, None)
    else:
        h = np.zeros(len(lab), np.int64)
        for c0 in range(0, len(lab), 4096):
            h[c0:c0 + 4096] = sl.G.paf(lab[c0:c0 + 4096][:, sl.orb].astype(np.int64))[:, reps] @ gw
        res = (lab, h)
    if len(_CACHE) > 20000: _CACHE.clear()
    _CACHE[key] = res
    return res


def solve_y(sl, Y, reps, gw, cap=1 << 16):
    labs = []; h = []
    for k in range(4):
        lab, hk = block_data(sl, Y[k], reps, gw, cap)
        if lab is None: return None, [None] * 4
        if len(lab) == 0: return [], [0] * 4
        labs.append(lab); h.append(hk)
    m = [len(x) for x in labs]
    # pairing with the smallest join cost: the smaller pair side is sorted, the larger streamed
    best = None
    for (a, b), (c, d) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        L_, R_ = m[a] * m[b], m[c] * m[d]
        if L_ > R_: (a, b), (c, d), L_, R_ = (c, d), (a, b), R_, L_
        if best is None or (L_, L_ + R_) < (best[0], best[0] + best[1]): best = (L_, R_, (a, b), (c, d))
    _, _, (a, b), (c, d) = best
    tw = np.int64(sl.T[reps] @ gw)
    HL = (h[a][:, None] + h[b][None, :]).ravel(); o = np.argsort(HL); HLs = HL[o]
    hits = []
    step = max(1, (1 << 22) // m[d])
    for i0 in range(0, m[c], step):
        HR = (tw - (h[c][i0:i0 + step][:, None] + h[d][None, :])).ravel()
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            ic, id_ = divmod(int(r), m[d]); ic += i0
            for p in o[lo[r]:hi[r]]:
                ia, ib = divmod(int(p), m[b])
                sel = {a: ia, b: ib, c: ic, d: id_}
                Xc = np.array([labs[k][sel[k]][sl.orb] for k in range(4)], np.int64)
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
