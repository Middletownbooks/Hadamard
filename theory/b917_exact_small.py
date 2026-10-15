r"""Exact core counts per stage-1 solution y in a small index-s, U = <2> slice of v = 7q, against the level Gaussian
of theory/b917_levels.py [M vs E].

Per block k: every 3 x s sign matrix with row sums sigma_k (and column 0 = alpha_k from y); its signature = PAF over
G = Z_7 x Z_q at the lag-class representatives (theory/b917_levels.lag_class_reps).  Complete 2+2 join: all pairs
(blocks 0, 1) against all pairs (blocks 2, 3) on the signature sum = 4v delta (64-bit linear hash, exact re-check of
every match with the full PAF over G).  The count is compared with H_E - C_all (energy exact, lattice Gaussian).
usage: python3 theory/b917_exact_small.py --q 31 --s 6 --ys 0,1,2,3"""
import sys, os, json, time, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b917_slice import Slice917
from b917_levels import LevelModel, stage1_ys, lag_class_reps, alpha_sigma
from q523_gauss import gauss_bits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def block_configs(S, alpha, sigma):
    """all 3 x (s+1) blocks with column 0 = alpha and class-row sums sigma"""
    s = S.s
    rows = []
    for i in range(3):
        n_plus = (s + int(sigma[i])) // 2
        opts = []
        for pos in itertools.combinations(range(s), n_plus):
            r = -np.ones(s, np.int8); r[list(pos)] = 1; opts.append(r)
        rows.append(np.array(opts, np.int8))
    n0, n1, n2 = len(rows[0]), len(rows[1]), len(rows[2])
    X = np.zeros((n0 * n1 * n2, 3, s + 1), np.int8)
    X[:, :, 0] = np.asarray(alpha, np.int8)
    i0, i1, i2 = np.meshgrid(np.arange(n0), np.arange(n1), np.arange(n2), indexing='ij')
    X[:, 0, 1:] = rows[0][i0.ravel()]; X[:, 1, 1:] = rows[1][i1.ravel()]; X[:, 2, 1:] = rows[2][i2.ravel()]
    return X


def signatures(S, X, reps, chunk=8192):
    hh = np.array([r[0] for r in reps]); zz = np.array([r[1] for r in reps])
    out = np.zeros((len(X), len(reps)), np.int64)
    for c0 in range(0, len(X), chunk):
        x = S.expand(X[c0:c0 + chunk]).astype(np.float64)
        F = np.fft.fft2(x)
        P = np.rint(np.fft.ifft2(np.abs(F) ** 2).real).astype(np.int64)
        out[c0:c0 + chunk] = P[:, hh, zz]
    return out


def count_y(S, Y, reps, T):
    al, sg = alpha_sigma(Y, S.f)
    Xs = [block_configs(S, al[k], sg[k]) for k in range(4)]
    Ps = [signatures(S, X, reps) for X in Xs]
    rng = np.random.default_rng(5)
    wv = rng.integers(1, 2 ** 62, size=len(reps), dtype=np.int64)
    h = [P @ wv for P in Ps]
    n = [len(x) for x in h]
    HL = (h[0][:, None] + h[1][None, :]).ravel()
    o = np.argsort(HL, kind='stable'); HLs = HL[o]; del HL
    tw = np.int64(T @ wv)
    count = 0; checked = 0
    step = max(1, (1 << 22) // n[3])
    for i0 in range(0, n[2], step):
        HR = (tw - (h[2][i0:i0 + step][:, None] + h[3][None, :])).ravel()
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            i2, i3 = divmod(int(r), n[3]); i2 += i0
            for p in o[lo[r]:hi[r]]:
                i0_, i1_ = divmod(int(p), n[1])
                checked += 1
                sig = Ps[0][i0_] + Ps[1][i1_] + Ps[2][i2] + Ps[3][i3]
                if (sig == T).all():
                    Xc = np.array([Xs[0][i0_], Xs[1][i1_], Xs[2][i2], Xs[3][i3]])
                    if S.is_core(Xc): count += 1
    return count, n, checked


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    q = int(opt('--q', '31')); s = int(opt('--s', '6')); N = int(opt('--N', '20000'))
    S = Slice917(q, s)
    Ys = stage1_ys(q, s)
    reps = lag_class_reps(S)
    T = np.zeros(len(reps), np.int64); T[0] = S.N4; T[1 + s // 2] = 0
    # the (1, 0) class is stage 1 (fixed by y): its target is 0 too (b = 0)
    out = os.path.join(ROOT, 'logs/b917/exact_small.jsonl')
    for yi in [int(x) for x in opt('--ys', '0').split(',')]:
        Y = Ys[yi]; t0 = time.time()
        cnt, n, checked = count_y(S, Y, reps, T)
        lm = LevelModel(S, Y)
        rec = dict(q=q, s=s, v=S.v, y_index=yi, y=Y.tolist(), block_list_sizes=n, exact_cores=cnt,
                   log2_exact=round(math.log2(cnt), 2) if cnt else None, H_E=round(lm.HE, 2))
        if lm.count_E:
            draw = lm.sampler(np.random.default_rng(9 + yi))
            from b917_levels import sample_paf
            V = sample_paf(S, draw, N, reps)
            g = gauss_bits(V, T)
            rec.update(C_all=g.get('bits'), log2E_gauss=round(lm.HE - g['bits'], 2) if g.get('bits') is not None else None)
            if cnt and rec.get('log2E_gauss') is not None:
                rec['exact_minus_gauss_bits'] = round(math.log2(cnt) - rec['log2E_gauss'], 2)
        rec['secs'] = round(time.time() - t0, 1)
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (int(3.5 * 2 ** 30), int(3.5 * 2 ** 30)))
    main()
