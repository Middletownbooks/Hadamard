r"""Level budget of the index-s, U = <2> slice of v = 7q (theory/b917_slice.py) per stage-1 solution y [M counts, E].

[P] Given y, block k is a 3 x s sign matrix (class columns) with row sums sigma_k (y(R_i) = alpha_i + f sigma_i) and the
column 0 fixed (alpha_k).  Its dual energy (lag 0 of the Z_q compression) is
    E_k = A_k^2 + f sum_j w_j^2 = A_k^2 + f (s + 24 n5 + 48 n7),   A_k = alpha_0 + 3 alpha_1 + 3 alpha_2,
(n5, n7 = numbers of columns with |w| = 5, 7), and sum_k E_k = 4v is the lag-0 dual equation: an exact condition on the
column types only.  Here:
  exact     number of configurations given y, per (n5, n7), by multinomials (python ints);
            the energy condition sum_k (n5 + 2 n7)_k = L* is then exact, and H_E = log2 #(x | y, energy);
  Gaussian  uniform sampling of x given (y, energy) (exact weights), the lattice Gaussian (theory/q523_gauss.gauss_bits)
            of sum_k PAF_G at the lag classes: the remaining dual classes, the cross classes, and both.
usage: python3 theory/b917_levels.py [--q 131 --s 26] [--ys 0,1,2] [--N 20000]"""
import sys, os, math, json, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b917_slice import Slice917
from q523_gauss import gauss_bits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cfg_json(q, s):
    return json.dumps(dict(v=7 * q, c=1, b=0, h=None, m=7 * q, Q=[], H1=[7], q=q, s=s, f=(q - 1) // s, mods=[7, q],
                           nQ=0, U=[[0, 3, 2]]))


def stage1_ys(q, s):
    """stage-1 solutions (up to block order and sign) as an array (n, 4, 7)"""
    import ab_sweep as W
    from ab_s1sd import cfg_from_json
    from ab_b0_stage1 import stage1
    sl = W.build(cfg_from_json(cfg_json(q, s)))
    sols, _ = stage1(sl)
    return np.array(sols, np.int64)


def alpha_sigma(Y, f):
    """Y (4, 7) on Z_7 -> alpha (4, 3), sigma (4, 3) per row R0, R1, R2"""
    Yr = np.asarray(Y)[:, [0, 1, 3]]
    al = np.where((Yr - 1) % (2 * f) == 0, 1, -1)
    sg = (Yr - al) // f
    assert ((al + f * sg) == Yr).all()
    return al, sg


def block_table(s, sig):
    """all (n5, n7, c, b, d) with weights: c = (c77, c7m, c5p, c5m) sym state counts, b/d = anti sums of x0/e;
    weight = C(s, n_s) * multinomial(sym) * C(n_a, .) C(n_a, .)  (count of 3 x s matrices with row sums sig)"""
    s0, s1, s2 = (int(x) for x in sig)
    if (s1 + s2) % 2 or (s1 - s2) % 2: return []
    E_sum, D_sum = (s1 + s2) // 2, (s1 - s2) // 2          # sum_sym e, sum_anti e
    out = []
    for ns in range(s + 1):
        na = s - ns
        if abs(D_sum) > na or (na + D_sum) % 2: continue
        cd = math.comb(na, (na + D_sum) // 2)
        for c77 in range(ns + 1):
            for c7m in range(ns + 1 - c77):
                for c5p in range(ns + 1 - c77 - c7m):
                    c5m = ns - c77 - c7m - c5p
                    if c77 + c5p - c7m - c5m != E_sum: continue
                    x0s = c77 + c5m - c7m - c5p
                    b = s0 - x0s
                    if abs(b) > na or (na + b) % 2: continue
                    w = math.comb(s, ns) * (math.factorial(ns) // (math.factorial(c77) * math.factorial(c7m) * math.factorial(c5p) * math.factorial(c5m))) \
                        * math.comb(na, (na + b) // 2) * cd
                    out.append((c5p + c5m, c77 + c7m, (c77, c7m, c5p, c5m), b, D_sum, w))
    return out


class LevelModel:
    def __init__(self, S, Y):
        self.S = S; s = S.s; f = S.f
        self.al, self.sg = alpha_sigma(Y, f)
        self.A = (self.al @ np.array([1, 3, 3])).astype(int)
        self.tabs = [block_table(s, self.sg[k]) for k in range(4)]
        self.H2 = sum(math.log2(math.comb(s, (s + int(x)) // 2)) for k in range(4) for x in self.sg[k])
        tot = [sum(t[-1] for t in tab) for tab in self.tabs]
        assert all(tot[k] == math.prod(math.comb(s, (s + int(x)) // 2) for x in self.sg[k]) for k in range(4))
        # energy: sum_k A_k^2 + f (4 s + 24 L) = 4v, L = sum_k (n5 + 2 n7)
        num = S.N4 - int((self.A ** 2).sum()) - 4 * s * f
        self.Lstar = num // (24 * f) if num % (24 * f) == 0 else None
        self.byL = []
        for tab in self.tabs:
            d = defaultdict(int)
            for (n5, n7, c, b, dd, w) in tab: d[n5 + 2 * n7] += w
            self.byL.append(dict(d))
        if self.Lstar is None:
            self.count_E = 0
        else:
            # convolution over 4 blocks
            acc = {0: 1}
            for k in range(4):
                new = defaultdict(int)
                for L1, c1 in acc.items():
                    for L2, c2 in self.byL[k].items():
                        if L1 + L2 <= self.Lstar: new[L1 + L2] += c1 * c2
                acc = new
            self.count_E = acc.get(self.Lstar, 0)
        self.HE = math.log2(self.count_E) if self.count_E else -math.inf

    def summary(self):
        return dict(alpha=self.al.tolist(), sigma=self.sg.tolist(), A=self.A.tolist(), Lstar=self.Lstar,
                    H2=round(self.H2, 2), HE=round(self.HE, 2), energy_bits=round(self.H2 - self.HE, 2))

    def sampler(self, rng):
        """uniform sampler of x (4, 3, s+1) given y and the dual energy"""
        S = self.S; s = S.s
        # split of L* over blocks, by exact weights
        Ls = [sorted(d) for d in self.byL]
        splits = []; wts = []
        for L0 in Ls[0]:
            for L1 in Ls[1]:
                for L2 in Ls[2]:
                    L3 = self.Lstar - L0 - L1 - L2
                    if L3 in self.byL[3]:
                        splits.append((L0, L1, L2, L3))
                        wts.append(self.byL[0][L0] * self.byL[1][L1] * self.byL[2][L2] * self.byL[3][L3])
        tw = sum(wts); P_split = np.array([float(w_ / tw) for w_ in wts])
        per = [defaultdict(list) for _ in range(4)]
        for k in range(4):
            for t in self.tabs[k]: per[k][t[0] + 2 * t[1]].append(t)
        def draw_block(k, L):
            opts = per[k][L]; ws = np.array([float(t[-1]) for t in opts]); ws /= ws.sum()
            n5, n7, (c77, c7m, c5p, c5m), b, d, _ = opts[rng.choice(len(opts), p=ws)]
            ns = n5 + n7; na = s - ns
            cols = []
            cols += [(1, 1, 1)] * c77 + [(-1, -1, -1)] * c7m + [(-1, 1, 1)] * c5p + [(1, -1, -1)] * c5m
            x0 = np.array([1] * ((na + b) // 2) + [-1] * ((na - b) // 2)); rng.shuffle(x0)
            e = np.array([1] * ((na + d) // 2) + [-1] * ((na - d) // 2)); rng.shuffle(e)
            cols += [(int(x0[i]), int(e[i]), -int(e[i])) for i in range(na)]
            cols = np.array(cols, np.int64); rng.shuffle(cols)
            X = np.zeros((3, s + 1), np.int64); X[:, 0] = self.al[k]; X[:, 1:] = cols.T
            return X
        def draw():
            L4 = splits[rng.choice(len(splits), p=P_split)]
            return np.array([draw_block(k, L4[k]) for k in range(4)])
        return draw


def lag_class_reps(S):
    """representatives (h, z) of the lag classes of G under U x M and negation: (0,0), (0, gamma^j) j < s/2, (1, 0),
    (1, gamma^j) j < s"""
    q, s, g = S.q, S.s, S.gamma
    dual = [(0, pow(g, j, q)) for j in range(s // 2)]
    cross = [(1, pow(g, j, q)) for j in range(s)]
    return [(0, 0)] + dual + [(1, 0)] + cross


def sample_paf(S, draw, N, reps, chunk=2000):
    V = np.zeros((N, len(reps)), np.int64); done = 0
    hh = np.array([r[0] for r in reps]); zz = np.array([r[1] for r in reps])
    while done < N:
        c = min(chunk, N - done)
        Xs = np.array([draw() for _ in range(c)])                     # (c, 4, 3, s+1)
        x = S.expand(Xs).astype(np.float64)                           # (c, 4, 7, q)
        F = np.fft.fft2(x)
        P = np.rint(np.fft.ifft2((np.abs(F) ** 2).sum(1)).real).astype(np.int64)   # (c, 7, q)
        V[done:done + c] = P[:, hh, zz]; done += c
    return V


def main():
    a = sys.argv[1:]
    opt = lambda n, d: a[a.index(n) + 1] if n in a else d
    q = int(opt('--q', '131')); s = int(opt('--s', '26')); N = int(opt('--N', '20000'))
    S = Slice917(q, s)
    t0 = time.time()
    Ys = stage1_ys(q, s)
    print(json.dumps(dict(v=S.v, q=q, s=s, f=S.f, stage1_solutions=len(Ys), secs=round(time.time() - t0, 1))), flush=True)
    sel = [int(x) for x in opt('--ys', '0,1,2').split(',')]
    reps = lag_class_reps(S)
    T = np.zeros(len(reps), np.int64); T[0] = S.N4
    h = s // 2
    idx_dual = list(range(0, 1 + h)); idx_cross = list(range(1 + h, len(reps)))
    out = os.path.join(ROOT, 'logs/b917/levels.jsonl')
    for yi in sel:
        Y = Ys[yi]; t1 = time.time()
        lm = LevelModel(S, Y)
        rec = dict(q=q, s=s, y_index=yi, y=Y.tolist(), **lm.summary())
        if lm.count_E:
            rng = np.random.default_rng(7 + yi)
            draw = lm.sampler(rng)
            # sanity: every sample has the right stage 1 and dual energy
            Xs = draw(); assert (S.ystage1(Xs)[:, :] == Y[:, [0, 1, 3]]).all()
            V = sample_paf(S, draw, N, reps)
            assert (V[:, 0] == S.N4).all(), "dual energy not exact"
            g_all = gauss_bits(V, T); g_d = gauss_bits(V[:, idx_dual], T[idx_dual]); g_c = gauss_bits(V[:, idx_cross], T[idx_cross])
            rec.update(C_all=g_all.get('bits'), rank_all=g_all.get('rank'), mahal_all=g_all.get('mahal'),
                       C_dual=g_d.get('bits'), rank_dual=g_d.get('rank'), mahal_dual=g_d.get('mahal'),
                       C_cross=g_c.get('bits'), rank_cross=g_c.get('rank'), mahal_cross=g_c.get('mahal'))
            if g_all.get('bits') is not None and np.isfinite(g_all['bits']):
                rec['log2E_per_y'] = round(lm.HE - g_all['bits'], 2)
            sd = V.std(0)
            rec['std_dual'] = [round(float(x), 1) for x in sd[idx_dual]]
            rec['std_cross'] = [round(float(x), 1) for x in sd[idx_cross]]
            rec['mean_cross'] = [round(float(x), 1) for x in V[:, idx_cross].mean(0)]
            rec['mean_dual'] = [round(float(x), 1) for x in V[:, idx_dual].mean(0)]
        rec['secs'] = round(time.time() - t1, 1)
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
