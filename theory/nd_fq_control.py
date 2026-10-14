r"""Exact control for the F_{p^2} uniform-cyclotomy slices at a small order (v = 55: Z_6 x F_9, s = 4) [M vs E].

1. Stage 1: every compressed solution y (theory/nd_s1_sd.py on the explicit structure of theory/nd_fq_split.py).
2. Stage 2 for each y, exhaustively: per block the configurations are the per-position patterns c in {+-1}^s with
   sum c = sigma(g) (all combinations), the signature is the class-system vector (for each class i the Z_F PAF of
   w_i = beta + p c_i, lags 0..F/2), and a 2+2 sorted join counts every configuration with sum_k = target
   (4v at lag 0, 0 elsewhere, for every class) - the exact number of stage-2 completions of y.
3. The same y through the estimators: theory/nd_fq.levels_classes (pattern-sampled joint Gaussian) and
   theory/nd_fq.class_bits (independent class systems).
4. One completion is turned into the full core on G = Z_F x F_q, checked by the full-group PAF (3-D FFT) and completed
   to a bordered Goethals-Seidel Hadamard matrix (independent check H H^T = N I in int64).
usage: python3 theory/nd_fq_control.py F p s [--maxy 20]"""
import sys, os, json, math, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nd_fq as FQ
import nd_fq_split as FS
import nd_s1_sd as S1

def stage1(F, p, s, max_sol=10 ** 6):
    C = FS.compressed(F, p, [(1, s)])
    C['MF'] = [1]
    P = S1.Problem(C['m'], 0, C['q'], explicit=C); P.C_explicit = C
    out, sols = P.solve(max_sol=max_sol, log=lambda x: None)
    return out, [np.array(x['y']) for x in sols]

def block_configs(F, p, s, beta_k, sigma_k):
    """all configurations of one block: (C (n, F, s) patterns, signature (n, s, F/2+1))"""
    pats = []
    for g in range(F):
        sg = int(sigma_k[g]); L = []
        for pos in itertools.combinations(range(s), (s + sg) // 2):
            c = -np.ones(s, np.int64); c[list(pos)] = 1; L.append(c)
        pats.append(np.array(L))
    idx = np.array(list(itertools.product(*[range(len(L)) for L in pats])), np.int64)
    Cc = np.stack([pats[g][idx[:, g]] for g in range(F)], 1)              # (n, F, s)
    W = beta_k[None, :, None] + p * Cc
    H = F // 2
    sig = np.stack([(W * np.roll(W, -h, axis=1)).sum(1) for h in range(H + 1)], 2)   # (n, s, H+1)
    return Cc, sig.reshape(len(idx), -1)

def stage2_exact(F, p, s, Y, max_left=6 * 10 ** 7):
    """exact number of stage-2 completions of y (2+2 join; left pairs stored as one 64-bit hash + index, every hash
    match verified on the exact signature); returns (count, alpha, first completion) or (None, ...) if too large"""
    alpha, sigma, beta, R = FQ.class_system(F, p, s, Y)
    v = F * p * p + 1; H = F // 2
    T = np.tile(np.array([4 * v] + [0] * H, np.int64), s)
    blocks = [block_configs(F, p, s, beta[k], sigma[k]) for k in range(4)]
    order = sorted(range(4), key=lambda k: len(blocks[k][1]))           # two smallest blocks on the left
    blocks_o = [blocks[k] for k in order]
    sigs = [b[1] for b in blocks_o]
    n0, n1 = len(sigs[0]), len(sigs[1])
    if n0 * n1 > max_left or len(sigs[2]) * len(sigs[3]) > 20 * max_left: return None, alpha, None
    rng = np.random.default_rng(12345)
    w = rng.integers(1, 2 ** 61, sigs[0].shape[1], dtype=np.int64)
    h0 = sigs[0] @ w; h1 = sigs[1] @ w                                   # int64 wrap-around hash (linear)
    left = (h0[:, None] + h1[None, :]).ravel()
    o = np.argsort(left); ls = left[o]
    hT = T @ w; h2 = sigs[2] @ w; h3 = sigs[3] @ w
    count = 0; first = None
    for c0 in range(len(sigs[2])):
        r = hT - h2[c0] - h3
        a = np.searchsorted(ls, r, 'left'); b = np.searchsorted(ls, r, 'right')
        for j in np.nonzero(b > a)[0]:
            for t in range(a[j], b[j]):
                i0, i1 = divmod(int(o[t]), n1)
                if ((sigs[0][i0] + sigs[1][i1] + sigs[2][c0] + sigs[3][j]) == T).all():
                    count += 1
                    if first is None: first = (i0, i1, c0, int(j))
    Cfirst = None
    if first is not None:
        Cfirst = [None] * 4
        for pos, k in enumerate(order): Cfirst[k] = blocks[k][0][first[pos]]
    return count, alpha, Cfirst

def nconfigs(F, p, s, Y):
    alpha, sigma, beta, R = FQ.class_system(F, p, s, Y)
    c = [math.prod(math.comb(s, (s + int(x)) // 2) for x in sigma[k]) for k in range(4)]
    return c[0] * c[1] + c[2] * c[3], c

def full_core(F, p, s, alpha, Cs):
    Sl = FQ.Slice(F, p, s)
    Xs = [Sl.x_of(alpha[k], Cs[k]) for k in range(4)]
    P = Sl.paf_full(Xs)
    m = F * p * p
    ok = P[0, 0, 0] == 4 * m and (P.ravel()[1:] == -4).all()
    rows = sorted(int(X.sum()) for X in Xs)
    return Xs, bool(ok), rows

def bordered_gs(Xs, F, p):
    """bordered Goethals-Seidel b = 1 completion over G = Z_F x Z_p x Z_p (abelian): group matrices A_k(g, h) =
    x_k(h - g), R(g, h) = [g + h = 0]; border rows/cols S x 1^T, T^T x 1, corner -S Sigma^T T / 4 over all S, T in H(4)
    up to signs; returns the matrix if some choice gives H H^T = N I (int64 check)"""
    G = [(a, b, c) for a in range(F) for b in range(p) for c in range(p)]
    n = len(G); idx = {g: i for i, g in enumerate(G)}
    def sub(g, h): return ((h[0] - g[0]) % F, (h[1] - g[1]) % p, (h[2] - g[2]) % p)
    D = np.array([[idx[sub(g, h)] for h in G] for g in G])
    flat = [np.array([X[g] for g in G]) for X in Xs]
    Mk = [f[D] for f in flat]
    neg = np.array([idx[((-g[0]) % F, (-g[1]) % p, (-g[2]) % p)] for g in G])
    Rm = np.zeros((n, n), np.int64); Rm[np.arange(n), neg] = 1
    Aa, Bb, Cc, Dd = Mk
    GS = np.block([[Aa, Bb @ Rm, Cc @ Rm, Dd @ Rm],
                   [-Bb @ Rm, Aa, Dd.T @ Rm, -Cc.T @ Rm],
                   [-Cc @ Rm, -Dd.T @ Rm, Aa, Bb.T @ Rm],
                   [-Dd @ Rm, Cc.T @ Rm, -Bb.T @ Rm, Aa]])
    N = 4 * n + 4
    Sig = np.array([[GS[I * n:(I + 1) * n, J * n:(J + 1) * n].sum(1)[0] for J in range(4)] for I in range(4)])
    H4 = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    H4 = [h for h in H4 if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]
    for Sm in H4:
        for Tm in H4:
            num = -(Sm @ Sig.T @ Tm)
            if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
            Hm = np.zeros((N, N), np.int64)
            Hm[:4, :4] = num // 4
            Hm[:4, 4:] = np.repeat(Sm, n, axis=1); Hm[4:, :4] = np.repeat(Tm.T, n, axis=0); Hm[4:, 4:] = GS
            if (Hm @ Hm.T == N * np.eye(N, dtype=np.int64)).all(): return Hm
    return None

if __name__ == "__main__":
    a = sys.argv[1:]; maxy = 20
    if "--maxy" in a: i = a.index("--maxy"); maxy = int(a[i + 1]); a = a[:i] + a[i + 2:]
    F, p, s = map(int, a)
    t0 = time.time()
    o, Ys = stage1(F, p, s, max_sol=4000)
    Ys = sorted(Ys, key=lambda Y: nconfigs(F, p, s, Y)[0])
    print(json.dumps(dict(F=F, p=p, s=s, v=F * p * p + 1, stage1=o)), flush=True)
    tot_exact = 0; tot_pred_pat = 0.0; tot_pred_ind = 0.0; built = False
    for i, Y in enumerate(Ys[:maxy]):
        cnt, alpha, Cf = stage2_exact(F, p, s, Y)
        if cnt is None:
            print(json.dumps(dict(y=i, note="too large for the exact join", configs=nconfigs(F, p, s, Y)[1])), flush=True); continue
        lv = FQ.levels_classes(F, p, s, Y.tolist(), N=20000)
        ind = FQ.class_bits(F, p, s, Y.tolist(), N=20000)
        pe = lv['H2'] - lv['all']['bits'] if lv['all'].get('in_lattice') else -math.inf
        rec = dict(y=i, exact_stage2=cnt, log2_exact=round(math.log2(cnt), 2) if cnt else None,
                   pattern_model=round(pe, 2), independent_model=ind['stage2_log2E_indep'], log2_S_Phi=ind['log2_S_Phi'],
                   H2=lv['H2'], C2=lv['all']['bits'])
        tot_exact += cnt; tot_pred_pat += 2 ** pe if pe > -1e9 else 0; tot_pred_ind += 2 ** ind['stage2_log2E_indep']
        if cnt and not built:
            Xs, ok, rows = full_core(F, p, s, alpha, Cf)
            rec.update(core_fullgroup_ok=ok, rows=rows)
            Hm = bordered_gs(Xs, F, p) if ok else None
            if Hm is not None:
                rec.update(hadamard_order=int(Hm.shape[0]), HHt_ok=True); built = True
                np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs', 'nd', f'H{Hm.shape[0]}_fq_F{F}_p{p}_s{s}.npy'), Hm.astype(np.int8))
        print(json.dumps(rec), flush=True)
    print(json.dumps(dict(summary=True, ys=min(len(Ys), maxy), exact_total=tot_exact,
                          pattern_model_total=round(tot_pred_pat, 2), independent_model_total=round(tot_pred_ind, 4),
                          secs=round(time.time() - t0, 1))), flush=True)
