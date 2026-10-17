r"""Half-shift templates beyond arcs: which S/D support patterns sigma make both halves of a b = 1 core feasible?

Setting (report/alpoge_methods_forensics.md §4.2).  m = 2h, h odd.  For a fixed support pattern sigma (S_k = {y in Z_h :
x_k(y) = x_k(y + h)}, D_k its complement) the b = 1 core equations split exactly into two independent halves on Z_h:
    S-half: c_k = +-1 on S_k,  sum_k PAF_{c_k}(t)  = -2   (t = 1..n, n = (h-1)/2),   sum_k |S_k| = 2h - 1,
    D-half: d_k = +-1 on D_k,  sum_k NPAF_{d_k}(t) =  0   (negacyclic),             sum_k |D_k| = 2h + 1.
This is the S/D split along the Z_2 axis of Z_2h = Z_2 x Z_h that our b = 1 solvers use (523, 647); the difference is
that sigma is fixed before the search, so the halves share no variable.

Lemma [P]: with V_t = sum_k |S_k cap (S_k - t)| and W_t = sum_k |D_k cap (D_k - t)| (cyclic), W_t = V_t + 2 for every
t != 0.  Proof: |D cap (D - t)| = |D| - |D cap (S - t)| = |D| - (|S| - |S cap (S - t)|) per block; sum over blocks and
use sum |D_k| - sum |S_k| = 2.  So one profile V decides both halves, and parity: every lag sum of a half has the parity
of its term count, so sigma is admissible only if every V_t is even (then W_t is even too), and V_t >= 2 (target -2).

Gaussian per half [E] (terms c(y)c(y+t) are pairwise uncorrelated, so Var = term count; lattice spacing 2):
    log2 E_S = sum|S_k| + sum_t log2(2 phi(-2; V_t)),   log2 E_D = sum|D_k| + sum_t log2(2 phi(0; V_t + 2)).
This script compares, per v: the published sigma (167, 179), plain anchored arcs (n, n+1, n, n; sizes must sum to 2h - 1), random sigma (parity
ignored), and sigma optimised by simulated annealing for log2 E_S + log2 E_D under the parity constraint.
usage: python3 forensics/template_opt.py [v1,v2,...] [--steps N] [--seeds K]"""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

LN2 = math.log(2.0)


def lagterm_tables(Vmax):
    """f_S(V) = log2(2 phi(-2; V)), f_D(V) = log2(2 phi(0; V + 2)) for even V >= 2; -inf otherwise"""
    fS = np.full(Vmax + 1, -np.inf); fD = np.full(Vmax + 1, -np.inf)
    for V in range(2, Vmax + 1, 2):
        fS[V] = math.log2(2.0) - 0.5 * math.log2(2 * math.pi * V) - (4.0 / (2 * V)) / LN2
        fD[V] = math.log2(2.0) - 0.5 * math.log2(2 * math.pi * (V + 2))
    return fS, fD


def profile(S):
    """V_t, t = 0..h-1: sum over blocks of |S_k cap (S_k - t)| (cyclic), via FFT"""
    F = np.fft.fft(S.astype(float), axis=1)
    return np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real).astype(np.int64)


def evaluate(S, fS, fD):
    h = S.shape[1]; n = (h - 1) // 2; V = profile(S)[1:n + 1]
    ok = bool(((V % 2 == 0) & (V >= 2)).all())
    eS = int(S.sum()) + float(fS[np.minimum(V, len(fS) - 1)].sum()) if ok else -np.inf
    eD = int((1 - S).sum()) + float(fD[np.minimum(V, len(fD) - 1)].sum()) if ok else -np.inf
    # parity-free Gaussian (for random sigma): treat every lag as lattice-2 with the right parity
    Vc = np.maximum(V + (V % 2), 2)
    eS0 = int(S.sum()) + float(fS[Vc].sum()); eD0 = int((1 - S).sum()) + float(fD[Vc].sum())
    return dict(ok=ok, logES=eS, logED=eD, logES_parityfree=eS0, logED_parityfree=eD0,
                bad_lags=int(((V % 2) | (V < 2)).sum()), Vmean=float(V.mean()))


def arcs(h, lengths):
    S = np.zeros((4, h), np.int8)
    for k, L in enumerate(lengths): S[k, h - L:] = 1
    return S


def anneal(h, steps, seed, fS, fD, S0=None, T0=2.0, T1=0.02, pen=6.0):
    """maximise sum_t [fS(V_t) + fD(V_t)] (sizes fixed: total |S| = 2h-1), parity as a penalty per bad lag"""
    rng = np.random.default_rng(seed); n = (h - 1) // 2
    if S0 is None:
        S = np.zeros((4, h), np.int8); sizes = [(2 * h - 1) // 4 + (1 if k < (2 * h - 1) % 4 else 0) for k in range(4)]
        for k in range(4): S[k, rng.choice(h, sizes[k], replace=False)] = 1
    else: S = S0.copy()
    V = profile(S)
    def score_lags(Vt):
        good = (Vt % 2 == 0) & (Vt >= 2)
        Vc = np.minimum(np.maximum(Vt + (Vt % 2), 2), len(fS) - 1)          # parity-free value, then a penalty per bad lag
        return float((fS[Vc] + fD[Vc]).sum() - pen * (~good).sum())
    cur = score_lags(V[1:n + 1]); best = (cur, S.copy())
    ts = np.arange(1, n + 1)
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / max(1, steps - 1))
        # move a cell y out of S_k and a cell z into S_j (j = k: a swap inside one block; j != k changes the sizes,
        # which matters for parity: sum_t V_t = sum_k C(|S_k|, 2) must be even)
        k = int(rng.integers(4)); j = k if rng.random() < 0.7 else int(rng.integers(4))
        ins = np.nonzero(S[k])[0]; outs = np.nonzero(S[j] == 0)[0]
        y = int(rng.choice(ins)); z = int(rng.choice(outs))
        if j == k and y == z: continue
        Sk = S[k].astype(np.int64)
        d1 = Sk[(y + ts) % h] + Sk[(y - ts) % h]                              # pairs lost by removing y from S_k
        Sj = S[j].astype(np.int64)
        if j == k: Sj[y] = 0
        d2 = Sj[(z + ts) % h] + Sj[(z - ts) % h]                              # pairs gained by adding z to S_j
        Vn = V[1:n + 1] - d1 + d2
        new = score_lags(Vn)
        if new >= cur or rng.random() < math.exp((new - cur) / T):
            S[k, y] = 0; S[j, z] = 1; V[1:n + 1] = Vn; V[h - n:] = Vn[::-1]; cur = new
            if cur > best[0]: best = (cur, S.copy())
    return best[1]


def published(v):
    from load_originals import load_any
    from halfshift_b1 import core_seqs
    name = {167: 'matrix_H668_v167_h_b1', 179: 'matrix_H716_v179_h_b1'}[v]
    H = load_any(name); m = v - 1; h = m // 2
    X = np.array([x for _, x in core_seqs(H, m)])
    return (X[:, :h] * X[:, h:] > 0).astype(np.int8)


def main():
    a = sys.argv[1:]
    vs = [int(x) for x in a[0].split(',')] if a and not a[0].startswith('--') else [167, 179, 251, 307, 419, 443, 523, 571, 647, 739]
    steps = int(a[a.index('--steps') + 1]) if '--steps' in a else 200000
    seeds = int(a[a.index('--seeds') + 1]) if '--seeds' in a else 2
    out = []
    for v in vs:
        h = (v - 1) // 2; n = (h - 1) // 2
        fS, fD = lagterm_tables(8 * h * h)
        rec = dict(v=v, h=h)
        if v in (167, 179): rec['published'] = evaluate(published(v), fS, fD)
        rec['arcs'] = evaluate(arcs(h, [n, n + 1, n, n]), fS, fD)            # sizes sum to 2h - 1
        rng = np.random.default_rng(5); R = []
        for _ in range(20):
            S = np.zeros((4, h), np.int8)
            for k, L in enumerate([n, n + 1, n, n]): S[k, rng.choice(h, L, replace=False)] = 1
            R.append(evaluate(S, fS, fD))
        rec['random_parityfree'] = dict(logES=round(float(np.mean([r['logES_parityfree'] for r in R])), 1),
                                        logED=round(float(np.mean([r['logED_parityfree'] for r in R])), 1),
                                        bad_lags=float(np.mean([r['bad_lags'] for r in R])))
        bestS, bestv, beststart = None, -np.inf, None
        starts = ('published', 'arcs', 'random') if v in (167, 179) else ('arcs', 'random')
        for sd in range(seeds):
            for start in starts:
                S0 = published(v) if start == 'published' else arcs(h, [n, n + 1, n, n]) if start == 'arcs' else None
                S = anneal(h, steps, 100 * sd + (start == 'arcs'), fS, fD, S0=S0)
                ev = evaluate(S, fS, fD)
                val = ev['logES'] + ev['logED'] if ev['ok'] else -np.inf
                rec.setdefault('runs', []).append(dict(start=start, seed=sd, parity_ok=ev['ok'], bad_lags=ev['bad_lags'],
                                                       logES_parityfree=round(ev['logES_parityfree'], 1),
                                                       logED_parityfree=round(ev['logED_parityfree'], 1)))
                if val > bestv: bestv, bestS, beststart = val, S, start
        ev = evaluate(bestS, fS, fD) if bestS is not None else None
        rec['annealed'] = ev; rec['annealed_from'] = beststart
        if bestS is not None:
            # how arc-like is the optimum: sign changes of each block's indicator around Z_h
            rec['annealed_sign_changes'] = [int((bestS[k] != np.roll(bestS[k], 1)).sum()) for k in range(4)]
            rec['annealed_sizes'] = [int(x) for x in bestS.sum(1)]
            np.save(os.path.join(HERE, '..', 'logs', 'gap', f'template_opt_sigma_v{v}.npy'), bestS)
        out.append(rec)
        f = lambda e: 'inf' if e is None else (f"{e['logES']:+.1f}/{e['logED']:+.1f}" if e['ok'] else f"parity-bad ({e['bad_lags']} lags)")
        print(json.dumps(dict(v=v, h=h, published=f(rec.get('published')), arcs=f(rec['arcs']),
                              random_parityfree=f"{rec['random_parityfree']['logES']:+.1f}/{rec['random_parityfree']['logED']:+.1f}",
                              annealed=f(ev), from_=rec['annealed_from'], sign_changes=rec.get('annealed_sign_changes'),
                              sizes=rec.get('annealed_sizes'))), flush=True)
    json.dump(out, open(os.path.join(HERE, '..', 'logs', 'gap', 'template_opt.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
