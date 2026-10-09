r"""Predicted CORE yield per designed gadget in a Paley + multiplier slice [E]: log2 E(stage 1) + log2 E(stage 2 | y),
with alpha sampled from the gadget's exact alpha lists and u uniform on the twisted orbits (same Gaussian estimator as
theory/paley_compfix.py).  Used to rank gadgets for the 573 route (report section 7.5); calibrated on 419, where the
published gadget and the from-scratch gadget are scored the same way.
usage: python3 theory/paley_coreyield.py v p n mu K Nm [n_gadgets] [seed]"""
import sys, os, json, math, itertools as it
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import paley_slicepipe as SP, paley_b1pipe as BP
from paley_compfix import gauss_logp

def core_yield(W, Jp, v, p, n, orbs, rng, N=3000):
    AL = SP.gadget_alpha_lists(W, Jp, p, n, orbs)
    if AL is None: return None
    m = n * p; TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m - 4 * (p - 1); half = n // 2
    bits1 = sum(math.log2(len(a)) for a in AL)
    As = np.stack([AL[k][rng.integers(len(AL[k]), size=N)] for k in range(4)], 1)          # (N, 4, n)
    Y = As + (p - 1) * W[None]
    S1 = BP.pafc(Y.reshape(-1, n)).reshape(N, 4, n).sum(1)[:, 1:half + 1]
    try: c1, _, _ = gauss_logp(S1, TJ[1:half + 1], 'S1')
    except AssertionError: return None
    # stage 2 at FIXED y (u is searched per y): average the per-y Gaussian yield over a few sampled y's
    s2 = [stage2_fixed_y(As[i], W, p, n, orbs, rng, N) for i in range(3)]
    s2 = [x for x in s2 if x is not None]
    if not s2: return None
    st2 = float(np.mean(s2))
    tw = [[o for o in orbs if W[k, o[0]] == 0] for k in range(4)]
    return dict(stage1=round(bits1 - c1, 1), stage2=round(st2, 1), core=round(bits1 - c1 + st2, 1),
                alpha_log2=[round(math.log2(len(a)), 1) for a in AL], u_log2=[len(t) for t in tw])

def stage2_fixed_y(A, W, p, n, orbs, rng, N=3000):
    half = n // 2; beta = A - W; tw = [[o for o in orbs if W[k, o[0]] == 0] for k in range(4)]
    U = np.zeros((N, 4, n), np.int64)
    for k in range(4):
        for o in tw[k]: U[:, k, o] = rng.choice((1, -1), N)[:, None]
    P = lambda M: np.rint(np.fft.ifft(np.abs(np.fft.fft(M, axis=-1)) ** 2, axis=-1).real).astype(np.int64)
    C = lambda f, g: np.rint(np.fft.ifft(np.conj(np.fft.fft(f, axis=-1)) * np.fft.fft(g, axis=-1), axis=-1).real).astype(np.int64)
    S2 = (P(beta[None]).sum(1) + p * P(U).sum(1))[:, 1:half + 1]
    X = (C(np.broadcast_to(beta, U.shape), U) - C(U, np.broadcast_to(beta, U.shape))).sum(1)[:, 1:half]
    try: c2, _, _ = gauss_logp(np.concatenate([S2, X], 1), np.zeros(S2.shape[1] + X.shape[1], np.int64), 'S2X')
    except AssertionError: return None
    return sum(len(t) for t in tw) - c2

def main(v, p, n, mu, K, Nm, n_g=200, seed=0, B=2, kmax=13):
    rng = np.random.default_rng(seed); orbs = SP.orbits(n, mu); r = len(orbs); osz = [len(o) for o in orbs]
    BG = SP.block_gadgets(n, orbs, kmax, B)
    BG = {sz: [g for g in lst if abs(g[1]) <= 1] for sz, lst in BG.items()}; BG = {sz: l for sz, l in BG.items() if l}
    sizes = sorted(BG); comps = [c for c in it.product(sizes, repeat=4) if sum(c) == K]
    rs = [2, 0, 0, 0]; got = []; tried = 0
    while len(got) < n_g and tried < 10 ** 6:
        tried += 1
        c = comps[rng.integers(len(comps))]
        pick = [BG[c[k]][rng.integers(len(BG[c[k]]))] for k in range(4)]
        W = np.array([pk[0] for pk in pick]); ts = [pk[1] for pk in pick]; s_ = sum(pk[2] for pk in pick)
        if not all(int(x) in (-2, -1, 0, 1) for x in s_[1:]): continue
        if int(sum(t * t for t in ts)) != K + int(s_[1:].sum()): continue
        num = 4 * v - sum((rs[k] - p * ts[k]) ** 2 for k in range(4))
        if num % p or not BP.sq4(num // p, [n - int((W[k] != 0).sum()) for k in range(4)]): continue
        cons = [(k, j) for k in range(4) for j in range(r) if W[k, orbs[j][0]] != 0]; Jp = None
        for _ in range(200):
            sel = [x for x in cons if rng.random() < Nm / max(K, 1)]
            if sum(osz[j] for (_, j) in sel) == Nm:
                Jp = np.zeros((4, n), bool)
                for (k, j) in sel: Jp[k, orbs[j]] = True
                break
        if Jp is None: continue
        y = core_yield(W, Jp, v, p, n, orbs, rng)
        if y is None: continue
        got.append(y)
    cores = sorted([g['core'] for g in got], reverse=True)
    out = dict(v=v, mu=mu, K=K, Nminus=Nm, gadgets=len(got), tried=tried,
               core_top5=cores[:5], core_p90=float(np.percentile(cores, 90)) if cores else None,
               stage2_median=float(np.median([g['stage2'] for g in got])) if got else None,
               best=max(got, key=lambda g: g['core']) if got else None)
    print(json.dumps(out), flush=True); return out

if __name__ == "__main__":
    a = list(map(int, sys.argv[1:]))
    main(*a[:6], *(a[6:8]))
