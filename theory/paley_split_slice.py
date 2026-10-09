r"""Stage split [E] of the b = 1 Paley-coset ansatz restricted to a MULTIPLIER SLICE: labels (alpha, type, sign) constant
on the orbits of g -> mu*g on a cyclic G' = Z_n.  Paley(p) on the Z_p factor is itself the orbit slice of the multiplier
group {mu = QR mod p, mu = 1 mod n} (report section 7.5), so "Paley + mult:mu" is the multiplier slice of order
(p-1)/2 * ord(mu): for 419 (p = 11, mu = 7 on Z38) that is the order-15 slice of the published core; for 573 (p = 11,
mu = 29 on Z52) the repo's d = 15 slice.

Same estimator and output as theory/paley_split_b1.py (stage 1 = (gadget, J/J', w, alpha) on orbits vs S1; stage 2 = u on
twisted orbits vs S2 and X, given one sampled y), with exact orbit-weighted counting:
  #(constant-orbit sets of total size K with J' part of size N-) from a generating function, alpha counted by an exact
  DP over orbit sizes with the required block sum.
usage: python3 theory/paley_split_slice.py v p n mu [N]"""
import sys, os, json, math
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_compfix import gauss_logp
from paley_wscreen import KN_options
from paley_split_b1 import paf, cc

def orbits(n, mu):
    seen, out = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * mu) % n
        seen.update(o); out.append(sorted(o))
    return out

def count_sets(sizes, K, Nm):
    """#(subset S of slots with sum size = K, subset J of S with sum size = Nm), slots = 4 copies of the orbits"""
    f = defaultdict(int); f[(0, 0)] = 1
    for s in sizes:
        g = defaultdict(int)
        for (a, c), val in f.items():
            g[(a, c)] += val
            if a + s <= K:
                g[(a + s, c)] += val                       # constant, J
                if c + s <= Nm: g[(a + s, c + s)] += val    # constant, J'
        f = g
    return f[(K, Nm)]

def alpha_count_and_sample(rng, osz, need):
    """# of +-1 vectors a on orbits (sizes osz) with sum osz*a = need, and one uniform sample (None if 0)"""
    # DP from the back: N[i][t] = #ways for orbits i.. to reach t
    tot = sum(osz); N = [defaultdict(int) for _ in range(len(osz) + 1)]; N[len(osz)][0] = 1
    for i in range(len(osz) - 1, -1, -1):
        for t, c in N[i + 1].items():
            N[i][t + osz[i]] += c; N[i][t - osz[i]] += c
    cnt = N[0].get(need, 0)
    if cnt == 0: return 0, None
    a = []; t = need
    for i in range(len(osz)):
        cp = N[i + 1].get(t - osz[i], 0); cm = N[i + 1].get(t + osz[i], 0)
        if rng.random() * (cp + cm) < cp: a.append(1); t -= osz[i]
        else: a.append(-1); t += osz[i]
    return cnt, a

def main(v, p, n, mu, N=20000, seed=1):
    m = n * p; assert m == v - 1
    rng = np.random.default_rng(seed)
    orbs = orbits(n, mu); r = len(orbs); osz = [len(o) for o in orbs]
    TJ = np.full(n, -4 * p, np.int64); TJ[0] = 4 * m - 4 * (p - 1)
    half = n // 2; rs = [2, 0, 0, 0]
    for (K, Nm, Np) in KN_options(v, 1, p, n):
        nsets = count_sets(osz * 4, K, Nm)
        if nsets == 0: print(json.dumps(dict(v=v, mu=mu, K=K, Nminus=Nm, note='no orbit-union of size K')), flush=True); continue
        A = []; W = []; la = []; tries = 0; n_reach = 0
        while len(A) < N and tries < 60 * N:
            tries += 1
            perm = rng.permutation(4 * r); tot = 0; const = []
            for s_ in perm:
                sz = osz[s_ % r]
                if tot + sz <= K and rng.random() < 0.5: const.append(s_); tot += sz
            if tot != K: continue
            # J' subset of total size Nm
            jp = []; tj = 0
            for s_ in rng.permutation(const):
                sz = osz[s_ % r]
                if tj + sz <= Nm and rng.random() < 0.7: jp.append(s_); tj += sz
            if tj != Nm: continue
            n_reach += 1
            a = np.zeros((4, n), np.int64); w = np.zeros((4, n), np.int64); ok = True; lg = 0.0
            isc = np.zeros(4 * r, bool); isc[const] = True; isj = np.zeros(4 * r, bool); isj[jp] = True
            for s_ in const:
                k, j = divmod(s_, r); sg = 1 if rng.random() < 0.5 else -1
                w[k, orbs[j]] = sg; a[k, orbs[j]] = -sg if isj[s_] else sg
            for k in range(4):
                tw = [j for j in range(r) if not isc[k * r + j]]
                need = rs[k] - (p - 1) * int(w[k].sum()) - int(a[k].sum())
                cnt, smp = alpha_count_and_sample(rng, [osz[j] for j in tw], need)
                if cnt == 0: ok = False; break
                lg += math.log2(cnt)
                for j, s_ in zip(tw, smp): a[k, orbs[j]] = s_
            if not ok: continue
            A.append(a); W.append(w); la.append(lg)
        if len(A) < 500:
            print(json.dumps(dict(v=v, mu=mu, K=K, Nminus=Nm, note='row sums / sizes rarely feasible', samples=len(A))), flush=True); continue
        A = np.array(A); W = np.array(W); la = np.array(la)
        Y = A + (p - 1) * W
        S1 = paf(Y).sum(1)[:, 1:half + 1]
        c1, r1, d1 = gauss_logp(S1, TJ[1:half + 1], 'S1')
        # H1 = log2 #(constant-orbit set, J' part) [exact] + #constant orbits (w signs) + log2(fraction of draws whose
        # row sums are feasible) + E log2 #alpha  (the set draws are not exactly uniform: [E])
        n_corb = float(np.mean([sum(1 for k in range(4) for j in range(r) if W[i][k, orbs[j][0]] != 0) for i in range(min(len(W), 2000))]))
        H1 = math.log2(nsets) + n_corb + math.log2(len(A) / n_reach) + float(np.mean(la))
        a0, w0 = A[0], W[0]; beta = a0 - w0; tw = (w0 == 0)
        U = np.zeros((N, 4, n), np.int64)
        for k in range(4):
            for j in range(r):
                if tw[k, orbs[j][0]]:
                    U[:, k, orbs[j]] = rng.choice((1, -1), N)[:, None]
        S2 = (paf(beta[None]).sum(1) + p * paf(U).sum(1))[:, 1:half + 1]
        if p % 4 == 1:     # chi even: the chi-coefficient is symmetric, CC(u, beta) + CC(beta, u), lags 0..n/2
            X = (cc(U, np.broadcast_to(beta, U.shape)) + cc(np.broadcast_to(beta, U.shape), U)).sum(1)[:, 0:half + 1]
        else:
            X = (cc(U, np.broadcast_to(beta, U.shape)) - cc(np.broadcast_to(beta, U.shape), U)).sum(1)[:, 1:half]
        c2, r2, d2 = gauss_logp(np.concatenate([S2, X], 1), np.zeros(S2.shape[1] + X.shape[1], np.int64), 'S2X')
        H2 = sum(1 for k in range(4) for j in range(r) if tw[k, orbs[j][0]])
        rec = dict(v=v, p=p, n=n, mu=mu, orbits=r, K=K, Nminus=Nm, Nplus=Np, H1=round(H1, 1), C1=round(c1, 1),
                   stage1_log2E=round(H1 - c1, 1), H2=H2, C2=round(c2, 1), stage2_log2E_per_y=round(H2 - c2, 1),
                   core_log2E=round(H1 + H2 - c1 - c2, 1), eqs1=r1, eqs2=r2, samples=len(A))
        print(json.dumps(rec), flush=True)

if __name__ == "__main__":
    v, p, n, mu = map(int, sys.argv[1:5]); N = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    main(v, p, n, mu, N)
