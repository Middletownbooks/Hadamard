r"""Cost of a b >= 3 Paley-coset option once the border-forced Q-compression and a w-gadget are fixed (p = 3 mod 4).

Facts used (report/paley_design_answer.md sections 2-4, all exact):
  * S1 at the Q-level forces sum_k |yhat_Q,k(xi)|^2 = 4b on Q^ \ 1: every block's compression is "constant + O(1)".
    We take the single-deviation family y_Q,k = c_k 1 + d_k delta_{q*} (the published 359 shape), sum d_k^2 = 4b,
    sum_k (c c_k + d_k)^2 = 4v, and the w-gadget balanced inside H'-cosets (so alpha_Q = y_Q).
  * Norm identity (p+1) K - 4 N_- = 4 n' - 4 b fixes (K, N_-) pairs.
  * Gadget family "one balanced pattern per block": block k carries a +-1 pattern of size K_k (even) with zero sum
    inside ONE H'-coset (so PAF_w lives on H'); alpha at the gadget positions is -w (J') or +w (J) as N_- dictates.
For each (compression, gadget) sampled, alpha is drawn uniformly with the prescribed coset sums and u uniformly on the
twisted positions; the Gaussian density of the full reduced system (and of stage 1 alone) at the target gives
log2 E for that fixed (compression, gadget).
usage: python3 theory/paley_option_cost.py v c p mods... [--Hgen a.b,..] [--gadgets G] [--N 60000]"""
import sys, os, json, math, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
from paley_ansatz import LABELS, LC
import paley_ls as LSm
from paley_compfix import gauss_logp

def cosets(P):
    Hs = sorted(P.H); n = P.n; rep = -np.ones(n, np.int64); q = 0; reps = []
    for g in range(n):
        if rep[g] < 0:
            for h in Hs: rep[P.add[g, h]] = q
            reps.append(g); q += 1
    return rep, reps

def compressions(v, b, c, hsz, maxn=6):
    """single-deviation compressions (c_k, d_k): y_Q = c_k on c-1 cosets, c_k + d_k on q*; alpha_Q parity = hsz mod 2,
    |values| <= hsz; returns list sorted by alpha entropy (largest first)"""
    vals = [x for x in range(-hsz, hsz + 1) if (x - hsz) % 2 == 0]
    ent = lambda x: math.log2(math.comb(hsz, (hsz - x) // 2))
    out = []
    for cs in it.combinations_with_replacement(vals, 4):
        for ds in it.product([d for d in range(-2 * hsz, 2 * hsz + 1, 2)], repeat=4):
            if sum(d * d for d in ds) != 4 * b: continue
            if any(abs(ck + dk) > hsz for ck, dk in zip(cs, ds)): continue
            if sum((c * ck + dk) ** 2 for ck, dk in zip(cs, ds)) != 4 * v: continue
            e = sum((c - 1) * ent(ck) + ent(ck + dk) for ck, dk in zip(cs, ds))
            out.append((round(e, 2), list(cs), list(ds)))
    out = sorted(set((e, tuple(cs), tuple(ds)) for e, cs, ds in out), reverse=True)
    return [(e, list(cs), list(ds)) for e, cs, ds in out[:maxn]]

def norm_pairs(v, b, p, n):
    return [(K, (p + 1) * K - 4 * n + 4 * b) for K in range(0, 4 * n + 1) if ((p + 1) * K - 4 * n + 4 * b) % 4 == 0
            and 0 <= ((p + 1) * K - 4 * n + 4 * b) // 4 <= K]

def run(v, c, p, mods, Hgen=None, n_gadgets=4, N=60000, seed=1, comp_index=0):
    b = c - 1; P = LSm.Prob(v, b, p, mods, Hgen); n = P.n; hsz = len(P.H)
    rep, reps = cosets(P); rng = np.random.default_rng(seed)
    comps = compressions(v, b, c, hsz)
    if not comps: return dict(v=v, c=c, p=p, mods=list(mods), feasible=False, why="no single-deviation compression")
    ent, cs, ds = comps[min(comp_index, len(comps) - 1)]
    KN = [(K, Nm // 4) for K, Nm in norm_pairs(v, b, p, n)]
    res = []
    def place_gadget():
        """random feasible gadget: K (norm identity), even part sizes per block, one coset per block, J'/J choice"""
        for _ in range(5000):
            K, Nm = KN[rng.integers(len(KN))]
            if K % 2: continue
            parts = rng.multinomial(K // 2, [0.25] * 4) * 2
            if (parts > hsz).any(): continue
            plan = []; Jm_left = Nm; ok = True
            for k in range(4):
                if parts[k] == 0: plan.append(None); continue
                good = None
                for q in rng.permutation(len(reps)):
                    pos_q = np.nonzero(rep == q)[0]; yq = cs[k] + (ds[k] if q == 0 else 0)
                    g = rng.choice(pos_q, parts[k], replace=False)
                    sg = np.array([1] * (parts[k] // 2) + [-1] * (parts[k] // 2)); rng.shuffle(sg)
                    for nJm in rng.permutation(min(Jm_left, parts[k]) + 1):
                        # J' positions: alpha = -w ; J: alpha = +w ; choose which positions are J' (first nJm after shuffle)
                        agad = np.where(np.arange(parts[k]) < nJm, -sg, sg)
                        need = yq - int(agad.sum()); nf = len(pos_q) - parts[k]
                        if (nf - need) % 2 == 0 and abs(need) <= nf:
                            good = (q, g, sg, agad, nJm); break
                    if good: break
                if good is None: ok = False; break
                Jm_left -= good[4]; plan.append(good)
            if ok and Jm_left == 0: return K, Nm, parts, plan
        return None
    for gi in range(n_gadgets):
        gd = place_gadget()
        if gd is None: res.append(dict(infeasible="no gadget placement found")); continue
        K, Nm, parts, plan = gd
        labs = np.zeros((N, 4, n), np.int64); abits = 0.0; ubits = 0
        for k in range(4):
            wv = np.zeros(n, np.int64); agad = np.zeros(n, np.int64)
            if plan[k] is not None:
                q_, g, sg, ag, _ = plan[k]; wv[g] = sg; agad[g] = ag
            for q in range(len(reps)):
                pos = np.nonzero(rep == q)[0]
                yq = cs[k] + (ds[k] if q == 0 else 0)
                fixed = [x for x in pos if wv[x] != 0]; free = [x for x in pos if wv[x] == 0]
                need = yq - (int(agad[fixed].sum()) if fixed else 0); nf = len(free)
                assert (nf - need) % 2 == 0 and abs(need) <= nf
                nminus = (nf - need) // 2; abits += math.log2(math.comb(nf, nminus))
                r = np.argsort(rng.random((N, nf)), axis=1)[:, :nminus]
                al = np.ones((N, nf), np.int64); np.put_along_axis(al, r, -1, axis=1)
                for j, x in enumerate(free):
                    u = rng.choice(np.array([1, -1]), N)
                    labs[:, k, x] = np.where(u == 1, np.where(al[:, j] == 1, LABELS.index((1, 1, 1)), LABELS.index((-1, 1, 1))),
                                             np.where(al[:, j] == 1, LABELS.index((1, 1, -1)), LABELS.index((-1, 1, -1))))
                    ubits += 1
                for x in fixed: labs[:, k, x] = LABELS.index((int(agad[x]), 0, int(wv[x])))
        R = np.zeros((N, n, 3), np.int64)
        for k in range(4):
            Ak = LC[labs[:, k]].astype(np.int64)
            for h in range(n):
                R[:, h, :] += np.einsum('ngi,ngj,ijc->nc', Ak, Ak[:, P.add[:, h]], P.C)
        T = P.T
        try:
            cb, eq, _ = gauss_logp(R.reshape(N, -1)[:, 1:], T.reshape(-1)[1:], 'full')
            Y = (R[:, :, 1] + R[:, :, 2]) // 2
            cb1, eq1, _ = gauss_logp((R[:, :, 0] + (p - 1) * Y)[:, 1:], (T[:, 0] + (p - 1) * T[:, 1])[1:], 's1')
        except AssertionError as e:
            res.append(dict(K=int(K), Nminus=int(Nm), parts=parts.tolist(), infeasible=str(e)[:80])); continue
        res.append(dict(K=int(K), Nminus=int(Nm), parts=parts.tolist(), alpha_bits=round(abits, 1), u_bits=int(ubits),
                        cbits=round(cb, 1), eqs=eq, log2E=round(abits + ubits - cb, 1),
                        stage1_cbits=round(cb1, 1), log2E_stage1=round(abits - cb1, 1),
                        stage2_deficit=round(ubits - (cb - cb1), 1)))
    return dict(v=v, c=c, b=b, p=p, mods=list(mods), n=n, H_size=hsz, KN_options=[list(t) for t in KN],
                compression=dict(alpha_entropy=ent, c=cs, d=ds), n_compressions_listed=len(comps), samples=res)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--Hgen": None, "--gadgets": "4", "--N": "60000", "--comp": "0", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, c, p = pos[:3]; mods = tuple(pos[3:])
    Hgen = None if opt["--Hgen"] is None else [tuple(map(int, e.split('.'))) for e in opt["--Hgen"].split(',')]
    print(json.dumps(run(v, c, p, mods, Hgen, int(opt["--gadgets"]), int(opt["--N"]), int(opt["--seed"]), int(opt["--comp"]))), flush=True)
