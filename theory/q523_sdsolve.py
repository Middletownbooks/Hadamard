r"""Exact stage-2 solver for the generalised-Paley slice by the Z_2 (S/D) split with a Z_3 class split and big-fibre
parents (report/paley_523_answer.md, section 3).  F = 2 F2 with 3 | F2 (F2 = 9 at 523 and 235).

Stage 2 given y is  sum_k P_k = T  in O_K[Z_F]  (theory/q523_lib.py; T = -4v at lag 0 in eta-coordinates).  With pairs
(g, g + F2), S(g) = Z(g) + Z(g + F2) and D(g) = (-1)^g (Z(g) - Z(g + F2)) on Z_F2 (g < F2):
    sum_k P(S_k) = T_S,  sum_k P(D_k) = T_D   on Z_F2   (T_S(h) = sum_{h' = h mod F2} T(h'), T_D the (-1)^h' twist),
exactly equivalent to stage 2.  Per pair the (c, c') choices give distinct S-values s with D-options; m(s) = #options.
Rows: the pairs g = r (mod 3).  Z3 level: the row sums phi_r of S; its equations are sum_k P(Phi_k) = T_S compressed to
Z_3.  Pipeline:
  1. per block and row, group the row configurations (one S-value per pair) by row sum; weight = prod m(s) (= number
     of full configurations), sampling weight = weight^beta;
  2. parent: one row sum per row and block, chosen among the heaviest groups, solving the Z3 equations exactly
     (2+2 join over candidate parents);
  3. S-step: per block a product list of weighted row samples inside the parent's fibre; exact 2+2 join on the S-level
     signature (lags 0..F2/2, eta-coordinates; numba, residue passes, every hit re-checked);
  4. D-check per S-solution: all D-options per block; Z3 class split on D (row sums of D), exact join on the D3
     equations, then on the full D-level signature;
  5. assembly Z -> (alpha, c) -> x, exact 2-D FFT check; cyclic CRT core file for bordered GS completion.
usage: python3 theory/q523_sdsolve.py F q s y.json [--K 64] [--beta 1.5] [--top 24] [--seed 0] [--planted] [--out f]"""
import sys, os, json, math, time, itertools as it
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, pz, expand_x, paf2, to_cyclic, target_P

# ------------------------------------------------------------------ signatures
@numba.njit(cache=True)
def seq_sig(X, T, nl):
    """X (N, n, s) int64 eta-coords of sequences on Z_n -> (N, nl * s) PAF signatures at lags 0..nl-1:
    P(h)_l = sum_ij T[i, j, l] sum_g X_i(g + h) X_j(g)"""
    N, n, s = X.shape
    out = np.zeros((N, nl * s), np.int64)
    cc = np.zeros((s, s), np.int64)
    for a in range(N):
        for h in range(nl):
            for i in range(s):
                for j in range(s): cc[i, j] = 0
            for g in range(n):
                g2 = (g + h) % n
                for i in range(s):
                    xi = X[a, g2, i]
                    if xi == 0: continue
                    for j in range(s): cc[i, j] += xi * X[a, g, j]
            for l in range(s):
                acc = 0
                for i in range(s):
                    for j in range(s): acc += T[i, j, l] * cc[i, j]
                out[a, h * s + l] = acc
    return out

# ------------------------------------------------------------------ 2+2 join (residue passes; pattern of theory/p571_hsolve.py)
@numba.njit(cache=True)
def join4(kA, kB, kC, kD, rA, rB, rC, rD, rT, P, bstart, border, dstart, dorder, kT, tbits, maxhits, c_lo, c_hi):
    size = 1 << tbits; mask = size - 1
    tkey = np.zeros(size, np.uint64); ta = np.zeros(size, np.int32); tb = np.zeros(size, np.int32)
    stamp = np.zeros(size, np.int32)
    hits = np.zeros((maxhits, 4), np.int64); nh = 0; overflow = 0; probes = 0
    cap = (size * 85) // 100
    for c in range(c_lo, c_hi):
        st = c + 1; filled = 0
        for a in range(kA.shape[0]):
            rb = (c - rA[a]) % P
            for t in range(bstart[rb], bstart[rb + 1]):
                b = border[t]; key = kA[a] + kB[b]
                if filled >= cap:
                    overflow += 1; continue
                slot = int(key >> np.uint64(64 - tbits)) & mask
                while stamp[slot] == st: slot = (slot + 1) & mask
                stamp[slot] = st; tkey[slot] = key; ta[slot] = a; tb[slot] = b; filled += 1
        probes += filled
        cc = (rT - c) % P
        for ci in range(kC.shape[0]):
            rd = (cc - rC[ci]) % P
            for t in range(dstart[rd], dstart[rd + 1]):
                d = dorder[t]; need = kT - kC[ci] - kD[d]
                slot = int(need >> np.uint64(64 - tbits)) & mask
                while stamp[slot] == st:
                    if tkey[slot] == need:
                        if nh < maxhits:
                            hits[nh, 0] = ta[slot]; hits[nh, 1] = tb[slot]; hits[nh, 2] = ci; hits[nh, 3] = d
                        nh += 1
                    slot = (slot + 1) & mask
                probes += 1
    return hits[:min(nh, maxhits)], nh, overflow, probes

def next_prime(x):
    x = max(2, int(x))
    while any(x % p == 0 for p in range(2, int(x ** 0.5) + 1)): x += 1
    return x

def csr_res(r, P):
    order = np.argsort(r, kind='stable').astype(np.int64)
    start = np.zeros(P + 1, np.int64); np.add.at(start, r + 1, 1)
    return np.cumsum(start), order

def join_all(sigs, target, seed=0, tbits=22, maxhits=1 << 20, max_passes=None):
    """exact 2+2 join of four (n_k, D) int64 signature arrays: all index quadruples with sum == target (verified).
    max_passes: stop after that many residue classes (a random subset of the pair space, for sampling)."""
    D = sigs[0].shape[1]; rng = np.random.default_rng(seed)
    wts = rng.integers(1, 2 ** 63, D, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    rho = rng.integers(1, 1 << 30, D, dtype=np.int64)
    hk = lambda S: (S.astype(np.int64).astype(np.uint64) * wts[None, :]).sum(1)
    keys = [hk(s) for s in sigs]; kT = hk(target[None, :])[0]
    nAB = len(keys[0]) * len(keys[1])
    P = next_prime(max(2, nAB / ((1 << tbits) * 0.5)))
    while True:
        res = lambda S: ((S.astype(np.int64) % P) * (rho % P)[None, :] % P).sum(1) % P
        rs = [res(s) for s in sigs]; rT = int(res(target[None, :])[0])
        bs, bo = csr_res(rs[1], P); ds, do = csr_res(rs[3], P)
        hi = P if max_passes is None else min(P, max_passes)
        hits, nh, ov, pr = join4(keys[0], keys[1], keys[2], keys[3], rs[0], rs[1], rs[2], rs[3], rT, P, bs, bo, ds, do,
                                 kT, tbits, maxhits, 0, hi)
        if ov == 0: break
        P = next_prime(2 * P)
    st = dict(probes=int(pr), passes=int(hi), P=int(P), hash_hits=int(nh), truncated=bool(nh > maxhits))
    if len(hits) == 0: return np.zeros((0, 4), np.int64), st
    tot = sigs[0][hits[:, 0]] + sigs[1][hits[:, 1]] + sigs[2][hits[:, 2]] + sigs[3][hits[:, 3]]
    ok = (tot == target[None, :]).all(1)
    return hits[ok], st

# ------------------------------------------------------------------ instance
class SDInst:
    def __init__(self, F, q, s, Y, T=None):
        assert F % 2 == 0 and (F // 2) % 3 == 0
        self.F, self.q, self.s, self.F2 = F, q, s, F // 2
        self.cy = Cyclo(q, s); self.v = F * q + 1; self.Y = np.array(Y, np.int64)
        self.T = target_P(F, self.v, s) if T is None else np.array(T, np.int64)
        F2 = self.F2; h = np.arange(F)
        self.TS = np.stack([self.T[h % F2 == a].sum(0) for a in range(F2)])
        self.TD = np.stack([(((-1) ** h)[:, None] * self.T)[h % F2 == a].sum(0) for a in range(F2)])
        self.nl = F2 // 2 + 1
        self.rows = [list(range(r, F2, 3)) for r in range(3)]
        # pair tables: per block, per pair g < F2: Svals (nS, s), mult (nS,), Dopts list of (m, s) arrays
        self.pairs = []
        for k in range(4):
            alph = [zcoords(self.cy.alphabet(int(y))) for y in Y[k]]
            blk = []
            for g in range(F2):
                sg = 1 if g % 2 == 0 else -1; byS = {}
                for a in alph[g]:
                    for b in alph[g + F2]:
                        byS.setdefault(tuple((a + b).tolist()), []).append(tuple((sg * (a - b)).tolist()))
                keys = sorted(byS)
                blk.append(dict(S=np.array(keys, np.int64), m=np.array([len(byS[t]) for t in keys], np.int64),
                                D=[np.array(byS[t], np.int64) for t in keys]))
            self.pairs.append(blk)
        # Z3 targets: compress T_S (Z_F2 -> Z_3)
        self.TS3 = np.stack([self.TS[np.arange(F2) % 3 == r].sum(0) for r in range(3)])
        self.TD3 = np.stack([self.TD[np.arange(F2) % 3 == r].sum(0) for r in range(3)])

    # -------------------------------------------------------------- rows
    def row_groups(self, k, r):
        """dict row-sum tuple -> (combos (n, 3) pair-S indices, weights (n,))"""
        prs = [self.pairs[k][g] for g in self.rows[r]]
        idx = np.array(list(it.product(*[range(len(p['S'])) for p in prs])), np.int64)
        sums = sum(prs[t]['S'][idx[:, t]] for t in range(len(prs)))
        w = np.prod(np.stack([prs[t]['m'][idx[:, t]] for t in range(len(prs))], 1), 1).astype(float)
        out = {}
        keys = [tuple(x) for x in sums.tolist()]
        order = sorted(range(len(keys)), key=lambda i: keys[i])
        for i in order: out.setdefault(keys[i], []).append(i)
        return {kk: (idx[v], w[v]) for kk, v in out.items()}

    def tsig3(self):
        return np.concatenate([self.TS3[0], self.TS3[1]])

    def tsigS(self):
        return self.TS[:self.nl].ravel()

    def tsigD(self):
        return self.TD[:self.nl].ravel()

    def tsigD3(self):
        return np.concatenate([self.TD3[0], self.TD3[1]])

def group_score(w, beta):
    """effective size 2^H of the beta-weighted distribution within a row group"""
    p = w ** beta; p = p / p.sum()
    return float(-(p * np.log2(p)).sum())

def choose_parent(inst, groups, beta=1.5, top=24, seed=0, want=1, log=None):
    """candidate parents per block: per row, `top` distinct row sums drawn from their beta-weighted marginal (typical
    parents: the Z3 lag-0 equation pins the row-sum energy, so only typical fibres occur); exact Z3 join; the solutions
    are ranked by the summed effective fibre sizes (group_score)"""
    T = inst.cy.T; rng = np.random.default_rng(seed)
    cands = []; sigs = []; scores = []
    for k in range(4):
        per_row = []
        for r in range(3):
            g = groups[k][r]; keys = list(g)
            pw = np.array([(g[kk][1] ** beta).sum() for kk in keys]); pw /= pw.sum()
            pick = rng.choice(len(keys), size=min(top, len(keys)), replace=False, p=pw)
            per_row.append([(group_score(g[keys[i]][1], beta), keys[i]) for i in pick])
        C = list(it.product(*per_row))
        rs = np.array([[np.array(c[r][1]) for r in range(3)] for c in C], np.int64)     # (n, 3, s)
        cands.append([tuple(c[r][1] for r in range(3)) for c in C])
        scores.append(np.array([sum(c[r][0] for r in range(3)) for c in C]))
        sigs.append(seq_sig(rs, T, 2))
    hits, st = join_all(sigs, inst.tsig3(), seed=seed, tbits=20)
    if log: log(dict(stage="Z3 parents", candidates=[len(c) for c in cands], solutions=int(len(hits)), **st))
    if len(hits) == 0: return []
    tot = sum(scores[k][hits[:, k]] for k in range(4))
    order = np.argsort(-tot)
    return [(tuple(cands[k][hits[i, k]] for k in range(4)), [float(scores[k][hits[i, k]]) for k in range(4)]) for i in order[:want]]

def fibre_list(inst, groups, k, parent_k, K, beta, rng):
    """product list for block k: per row up to K configurations sampled without replacement prop. to w^beta
    (Gumbel top-K); returns S sequences (n, F2, s) and per-pair S-index arrays (n, F2)"""
    sel = []
    for r in range(3):
        idx, w = groups[k][r][parent_k[r]]
        if len(w) <= K: pick = np.arange(len(w))
        else:
            gk = beta * np.log(w) + rng.gumbel(size=len(w)); pick = np.argsort(-gk)[:K]
        sel.append(idx[pick])
    n = int(np.prod([len(x) for x in sel]))
    P = np.zeros((n, inst.F2), np.int64)
    grid = np.array(list(it.product(*[range(len(x)) for x in sel])), np.int64)
    for r in range(3):
        for t, g in enumerate(inst.rows[r]): P[:, g] = sel[r][grid[:, r], t]
    S = np.zeros((n, inst.F2, inst.s), np.int64)
    for g in range(inst.F2): S[:, g] = inst.pairs[k][g]['S'][P[:, g]]
    return S, P

# ------------------------------------------------------------------ D-check
def d_options(inst, k, pidx):
    """all D sequences for block k given S pair-indices pidx (F2,): (n, F2, s)"""
    opts = [inst.pairs[k][g]['D'][pidx[g]] for g in range(inst.F2)]
    n = int(np.prod([len(o) for o in opts]))
    grid = np.array(list(it.product(*[range(len(o)) for o in opts])), np.int64) if n <= 2_000_000 else None
    if grid is None: raise RuntimeError("too many D options")
    return np.stack([opts[g][grid[:, g]] for g in range(inst.F2)], 1)

def d_check(inst, pidx4, seed=0):
    """D-level for one S-solution: Z3 class split on D, then full D-level join.  Returns list of 4 D sequences."""
    T = inst.cy.T; Dl = [d_options(inst, k, pidx4[k]) for k in range(4)]
    # Z3 class split: row sums of D
    D3 = [np.stack([d[:, inst.rows[r]].sum(1) for r in range(3)], 1) for d in Dl]
    keys3 = []; uniq = []
    for k in range(4):
        u, inv = np.unique(D3[k].reshape(len(D3[k]), -1), axis=0, return_inverse=True)
        uniq.append((u.reshape(-1, 3, inst.s), inv.ravel()))
    sig3 = [seq_sig(u, T, 2) for (u, _) in uniq]
    h3, _ = join_all(sig3, inst.tsigD3(), seed=seed, tbits=18)
    out = []
    for q4 in h3:
        sub = [Dl[k][uniq[k][1] == q4[k]] for k in range(4)]
        sg = [seq_sig(sub[k], T, inst.nl) for k in range(4)]
        hf, _ = join_all(sg, inst.tsigD(), seed=seed + 1, tbits=16)
        for h in hf: out.append([sub[k][h[k]] for k in range(4)])
    return out, dict(d_options_log2=[round(math.log2(len(d)), 2) for d in Dl], d3_solutions=int(len(h3)))

# ------------------------------------------------------------------ assembly
def assemble(inst, S4, D4):
    """S, D sequences (4, F2, s) -> Z (4, F, s) -> (alpha, c) configs (4, F, 1 + s); raises if not a valid slice point"""
    F, F2, s = inst.F, inst.F2, inst.s
    Z = np.zeros((4, F, s), np.int64)
    for k in range(4):
        for g in range(F2):
            sg = 1 if g % 2 == 0 else -1
            Dp = sg * D4[k][g]; Z[k, g] = (S4[k][g] + Dp) // 2; Z[k, g + F2] = (S4[k][g] - Dp) // 2
    AC = np.zeros((4, F, 1 + s), np.int64)
    for k in range(4):
        for g in range(F):
            al = inst.cy.alphabet(int(inst.Y[k, g])); zc = zcoords(al)
            j = np.nonzero((zc == Z[k, g]).all(1))[0]
            assert len(j) == 1, "assembled coset not in the alphabet"
            AC[k, g] = al[j[0]]
    return Z, AC

def verify(inst, AC):
    X = expand_x(AC, inst.cy); Q = paf2(X)
    m = inst.F * inst.q
    ok = bool(Q[0, 0] == 4 * m and (Q.ravel()[1:] == -4).all())
    return ok, X

def write_core(inst, X, path, note=""):
    Xc = to_cyclic(X); m = inst.F * inst.q
    with open(path, 'w') as fh:
        fh.write(f"v={m + 1} N={4 * (m + 1)} b=1 m={m} cyclic Z{m} (CRT of Z{inst.F} x Z{inst.q}); generalised Paley slice s={inst.s}; theory/q523_sdsolve.py {note}\n")
        for nm, row in zip('ABCD', Xc): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')

# ------------------------------------------------------------------ driver
def run(inst, K=64, beta=1.5, top=24, seed=0, n_parents=1, max_s=None, log=print, out=None, max_passes=None):
    t0 = time.time(); rng = np.random.default_rng(seed)
    groups = [[inst.row_groups(k, r) for r in range(3)] for k in range(4)]
    parents = choose_parent(inst, groups, beta, top, seed, want=n_parents, log=log)
    found = []
    for pi, (par, sc) in enumerate(parents):
        lists = [fibre_list(inst, groups, k, par[k], K, beta, rng) for k in range(4)]
        sigs = [seq_sig(L[0], inst.cy.T, inst.nl) for L in lists]
        t1 = time.time()
        hits, st = join_all(sigs, inst.tsigS(), seed=seed + pi, tbits=24, max_passes=max_passes)
        rec = dict(stage="S-step", parent=pi, parent_scores=[round(x, 2) for x in sc], list_log2=[round(math.log2(len(L[0])), 2) for L in lists],
                   s_solutions=int(len(hits)), secs=round(time.time() - t1, 1), **st)
        log(rec)
        nd = 0; dopt = []
        for hq in (hits if max_s is None else hits[:max_s]):
            pidx4 = [lists[k][1][hq[k]] for k in range(4)]
            S4 = [lists[k][0][hq[k]] for k in range(4)]
            sols, info = d_check(inst, pidx4, seed=seed)
            nd += info['d3_solutions']; dopt.append(sum(info['d_options_log2']))
            for D4 in sols:
                Z, AC = assemble(inst, S4, D4)
                ok, X = verify(inst, AC)
                found.append(dict(ok=ok, AC=AC))
                log(dict(stage="CORE", ok=ok, secs=round(time.time() - t0, 1)))
                if ok and out: write_core(inst, X, out, f"seed={seed}")
                if ok: return found
        log(dict(stage="D-check", parent=pi, checked=int(len(hits) if max_s is None else min(len(hits), max_s)),
                 d3_solutions=nd, mean_hD=round(float(np.mean(dopt)), 2) if dopt else None, secs=round(time.time() - t0, 1)))
    return found

def planted_target(inst, rng):
    """random configuration given y; returns its P (F, s) and the configuration"""
    F = inst.F; AC = np.zeros((4, F, 1 + inst.s), np.int64)
    for k in range(4):
        for g in range(F):
            al = inst.cy.alphabet(int(inst.Y[k, g])); AC[k, g] = al[rng.integers(len(al))]
    P = pz(zcoords(AC), inst.cy.T).sum(0)
    return P, AC

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--K": "64", "--beta": "1.5", "--top": "24", "--seed": "0", "--out": "", "--parents": "1",
                             "--max-s": "", "--passes": ""}
    pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y)
    if '--planted' in flags:
        P, AC0 = planted_target(inst, np.random.default_rng(int(opt["--seed"]) + 999))
        inst = SDInst(F, q, s, Y, T=P)
        print(json.dumps(dict(planted=True)), flush=True)
    log = lambda d: print(json.dumps(d), flush=True)
    run(inst, int(opt["--K"]), float(opt["--beta"]), int(opt["--top"]), int(opt["--seed"]), int(opt["--parents"]),
        int(opt["--max-s"]) if opt["--max-s"] else None, log, opt["--out"] or None,
        int(opt["--passes"]) if opt["--passes"] else None)
