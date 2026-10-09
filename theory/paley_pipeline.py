r"""General exact two-stage pipeline for b >= 3 Paley-coset cores with Q = G'/H' an elementary abelian 2-group
(report/paley_design_answer.md sections 2, 7): fixed Q-compression (border-forced), sampled gadgets, exact stage 1,
exact stage 2, then completion to a full Hadamard matrix with a border found for the compression (c = 4: every
single-bump compression has one, scripts/bge3_border_check.py border_c4).

Stage 1: all alpha with the prescribed coset sums; alpha = +w on J positions, -w on J' positions; exact 4-list join
         of sum_k PAF(y_k) = T_J (linear 64-bit hash, hits re-checked; small lists here).
Stage 2: u on the twisted positions.  Exact: group each block's u by its coset sums u_Q; the Q-level part of S2
         (characters of Q are real) is sum_k [betahat_Q,k(xi)^2 + p uhat_Q,k(xi)^2] = 4v for every xi, a condition on
         the groups only; compatible group 4-tuples are then joined on the full lag-domain signature
         (sum PAF(beta) + p PAF(u), X) = (4v delta, 0).
usage: python3 theory/paley_pipeline.py run179 [max_gadgets] [seed]
       python3 theory/paley_pipeline.py check479          (rebuild + verify H(1916) from the 479 core)"""
import sys, os, json, time, itertools as it, math
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
import paley_ls as LSm
from paley_ansatz import LABELS

class Pipe:
    def __init__(s, v, b, p, mods, Hgen):
        s.v, s.b, s.p = v, b, p
        s.P = LSm.Prob(v, b, p, mods, Hgen); P = s.P; n = P.n
        Hs = sorted(P.H); rep = -np.ones(n, np.int64); reps = []; q = 0
        for g in range(n):
            if rep[g] < 0:
                for h in Hs: rep[P.add[g, h]] = q
                reps.append(g); q += 1
        s.rep, s.reps, s.c = rep, reps, q
        # characters of Q: Q must be elementary abelian 2-group; build from coset representatives
        s.Qadd = np.array([[rep[P.add[reps[a], reps[b_]]] for b_ in range(q)] for a in range(q)])
        assert all(s.Qadd[a, a] == 0 for a in range(q)), "Q must be elementary abelian 2-group"
        # characters: functions chi: Q -> +-1 with chi(a+b) = chi(a) chi(b): brute force over generators
        chars = []
        for signs in it.product((1, -1), repeat=q):
            ch = np.array(signs)
            if ch[0] == 1 and all(ch[s.Qadd[a, b_]] == ch[a] * ch[b_] for a in range(q) for b_ in range(q)): chars.append(ch)
        s.QCH = np.array(chars); assert len(chars) == q
        s.TJ = (P.T[:, 0] + (p - 1) * P.T[:, 1]).astype(np.int64)
        s.neg = np.array([int(P.G.idx(-P.G.co[g])) for g in range(n)])

    def paf_rows(s, Y):
        return np.stack([(Y * Y[:, s.P.add[:, h]]).sum(1) for h in range(s.P.n)], 1)

    def stage1(s, alphaQ, W, Jm):
        P = s.P; n = P.n; p = s.p; lists = []
        for k in range(4):
            choices = []
            for q in range(s.c):
                pos = list(np.nonzero(s.rep == q)[0])
                fixed = {g: (int(W[k, g]) if Jm[k, g] else -int(W[k, g])) for g in pos if W[k, g] != 0}
                free = [g for g in pos if g not in fixed]
                need = int(alphaQ[k][q]) - sum(fixed.values())
                opts = [dict(zip(free, vals)) | fixed for vals in it.product((1, -1), repeat=len(free)) if sum(vals) == need]
                if not opts: return []
                choices.append(opts)
            A = []
            for combo in it.product(*choices):
                a = np.zeros(n, np.int64)
                for d in combo:
                    for g, val in d.items(): a[g] = val
                A.append(a)
            A = np.array(A); lists.append((A, s.paf_rows(A + (p - 1) * W[k][None, :])))
        # exact 2+2 join
        T = s.TJ
        S01 = (lists[0][1][:, None, :] + lists[1][1][None, :, :]).reshape(-1, n)
        H = defaultdict(list)
        for i, r in enumerate(S01): H[r.tobytes()].append(i)
        out = []
        S23 = (T[None, None, :] - lists[2][1][:, None, :] - lists[3][1][None, :, :]).reshape(-1, n)
        for j, r in enumerate(S23):
            for i in H.get(r.tobytes(), []):
                a, b_ = divmod(i, len(lists[1][0])); c_, d = divmod(j, len(lists[3][0]))
                out.append(np.array([lists[0][0][a], lists[1][0][b_], lists[2][0][c_], lists[3][0][d]]))
        return out

    def stage2(s, A, W, max_sols=2):
        P = s.P; n = P.n; p = s.p; beta = A - W
        T = np.concatenate([np.eye(1, n, 0, dtype=np.int64)[0] * 4 * s.v, np.zeros(n, np.int64)])
        betaQ = np.array([[beta[k][s.rep == q].sum() for q in range(s.c)] for k in range(4)])
        tgt = 4 * s.v - ((s.QCH @ betaQ.T) ** 2).sum(1)
        if (tgt % p).any(): return [], 0
        tgt //= p                                                   # sum_k uhat_Q,k(xi)^2 = tgt(xi)
        blocks = []
        for k in range(4):
            tw = np.nonzero(W[k] == 0)[0]; m = len(tw)
            idx = np.arange(2 ** m, dtype=np.int64)
            Ub = np.zeros((2 ** m, n), np.int64); Ub[:, tw] = 1 - 2 * ((idx[:, None] >> np.arange(m)) & 1)
            Pu = s.paf_rows(Ub); Pb = s.paf_rows(beta[k][None, :])[0]
            Xk = np.stack([(beta[k][None, :] * Ub[:, P.add[:, h]]).sum(1) - (Ub * beta[k][P.add[:, h]][None, :]).sum(1) for h in range(n)], 1)
            sig = np.concatenate([Pb[None, :] + p * Pu, Xk], 1)
            uQ = np.stack([Ub[:, s.rep == q].sum(1) for q in range(s.c)], 1)
            spec = (uQ @ s.QCH.T) ** 2
            keys = (spec * (1 << (11 * np.arange(s.c, dtype=np.int64)))[None, :]).sum(1)
            gk, inv = np.unique(keys, return_inverse=True)
            groups = [np.nonzero(inv == g)[0] for g in range(len(gk))]
            blocks.append((Ub, sig, spec[[g[0] for g in groups]], groups))
        # Q-level compatible group 4-tuples
        def pairs(i, j):
            d = defaultdict(list)
            for a, sa in enumerate(blocks[i][2]):
                S = sa[None, :] + blocks[j][2]
                ok = np.nonzero((S <= tgt[None, :]).all(1))[0]
                for b_ in ok: d[S[b_].tobytes()].append((a, int(b_)))
            return d
        p01 = pairs(0, 1); p23 = pairs(2, 3)
        wts = np.random.default_rng(5).integers(1, 2 ** 62, 2 * n, dtype=np.int64).astype(np.uint64)
        key = lambda M: (M.astype(np.uint64) * wts[None, :]).sum(1)
        tkey = (T.astype(np.uint64) * wts).sum()
        sols = []; work = 0
        for kb, lst in p01.items():
            need = (tgt - np.frombuffer(kb, np.int64)).astype(np.int64).tobytes()
            if need not in p23: continue
            # union over (A,B) group pairs and (C,D) group pairs with these Q-level sums
            idA = []; idB = []
            K01 = []
            for (ga, gb) in lst:
                ia = blocks[0][3][ga]; ib = blocks[1][3][gb]
                kk = (key(blocks[0][1][ia])[:, None] + key(blocks[1][1][ib])[None, :]).reshape(-1)
                K01.append((kk, ia, ib))
            allk = np.concatenate([x[0] for x in K01]); order = np.argsort(allk); Ks = allk[order]
            offs = np.cumsum([0] + [len(x[0]) for x in K01])
            for (gc, gd) in p23[need]:
                ic = blocks[2][3][gc]; idd = blocks[3][3][gd]
                needk = (tkey - key(blocks[2][1][ic])[:, None] - key(blocks[3][1][idd])[None, :]).reshape(-1)
                lo = np.searchsorted(Ks, needk, 'left'); hi = np.searchsorted(Ks, needk, 'right')
                work += len(needk)
                for jj in np.nonzero(hi > lo)[0]:
                    c_, d = divmod(int(jj), len(idd))
                    for h in order[lo[jj]:hi[jj]]:
                        seg = int(np.searchsorted(offs, h, 'right') - 1); loc = int(h - offs[seg])
                        kk, ia, ib = K01[seg]; a, b_ = divmod(loc, len(ib))
                        rows = blocks[0][1][ia[a]] + blocks[1][1][ib[b_]] + blocks[2][1][ic[c_]] + blocks[3][1][idd[d]]
                        if not (rows == T).all(): continue
                        sols.append(np.array([blocks[0][0][ia[a]], blocks[1][0][ib[b_]], blocks[2][0][ic[c_]], blocks[3][0][idd[d]]]))
                        if len(sols) >= max_sols: return sols, work
        return sols, work

    def labels(s, A, U, W):
        lab = np.zeros(A.shape, np.int64)
        for k in range(4):
            for g in range(A.shape[1]):
                lab[k, g] = LABELS.index((int(A[k, g]), 0, int(W[k, g]))) if W[k, g] else LABELS.index((int(A[k, g]), 1, int(U[k, g])))
        return lab

def full_matrix_c4(F, mods, qidx, b):
    """c = 4 bordered matrix: core = 347-template GS layout (scripts/bge3_template347.py, trivial rho/gam, signs) on the
    group with radices `mods`, border (K, Kc) from scripts/bge3_border_check.border_c4 for the Q = Z2^2 compression
    (qidx(g) = Q-coset of g, ravel index of Z2^2).  Corner recomputed from orthogonality: C = -Kbig Core^T Kcbig/(16 h).
    Returns (H, info) with an exact H H^T = N I check."""
    cwd = os.getcwd(); os.chdir(ROOT)
    from bge3_border_check import border_c4, T as TPL
    os.chdir(cwd)
    LAYOUT = json.loads(str(TPL['layout']))
    F = np.asarray(F); m = F.shape[1]; c = 4; h = m // c; off = 4 * b; N = 4 * m + off
    co = np.array(np.unravel_index(np.arange(m), mods)).T
    idx = lambda X: np.ravel_multi_index(tuple(np.moveaxis(np.asarray(X) % np.array(mods), -1, 0)), mods)
    sig = {'h-g': idx(co[None, :, :] - co[:, None, :]), 'g+h': idx(co[:, None, :] + co[None, :, :]),
           'g-h': idx(co[:, None, :] - co[None, :, :]), '-g-h': idx(-co[:, None, :] - co[None, :, :])}
    Core = np.zeros((4 * m, 4 * m), np.int64)
    for key, L in LAYOUT.items():
        I_, J_ = int(key[0]), int(key[1]); sg = L['rho_sign'] * L['gam_sign']
        Core[I_ * m:(I_ + 1) * m, J_ * m:(J_ + 1) * m] = sg * F[L['k']][sig[L['sigma']]]
    xbar = np.array([[int(F[k][qidx == q].sum()) for q in range(4)] for k in range(4)])
    res = border_c4(xbar)
    if res is None: return None, dict(border=False, xbar=xbar.tolist())
    K, Kc, _ = res
    Kbig = np.concatenate([K[:, J * 4 + qidx] for J in range(4)], axis=1)
    Kcbig = np.concatenate([Kc[:, J * 4 + qidx] for J in range(4)], axis=1).T
    num = -(Kbig @ Core.T @ Kcbig)
    if (num % (16 * h)).any(): return None, dict(border=True, corner_integral=False)
    C = num // (16 * h)
    if not (np.abs(C) == 1).all(): return None, dict(border=True, corner_pm1=False)
    H = np.zeros((N, N), np.int64); H[:off, :off] = C; H[:off, off:] = Kbig; H[off:, :off] = Kcbig; H[off:, off:] = Core
    return H, dict(border=True, corner_pm1=True, N=N, hadamard=bool(np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64))))

def check479():
    """completion check: rebuild H(1916) from the published 479 core (b = 3, G' x Z_p = Z2 x Z2 x Z119, Q = Z2^2) with a
    c = 4 border found automatically by border_c4, and verify H H^T = N I exactly"""
    from hadlib import Core as HCore
    core = HCore(os.path.join(ROOT, 'data/core_v479_b3.txt'))
    mods = core.mods; m = core.m
    co = np.array(np.unravel_index(np.arange(m), mods)).T
    qidx = co[:, 0] * 2 + co[:, 1]                       # H = Z119 (last factor), Q = Z2 x Z2
    t0 = time.time()
    H, info = full_matrix_c4(core.X, mods, qidx, 3)
    info['secs'] = round(time.time() - t0, 1)
    print(json.dumps(info))
    return info

def run179(max_gadgets=200, seed=0):
    v, b, p = 179, 3, 11
    pipe = Pipe(v, b, p, (2, 8), [(0, 2)])
    P = pipe.P; n = P.n; rng = np.random.default_rng(seed)
    t = int(P.G.idx(np.array([0, 2])))
    # compression: c4 single-bump family; choose constants c_k and bumps d_k (three +-2 at one coset each, one block flat)
    comps = []
    for cs in it.product(range(-4, 5, 2), repeat=4):
        for flat in range(4):
            for sg in it.product((2, -2), repeat=3):
                ds = list(sg); ds.insert(flat, 0)
                for qs in it.product(range(4), repeat=3):
                    qq = list(qs); qq.insert(flat, 0)
                    xq = [[cs[k]] * 4 for k in range(4)]
                    ok = True
                    for k in range(4):
                        xq[k][qq[k]] += ds[k]
                        if any(abs(x) > 4 for x in xq[k]): ok = False
                    if not ok: continue
                    if sum(sum(x) ** 2 for x in xq) != 4 * v: continue
                    comps.append(xq)
    if not comps:                                     # report section 7.2: no balanced c = 4 compression exists for 179
        print(json.dumps(dict(v=v, b=b, p=p, balanced_single_bump_compressions=0, note='179 has no admissible compression in this family')))
        return []
    rng.shuffle(comps)
    pairs = []
    for c0 in range(n):
        for d in (t, 2 * t % n if False else int(P.add[t, t])):
            c1 = int(P.add[c0, d]); pr = tuple(sorted((c0, c1)))
            if pr not in pairs: pairs.append(pr)
    t0 = time.time(); n1 = 0; found = []
    log = open(os.path.join(ROOT, 'theory', 'logs', 'paley_pipeline_179.log'), 'w')
    for gi in range(max_gadgets):
        xq = comps[gi % len(comps)]
        blocks3 = rng.choice(4, 3, replace=False)
        W = np.zeros((4, n), np.int64); Jm = np.zeros((4, n), bool); gpos = []
        for k in blocks3:
            c0, c1 = pairs[rng.integers(len(pairs))]; s_ = int(rng.choice([1, -1]))
            W[k, c0] = s_; W[k, c1] = -s_; gpos += [(k, c0), (k, c1)]
        kj, gj = gpos[rng.integers(len(gpos))]; Jm[kj, gj] = True        # N_+ = 1, N_- = 5
        alphaQ = [np.array(x) for x in xq]
        # alphaQ is per coset in pipe.rep order: cosets indexed by rep (q = 0..3); xq built in that order
        S1 = pipe.stage1(alphaQ, W, Jm); n1 += len(S1)
        for A in S1[:4]:
            Us, wk = pipe.stage2(A, W)
            rec = dict(gadget=gi, stage1=len(S1), stage2_found=len(Us), work=wk, secs=round(time.time() - t0, 1))
            log.write(json.dumps(rec) + '\n'); log.flush()
            for U in Us:
                lab = pipe.labels(A, U, W)
                ok, X = LSm.verify(P, lab)
                F = np.array([x.reshape(-1) for x in X])                      # on Z2 x Z8 x Z11
                co = np.array(np.unravel_index(np.arange(F.shape[1]), (2, 8, 11))).T
                qidx = co[:, 0] * 2 + (co[:, 1] % 2)
                Hm, info = full_matrix_c4(F, (2, 8, 11), qidx, b)
                found.append(dict(gadget=gi, compression=[list(map(int, x)) for x in xq], core_ok=bool(ok), matrix=info, lab=lab.tolist()))
                log.write(json.dumps(dict(CORE=True, gadget=gi, core_ok=bool(ok), **info)) + '\n'); log.flush()
                if info.get('hadamard'):
                    with open(os.path.join(ROOT, 'data', f'core_v179_b3_paley_{len(found)}.txt'), 'w') as fh:
                        fh.write(f'v=179 N=716 b=3 m=176 group=2x8x11 (mixed radix, last index fastest) Paley-coset ansatz p=11, '
                                 f'found by theory/paley_pipeline.py (exact two-stage), border: bge3_border_check.border_c4\n')
                        for nm, row in zip('ABCD', F): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')
            if found: break
        if found: break
    summary = dict(gadgets_tried=gi + 1, stage1_solutions=n1, cores=len(found), secs=round(time.time() - t0, 1))
    log.write(json.dumps(summary) + '\n'); log.close()
    json.dump(dict(summary=summary, cores=found), open(os.path.join(ROOT, 'theory', 'logs', 'paley_pipeline_179.json'), 'w'))
    print(json.dumps(summary))

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "run179": run179(int(sys.argv[2]) if len(sys.argv) > 2 else 200, int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    elif cmd == "check479": check479()
