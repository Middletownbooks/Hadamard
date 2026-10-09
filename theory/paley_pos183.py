r"""Positive control for the exact two-stage design: Paley-coset cores at v = 183 (c = 8, b = 7, p = 11, G' = Z2 x Z2 x Z4,
H' = <(0,0,2)>) whose Q-compression is accepted by the PUBLISHED 359 border, completed to H(732) and verified.

Design (report/paley_design_answer.md section 6):
  * compression (x-coordinates, cosets q = 4a + 2b + (e mod 2), q* = 7): A = (-2)^7 0, B = 2^7 0, C = (-2)^7 0,
    D = (-2)^7 2  -- the only single-deviation compressions the 359 border accepts (theory/paley_c8family.py test);
  * gadget: norm identity 3K - N_- = n' - b = 9 forces K = 4, N_- = 3 (K = 3 would need an unbalanced coset);
    two balanced pairs {c, c + t} (t = (0,0,2)) in two different blocks (a same-block pair of pairs leaves +-200 off H');
  * stage 1: all alpha with the compression, exact check of sum_k PAF(y_k) = T_J;
  * stage 2: u on the 60 twisted positions, exact 4-list join grouped by the Q-compression of u (meet-in-the-middle with
    a linear 64-bit hash, hits re-checked), conditions S2: sum PAF(beta) + p PAF(u) = 4v delta and X = 0;
  * every core found is expanded (FFT check of sum PAF) and spliced into the 359 border (H H^T = 732 I checked exactly).
usage: python3 theory/paley_pos183.py [max_cores]"""
import sys, os, json, time, itertools as it, math
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
import paley_ls as LSm
from paley_ansatz import LABELS
from paley_c8family import C8Family

V, B, P_, MODS, HGEN = 183, 7, 11, (2, 2, 4), [(0, 0, 2)]
XQ = [[-2] * 7 + [0], [2] * 7 + [0], [-2] * 7 + [0], [-2] * 7 + [2]]

def paf_vec(P, f):
    return np.array([(f * f[P.add[:, h]]).sum() for h in range(P.n)])

def cc_vec(P, f, g):
    return np.array([(f * g[P.add[:, h]]).sum() for h in range(P.n)])

def setup():
    P = LSm.Prob(V, B, P_, MODS, HGEN)
    co = P.G.co; cos = co[:, 0] * 4 + co[:, 1] * 2 + (co[:, 2] % 2)
    t = int(P.G.idx(np.array([0, 0, 2])))
    return P, cos, t

def gadgets(P, cos, t):
    """two balanced pairs in two different blocks, N_- = 3 (one J position), w = (s, -s) on {c, c + t}"""
    out = []
    pairs = sorted({tuple(sorted((g, int(P.add[g, t])))) for g in range(P.n)})
    for kb in it.combinations(range(4), 2):
        for pr in it.product(pairs, repeat=2):
            for sg in it.product((1, -1), repeat=2):
                for jpos in range(4):                                     # which of the 4 gadget positions is J
                    W = np.zeros((4, P.n), np.int64); J = np.zeros((4, P.n), bool)
                    pos = []
                    for (k, (c0, c1), s) in zip(kb, pr, sg):
                        W[k, c0] = s; W[k, c1] = -s; pos += [(k, c0), (k, c1)]
                    k_, g_ = pos[jpos]; J[k_, g_] = True
                    out.append((W, J))
    return out

def stage1(P, cos, W, J):
    """all alpha (4 x n) with coset sums XQ, alpha = +w on J, -w on the other gadget positions; exact S1"""
    p = P.p; n = P.n
    TJ = (P.T[:, 0] + (p - 1) * P.T[:, 1]).astype(np.int64)
    per_block = []
    for k in range(4):
        choices = []
        for q in range(8):
            pos = list(np.nonzero(cos == q)[0]); fixed = {g: (W[k, g] if J[k, g] else -W[k, g]) for g in pos if W[k, g] != 0}
            free = [g for g in pos if g not in fixed]
            need = XQ[k][q] - sum(fixed.values())
            opts = []
            for vals in it.product((1, -1), repeat=len(free)):
                if sum(vals) == need: opts.append(dict(zip(free, vals)) | fixed)
            if not opts: return []
            choices.append(opts)
        alphas = []
        for combo in it.product(*choices):
            a = np.zeros(n, np.int64)
            for d in combo:
                for g, val in d.items(): a[g] = val
            alphas.append(a)
        per_block.append(alphas)
    sols = []
    for combo in it.product(*per_block):
        Y = [combo[k] + (p - 1) * W[k] for k in range(4)]
        if (sum(paf_vec(P, y) for y in Y) == TJ).all(): sols.append(np.array(combo))
    return sols

QCO = np.array([[q // 4, (q // 2) % 2, q % 2] for q in range(8)])
QCH = np.array([[(-1) ** int((xi * qc).sum()) for qc in QCO] for xi in QCO])          # characters of Q = Z2^3

def q_level_solutions(P, A, W, cos):
    """all (u_Q,A..u_Q,D): u_Q,k(q) in {-2,0,2} (two twisted positions) or {0} (gadget coset), sum_q u_Q,k = 0,
    sum_k uhat_Q,k(xi)^2 = (4v - sum_k betahat_Q,k(xi)^2) / p for every xi (exact Q-level part of S2)"""
    p = P.p; beta = A - W
    betaQ = np.array([[beta[k][cos == q].sum() for q in range(8)] for k in range(4)])
    tgt = 4 * V - ((QCH @ betaQ.T) ** 2).sum(1)
    if (tgt % p).any(): return []
    tgt //= p
    per = []
    for k in range(4):
        ntw = [int(((cos == q) & (W[k] == 0)).sum()) for q in range(8)]
        vecs = np.array(list(it.product(*[list(range(-m, m + 1, 2)) for m in ntw])), np.int64)
        spec = (vecs @ QCH.T) ** 2
        keep = spec[:, 0] == tgt[0] if False else (vecs.sum(1) ** 2 * 0 == 0)
        per.append((vecs, spec))
    # join on the full spectrum vector (8 values incl. xi = 1): pairs (0,1) and (2,3)
    def table(i, j):
        T = defaultdict(list)
        Si, Sj = per[i][1], per[j][1]
        for a in range(len(Si)):
            S = Si[a][None, :] + Sj
            ok = (S <= tgt[None, :]).all(1)
            for b_ in np.nonzero(ok)[0]: T[S[b_].tobytes()].append((a, int(b_)))
        return T
    T01 = table(0, 1); T23 = table(2, 3); idx = []
    for kb, lst in T01.items():
        need = (tgt - np.frombuffer(kb, np.int64)).astype(np.int64).tobytes()
        if need not in T23: continue
        L1 = np.array(lst, np.int64); L2 = np.array(T23[need], np.int64)
        I = np.repeat(np.arange(len(L1)), len(L2)); J = np.tile(np.arange(len(L2)), len(L1))
        idx.append(np.stack([L1[I, 0], L1[I, 1], L2[J, 0], L2[J, 1]], 1))
    if not idx: return np.zeros((0, 4, 8), np.int64)
    idx = np.concatenate(idx)
    return np.stack([per[k][0][idx[:, k]] for k in range(4)], 1)          # (N, 4, 8)

def stage2(P, A, W, cos, max_sols=4, QS=None):
    """exact.  For |H'| = 2 the S2/X conditions split by the restriction lambda = psi(t) = +-1:
         lambda = +1 (Q-level): only the coset sums u_Q;   lambda = -1: only d(q) = u(s_q) - u(s_q + t),
       and d(q) = +-2 exactly on the twisted cosets with u_Q(q) = 0 (one free sign each), else d = 0.
       So: Q-level solutions -> group by zero pattern -> per pattern an exact 2+2 join on the lambda = -1 conditions,
       written with integers: |bhat + i sqrt(p) uhat|^2 = A + sqrt(p) B, A = |bhat|^2 + p |uhat|^2, B = 2 Im(conj(bhat) uhat)
       (bhat, uhat in Z[i] since the lambda = -1 characters take values in {+-1, +-i})."""
    p = P.p; n = P.n; beta = A - W
    co = P.G.co; t = int(P.G.idx(np.array([0, 0, 2])))
    reps = [int(np.nonzero(cos == q)[0][0]) for q in range(8)]
    for q in range(8): assert int(P.add[reps[q], t]) in set(np.nonzero(cos == q)[0].tolist())
    # lambda = -1 characters of Z2 x Z2 x Z4: psi(a,b,e) = (-1)^(c1 a + c2 b) i^(c3 e), c3 odd
    chars = [(c1, c2, c3) for c1 in (0, 1) for c2 in (0, 1) for c3 in (1, 3)]
    def chi_val(c, g):
        a, b, e = co[g]; return ((-1) ** (c[0] * a + c[1] * b)) * (1j ** ((c[2] * e) % 4))
    Bh = np.array([[sum(beta[k][g] * chi_val(c, g) for g in range(n)) for c in chars] for k in range(4)])   # (4, 8) Z[i]
    psi_s = np.array([[chi_val(c, reps[q]) for c in chars] for q in range(8)])                           # (8 cosets, 8 chars)
    if QS is None: QS = q_level_solutions(P, A, W, cos)
    groups = {}
    if len(QS):
        tw = np.array([[W[k, reps[q]] == 0 for q in range(8)] for k in range(4)])          # (4, 8) coset has twisted pos
        zm = ((QS == 0) & tw[None, :, :])                                                  # (N, 4, 8)
        keys = (zm.reshape(len(QS), 32) * (1 << np.arange(32, dtype=np.int64))[None, :]).sum(1)
        uk, first = np.unique(keys, return_index=True)
        for kv, i0 in zip(uk.tolist(), first.tolist()):
            Z = tuple(tuple(int(q) for q in range(8) if (kv >> (8 * k + q)) & 1) for k in range(4))
            groups[Z] = [QS[i0]]
    sols = []; work = 0
    def block_sigs(k, Zk):
        m = len(Zk); S = 1 - 2 * ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)                 # (2^m, m)
        uh = 2 * (S @ psi_s[list(Zk), :]) if m else np.zeros((1, 8), complex)                          # (2^m, 8)
        b = Bh[k][None, :]
        Aint = np.rint(np.abs(b) ** 2 + p * np.abs(uh) ** 2).astype(np.int64)
        Bint = np.rint(2 * (np.conj(b) * uh).imag).astype(np.int64)
        return S, np.concatenate([Aint, Bint], 1)
    Tl = np.concatenate([np.full(8, 4 * V, np.int64), np.zeros(8, np.int64)])
    for Z, members in groups.items():
        sg = [block_sigs(k, Z[k]) for k in range(4)]
        S01 = (sg[0][1][:, None, :] + sg[1][1][None, :, :]).reshape(-1, 16)
        H = defaultdict(list)
        for i, r in enumerate(S01): H[r.tobytes()].append(i)
        S23 = (Tl[None, None, :] - sg[2][1][:, None, :] - sg[3][1][None, :, :]).reshape(-1, 16)
        work += len(S01) + len(S23)
        for j, r in enumerate(S23):
            for i in H.get(r.tobytes(), []):
                a, b_ = divmod(i, len(sg[1][0])); c, d = divmod(j, len(sg[3][0]))
                sig_choice = [sg[0][0][a], sg[1][0][b_], sg[2][0][c], sg[3][0][d]]
                uq = members[0]
                U = np.zeros((4, n), np.int64)
                for k in range(4):
                    for q in range(8):
                        pos = [g for g in np.nonzero(cos == q)[0] if W[k, g] == 0]
                        if not pos: continue
                        s0, s1 = reps[q], int(P.add[reps[q], t])
                        if uq[k][q] != 0: U[k, s0] = U[k, s1] = uq[k][q] // 2
                        else:
                            sgn = sig_choice[k][list(Z[k]).index(q)]; U[k, s0] = sgn; U[k, s1] = -sgn
                sols.append(U)
                if len(sols) >= max_sols: return sols, work, len(QS), len(groups)
    return sols, work, len(QS), len(groups)

def labels_from(A, U, W):
    lab = np.zeros(A.shape, np.int64)
    for k in range(4):
        for g in range(A.shape[1]):
            lab[k, g] = LABELS.index((int(A[k, g]), 0, int(W[k, g]))) if W[k, g] else LABELS.index((int(A[k, g]), 1, int(U[k, g])))
    return lab

def main(max_cores=3):
    P, cos, t = setup(); t0 = time.time()
    S = C8Family(); mods, n_full, co_full, sig, chars, cos_full = S.group(2)
    lam = np.where((co_full[:, 0] + co_full[:, 1]) % 2 == 0, 1, -1)
    gs = gadgets(P, cos, t)
    n_s1 = 0; found = []; s2work = 0; nQ_tot = [0]; Qcache = {}; seen_orbits = set()
    for gi, (W, J) in enumerate(gs):
        for A in stage1(P, cos, W, J):
            # per-block symmetries preserving compression, PAFs and X: translation by t in H', reflection g -> -g
            # (with u -> -u, the NQR multiplier on Z_p, so X_k is unchanged); process one stage-1 solution per orbit
            key = []
            for k in range(4):
                imgs = []
                for tr in (0, t):
                    for rf in (False, True):
                        perm = np.array([int(P.add[int(P.G.idx(-P.G.co[g])) if rf else g, tr]) for g in range(P.n)])
                        a2 = np.zeros(P.n, np.int64); w2 = np.zeros(P.n, np.int64)
                        a2[perm] = A[k]; w2[perm] = W[k]
                        imgs.append(a2.tobytes() + w2.tobytes())
                key.append(min(imgs))
            key = tuple(key)
            if key in seen_orbits: continue
            seen_orbits.add(key)
            n_s1 += 1
            hk = tuple(tuple(np.nonzero(W[k] != 0)[0].tolist()) for k in range(4))
            if hk not in Qcache: Qcache[hk] = q_level_solutions(P, A, W, cos)
            Us, wk, nq, ng = stage2(P, A, W, cos, QS=Qcache[hk]); s2work += wk; nQ_tot[0] += nq
            print(json.dumps(dict(gadget=gi, stage1=n_s1, Q_solutions=nq, zero_patterns=ng, found=len(Us), secs=round(time.time() - t0, 1))), flush=True)
            for U in Us:
                lab = labels_from(A, U, W)
                ok, X = LSm.verify(P, lab)
                Xf = np.array([x.reshape(-1) for x in X])                  # x-coordinates on Z2 x Z2 x Z4 x Z11
                H, info = S.build(lam[None, :] * Xf, 2)
                rec = dict(core_ok=bool(ok), **info)
                found.append(dict(gadget=gi, rec=rec, lab=lab.tolist()))
                print(json.dumps(dict(gadget=gi, stage1_solutions_so_far=n_s1, **rec, secs=round(time.time() - t0, 1))), flush=True)
                if rec.get('hadamard'):
                    np.save(os.path.join(ROOT, 'theory', 'logs', f'H732_paley_pos183_{len(found)}.npy'), H.astype(np.int8))
                if len(found) >= max_cores: break
            if len(found) >= max_cores: break
        if len(found) >= max_cores: break
    summary = dict(gadgets=len(gs), gadgets_scanned=gi + 1, stage1_solutions=n_s1, Q_level_solutions=nQ_tot[0], cores=len(found),
                   hadamard=sum(1 for f in found if f['rec'].get('hadamard')), stage2_hash_work=s2work, secs=round(time.time() - t0, 1))
    print(json.dumps(summary), flush=True)
    json.dump(dict(summary=summary, cores=found), open(os.path.join(ROOT, 'theory', 'logs', 'paley_pos183.json'), 'w'))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
