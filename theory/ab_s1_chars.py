r"""Exact stage 1 by the character decomposition over Q = Z_2^2 (the "coset-wise meet in the middle") [P, M].

Setting: a c = 4 (b = 3) staged slice with G' = Q x H' (Q = the first two axes, U acting on H' only), stage-1 target
T1 on G'.  For a character psi of Q put Z_k[psi] = sum_q psi(q) y_k(q, .) (a U-invariant vector on H').  [P]
    sum_k PAF_{G'}(y_k) = T1   <=>   for every psi:  sum_k PAF_{H'}(Z_k[psi]) = T^psi,  T^psi(c) = sum_q' psi(q') T1(q', c),
because PAF_{G'}(y)(q', c) Fourier-transformed over q' is PAF_{H'}(Z[psi])(c) (real characters, Q = Z_2^2).
Each block is a choice, per H'-orbit o, of a letter column (y(q, o))_q; its four character values are the Hadamard
transform of the column.  Characters psi_0 (trivial), psi_1, psi_2, psi_3.  Algorithm:
  1. per character: all block vectors Z[psi] (values = signed sums of 4 letters) with energy <= T^psi(0) and the
     exact weighted sum (psi_0: a block row sum of an admissible row quadruple, taken >= 0; psi != psi_0: 0 or +-2,
     since |zhat(psi)|^2 <= 4b = 12 and zhat is even); the 4-block system is solved on the distinct PAF signatures
     (multiset quads, energy-split hashed join: in energy order the two lighter signatures carry at most half the
     energy, so the light pairs are tabulated by hash and the other pairs streamed; every quad verified exactly);
  2. a block vector is kept only if its signature occurs in a solution of its character (all four characters);
  3. the blocks: per kept psi_0 vector a0, every (a1, a2, a3) of kept vectors whose (block, orbit) columns of
     character values are all transforms of letter columns (pairwise tables, then pairs, then triples);
  4. per psi_0 signature quad, a hashed 2+2 join of the four block lists on the psi_1..psi_3 PAFs against
     (T^psi1, T^psi2, T^psi3) (the psi_0 part holds by construction);
  5. every result is converted back to y and verified exactly (sum_k PAF_{G'} = T1).
Block signs are fixed by a non-negative row sum and the output is canonical up to block order.
Writes the solutions to OUT as raw int8 rows (4 x L1 orbit values each, canonical; numpy.fromfile(OUT, int8)
.reshape(-1, 4, L1)), every one verified; OUT.ckpt is the resume point (restart-safe).
usage: python3 theory/ab_s1_chars.py CFG.json OUT.jsonl [--plant SEED] [--maxB N] [--probe]
  --plant SEED: planted control (a random realistic quadruple with its own target; reports, per character, whether
                its vectors and signature quad were found, and whether the quadruple itself was recovered)
  --probe:      stop after the character systems; estimate the block count from 64 sampled psi_0 vectors
  --maxB N:     give up (undecided) when the block lists exceed N entries (default 4e7)"""
import sys, os, json, time, math, itertools
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD

CHARS = [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]      # psi(q), q index = 2 a + b


def structure(cfg, sl):
    Gp = sl.Gp; borb = np.asarray(sl.borb); co = Gp.co; n = Gp.n
    assert tuple(Gp.mods[:2]) == (2, 2)
    qix = co[:, 0] * 2 + co[:, 1]
    base = {}
    for x in range(n):
        if qix[x] == 0: base[tuple(co[x, 2:].tolist())] = int(borb[x])
    olab = sorted(set(base.values())); oid = {g: i for i, g in enumerate(olab)}; no = len(olab)
    hor = np.array([oid[base[tuple(co[x, 2:].tolist())]] for x in range(n)])
    j_of = {}
    for x in range(n): j_of[(int(qix[x]), int(hor[x]))] = int(borb[x])
    assert len(j_of) == len(sl.breps) == 4 * no
    H0 = np.nonzero(qix == 0)[0]
    wo = np.array([int((hor[H0] == o).sum()) for o in range(no)], np.int64)
    Up = [Gp.act(u) for u in sl.Sp]
    reps, _ = AGs.lag_reps(Gp, Up)
    crep = [int(c) for c in reps if qix[c] == 0]
    assert crep[0] == 0
    NH = np.zeros((len(crep), no, no), np.int64)
    for ci, c in enumerate(crep):
        xc = Gp.add(H0, np.full(len(H0), c)); np.add.at(NH[ci], (hor[H0], hor[xc]), 1)
    T1 = np.asarray(sl.T1, np.int64)
    Tpsi = []
    for ch in CHARS:
        row = []
        for c in crep:
            hc = tuple(co[c, 2:].tolist())
            row.append(int(sum(ch[int(qix[x])] * int(T1[x]) for x in range(n) if tuple(co[x, 2:].tolist()) == hc)))
        Tpsi.append(np.array(row, np.int64))
    lets = [[sorted(int(v) for v in sl.letters[j_of[(q, o)]]) for q in range(4)] for o in range(no)]
    return dict(no=no, wo=wo, NH=NH, crep=crep, Tpsi=Tpsi, lets=lets, j_of=j_of)


def columns(lets_o):
    """valid columns of one orbit: map transform tuple (z0, z1, z2, z3) -> column (y_q)_q"""
    out = {}
    for col in itertools.product(*lets_o):
        z = tuple(int(sum(ch[q] * col[q] for q in range(4))) for ch in CHARS)
        out[z] = col
    return out


def block_vectors(st, pi, wsums):
    no, wo = st['no'], st['wo']
    T0 = int(st['Tpsi'][pi][0])
    cur = np.zeros((1, 0), np.int64); E = np.zeros(1, np.int64)
    for o in range(no):
        V = np.array(sorted({z[pi] for z in st['cols'][o]}), np.int64)
        En = (E[:, None] + wo[o] * V[None, :] ** 2).ravel(); idx = np.nonzero(En <= T0)[0]
        cur = np.concatenate([cur[idx // len(V)], V[idx % len(V)][:, None]], 1); E = En[idx]
    ws = cur @ wo
    return cur[np.isin(ws, sorted(wsums))]


def pafs(st, Z):
    P = np.zeros((len(Z), len(st['crep'])), np.int64)
    for ci in range(len(st['crep'])): P[:, ci] = np.einsum('ni,ij,nj->n', Z, st['NH'][ci], Z)
    return P


def expand_ranges(lo, cnt):
    """(owner, value) for the concatenation of the integer ranges [lo_i, lo_i + cnt_i)"""
    cnt = np.asarray(cnt, np.int64); tot = int(cnt.sum())
    if tot == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64)
    rep = np.repeat(np.arange(len(cnt)), cnt)
    return rep, np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt) + np.repeat(np.asarray(lo, np.int64), cnt)


def sig_solve(P, T, chunk=1 << 22):
    """multiset quads of PAF signatures summing to T.  Returns (Pu, inv, Q): the distinct signatures sorted by energy
    (column 0), the signature index of every vector, and the quads i <= j <= k <= l (exactly verified).  In index
    order E_i + E_j <= T0 / 2, so the light pairs (i, j) are tabulated by hash and the pairs (k, l), k >= j, streamed."""
    P = np.asarray(P, np.int64); T = np.asarray(T, np.int64)
    if len(P) == 0: return P, np.zeros(0, np.int64), np.zeros((0, 4), np.int64)
    Pu, inv = np.unique(P, axis=0, return_inverse=True); inv = np.asarray(inv).ravel()
    o = np.argsort(Pu[:, 0], kind='stable'); r = np.empty_like(o); r[o] = np.arange(len(o)); Pu = Pu[o]; inv = r[inv]
    nu = len(Pu); E = Pu[:, 0]; T0 = int(T[0]); ar = np.arange(nu)
    Wv = np.random.default_rng(7).integers(1, 2 ** 61, Pu.shape[1]).astype(np.int64)
    h = Pu @ Wv; th = np.int64(T @ Wv)
    jmax = np.searchsorted(E, (T0 - 2 * E) // 2, 'right')
    LI, LJ = expand_ranges(ar, np.maximum(jmax - ar, 0))
    if len(LI) == 0: return Pu, inv, np.zeros((0, 4), np.int64)
    LI = LI.astype(np.int32); LJ = LJ.astype(np.int32)
    hL = h[LI] + h[LJ]; oL = np.argsort(hL, kind='stable'); hLs = hL[oL]; del hL
    EL = E[LI] + E[LJ]; eLmin, eLmax = int(EL.min()), int(EL.max()); del EL
    lo = np.maximum(ar, np.searchsorted(E, T0 - eLmax - E, 'left'))
    hi = np.searchsorted(E, T0 - eLmin - E, 'right'); cnt = np.maximum(hi - lo, 0)
    out = []; k0 = 0; cum = np.cumsum(cnt)
    while k0 < nu:
        base = int(cum[k0 - 1]) if k0 else 0
        k1 = max(k0 + 1, int(np.searchsorted(cum, base + chunk, 'right')))
        rep, LL = expand_ranges(lo[k0:k1], cnt[k0:k1]); KK = rep + k0
        if len(KK):
            need = th - (h[KK] + h[LL])
            a = np.searchsorted(hLs, need, 'left'); b = np.searchsorted(hLs, need, 'right'); m = b > a
            if m.any():
                KK, LL, a, b = KK[m], LL[m], a[m], b[m]
                rep2, p = expand_ranges(a, b - a); li = oL[p]
                I, J, K, L = LI[li].astype(np.int64), LJ[li].astype(np.int64), KK[rep2], LL[rep2]
                ok = J <= K; I, J, K, L = I[ok], J[ok], K[ok], L[ok]
                ok = ((Pu[I] + Pu[J] + Pu[K] + Pu[L]) == T).all(1)
                out.append(np.stack([I[ok], J[ok], K[ok], L[ok]], 1))
        k0 = k1
    Q = np.concatenate(out) if out else np.zeros((0, 4), np.int64)
    return Pu, inv, Q


def nvec_quads(Q, fib):
    """number of vector multisets over the signature quads (fib[u] = number of vectors with signature u)"""
    tot = 0
    for q in Q:
        c = 1
        for u, m in Counter(int(x) for x in q).items(): c *= math.comb(int(fib[u]) + m - 1, m)
        tot += c
    return tot


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); OUT = a[1]
    maxB = int(float(a[a.index('--maxB') + 1])) if '--maxB' in a else 40_000_000
    sl = W.build(cfg); t0 = time.time()
    plant = int(a[a.index('--plant') + 1]) if '--plant' in a else None
    if plant is not None:                       # planted control: a random slice quadruple with its own target
        rng = np.random.default_rng(plant)
        while True:                             # realistic energies: smallest letters, plus 1-2 larger ones per block
            Yo = []
            for _ in range(4):
                row = [int(rng.choice([v for v in sl.letters[j] if abs(v) == min(abs(x) for x in sl.letters[j])])) for j in range(len(sl.breps))]
                for j in rng.choice(len(sl.breps), size=int(rng.integers(1, 3)), replace=False):
                    big = [v for v in sl.letters[j] if abs(v) > min(abs(x) for x in sl.letters[j])]
                    row[j] = int(rng.choice(big))
                Yo.append(row)
            Yo = np.array(Yo)
            full = Yo[:, np.asarray(sl.borb)]
            rs = full.sum(1)
            if (rs != 0).all() and int(sl.Gp.paf(full).sum(0)[0]) <= 1.3 * int(sl.T1[0]): break
        Yo = np.array([r if rs[k] > 0 else -r for k, r in enumerate(Yo)]); full = Yo[:, np.asarray(sl.borb)]
        sl.T1 = sl.Gp.paf(full).sum(0)
        planted_key = tuple(sorted(tuple(int(x) for x in r) for r in Yo))
    st = structure(cfg, sl); no = st['no']
    st['cols'] = [columns(st['lets'][o]) for o in range(no)]
    if plant is None:
        rq, A = SD.row_quads_of(cfg, sl)
        rows = sorted({abs(r) for q in rq for r in q}); zh = {pi: [-2, 0, 2] for pi in (1, 2, 3)}
    else:
        rq = [tuple(sorted(int(x) for x in full.sum(1)))]; rows = sorted(set(rq[0]))
        qixp = sl.Gp.co[:, 0] * 2 + sl.Gp.co[:, 1]
        zh = {pi: sorted({s_ * int((np.array([CHARS[pi][int(t)] for t in qixp]) * full[k]).sum()) for k in range(4) for s_ in (1, -1)})
              for pi in (1, 2, 3)}
        Zpl = {pi: np.array([[sum(CHARS[pi][qq] * int(Yo[k][st['j_of'][(qq, o)]]) for qq in range(4)) for o in range(no)]
                             for k in range(4)], np.int64) for pi in range(4)}
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    log(stage='setup', n_o=no, weights=st['wo'].tolist(), Tpsi=[t.tolist() for t in st['Tpsi']], row_quads=rq)
    V, SG, plant_rep = {}, {}, {}
    for pi in (1, 2, 3, 0):                                      # nontrivial characters first: any empty system decides
        ws = rows if pi == 0 else zh[pi]                        # psi_0: non-negative row sums fix the block signs
        Z = block_vectors(st, pi, ws); P = pafs(st, Z)
        Pu, inv, Q = sig_solve(P, st['Tpsi'][pi])
        fib = np.bincount(inv, minlength=len(Pu))
        V[pi], SG[pi] = Z, (Pu, inv, Q)
        log(stage='character', psi=pi, vectors=len(Z), signatures=len(Pu), signature_quads=len(Q),
            solutions=nvec_quads(Q, fib) if len(Q) < 200_000 else 'not counted (> 2e5 signature quads)')
        if plant is not None:
            rowsZ = {tuple(r): i for i, r in enumerate(Z.tolist())}
            idx = [rowsZ.get(tuple(r)) for r in Zpl[pi].tolist()]
            sq = tuple(sorted(int(inv[i]) for i in idx)) if None not in idx else None
            plant_rep[pi] = dict(vectors_found=None not in idx,
                                 signature_quad_found=sq is not None and bool((np.sort(Q, 1) == np.array(sq)).all(1).any()))
            log(stage='planted-character', psi=pi, **plant_rep[pi])
        if not len(Q): break
    nfound, nbad, planted_found = 0, 0, False
    if all(pi in SG and len(SG[pi][2]) for pi in range(4)):
        for pi in range(4):                                      # keep only vectors whose signature occurs in a solution
            Pu, inv, Q = SG[pi]; keep = np.isin(inv, np.unique(Q))
            V[pi] = V[pi][keep]; SG[pi] = (Pu, inv[keep], Q)
        log(stage='kept', vectors=[len(V[pi]) for pi in range(4)])
        vidx = {pi: {int(x): k for k, x in enumerate(sorted({int(z[pi]) for o in range(no) for z in st['cols'][o]}))} for pi in range(4)}
        Ix = {pi: np.vectorize(lambda x, pi=pi: vidx[pi][int(x)])(V[pi]).astype(np.int64) for pi in range(4)}
        valid = []
        for o in range(no):
            Tb = np.zeros(tuple(len(vidx[pi]) for pi in range(4)), bool)
            for z in st['cols'][o]: Tb[vidx[0][z[0]], vidx[1][z[1]], vidx[2][z[2]], vidx[3][z[3]]] = True
            valid.append(Tb)
        T0i = {i: [valid[o].any(axis=tuple(x for x in (1, 2, 3) if x != i)) for o in range(no)] for i in (1, 2, 3)}
        V012 = [valid[o].any(axis=3) for o in range(no)]
        Pc = {pi: pafs(st, V[pi]) for pi in (1, 2, 3)}
        Ttar = np.concatenate([st['Tpsi'][pi] for pi in (1, 2, 3)])
        rngW = np.random.default_rng(11)
        Wc = {pi: rngW.integers(1, 2 ** 61, Pc[pi].shape[1]).astype(np.int64) for pi in (1, 2, 3)}
        Hc = {pi: Pc[pi] @ Wc[pi] for pi in (1, 2, 3)}
        th = np.array([st['Tpsi'][pi] @ Wc[pi] for pi in (1, 2, 3)], np.int64).sum()

        def blocks_of(a0, cap=1 << 25):
            z0 = Ix[0][a0]; c = {}
            for i in (1, 2, 3):
                m = np.ones(len(V[i]), bool)
                for o in range(no): m &= T0i[i][o][z0[o], Ix[i][:, o]]
                c[i] = np.nonzero(m)[0]
            if not (len(c[1]) and len(c[2]) and len(c[3])): return None
            ok12 = np.ones((len(c[1]), len(c[2])), bool)
            for o in range(no): ok12 &= V012[o][z0[o]][Ix[1][c[1], o][:, None], Ix[2][c[2], o][None, :]]
            p1, p2 = np.nonzero(ok12)
            if not len(p1): return None
            A1, A2, A3 = [], [], []
            step = max(1, cap // len(c[3]))
            for s0 in range(0, len(p1), step):
                q1, q2 = c[1][p1[s0:s0 + step]], c[2][p2[s0:s0 + step]]
                ok = np.ones((len(q1), len(c[3])), bool)
                for o in range(no): ok &= valid[o][z0[o]][Ix[1][q1, o][:, None], Ix[2][q2, o][:, None], Ix[3][c[3], o][None, :]]
                t, t3 = np.nonzero(ok)
                A1.append(q1[t]); A2.append(q2[t]); A3.append(c[3][t3])
            A1, A2, A3 = np.concatenate(A1), np.concatenate(A2), np.concatenate(A3)
            return (A1, A2, A3) if len(A1) else None

        Pu0, inv0, Q0 = SG[0]
        if '--probe' in a:
            smp = np.random.default_rng(3).choice(len(V[0]), size=min(64, len(V[0])), replace=False); tb = time.time()
            sz = []
            for x in smp:
                r = blocks_of(int(x)); sz.append(0 if r is None else len(r[0]))
            log(stage='probe', sampled=len(smp), mean_blocks=float(np.mean(sz)), max_blocks=int(max(sz)),
                est_total_blocks=float(np.mean(sz) * len(V[0])), est_secs=round((time.time() - tb) / len(smp) * len(V[0]), 1))
            return
        Bg = defaultdict(list); nB = 0
        for a0 in range(len(V[0])):
            if a0 % 2000 == 0: log(stage='blocks-progress', done=a0, of=len(V[0]), blocks=nB)
            r = blocks_of(a0)
            if r is None: continue
            Bg[int(inv0[a0])].append((np.full(len(r[0]), a0, np.int32), r[0].astype(np.int32), r[1].astype(np.int32), r[2].astype(np.int32)))
            nB += len(r[0])
            if nB > maxB:
                log(stage='done', undecided=f'block lists over the cap ({maxB})', blocks=nB, psi0_vectors_done=a0)
                return
        # every block's letters in orbit-representative order j (columns from per-orbit lookup tables)
        colLUT = []
        for o in range(no):
            Lt = np.zeros(tuple(len(vidx[pi]) for pi in range(4)) + (4,), np.int8)
            for z, col in st['cols'][o].items(): Lt[vidx[0][z[0]], vidx[1][z[1]], vidx[2][z[2]], vidx[3][z[3]]] = col
            colLUT.append(Lt)
        jpos = np.array([[st['j_of'][(qq, o)] for qq in range(4)] for o in range(no)])
        L1 = 4 * no; RB = 4 * L1; wj = np.zeros(L1, np.int64)
        for o in range(no): wj[jpos[o]] = st['wo'][o]
        B = {}
        for u, lst in Bg.items():
            b = [np.concatenate([x[t] for x in lst]) for t in range(4)]
            R = np.zeros((len(b[0]), L1), np.int8)
            for o in range(no): R[:, jpos[o]] = colLUT[o][Ix[0][b[0], o], Ix[1][b[1], o], Ix[2][b[2], o], Ix[3][b[3], o]]
            B[u] = (b, Hc[1][b[1]] + Hc[2][b[2]] + Hc[3][b[3]], R)
        del Bg
        log(stage='blocks', psi0_vectors=len(V[0]), signatures_with_blocks=len(B), column_valid_blocks=nB)

        def canon(Y):
            """block order, and the sign of zero-row-sum blocks (letters are odd, so y_0 != 0), normalized"""
            rs = (Y.astype(np.int64) * wj).sum(2)
            Y = np.where(((rs == 0) & (Y[:, :, 0] < 0))[:, :, None], -Y, Y)
            Yv = np.ascontiguousarray(Y).view(np.dtype((np.void, Y.shape[2]))).reshape(len(Y), 4)
            return np.take_along_axis(Y, np.argsort(Yv, axis=1)[:, :, None], axis=1)

        def verify(Y, chunk=20000):
            ok = np.ones(len(Y), bool); borb = np.asarray(sl.borb); T1 = np.asarray(sl.T1, np.int64)
            for s0 in range(0, len(Y), chunk):
                full = Y[s0:s0 + chunk][:, :, borb].astype(np.int64)
                ok[s0:s0 + chunk] = (sl.Gp.paf(full).sum(1) == T1[None, :]).all(1)
            return ok

        # OUT holds the solutions as raw int8 rows (4 x L1 letters each, canonical); OUT.ckpt the resume point
        ck = OUT + '.ckpt'; start = 0
        if os.path.exists(ck) and os.path.exists(OUT):
            c = json.load(open(ck)); start, nfound, nbad = c['next_quad'], c['rows'], c['unverified']
            with open(OUT, 'r+b') as f: f.truncate(nfound * RB)
        elif os.path.exists(OUT): os.remove(OUT)
        cache = {}; csize = 0; pairs_done = 0; t_last = 0.0; over = 0
        splits = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
        fo = open(OUT, 'ab')
        for si in range(start, len(Q0)):
            if time.time() - t_last > 60:
                log(stage='progress', psi0_quads_done=si, of=len(Q0), stage1_solutions=nfound, pairs=pairs_done); t_last = time.time()
            u = [int(x) for x in Q0[si]]
            if all(x in B for x in u):
                sz = [len(B[x][1]) for x in u]
                (l0, l1), (r0, r1) = min(splits, key=lambda s: sz[s[0][0]] * sz[s[0][1]] + sz[s[1][0]] * sz[s[1][1]])
                if sz[l0] * sz[l1] < sz[r0] * sz[r1]: (l0, l1), (r0, r1) = (r0, r1), (l0, l1)
                if sz[r0] * sz[r1] > 8e7: over += 1; continue
                key = (u[l0], u[l1])
                if key not in cache:
                    hL = (B[u[l0]][1][:, None] + B[u[l1]][1][None, :]).ravel(); oo = np.argsort(hL)
                    if csize + len(hL) > 3e7: cache.clear(); csize = 0
                    cache[key] = (hL[oo], oo); csize += len(hL); pairs_done += len(hL)
                hLs, oL = cache[key]; n1 = sz[l1]
                hR = (th - (B[u[r0]][1][:, None] + B[u[r1]][1][None, :])).ravel(); pairs_done += len(hR)
                lo = np.searchsorted(hLs, hR, 'left'); hi = np.searchsorted(hLs, hR, 'right')
                m = np.nonzero(hi > lo)[0]
                if len(m):
                    rep, pos = expand_ranges(lo[m], hi[m] - lo[m])
                    i2, i3 = np.divmod(m[rep], sz[r1]); i0, i1 = np.divmod(oL[pos], n1)
                    picks = [(u[l0], i0), (u[l1], i1), (u[r0], i2), (u[r1], i3)]
                    tot = sum(np.concatenate([Pc[pi][B[g][0][pi][ii]] for pi in (1, 2, 3)], 1) for g, ii in picks)
                    ok = (tot == Ttar).all(1)
                    if ok.any():
                        Y = np.stack([B[g][2][ii[ok]] for g, ii in picks], 1)
                        Y = np.unique(canon(Y).reshape(-1, RB), axis=0).reshape(-1, 4, L1)
                        nbad += int((~verify(Y)).sum())
                        fo.write(np.ascontiguousarray(Y, np.int8).tobytes()); nfound += len(Y)
            if si % 200 == 199 or si == len(Q0) - 1:
                fo.flush(); json.dump(dict(next_quad=si + 1, rows=nfound, unverified=nbad), open(ck, 'w'))
        fo.close()
        if over: log(stage='done', undecided=f'{over} psi_0 quads over the pair cap', stage1_solutions_so_far=nfound); return
        if plant is not None and nfound:
            Yp = canon(np.array(Yo, np.int8)[None]).reshape(RB)
            Yall = np.fromfile(OUT, np.int8).reshape(-1, RB)
            planted_found = bool((Yall == Yp[None, :]).all(1).any())
    extra = {}
    if plant is not None: extra = dict(planted_found=planted_found, planted_characters=plant_rep)
    log(stage='done', stage1_solutions=nfound, all_verified=nbad == 0, **extra)


if __name__ == "__main__":
    main()
