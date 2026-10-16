r"""Exact stage 1 by the character decomposition over Q = Z_2^2 (the "coset-wise meet in the middle") [P, M].

Setting: a c = 4 (b = 3) staged slice with G' = Q x H' (Q = the first two axes, U acting on H' only), stage-1 target
T1 on G'.  For a character psi of Q put Z_k[psi] = sum_q psi(q) y_k(q, .) (a U-invariant vector on H').  [P]
    sum_k PAF_{G'}(y_k) = T1   <=>   for every psi:  sum_k PAF_{H'}(Z_k[psi]) = T^psi,  T^psi(c) = sum_q' psi(q') T1(q', c),
because PAF_{G'}(y)(q', c) Fourier-transformed over q' is PAF_{H'}(Z[psi])(c) (real characters, Q = Z_2^2).
Each block is a choice, per H'-orbit o, of a letter column (y(q, o))_q; its four character values are the Hadamard
transform of the column.  Algorithm:
  1. per character: all block vectors Z[psi] (values = signed sums of 4 letters) with energy <= T^psi(0) and the
     exact weighted sum (psi = 1: +-a block row sum of an admissible row quadruple; psi != 1: 0 or +-2, since
     |zhat(psi)|^2 <= 4b = 12 and zhat is even); then the 4-block system by a hashed 2+2 join (unordered quads);
  2. the cross-character join: a stage-1 solution is a choice of one solution per character and an assignment of
     their vectors to the 4 blocks such that every (block, orbit) column of character values is the transform of a
     letter column.  Joined psi_1 -> psi_2 -> psi_3 on top of each psi_0 solution, with exact column tests;
  3. every result is converted back to y and verified exactly (sum_k PAF_{G'} = T1).
Block signs are fixed by a positive row sum (psi_0 weighted sum > 0; row sums of c = 4 slices are non-zero) and the
output is canonical up to block order.  Writes JSON lines {"y": 4 x L1 orbit values, "verified": true} to OUT.
usage: python3 theory/ab_s1_chars.py CFG.json OUT.jsonl"""
import sys, os, json, time, itertools
from collections import defaultdict
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
    vals = sorted({z[pi] for o in range(no) for z in st['cols'][o]})
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


def solve_sys(P, T, chunk=1 << 24):
    """unordered quads a <= b <= c <= d with P[a] + P[b] + P[c] + P[d] == T (hashed 2+2, verified)"""
    nv = len(P); Wv = np.random.default_rng(7).integers(1, 2 ** 61, P.shape[1]).astype(np.int64)
    h = P @ Wv; th = np.int64(int(np.asarray(T, np.int64) @ Wv)); E = P[:, 0]; T0 = int(T[0])
    I, J = np.triu_indices(nv); keep = E[I] + E[J] <= T0; I, J = I[keep], J[keep]
    HP = h[I] + h[J]; o = np.argsort(HP); HPs = HP[o]
    sols = set()
    for s0 in range(0, len(HP), chunk):
        need = th - HP[s0:s0 + chunk]
        lo = np.searchsorted(HPs, need, 'left'); hi = np.searchsorted(HPs, need, 'right')
        for p in np.nonzero(hi > lo)[0]:
            a, b = int(I[s0 + p]), int(J[s0 + p])
            for r in o[lo[p]:hi[p]]:
                c, d = int(I[r]), int(J[r]); q = tuple(sorted((a, b, c, d)))
                if q not in sols and (P[a] + P[b] + P[c] + P[d] == T).all(): sols.add(q)
    return sorted(sols)


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); OUT = a[1]
    sl = W.build(cfg); t0 = time.time()
    st = structure(cfg, sl); no = st['no']
    st['cols'] = [columns(st['lets'][o]) for o in range(no)]
    rq, A = SD.row_quads_of(cfg, sl)
    rows = sorted({abs(r) for q in rq for r in q})
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    log(stage='setup', n_o=no, weights=st['wo'].tolist(), Tpsi=[t.tolist() for t in st['Tpsi']], row_quads=rq)
    V, S = {}, {}
    for pi in range(4):
        ws = rows if pi == 0 else [-2, 0, 2]                    # psi_0: positive row sums fix the block signs
        Z = block_vectors(st, pi, ws); P = pafs(st, Z)
        V[pi], S[pi] = Z, solve_sys(P, st['Tpsi'][pi])
        log(stage='character', psi=pi, vectors=len(Z), solutions=len(S[pi]))
        if not S[pi]: break
    res = []
    if all(S.get(pi) for pi in range(4)):
        # value indices and pair-compatibility matrices C[i][a, b]: psi_0 vector a and psi_i vector b column-compatible
        vidx = {pi: {int(x): k for k, x in enumerate(sorted({int(z[pi]) for o in range(no) for z in st['cols'][o]}))} for pi in range(4)}
        def idxmat(pi): return np.vectorize(lambda x: vidx[pi][int(x)])(V[pi]) if len(V[pi]) else np.zeros((0, no), np.int64)
        I0 = idxmat(0); C = {}
        for i in (1, 2, 3):
            Ii = idxmat(i); M = np.ones((len(V[0]), len(V[i])), bool)
            for o in range(no):
                Tb = np.zeros((len(vidx[0]), len(vidx[i])), bool)
                for z in st['cols'][o]: Tb[vidx[0][z[0]], vidx[i][z[i]]] = True
                M &= Tb[I0[:, o][:, None], Ii[:, o][None, :]]
            C[i] = M
        log(stage='compat', density={i: round(float(C[i].mean()), 5) for i in C})
        T012 = [{(z[0], z[1], z[2]) for z in st['cols'][o]} for o in range(no)]
        def index(sols):
            ix = defaultdict(list)
            for si, q in enumerate(sols):
                for v_ in set(q): ix[v_].append(si)
            return ix
        ix = {pi: index(S[pi]) for pi in (1, 2, 3)}
        def assignments(p0, pi):
            out = set(); Ci = C[pi]
            for b in np.nonzero(Ci[p0[0]])[0].tolist():
                for si in ix[pi].get(b, []):
                    for perm in set(itertools.permutations(S[pi][si])):
                        if perm[0] == b and all(Ci[p0[k], perm[k]] for k in range(1, 4)): out.add(perm)
            return sorted(out)
        nfound = 0; seen = set()
        if os.path.exists(OUT):                                  # restart-safe: keep what is already written
            for l in open(OUT):
                try: seen.add(tuple(tuple(r) for r in json.loads(l)['y'])); res.append(json.loads(l))
                except Exception: pass
            nfound = len(seen)
        for si0, q0 in enumerate(S[0]):
            if si0 % 20000 == 0: log(stage='progress', psi0_done=si0, of=len(S[0]), stage1_solutions=nfound)
            p0 = q0
            A1 = assignments(p0, 1)
            if not A1: continue
            A2 = assignments(p0, 2)
            if not A2: continue
            A3 = assignments(p0, 3)
            if not A3: continue
            z0s = [V[0][i] for i in p0]
            for p1 in A1:
                z1s = [V[1][i] for i in p1]
                for p2 in A2:
                    z2s = [V[2][i] for i in p2]
                    if not all((int(z0s[k][o]), int(z1s[k][o]), int(z2s[k][o])) in T012[o] for k in range(4) for o in range(no)):
                        continue
                    for p3 in A3:
                        z3s = [V[3][i] for i in p3]
                        Y = []; ok = True
                        for k in range(4):
                            row = [None] * (4 * no)
                            for o in range(no):
                                col = st['cols'][o].get((int(z0s[k][o]), int(z1s[k][o]), int(z2s[k][o]), int(z3s[k][o])))
                                if col is None: ok = False; break
                                for q in range(4): row[st['j_of'][(q, o)]] = int(col[q])
                            if not ok: break
                            Y.append(row)
                        if not ok: continue
                        key = tuple(sorted(tuple(r) for r in Y))
                        if key in seen: continue
                        seen.add(key)
                        full = np.array([[Y[k][int(sl.borb[x])] for x in range(sl.Gp.n)] for k in range(4)], np.int64)
                        ver = bool((sl.Gp.paf(full).sum(0) == sl.T1).all())
                        res.append(dict(y=[list(r) for r in key], verified=ver))
                        with open(OUT, 'a') as f: f.write(json.dumps(res[-1]) + "\n")
                        nfound += 1
    log(stage='done', stage1_solutions=len(res), all_verified=all(r['verified'] for r in res))


if __name__ == "__main__":
    main()
