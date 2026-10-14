r"""647 x183 (and x353): doubly-staged slice, cross level given BOTH compressions [M for counts, E for the Gaussian].

Slice: x_k on Z_646 = Z_2 x Z_17 x Z_19 constant on the orbits of M = <mu>.  Compressions: y_B along Z_19 (onto Z_34)
and y_A along Z_17 (onto Z_38), both exact stage-1 solutions (theory/nd_s1_sd.py, logs/nd/s1sols/).  They must agree
on the Z_2 level (per block: the sums over the two Z_2-cosets).  Given an aligned pair (block order, negation of
row-0 blocks, translation by the M-fixed element of order 2), every block's labels must reproduce both compressions:
per block and Z_2-coset e, a 0/1-free linear system on the orbit labels, enumerated exactly by DP (counts and lists).
Then:
  H = log2 #(label configurations consistent with y_A and y_B)  (exact),
  C = lattice Gaussian of the full core equations (one lag per <M,-1>-orbit of Z_646) under uniform sampling,
  E[completions | y_A, y_B] = 2^(H - C), and if the lists are small enough an exact 2+2 join over the blocks finds
  every core with these compressions (each verified by the cyclic PAF on Z_646).
usage: python3 theory/nd_647_joint.py MU [--pairs 3] [--N 20000] [--join] [--maxlist 2000000]"""
import sys, os, json, math, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import mgroup, lag_reps
from q523_gauss import gauss_bits

m = 646

def structure(mu):
    Mg = mgroup(m, mu)
    orb = -np.ones(m, np.int64); reps = []
    for g in range(m):
        if orb[g] >= 0: continue
        i = len(reps)
        for a in Mg: orb[g * a % m] = i
        reps.append(g)
    return Mg, orb, len(reps)

def compress(x, K):
    """x on Z_646 -> y on Z_{646/K}: y(b) = sum_{g = b mod F} x(g)"""
    F = m // K; return np.array([x[b::F].sum() for b in range(F)])

def z2data(Y, F):
    return [(int(a), int(b)) for a, b in zip(np.array(Y)[:, 0::2].sum(1), np.array(Y)[:, 1::2].sum(1))]

def align(YA, YB):
    """find block order / negation / e-swap of YA so that its Z_2 data equal YB's block by block; returns aligned YA or
    None.  e-swap = translation by 19 on Z_38 (order-2 element, fixed by every multiplier)."""
    FA, FB = 38, 34
    zb = z2data(YB, FB)
    variants = []
    for k in range(4):
        y = np.array(YA[k]); vs = []
        for neg in (1, -1):
            for sh in (0, 19):
                yy = neg * np.roll(y, -sh); d = (int(yy[0::2].sum()), int(yy[1::2].sum()))
                vs.append((d, yy))
        variants.append(vs)
    for perm in itertools.permutations(range(4)):
        out = []
        for kb in range(4):
            cand = [yy for d, yy in variants[perm[kb]] if d == zb[kb]]
            if not cand: break
            out.append(cand[0])
        if len(out) == 4: return [o.tolist() for o in out]
    return None

def block_lists(Mg, orb, nor, ya, yb, maxlist):
    """all label vectors (one +-1 per M-orbit) of one block with compressions ya (Z_38) and yb (Z_34): DP over orbits
    with the residual constraint vector as state; returns (count, list or None)"""
    # coefficient of orbit o in each constraint: yb rows = Z_34 positions, ya rows = Z_38 positions
    pts = defaultdict(list)
    for g in range(m): pts[int(orb[g])].append(g)
    rows_b = list(range(34)); rows_a = list(range(38))
    coef = []
    for o in range(nor):
        cb = np.zeros(34, np.int64); ca = np.zeros(38, np.int64)
        for g in pts[o]: cb[g % 34] += 1; ca[g % 38] += 1
        coef.append(np.concatenate([cb, ca]))
    coef = np.array(coef); target = np.concatenate([np.array(yb), np.array(ya)])
    # reduce to independent rows: keep the distinct coefficient rows (positions in the same M-orbit repeat)
    _, keep = np.unique(coef.T, axis=0, return_index=True); keep = np.sort(keep)
    C = coef[:, keep]; T = target[keep]
    # order orbits so that constraints close early (by first constraint index they touch, then by last)
    order = sorted(range(nor), key=lambda o: (np.nonzero(C[o])[0].max(), np.nonzero(C[o])[0].min()))
    last = {}
    for pos, o in enumerate(order):
        for r in np.nonzero(C[o])[0]: last[int(r)] = pos
    maxrem = np.zeros((nor + 1, C.shape[1]), np.int64)
    for pos in range(nor - 1, -1, -1): maxrem[pos] = maxrem[pos + 1] + np.abs(C[order[pos]])
    # DP: state = residual vector (tuple); forward with counts and back-pointers for listing
    states = [{tuple(T.tolist()): 1}]
    for pos, o in enumerate(order):
        nxt = defaultdict(int); co = C[o]
        for st, c in states[-1].items():
            sv = np.array(st)
            for lab in (1, -1):
                ns = sv - lab * co
                if (np.abs(ns) > maxrem[pos + 1]).any(): continue
                nxt[tuple(ns.tolist())] += c
        states.append(dict(nxt))
    zero = tuple([0] * C.shape[1]); cnt = states[-1].get(zero, 0)
    if cnt == 0 or cnt > maxlist: return cnt, None
    # list all solutions by backtracking through the DP layers
    sols = []
    def back(pos, st, lab):
        if pos == 0:
            sols.append(lab[::-1].copy()); return
        o = order[pos - 1]; co = C[o]
        for l in (1, -1):
            prev = tuple((np.array(st) + l * co).tolist())
            if prev in states[pos - 1]:
                lab.append(l); back(pos - 1, prev, lab); lab.pop()
    back(nor, zero, [])
    L = np.zeros((len(sols), nor), np.int8)
    for i, s_ in enumerate(sols):
        for pos, o in enumerate(order): L[i, o] = s_[pos]
    return cnt, L

def main(mu, npairs=3, N=20000, do_join=False, maxlist=2000000):
    Mg, orb, nor = structure(mu)
    A = json.load(open(f'logs/nd/s1sols/cyc_646_{mu}_17.json'))['solutions']
    B = json.load(open(f'logs/nd/s1sols/cyc_646_{mu}_19.json'))['solutions']
    lr = lag_reps(m, Mg); T = np.array([4 * m] + [-4] * (len(lr) - 1), np.int64)
    done = 0; rng = np.random.default_rng(1)
    for ib, xb in enumerate(B):
        for ia, xa in enumerate(A):
            YA = align(xa['y'], xb['y'])
            if YA is None: continue
            t0 = time.time()
            res = [block_lists(Mg, orb, nor, YA[k], xb['y'][k], maxlist) for k in range(4)]
            cnts = [r[0] for r in res]
            rec = dict(mu=mu, yA=ia, yB=ib, block_counts=cnts)
            if all(c > 0 for c in cnts):
                H = sum(math.log2(c) for c in cnts); rec['H'] = round(H, 2)
                if all(r[1] is not None for r in res):
                    Ls = [r[1] for r in res]
                    V = np.zeros((N, len(lr)), np.int64)
                    for st in range(0, N, 1000):
                        c = min(1000, N - st); acc = 0
                        for k in range(4):
                            lab = Ls[k][rng.integers(len(Ls[k]), size=c)]
                            x = lab[:, orb].astype(np.float64); fx = np.fft.rfft(x, axis=1)
                            acc = acc + np.fft.irfft(fx * fx.conj(), n=m, axis=1)
                        V[st:st + c] = np.rint(acc[:, lr]).astype(np.int64)
                    gb = gauss_bits(V, T)
                    rec.update(C=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'),
                               log2E=round(H - gb['bits'], 2) if gb.get('in_lattice') else None, list_sizes=[len(l) for l in Ls])
                    if do_join:
                        rec.update(join_core(Ls, orb, lr, T))
            rec['secs'] = round(time.time() - t0, 1)
            print(json.dumps(rec), flush=True)
            done += 1
            if done >= npairs: return
    print(json.dumps(dict(note="no more aligned pairs", done=done)))

def join_core(Ls, orb, lr, T):
    """exact 2+2 join of the four block lists on the full lag signature; every hit verified on Z_646"""
    sig = []
    for L in Ls:
        S = np.zeros((len(L), len(lr)), np.int64)
        for st in range(0, len(L), 4000):
            x = L[st:st + 4000][:, orb].astype(np.float64); fx = np.fft.rfft(x, axis=1)
            S[st:st + 4000] = np.rint(np.fft.irfft(fx * fx.conj(), n=m, axis=1)[:, lr]).astype(np.int64)
        sig.append(S)
    w = np.random.default_rng(7).integers(1, 2 ** 61, len(lr), dtype=np.int64)
    h = [S @ w for S in sig]
    order = sorted(range(4), key=lambda k: len(Ls[k]))
    a, b, c, d = order
    left = (h[a][:, None] + h[b][None, :]).ravel(); o = np.argsort(left); ls = left[o]
    hT = T @ w; hits = []
    for i2 in range(len(h[c])):
        r = hT - h[c][i2] - h[d]
        lo = np.searchsorted(ls, r, 'left'); hi = np.searchsorted(ls, r, 'right')
        for j in np.nonzero(hi > lo)[0]:
            for t in range(lo[j], hi[j]):
                i0, i1 = divmod(int(o[t]), len(h[b]))
                if ((sig[a][i0] + sig[b][i1] + sig[c][i2] + sig[d][j]) == T).all():
                    lab = {a: Ls[a][i0], b: Ls[b][i1], c: Ls[c][i2], d: Ls[d][j]}
                    X = np.array([lab[k][orb] for k in range(4)], np.int64)
                    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X)
                    ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
                    hits.append(dict(core=X.tolist(), verified=ok, rows=X.sum(1).tolist()))
                    if len(hits) >= 5: return dict(cores=hits)
    return dict(cores=hits)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--pairs": "3", "--N": "20000", "--maxlist": "2000000"}; flags = set(); pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(int(a[i])); i += 1
    main(pos[0], int(opt["--pairs"]), int(opt["--N"]), '--join' in flags, int(opt["--maxlist"]))

def ya_options(Mg, orb, nor, yb_k):
    """for one block with fixed y_B (Z_34 compression), all achievable y_A vectors (Z_38 orbit values) with the number
    of label configurations giving each: DP over orbits, state = (residual of the Z_34 rows, partial Z_38 sums)"""
    pts = defaultdict(list)
    for g in range(m): pts[int(orb[g])].append(g)
    cb = np.zeros((nor, 34), np.int64); ca = np.zeros((nor, 38), np.int64)
    for o in range(nor):
        for g in pts[o]: cb[o, g % 34] += 1; ca[o, g % 38] += 1
    _, kb = np.unique(cb.T, axis=0, return_index=True); kb = np.sort(kb)
    _, ka = np.unique(ca.T, axis=0, return_index=True); ka = np.sort(ka)
    CB = cb[:, kb]; CA = ca[:, ka]; TB = np.array(yb_k)[kb]
    order = sorted(range(nor), key=lambda o: (np.nonzero(CB[o])[0].max() if CB[o].any() else -1))
    maxremB = np.zeros((nor + 1, CB.shape[1]), np.int64)
    for pos in range(nor - 1, -1, -1): maxremB[pos] = maxremB[pos + 1] + np.abs(CB[order[pos]])
    states = {(tuple(TB.tolist()), tuple([0] * CA.shape[1])): 1}
    for pos, o in enumerate(order):
        nxt = defaultdict(int)
        for (rb, pa), c in states.items():
            rbv = np.array(rb); pav = np.array(pa)
            for lab in (1, -1):
                nr = rbv - lab * CB[o]
                if (np.abs(nr) > maxremB[pos + 1]).any(): continue
                nxt[(tuple(nr.tolist()), tuple((pav + lab * CA[o]).tolist()))] += c
        states = nxt
    zero = tuple([0] * CB.shape[1])
    out = {pa: c for (rb, pa), c in states.items() if rb == zero}
    return ka, out

def joint_from_yB(mu, ib, maxpairs=4 * 10 ** 7):
    """joint stage 1 for a fixed y_B: per block the achievable y_A (Z_38 orbit values) with multiplicities, then an
    exact 2+2 join on the stage-1A equations (sum_k PAF_{Z_38}(y_A) = (2520, -68, ...), rows equal to y_B's rows).
    Returns the list of (y_A per block as Z_38 sequences, multiplicities)"""
    Mg, orb, nor = structure(mu)
    B = json.load(open(f'logs/nd/s1sols/cyc_646_{mu}_19.json'))['solutions']
    yb = B[ib]['y']
    # Z_38 orbit structure (M acting on Z_38) for expanding orbit values to sequences
    M38 = sorted({a % 38 for a in Mg})
    lags = lag_reps(38, M38); T = np.array([2520] + [-68] * (len(lags) - 1), np.int64)
    blocks = []
    for k in range(4):
        ka, opts = ya_options(Mg, orb, nor, yb[k])
        # ka indexes Z_38 positions that are representatives of distinct coefficient columns = orbit reps
        reps = list(ka)
        # map every Z_38 position to its representative's index
        rep_of = {}
        for b in range(38):
            for idx, r in enumerate(reps):
                if any((r * a) % 38 == b for a in M38): rep_of[b] = idx; break
        seqs = []; mult = []
        for pa, c in opts.items():
            y = np.array([pa[rep_of[b]] for b in range(38)])
            seqs.append(y); mult.append(c)
        blocks.append((np.array(seqs).reshape(-1, 38), np.array(mult, dtype=object)))
    sizes = [len(b[0]) for b in blocks]
    sig = [np.stack([(S * np.roll(S, -h, axis=1)).sum(1) for h in lags], 1) if len(S) else np.zeros((0, len(lags)), np.int64)
           for S, _ in blocks]
    rec = dict(mu=mu, yB=ib, ya_options=sizes, ya_mult_log2=[round(math.log2(sum(b[1])), 1) if len(b[1]) else None for b in blocks])
    print(json.dumps(dict(stage='ya_options', yB=ib, sizes=sizes)), flush=True)
    if min(sizes) == 0: rec['joint'] = 0; return rec, []
    w = np.random.default_rng(3).integers(1, 2 ** 61, len(lags), dtype=np.int64)
    h = [s_ @ w for s_ in sig]
    # memory-bounded join: sort the smaller pair side once (if it fits), stream the other side
    order = sorted(range(4), key=lambda k: sizes[k]); a0, a1, a2, a3 = order
    if sizes[a0] * sizes[a1] > 3 * 10 ** 8:
        rec['joint'] = None; rec['note'] = 'left side too large'; return rec, []
    left = (h[a0][:, None] + h[a1][None, :]).ravel(); o = np.argsort(left); ls = left[o]
    hT = T @ w; sols = []
    for i2 in range(sizes[a2]):
        r = hT - h[a2][i2] - h[a3]
        lo = np.searchsorted(ls, r, 'left'); hi = np.searchsorted(ls, r, 'right')
        for j in np.nonzero(hi > lo)[0]:
            for t in range(lo[j], hi[j]):
                i0, i1 = divmod(int(o[t]), sizes[a1])
                if ((sig[a0][i0] + sig[a1][i1] + sig[a2][i2] + sig[a3][j]) == T).all():
                    idx = {a0: i0, a1: i1, a2: i2, a3: int(j)}
                    YA = [blocks[k][0][idx[k]] for k in range(4)]
                    mult = [blocks[k][1][idx[k]] for k in range(4)]
                    sols.append(([y.tolist() for y in YA], [int(x) for x in mult]))
                    if len(sols) >= 2000: break
            if len(sols) >= 2000: break
        if len(sols) >= 2000: break
    rec['joint'] = len(sols)
    if sols: rec['log2_mult_first'] = round(sum(math.log2(x) for x in sols[0][1]), 1)
    return rec, sols

def cross_level(mu, YA, YB, N=20000, do_join=True, maxlist=3 * 10 ** 6, seed=1):
    """given a joint stage-1 solution (y_A on Z_38, y_B on Z_34, blocks aligned): per-block label lists consistent with
    both, H (exact), the lattice Gaussian C of the full equations, and an exact 2+2 join for cores if the lists fit"""
    Mg, orb, nor = structure(mu)
    lr = lag_reps(m, Mg); T = np.array([4 * m] + [-4] * (len(lr) - 1), np.int64)
    res = [block_lists(Mg, orb, nor, YA[k], YB[k], maxlist) for k in range(4)]
    cnts = [r[0] for r in res]; rec = dict(block_counts=cnts)
    if min(cnts) == 0: return rec
    H = sum(math.log2(c) for c in cnts); rec['H'] = round(H, 2)
    if any(r[1] is None for r in res): rec['note'] = 'lists too large to enumerate'; return rec
    Ls = [r[1] for r in res]; rng = np.random.default_rng(seed)
    V = np.zeros((N, len(lr)), np.int64)
    for st in range(0, N, 1000):
        c = min(1000, N - st); acc = 0
        for k in range(4):
            lab = Ls[k][rng.integers(len(Ls[k]), size=c)]
            x = lab[:, orb].astype(np.float64); fx = np.fft.rfft(x, axis=1)
            acc = acc + np.fft.irfft(fx * fx.conj(), n=m, axis=1)
        V[st:st + c] = np.rint(acc[:, lr]).astype(np.int64)
    gb = gauss_bits(V, T)
    rec.update(C=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'),
               log2E=round(H - gb['bits'], 2) if gb.get('in_lattice') else None)
    if do_join:
        p2 = sorted(len(L) for L in Ls)
        rec['join_left_pairs_log2'] = round(math.log2(p2[0] * p2[1]), 1)
        if p2[0] * p2[1] <= 2 * 10 ** 8:
            rec.update(join_core(Ls, orb, lr, T))
    return rec
