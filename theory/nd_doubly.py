r"""Doubly-staged slices, general form (647 x183 and its small controls) [M for every exact step, E for Gaussians].

Slice: x_k on Z_m constant on the orbits of M = <mu>; qA, qB are coprime Hall factors of m (647: m = 2*17*19 with
qA = 17, qB = 19; 971: m = 970 with qB = 97 (stage 1 on Z_10) and qA = 10 (the dual compression onto Z_97)).
Compressions: y_B along Z_qB (onto Z_{m/qB}) and y_A along Z_qA (onto Z_{m/qA}).  Pipeline, all exact:
  1. stage 1B: every compressed solution y_B (theory/nd_s1_sd.py, F = m/qB = 2 qA, n = qA odd);
  2. joint stage 1: for one y_B, per block every achievable y_A (orbit values on Z_{m/qA}) with its label multiplicity
     (DP over the M-orbits of Z_m, state = residual of the y_B rows and partial y_A sums), then an exact 2+2 join on the
     stage-1A equations (sum_k PAF(y_A) = (4m - 4(qA - 1), -4 qA, ...), rows = y_B's rows);
  3. cross level: per block every label vector consistent with both compressions (DP), and an exact 2+2 join of the
     four lists on the full core equations (one lag per <M, -1>-orbit of Z_m); every hit verified by the cyclic PAF;
  4. completion: bordered Goethals-Seidel over Z_m (circulants), H H^T = N I checked in int64.
usage: python3 theory/nd_doubly.py m mu qA qB [--yB N] [--maxjoint 50] [--out FILE]"""
import sys, os, json, math, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import mgroup, lag_reps, Problem
from q523_gauss import gauss_bits

class Doubly:
    def __init__(self, m, mu, qA, qB):
        self.m, self.mu, self.qA, self.qB = m, mu, qA, qB
        assert m % qA == 0 and m % qB == 0 and math.gcd(qA, m // qA) == 1 and math.gcd(qB, m // qB) == 1
        self.FA, self.FB = m // qA, m // qB
        Mg = mgroup(m, mu); self.Mg = Mg
        orb = -np.ones(m, np.int64); reps = []
        for g in range(m):
            if orb[g] >= 0: continue
            i = len(reps)
            for a in Mg: orb[g * a % m] = i
            reps.append(g)
        self.orb = orb; self.nor = len(reps)
        pts = defaultdict(list)
        for g in range(m): pts[int(orb[g])].append(g)
        self.cA = np.zeros((self.nor, self.FA), np.int64); self.cB = np.zeros((self.nor, self.FB), np.int64)
        for o in range(self.nor):
            for g in pts[o]: self.cA[o, g % self.FA] += 1; self.cB[o, g % self.FB] += 1
        _, kA = np.unique(self.cA.T, axis=0, return_index=True); self.kA = np.sort(kA)
        _, kB = np.unique(self.cB.T, axis=0, return_index=True); self.kB = np.sort(kB)
        self.lr = lag_reps(m, Mg)
        self.T = np.array([4 * m] + [-4] * (len(self.lr) - 1), np.int64)
        MA = sorted({a % self.FA for a in Mg}); self.MA = MA
        self.lagsA = lag_reps(self.FA, MA)
        self.TA = np.array([4 * m - 4 * (qA - 1)] + [-4 * qA] * (len(self.lagsA) - 1), np.int64)
        # every Z_FA position -> index of its representative column in kA
        self.repA = np.array([next(i for i, r in enumerate(self.kA) if any((r * a) % self.FA == b for a in MA))
                              for b in range(self.FA)])

    # ---------------- stage 1B ----------------
    def stage1B(self, max_sol=200):
        P = Problem(self.m, self.mu, self.qB)
        out, sols = P.solve(max_sol=max_sol, log=lambda x: None)
        return out, [s_['y'] for s_ in sols]

    # ---------------- joint stage 1 ----------------
    def ya_options(self, yb_k):
        CB = self.cB[:, self.kB]; CA = self.cA[:, self.kA]; TB = np.array(yb_k)[self.kB]
        nor = self.nor
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
        return {pa: c for (rb, pa), c in states.items() if rb == zero}

    def joint(self, yb, maxsol=2000, maxleft=3 * 10 ** 8, TA=None):
        TA = self.TA if TA is None else TA
        blocks = []
        for k in range(4):
            opts = self.ya_options(yb[k])
            seqs = np.array([[pa[self.repA[b]] for b in range(self.FA)] for pa in opts], np.int64).reshape(-1, self.FA)
            mult = np.array(list(opts.values()), dtype=object)
            blocks.append((seqs, mult))
        sizes = [len(b[0]) for b in blocks]
        rec = dict(ya_options=sizes)
        if min(sizes) == 0: rec['joint'] = 0; return rec, []
        # rows of y_A must equal rows of y_B (same block): filter
        rows = [int(np.sum(yb[k])) for k in range(4)]
        for k in range(4):
            keep = blocks[k][0].sum(1) == rows[k]
            blocks[k] = (blocks[k][0][keep], blocks[k][1][keep])
        sizes = [len(b[0]) for b in blocks]; rec['ya_options_rowok'] = sizes
        if min(sizes) == 0: rec['joint'] = 0; return rec, []
        sig = [np.stack([(S * np.roll(S, -h, axis=1)).sum(1) for h in self.lagsA], 1) for S, _ in blocks]
        w = np.random.default_rng(3).integers(1, 2 ** 61, len(self.lagsA), dtype=np.int64)
        h = [s_ @ w for s_ in sig]
        order = sorted(range(4), key=lambda k: sizes[k]); a0, a1, a2, a3 = order
        if sizes[a0] * sizes[a1] > maxleft: rec['joint'] = None; rec['note'] = 'left side too large'; return rec, []
        with np.errstate(over='ignore'):
            left = (h[a0][:, None] + h[a1][None, :]).ravel()
        o = np.argsort(left); ls = left[o]; hT = TA @ w; sols = []
        for i2 in range(sizes[a2]):
            with np.errstate(over='ignore'):
                r = hT - h[a2][i2] - h[a3]
            lo = np.searchsorted(ls, r, 'left'); hi = np.searchsorted(ls, r, 'right')
            for j in np.nonzero(hi > lo)[0]:
                for t in range(lo[j], hi[j]):
                    i0, i1 = divmod(int(o[t]), sizes[a1])
                    if ((sig[a0][i0] + sig[a1][i1] + sig[a2][i2] + sig[a3][j]) == TA).all():
                        idx = {a0: i0, a1: i1, a2: i2, a3: int(j)}
                        sols.append(([blocks[k][0][idx[k]].tolist() for k in range(4)], [int(blocks[k][1][idx[k]]) for k in range(4)]))
                        if len(sols) >= maxsol: break
                if len(sols) >= maxsol: break
            if len(sols) >= maxsol: break
        rec['joint'] = len(sols)
        if sols: rec['H_first'] = round(sum(math.log2(x) for x in sols[0][1]), 2)
        return rec, sols

    # ---------------- cross level ----------------
    def block_lists(self, ya, yb, maxlist):
        C = np.concatenate([self.cB[:, self.kB], self.cA[:, self.kA]], 1)
        Tt = np.concatenate([np.array(yb)[self.kB], np.array(ya)[self.kA]])
        nor = self.nor
        order = sorted(range(nor), key=lambda o: (np.nonzero(C[o])[0].max(), np.nonzero(C[o])[0].min()))
        maxrem = np.zeros((nor + 1, C.shape[1]), np.int64)
        for pos in range(nor - 1, -1, -1): maxrem[pos] = maxrem[pos + 1] + np.abs(C[order[pos]])
        states = [{tuple(Tt.tolist()): 1}]
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
        sols = []
        def back(pos, st, lab):
            if pos == 0: sols.append(lab[::-1].copy()); return
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

    def cross(self, YA, YB, N=20000, maxlist=3 * 10 ** 6, do_join=True, maxleft=2 * 10 ** 8, seed=1, T=None, verify=True, maxcores=3):
        m = self.m; orb = self.orb; lr = self.lr; T = self.T if T is None else T
        res = [self.block_lists(YA[k], YB[k], maxlist) for k in range(4)]
        cnts = [r[0] for r in res]; rec = dict(block_counts=cnts)
        if min(cnts) == 0: return rec, []
        rec['H'] = round(sum(math.log2(c) for c in cnts), 2)
        if any(r[1] is None for r in res): rec['note'] = 'lists too large'; return rec, []
        Ls = [r[1] for r in res]
        sig = []
        for L in Ls:
            S = np.zeros((len(L), len(lr)), np.int64)
            for st in range(0, len(L), 4000):
                x = L[st:st + 4000][:, orb].astype(np.float64); fx = np.fft.rfft(x, axis=1)
                S[st:st + 4000] = np.rint(np.fft.irfft(fx * fx.conj(), n=m, axis=1)[:, lr]).astype(np.int64)
            sig.append(S)
        rng = np.random.default_rng(seed)
        V = sum(sig[k][rng.integers(len(sig[k]), size=N)] for k in range(4))
        gb = gauss_bits(V, T)
        rec.update(C=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'),
                   log2E=round(rec['H'] - gb['bits'], 2) if gb.get('in_lattice') else None)
        cores = []
        if do_join:
            sizes = [len(L) for L in Ls]
            order = sorted(range(4), key=lambda k: sizes[k]); a0, a1, a2, a3 = order
            rec['left_pairs_log2'] = round(math.log2(sizes[a0] * sizes[a1]), 1)
            if sizes[a0] * sizes[a1] <= maxleft:
                w = np.random.default_rng(7).integers(1, 2 ** 61, len(lr), dtype=np.int64)
                h = [S @ w for S in sig]
                with np.errstate(over='ignore'):
                    left = (h[a0][:, None] + h[a1][None, :]).ravel()
                o = np.argsort(left); ls = left[o]; hT = T @ w
                for i2 in range(sizes[a2]):
                    with np.errstate(over='ignore'):
                        r = hT - h[a2][i2] - h[a3]
                    lo = np.searchsorted(ls, r, 'left'); hi = np.searchsorted(ls, r, 'right')
                    for j in np.nonzero(hi > lo)[0]:
                        for t in range(lo[j], hi[j]):
                            i0, i1 = divmod(int(o[t]), sizes[a1])
                            if ((sig[a0][i0] + sig[a1][i1] + sig[a2][i2] + sig[a3][j]) == T).all():
                                idx = {a0: i0, a1: i1, a2: i2, a3: int(j)}
                                X = np.array([Ls[k][idx[k]][orb] for k in range(4)], np.int64)
                                P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X)
                                if not verify or (P[0] == 4 * m and (P[1:] == -4).all()): cores.append(X)
                                if len(cores) >= maxcores: break
                        if len(cores) >= maxcores: break
                    if len(cores) >= maxcores: break
                rec['cores'] = len(cores)
        return rec, cores

def planted(D, seed=0):
    """plant a random label configuration; its own compressions and PAF sums are the targets; check that the joint step
    returns its y_A and the cross level returns the planted labels"""
    rng = np.random.default_rng(seed); m = D.m
    lab = rng.choice([-1, 1], (4, D.nor)); X = lab[:, D.orb]
    yb = [[int(X[k][b::D.FB].sum()) for b in range(D.FB)] for k in range(4)]
    ya = [[int(X[k][b::D.FA].sum()) for b in range(D.FA)] for k in range(4)]
    TA = np.array([sum(int((np.array(ya[k]) * np.roll(np.array(ya[k]), -h)).sum()) for k in range(4)) for h in D.lagsA], np.int64)
    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X.astype(float))
    T = P[D.lr]
    rj, sols = D.joint(yb, TA=TA, maxsol=10 ** 6)
    found_A = any(all(list(s_[0][k]) == ya[k] for k in range(4)) for s_ in sols)
    rc, cores = D.cross(ya, yb, T=T, verify=False, N=2000, maxcores=10 ** 6)
    found_X = any((c == X).all() for c in cores) if cores else False
    return dict(seed=seed, ya_options=rj.get('ya_options'), joint=rj.get('joint'), planted_yA_found=found_A,
                block_counts=rc.get('block_counts'), cross_hits=len(cores), planted_core_found=found_X)

def bordered_gs_cyclic(X):
    """bordered GS completion of a b = 1 core of circulants (pattern of scripts/verify_h2292.py); H or None"""
    m = X.shape[1]; N = 4 * m + 4
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])
    R = np.fliplr(np.eye(m, dtype=np.int64))
    A, B, C, D = [circ(x).astype(np.int64) for x in X]
    core = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])
    Sig = np.array([[core[I * m:(I + 1) * m, J * m:(J + 1) * m].sum(1)[0] for J in range(4)] for I in range(4)])
    H4 = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    H4 = [h for h in H4 if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]
    for S in H4:
        for T_ in H4:
            num = -(S @ Sig.T @ T_)
            if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
            H = np.zeros((N, N), np.int64)
            H[:4, :4] = num // 4
            H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(T_.T, m, axis=0); H[4:, 4:] = core
            if np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64)): return H
    return None

def run(m, mu, qA, qB, nyB=20, maxjoint=50, out=""):
    D = Doubly(m, mu, qA, qB)
    t0 = time.time()
    o1, YBs = D.stage1B(max_sol=max(nyB, 1))
    print(json.dumps(dict(stage="1B", m=m, mu=mu, qA=qA, qB=qB, orbits=D.nor, result=o1, secs=round(time.time() - t0, 1))), flush=True)
    for ib, yb in enumerate(YBs[:nyB]):
        rj, sols = D.joint(yb)
        print(json.dumps(dict(stage="joint", yB=ib, **rj, secs=round(time.time() - t0, 1))), flush=True)
        for js, (YA, mult) in enumerate(sols[:maxjoint]):
            rc, cores = D.cross(YA, yb)
            print(json.dumps(dict(stage="cross", yB=ib, joint=js, **rc, secs=round(time.time() - t0, 1))), flush=True)
            if cores:
                X = cores[0]; H = bordered_gs_cyclic(X)
                rec = dict(stage="core", v=m + 1, yB=ib, joint=js, hadamard=H is not None, N=int(H.shape[0]) if H is not None else None,
                           core=X.tolist())
                print(json.dumps({k: rec[k] for k in rec if k != 'core'}), flush=True)
                if out: json.dump(rec, open(out, "w"))
                return rec
    return None

if __name__ == "__main__":
    if "--plant" in sys.argv:
        i = sys.argv.index("--plant"); n = int(sys.argv[i + 1]); args = [int(x) for x in sys.argv[1:i]]
        D = Doubly(*args)
        for sd in range(n): print(json.dumps(planted(D, sd)), flush=True)
        sys.exit()
    a = sys.argv[1:]; opt = {"--yB": "20", "--maxjoint": "50", "--out": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    run(*pos, nyB=int(opt["--yB"]), maxjoint=int(opt["--maxjoint"]), out=opt["--out"])
