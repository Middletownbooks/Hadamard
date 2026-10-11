r"""Stage 2 by a two-level (Wagner) tree on the S/D split: the first merge is on the Z3 level, done at the level of
parents (row sums), the second on the fine S-signature; then a D-lag-0 filter and the D-check
(report/paley_523_answer.md, section 3).

Per block k, the S-configurations (one S-value per pair g < F2, F2 = F/2) are grouped by their Z3 value phi = (row sums
of the rows g = r mod 3).  The Z3 signature Z3sig(phi) (the Z_3-compressed S-level PAF; 8 integers) is additive over
blocks, so for a middle value c every left parent pair (phi_0, phi_1) with Z3sig(phi_0) + Z3sig(phi_1) = c combines
with every right parent pair with Z3sig(phi_2) + Z3sig(phi_3) = T3 - c.  The S-solutions are the quadruples of
configurations from these fibres whose full S-signatures (S-level PAF on Z_F2, lags 0..F2/2) sum to T_S: a 2+2 join
whose sides are unions of fibre products over many parents (so a side is not limited to one parent's fibres).
Fibres hold up to K configurations per row (sampled without replacement with weight prod m(s)^beta, m(s) = number of
D-options of the pair value s).  Join: residue passes (additive hash mod P), cache-sized open-addressing table,
numba; every hit is re-checked exactly.
D-lag-0 filter: the D-level lag-0 value sum_g |D(g)|^2 depends on S except on pairs whose options differ in norm; an
S-solution is kept only if the D-level lag-0 target is reachable (exact subset-sum over those pairs).
D-check: theory/q523_sdsolve.d_check (Z3 class split on D, then the full D-level join); assembly and 2-D FFT check.
usage: python3 theory/q523_wagner.py F q s y.json [--beta 1.5] [--rows 24] [--K 24] [--logN 30] [--ncent 4]
          [--seed 0] [--table-bits 20] [--out core.txt] [--planted]"""
import sys, os, json, math, time, itertools as it
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import zcoords, pz
from q523_sdsolve import SDInst, seq_sig, next_prime, d_check, assemble, verify, write_core, planted_target, join_all
from q523_dcheck_nb import DCheckNB as DChecker

# ------------------------------------------------------------------ candidates and fibres
def block_candidates(inst, k, beta, nrow, K, rng):
    """per row: the nrow row sums with the largest beta-mass, each with up to K configurations (Gumbel top-K);
    returns rows = [list of (rowsum tuple, configs (n, 3) pair-S indices, log2 beta-mass)] for r = 0, 1, 2"""
    rows = []
    for r in range(3):
        g = inst.row_groups(k, r)
        keys = list(g); mass = np.array([(g[kk][1] ** beta).sum() for kk in keys])
        top = np.argsort(-mass)[:nrow]
        lst = []
        for i in top:
            idx, w = g[keys[i]]
            if len(w) > K:
                gk = beta * np.log(w) + rng.gumbel(size=len(w)); pick = np.argsort(-gk)[:K]; idx = idx[pick]
            lst.append((keys[i], idx, math.log2(mass[i])))
        rows.append(lst)
    return rows

def z3sig_of(inst, rowsums):
    """rowsums (n, 3, s) -> Z3 signatures (n, 8)"""
    return seq_sig(np.asarray(rowsums, np.int64), inst.cy.T, 2)

def cand_table(inst, rows):
    """all row-sum triples: arrays (n, 3) of row-candidate indices, Z3 signatures, log2 fibre sizes, log2 mass"""
    n0, n1, n2 = (len(rows[r]) for r in range(3))
    grid = np.array(list(it.product(range(n0), range(n1), range(n2))), np.int64)
    rs = np.stack([np.array([rows[r][grid[i, r]][0] for i in range(len(grid))]) for r in range(3)], 1)
    sig = z3sig_of(inst, rs)
    lsz = sum(np.array([math.log2(len(rows[r][j][1])) for j in range(len(rows[r]))])[grid[:, r]] for r in range(3))
    lms = sum(np.array([rows[r][j][2] for j in range(len(rows[r]))])[grid[:, r]] for r in range(3))
    return grid, sig, lsz, lms

class SigHash:
    """linear 64-bit hash of Z3 signatures, sorted for vectorised lookups (hits are re-checked exactly)"""
    def __init__(self, sig, w):
        self.sig = sig; self.w = w
        self.k = (sig.astype(np.uint64) * w[None, :]).sum(1); self.o = np.argsort(self.k, kind='stable'); self.ks = self.k[self.o]
    def key(self, row): return (np.asarray(row, np.int64).astype(np.uint64) * self.w).sum()

def parent_pairs(sigA, kA, hB, c):
    """all (a, b) with sigA[a] + sigB[b] == c (vectorised: needed key = key(c) - kA, linear hash)"""
    need = hB.key(c) - kA
    lo = np.searchsorted(hB.ks, need, 'left'); hi = np.searchsorted(hB.ks, need, 'right')
    out = []
    for a in np.nonzero(hi > lo)[0]:
        for t in range(lo[a], hi[a]):
            b = hB.o[t]
            if (sigA[a] + hB.sig[b] == c).all(): out.append((int(a), int(b)))
    return out

def choose_middle(inst, C, rng, ncand=200):
    """pick c maximising min(log2 left configuration pairs, log2 right ones) over c values of random left parent pairs
    drawn by beta-mass (the right side must reach T3 - c)"""
    (g0, s0, l0, m0), (g1, s1, l1, m1), (g2, s2, l2, m2), (g3, s3, l3, m3) = C
    w = np.random.default_rng(12345).integers(1, 2 ** 62, s0.shape[1], dtype=np.int64).astype(np.uint64)
    h1, h3 = SigHash(s1, w), SigHash(s3, w)
    k0 = (s0.astype(np.uint64) * w[None, :]).sum(1); k2 = (s2.astype(np.uint64) * w[None, :]).sum(1)
    T3 = inst.tsig3()
    p0 = np.exp2(m0 - m0.max()); p0 /= p0.sum(); p1 = np.exp2(m1 - m1.max()); p1 /= p1.sum()
    best = None; tried = set()
    for _ in range(ncand):
        a = rng.choice(len(s0), p=p0); b = rng.choice(len(s1), p=p1); c = s0[a] + s1[b]
        kk = c.tobytes()
        if kk in tried: continue
        tried.add(kk)
        R = parent_pairs(s2, k2, h3, T3 - c)
        if not R: continue
        L = parent_pairs(s0, k0, h1, c)
        wl = math.log2(sum(2.0 ** (l0[i] + l1[j]) for i, j in L)); wr = math.log2(sum(2.0 ** (l2[i] + l3[j]) for i, j in R))
        sc = min(wl, wr)
        if best is None or sc > best[0]: best = (sc, c, L, R, wl, wr)
    return best

def materialise(inst, k, rows, grid_row, K):
    """S-configurations of one fibre: (n, F2) pair-S indices and S sequences (n, F2, s)"""
    sel = [rows[r][grid_row[r]][1] for r in range(3)]
    n = int(np.prod([len(x) for x in sel]))
    pidx = np.zeros((n, inst.F2), np.int64)
    g3 = np.array(list(it.product(*[range(len(x)) for x in sel])), np.int64)
    for r in range(3):
        for t, g in enumerate(inst.rows[r]): pidx[:, g] = sel[r][g3[:, r], t]
    S = np.zeros((n, inst.F2, inst.s), np.int64)
    for g in range(inst.F2): S[:, g] = inst.pairs[k][g]['S'][pidx[:, g]]
    return pidx, S

# ------------------------------------------------------------------ multi-parent join kernel
@numba.njit(cache=True)
def mjoin(keyA, resA, offA, keyB, rstB, offB, keyC, resC, offC, keyD, rstD, offD,
          lpA, lpB, rpC, rpD, P, rT, kT, tbits, maxhits, p_lo, p_hi):
    """left pairs: for j, a in fibre lpA[j] of A, b in fibre lpB[j] of B with (resA[a] + resB[b]) % P == p;
    B fibres are CSR'd by residue: rstB[f, r] .. rstB[f, r + 1] index into keyB order (already sorted per fibre).
    right pairs similarly with residue class (rT - p) % P.  Returns hits (global indices a, b, c, d)."""
    size = 1 << tbits; mask = size - 1
    tkey = np.zeros(size, np.uint64); ta = np.zeros(size, np.int64); tb = np.zeros(size, np.int64)
    stamp = np.zeros(size, np.int32)
    hits = np.zeros((maxhits, 4), np.int64); nh = 0; overflow = 0; probes = 0
    cap = (size * 80) // 100
    for p in range(p_lo, p_hi):
        st = p + 1; filled = 0
        for j in range(lpA.shape[0]):
            fa = lpA[j]; fb = lpB[j]
            for a in range(offA[fa], offA[fa + 1]):
                rb = (p - resA[a]) % P
                for t in range(rstB[fb, rb], rstB[fb, rb + 1]):
                    key = keyA[a] + keyB[t]
                    if filled >= cap:
                        overflow += 1; continue
                    slot = int(key >> np.uint64(64 - tbits)) & mask
                    while stamp[slot] == st: slot = (slot + 1) & mask
                    stamp[slot] = st; tkey[slot] = key; ta[slot] = a; tb[slot] = t; filled += 1
        probes += filled
        q = (rT - p) % P
        for j in range(rpC.shape[0]):
            fc = rpC[j]; fd = rpD[j]
            for c in range(offC[fc], offC[fc + 1]):
                rd = (q - resC[c]) % P
                for t in range(rstD[fd, rd], rstD[fd, rd + 1]):
                    need = kT - keyC[c] - keyD[t]
                    slot = int(need >> np.uint64(64 - tbits)) & mask
                    while stamp[slot] == st:
                        if tkey[slot] == need:
                            if nh < maxhits:
                                hits[nh, 0] = ta[slot]; hits[nh, 1] = tb[slot]; hits[nh, 2] = c; hits[nh, 3] = t
                            nh += 1
                        slot = (slot + 1) & mask
                    probes += 1
    return hits[:min(nh, maxhits)], nh, overflow, probes

class Side:
    """concatenated fibres of one block for the join: keys/residues sorted per fibre by residue, CSR by residue;
    configurations kept as uint8 pair-S indices (S is rebuilt from them on demand)"""
    def __init__(self, fibres, wts, rho, P):
        keys = []; res = []; offs = [0]; rst = []; pid = []
        for pidx, sig in fibres:
            k = (sig.astype(np.int64).astype(np.uint64) * wts[None, :]).sum(1)
            r = ((sig.astype(np.int64) % P) * (rho % P)[None, :] % P).sum(1) % P
            o = np.argsort(r, kind='stable')
            pid.append(pidx[o].astype(np.uint8)); keys.append(k[o]); res.append(r[o])
            cnt = np.bincount(r, minlength=P); st = np.zeros(P + 1, np.int64); st[1:] = np.cumsum(cnt)
            rst.append(st + offs[-1]); offs.append(offs[-1] + len(k))
        self.key = np.concatenate(keys); self.res = np.concatenate(res).astype(np.int64)
        self.off = np.array(offs, np.int64); self.rst = np.stack(rst) if rst else np.zeros((0, P + 1), np.int64)
        self.pidx_all = np.concatenate(pid)

def S_of(inst, k, pidx):
    """S sequences (n, F2, s) of block k from pair-S indices (n, F2)"""
    pidx = np.atleast_2d(pidx).astype(np.int64)
    S = np.zeros((pidx.shape[0], inst.F2, inst.s), np.int64)
    for g in range(inst.F2): S[:, g] = inst.pairs[k][g]['S'][pidx[:, g]]
    return S

# ------------------------------------------------------------------ D-lag-0 filter
def dlag0_filter(inst, pidx4):
    """True if the D-level lag-0 target is reachable given S (exact: per pair the set of |D|^2 over its options)"""
    T = inst.cy.T; base = np.zeros(inst.s, np.int64); var = []
    for k in range(4):
        for g in range(inst.F2):
            Ds = inst.pairs[k][g]['D'][pidx4[k][g]]
            nr = np.unique(seq_sig(Ds[:, None, :], T, 1), axis=0)
            if len(nr) == 1: base += nr[0]
            else: var.append(nr)
    reach = {tuple(base)}
    for nr in var:
        reach = {tuple(np.array(x) + v) for x in reach for v in nr}
        if len(reach) > 200000: return True
    return tuple(inst.TD[0]) in reach

# ------------------------------------------------------------------ driver
def dcheck_list(inst, dchk, P4, start, out, log, stopfile, ckpt, seed, beta, t0, stats, found):
    """D-check the S-solutions P4[start:] (array (n, 4, F2) of pair indices); progress is written to ckpt + '.json'
    every 2000 checks so that a restarted worker resumes there (resume_ckpt).  Returns (checked, d3_solutions, core)."""
    nf = 0; nd3 = 0
    for i in range(start, len(P4)):
        p4 = [P4[i, k].astype(np.int64) for k in range(4)]
        nf += 1
        dsols, info = dchk.check(p4)
        nd3 += info['d3_solutions']
        for D4 in dsols:
            S4 = [S_of(inst, k, p4[k])[0] for k in range(4)]
            Z, AC = assemble(inst, S4, D4); ok, X = verify(inst, AC)
            stats["cores"] += int(ok); found.append(AC)
            log(dict(stage="CORE", ok=ok, secs=round(time.time() - t0, 1)))
            if ok and out: write_core(inst, X, out, f"q523_wagner seed={seed} beta={beta}")
            if ok and stopfile: open(stopfile, 'w').write(json.dumps(dict(seed=seed, out=out)))
            if ok: return nf, nd3, True
        if ckpt and (i + 1) % 2000 == 0:
            with open(ckpt + '.json', 'w') as fh: json.dump(dict(done=i + 1, n=len(P4)), fh)
    return nf, nd3, False

def resume_ckpt(inst, path, dchk, out, log, stopfile, seed=0, beta=1.0):
    """finish the D-checks of a checkpointed S-join (path = .npy of pair indices, progress in path[:-4] + '.json')"""
    t0 = time.time(); P4 = np.load(path); base = path[:-4]; start = 0
    if os.path.exists(base + '.json'):
        try: start = json.load(open(base + '.json'))['done']
        except Exception: start = 0
    stats = dict(cores=0); found = []
    nf, nd3, core = dcheck_list(inst, dchk, P4, start, out, log, stopfile, base, seed, beta, t0, stats, found)
    log(dict(stage="D-check (resumed)", file=os.path.basename(path), start=start, checked=nf, d3_solutions=nd3,
             secs=round(time.time() - t0, 1)))
    for f in (path, base + '.json'):
        if os.path.exists(f): os.remove(f)
    return found

def run(inst, beta=1.5, nrow=24, K=24, logN=30.0, ncent=4, seed=0, tbits=20, out=None, log=print, max_dcheck=None,
        deadline=None, stopfile=None, predict=True, ckpt=None):
    t0 = time.time(); rng = np.random.default_rng(seed)
    rows = [block_candidates(inst, k, beta, nrow, K, rng) for k in range(4)]
    C = [cand_table(inst, rows[k]) for k in range(4)]
    log(dict(stage="candidates", per_block=[len(c[0]) for c in C], secs=round(time.time() - t0, 1)))
    T = inst.cy.T; tS = inst.tsigS(); found = []
    D = tS.shape[0]
    wts = rng.integers(1, 2 ** 63, D, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    rho = rng.integers(1, 1 << 30, D, dtype=np.int64)
    stats = dict(s_solutions=0, passed_filter=0, d3_solutions=0, cores=0, pairs_log2=[]); dchk = DChecker(inst, seed)
    for ci in range(ncent):
        if deadline is not None and time.time() > deadline: break
        if stopfile is not None and os.path.exists(stopfile): break
        best = choose_middle(inst, C, rng)
        if best is None:
            log(dict(stage="middle", note="no middle value with both sides non-empty")); continue
        sc, c, L, R, wl, wr = best
        # keep the heaviest parent pairs until each side has about 2^logN configuration pairs
        def trim(Pp, la, lb):
            Pp = sorted(Pp, key=lambda ij: -(la[ij[0]] + lb[ij[1]])); keep = []; tot = 0.0
            for i, j in Pp:
                keep.append((i, j)); tot += 2 ** (la[i] + lb[j])
                if math.log2(tot) >= logN: break
            return keep, math.log2(tot)
        Lk, wL = trim(L, C[0][2], C[1][2]); Rk, wR = trim(R, C[2][2], C[3][2])
        # materialise the fibres used
        need = [sorted({p[0] for p in Lk}), sorted({p[1] for p in Lk}), sorted({p[0] for p in Rk}), sorted({p[1] for p in Rk})]
        P = next_prime(max(2, 2 ** wL / ((1 << tbits) * 0.6)))
        sides = []; fid = []
        for k in range(4):
            fibres = []
            for ci_ in need[k]:
                pidx, S = materialise(inst, k, rows[k], C[k][0][ci_], K)
                fibres.append((pidx, seq_sig(S, T, inst.nl))); del S
            sides.append(Side(fibres, wts, rho, P)); fid.append({x: i for i, x in enumerate(need[k])})
        lpA = np.array([fid[0][i] for i, j in Lk], np.int64); lpB = np.array([fid[1][j] for i, j in Lk], np.int64)
        rpC = np.array([fid[2][i] for i, j in Rk], np.int64); rpD = np.array([fid[3][j] for i, j in Rk], np.int64)
        kT = (tS.astype(np.int64).astype(np.uint64) * wts).sum(); rT = int(((tS % P) * (rho % P) % P).sum() % P)
        A, B, Cc, Dd = sides
        # Gaussian prediction for this join: random left pairs (parent pair chosen by its share of pairs) + right pairs
        from q523_gauss import gauss_bits
        def draw(Pp, SA, SB, ka, kb, n):
            sz = np.array([(SA.off[i + 1] - SA.off[i]) * (SB.off[j + 1] - SB.off[j]) for i, j in Pp], float)
            jj = rng.choice(len(Pp), size=n, p=sz / sz.sum())
            ia = np.array([rng.integers(SA.off[Pp[j][0]], SA.off[Pp[j][0] + 1]) for j in jj])
            ib = np.array([rng.integers(SB.off[Pp[j][1]], SB.off[Pp[j][1] + 1]) for j in jj])
            return seq_sig(S_of(inst, ka, SA.pidx_all[ia]), T, inst.nl) + seq_sig(S_of(inst, kb, SB.pidx_all[ib]), T, inst.nl)
        if predict:
            Lp = list(zip(lpA, lpB)); Rp = list(zip(rpC, rpD))
            V = draw(Lp, A, B, 0, 1, 20000) + draw(Rp, Cc, Dd, 2, 3, 20000)
            g = gauss_bits(V, tS)
            pred = wL + wR - g['bits']
            log(dict(stage="prediction", centre=ci, fine_bits=g['bits'], mahal=g.get('mahal'), predicted_s_solutions_log2=round(pred, 2)))
        t1 = time.time()
        hits, nh, ov, pr = mjoin(A.key, A.res, A.off, B.key, B.rst, B.off, Cc.key, Cc.res, Cc.off, Dd.key, Dd.rst, Dd.off,
                                 lpA, lpB, rpC, rpD, P, rT, kT, tbits, 1 << 22, 0, P)
        # exact re-check
        sols = []
        for a, b, cc, dd in hits:
            p4 = [A.pidx_all[a].astype(np.int64), B.pidx_all[b].astype(np.int64), Cc.pidx_all[cc].astype(np.int64), Dd.pidx_all[dd].astype(np.int64)]
            S4 = [S_of(inst, k, p4[k])[0] for k in range(4)]
            tot = sum(seq_sig(x[None], T, inst.nl)[0] for x in S4)
            if (tot == tS).all(): sols.append((p4, S4))
        stats["s_solutions"] += len(sols); stats["pairs_log2"].append(round(math.log2(max(1, pr)), 2))
        log(dict(stage="S-join", centre=ci, parents=(len(Lk), len(Rk)), side_log2=(round(wL, 2), round(wR, 2)), P=P,
                 hash_hits=int(nh), overflow=int(ov), probes_log2=round(math.log2(max(1, pr)), 2), s_solutions=len(sols),
                 secs=round(time.time() - t1, 1)))
        t2 = time.time()
        use = sols if max_dcheck is None else sols[:max_dcheck]
        P4 = np.array([[p4[k].astype(np.uint8) for k in range(4)] for p4, S4 in use], np.uint8).reshape(len(use), 4, inst.F2)
        cbase = None
        if ckpt and len(use):
            cbase = f"{ckpt}_c{ci}"; np.save(cbase + '.npy', P4)
            with open(cbase + '.json', 'w') as fh: json.dump(dict(done=0, n=len(P4)), fh)
        nf, nd3, core = dcheck_list(inst, dchk, P4, 0, out, log, stopfile, cbase, seed, beta, t0, stats, found)
        if cbase:
            for f in (cbase + '.npy', cbase + '.json'):
                if os.path.exists(f): os.remove(f)
        if core: return found, stats
        stats["passed_filter"] += nf; stats["d3_solutions"] += nd3
        log(dict(stage="D-check", centre=ci, checked=nf, d3_solutions=nd3, secs=round(time.time() - t2, 1),
                 total_s=stats["s_solutions"], total_secs=round(time.time() - t0, 1)))
    return found, stats

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--beta": "1.5", "--rows": "24", "--K": "24", "--logN": "30", "--ncent": "4", "--seed": "0",
                             "--table-bits": "20", "--out": "", "--max-dcheck": ""}
    pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y)
    if '--planted' in flags:
        Pt, AC0 = planted_target(inst, np.random.default_rng(int(opt["--seed"]) + 999)); inst = SDInst(F, q, s, Y, T=Pt)
        print(json.dumps(dict(planted=True)), flush=True)
    log = lambda d: print(json.dumps(d), flush=True)
    found, st = run(inst, float(opt["--beta"]), int(opt["--rows"]), int(opt["--K"]), float(opt["--logN"]), int(opt["--ncent"]),
                    int(opt["--seed"]), int(opt["--table-bits"]), opt["--out"] or None, log,
                    int(opt["--max-dcheck"]) if opt["--max-dcheck"] else None)
    log(dict(stage="done", **{k: v for k, v in st.items()}))
