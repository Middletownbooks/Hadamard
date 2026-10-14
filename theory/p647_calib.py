r"""Calibration of the 647 x353 cross-level Gaussian against exact counts [M vs E].

For a joint solution j of y_B #0 (theory/p647_common.py):
  1. Gaussian predictions (theory/q523_gauss.py, lattice-exact):
       C_S^u  = S-level bits for sigma drawn uniformly from the distinct sigma of each block   -> exact S-join count
       C_S^p  = S-level bits for (L0, L1) pairs drawn uniformly (sigma weighted by multiplicity) -> sum over
                S-solutions of the D-candidates (prod_k mult(sigma_k))
       C      = S and D together, pairs uniform -> expected cores 2^(H - C)
  2. the exact S-join (theory/p647_join.join4p), restart-safe per left energy EL (each finished EL appends a line to
     logs/p647/calib/j{j}_els.jsonl with its pair counts, time and S-solutions);
  3. the exact D-check of every S-solution (theory/p647_cross.dcheck) -> cores (each verified by the full PAF);
  4. the D-tail: over every D-candidate of the S-solutions, the exact number hitting the D-target on subsets of the
     cross lags, against the Gaussian fitted to a sample of the same candidates' residuals (as theory/q523_dtail_check).
Output: one JSON line per joint solution in logs/p647/calib/calib.jsonl.
usage: python3 theory/p647_calib.py J [J ...] [--N 20000] [--tail-frac 1.0]"""
import sys, os, json, math, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from p647_join import join4p
from nd_cross import orbit_maps
from q523_gauss import gauss_bits

CAL = os.path.join(PC.LOG, "calib")

def gauss_preds(X, blocks, N=20000, seed=1):
    rng = np.random.default_rng(seed)
    TE = 2 * (X.TS + X.TD)
    VSu = 0; VSp = 0; VF = 0
    for B in blocks:
        su = rng.integers(B.N, size=N)
        VSu = VSu + B.Ssig[su]
        i = rng.integers(B.n0, size=N); j = rng.integers(B.n1, size=N)
        codes = ((B.A0[i].astype(np.int64) + 1) // 2) @ B.W + ((B.A1[j].astype(np.int64) + 1) // 2) @ B.W
        sp = np.searchsorted(B.U, codes)
        VSp = VSp + B.Ssig[sp]
        VF = VF + np.concatenate([B.Ssig[sp], B.P0[i] + B.P1[j]], 1)
    gu = gauss_bits(VSu, X.TS); gp = gauss_bits(VSp, X.TS)
    gF = gauss_bits(VF, np.concatenate([X.TS, TE]))
    Hs = sum(math.log2(B.N) for B in blocks); H = sum(math.log2(B.n0 * B.n1) for B in blocks)
    return dict(Hs=round(Hs, 3), H=round(H, 3), CSu=gu.get('bits'), CSu_rank=gu.get('rank'), CSu_mahal=gu.get('mahal'),
                CSp=gp.get('bits'), CSp_rank=gp.get('rank'), C=gF.get('bits'), C_rank=gF.get('rank'), C_mahal=gF.get('mahal'),
                pred_S_log2=round(Hs - gu['bits'], 3), pred_Spairs_log2=round(H - gp['bits'], 3),
                pred_cores_log2=round(H - gF['bits'], 3), CDgS=round(gF['bits'] - gp['bits'], 3))

@njit(cache=True)
def _dtail(sols, pr0, bd0, n10, pr1, bd1, n11, pr2, bd2, n12, pr3, bd3, n13,
           P00, P10, P01, P11, P02, P12, P03, P13, TE, subs, nsub, cnt, smp, every, nsmp_max):
    """enumerate every D-candidate of every S-solution; cnt[q] += 1 when the residual vanishes on subset q (subs[q, :]
    is a 0/1 mask over the lags); every 'every'-th candidate's residual goes to smp"""
    L = TE.shape[0]; ns = 0; tot = 0
    r0 = np.empty(L, np.int64); r1 = np.empty(L, np.int64); r2 = np.empty(L, np.int64); r3 = np.empty(L, np.int64)
    for s in range(sols.shape[0]):
        for x0 in range(bd0[sols[s, 0]], bd0[sols[s, 0] + 1]):
            p = pr0[x0]
            for t in range(L): r0[t] = P00[p // n10, t] + P10[p % n10, t] - TE[t]
            for x1 in range(bd1[sols[s, 1]], bd1[sols[s, 1] + 1]):
                p = pr1[x1]
                for t in range(L): r1[t] = r0[t] + P01[p // n11, t] + P11[p % n11, t]
                for x2 in range(bd2[sols[s, 2]], bd2[sols[s, 2] + 1]):
                    p = pr2[x2]
                    for t in range(L): r2[t] = r1[t] + P02[p // n12, t] + P12[p % n12, t]
                    for x3 in range(bd3[sols[s, 3]], bd3[sols[s, 3] + 1]):
                        p = pr3[x3]
                        for t in range(L): r3[t] = r2[t] + P03[p // n13, t] + P13[p % n13, t]
                        for q in range(nsub):
                            z = True
                            for t in range(L):
                                if subs[q, t] and r3[t] != 0:
                                    z = False; break
                            if z: cnt[q] += 1
                        if tot % every == 0 and ns < nsmp_max:
                            for t in range(L): smp[ns, t] = r3[t]
                            ns += 1
                        tot += 1
    return tot, ns

def dtail(X, blocks, sols, subsets, nsmp=200000):
    sols = np.ascontiguousarray(np.asarray(sols, np.int64).reshape(-1, 4))
    B = blocks; TE = 2 * (X.TS + X.TD)
    cand = (B[0].mult[sols[:, 0]] * B[1].mult[sols[:, 1]] * B[2].mult[sols[:, 2]] * B[3].mult[sols[:, 3]]).astype(np.int64)
    total = int(cand.sum()); every = max(1, total // nsmp)
    subs = np.zeros((len(subsets), len(TE)), np.bool_)
    for q, sb in enumerate(subsets): subs[q, sb] = True
    cnt = np.zeros(len(subsets), np.int64); smp = np.zeros((nsmp + 1, len(TE)), np.int64)
    tot, ns = _dtail(sols, B[0].pairs, B[0].bd, B[0].n1, B[1].pairs, B[1].bd, B[1].n1, B[2].pairs, B[2].bd, B[2].n1,
                     B[3].pairs, B[3].bd, B[3].n1, B[0].P0, B[0].P1, B[1].P0, B[1].P1, B[2].P0, B[2].P1, B[3].P0, B[3].P1,
                     TE, subs, len(subsets), cnt, smp, every, nsmp + 1)
    assert tot == total
    V = smp[:ns]
    out = []
    for q, sb in enumerate(subsets):
        g = gauss_bits(V[:, sb], np.zeros(len(sb), np.int64))
        exact = (math.log2(total / cnt[q]) if cnt[q] else None)
        out.append(dict(lags=[int(X.lags[t]) for t in sb], hits=int(cnt[q]), exact_bits=round(exact, 2) if exact else None,
                        gauss_bits=g.get('bits'), rank=g.get('rank')))
    gall = gauss_bits(V, np.zeros(len(TE), np.int64))
    return dict(candidates=total, candidates_log2=round(math.log2(total), 3) if total else None, sample=int(ns),
                tail=out, gauss_full_bits=gall.get('bits'), gauss_full_rank=gall.get('rank'), gauss_full_mahal=gall.get('mahal'))

# cross lags (a != 0 mod 17 and b != 0 mod 19) in X.lags order
def cross_lag_positions(X):
    return [t for t, h in enumerate(X.lags) if h % 17 and h % 19]

def calibrate(j, X, ys, sols, maps, NH, oi, N=20000, log=print, tail=True):
    os.makedirs(CAL, exist_ok=True)
    yb = ys[0]; ya = X.J.ya_full(sols[j]['layers'])
    t0 = time.time()
    blocks = P.build(X, yb, ya, maps, NH, oi)
    rec = dict(j=j, mult_log2=sols[j]['mult_log2'], N=[B.N for B in blocks], n=[(B.n0, B.n1) for B in blocks],
               build_secs=round(time.time() - t0, 1))
    rec.update(gauss_preds(X, blocks, N=N))
    log(json.dumps(dict(stage="pred", **rec)))
    # exact S-join, restart-safe per EL
    elp = os.path.join(CAL, f"j{j}_els.jsonl")
    done = {}; Ssols = []
    if os.path.exists(elp):
        for line in open(elp):
            try:
                r = json.loads(line); done[r['EL']] = r; Ssols.extend(map(tuple, r['sols']))
            except Exception: pass
    S = [B.Ssig for B in blocks]
    _, st0 = join4p(S, X.TS, stats_only=True)
    todo = [EL for EL in range(int(X.TS[0]) + 1) if EL not in done]
    info = {}
    fo = open(elp, "a")
    def el_cb(EL, d): info[EL] = d
    def done_cb(EL, found):
        d = info.get(EL, dict(left=0, right=0, key_matches=0, nsol=0, secs=0.0))
        fo.write(json.dumps(dict(EL=EL, **d, sols_list_len=len(found), sols=[list(t) for t in found])) + "\n"); fo.flush()
    t1 = time.time()
    newsols, st = join4p(S, X.TS, order=st0['order'], els=todo, el_cb=el_cb, done_cb=done_cb)
    fo.close()
    Ssols = sorted(set(Ssols) | set(newsols))
    rec.update(S_solutions=len(Ssols), S_log2=round(math.log2(len(Ssols)), 3) if Ssols else None,
               join_order=st0['order'], work_log2=st0['work_log2'], join_secs_this_run=round(time.time() - t1, 1),
               join_ns_per_pair=st.get('ns_per_pair'), resumed_els=len(done))
    # D-check
    t2 = time.time()
    TE = 2 * (X.TS + X.TD)
    hits, cand = P.dcheck(blocks, Ssols, TE) if Ssols else (np.zeros((0, 5), np.int64), np.zeros(0, np.int64))
    cores = []
    for h in hits:
        xs = P.assemble(X, blocks, h[1:])
        ok = P.is_core(xs)
        cores.append(dict(sol=int(h[0]), pairs=[int(v) for v in h[1:]], verified=ok))
        if ok:
            path = os.path.join(PC.LOG, f"CORE_calib_j{j}_{len(cores)}.json")
            json.dump(dict(j=j, yb=[list(map(int, y)) for y in yb], ya=ya, core=xs.tolist()), open(path, "w"))
            log(json.dumps(dict(stage="CORE", j=j, path=path)))
    rec.update(S_pairs=int(cand.sum()), S_pairs_log2=round(math.log2(cand.sum()), 3) if len(cand) else None,
               D_hits=len(hits), cores=sum(c['verified'] for c in cores), dcheck_secs=round(time.time() - t2, 1))
    np.save(os.path.join(CAL, f"j{j}_ssols.npy"), np.asarray(Ssols, np.int64).reshape(-1, 4))
    log(json.dumps(dict(stage="exact", **{k: rec[k] for k in ('j', 'S_solutions', 'S_log2', 'pred_S_log2', 'S_pairs_log2',
                                                             'pred_Spairs_log2', 'D_hits', 'cores', 'pred_cores_log2')})))
    if tail and Ssols:
        cl = cross_lag_positions(X)
        subsets = [[t] for t in cl[:4]] + [cl[:2], cl[2:4], cl[:3], cl[:4], cl[4:8], cl[:6]]
        t3 = time.time()
        rec['dtail'] = dtail(X, blocks, Ssols, subsets)
        rec['dtail_secs'] = round(time.time() - t3, 1)
    rec['secs'] = round(time.time() - t0, 1)
    with open(os.path.join(CAL, "calib.jsonl"), "a") as f: f.write(json.dumps(rec) + "\n")
    log(json.dumps(dict(stage="done", j=j, secs=rec['secs'])))
    return rec

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "20000"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    for j in pos:
        calibrate(j, X, ys, sols, maps, NH, oi, N=int(opt["--N"]), log=lambda s: print(s, flush=True))

# ---------------------------------------------------------------- energy-stratified Gaussian (per left-energy class)
def strat_preds(X, blocks, order, N=3000, seed=3, weighted=False, min_frac=1e-4):
    """per left-energy class EL: sample N left pairs (a, b) with e_a + e_b = EL and N right pairs (c, d) with
    e_c + e_d = T0 - EL, exactly in proportion to the class (uniform over distinct sigma, or proportional to the
    multiplicity if weighted); fit the lattice Gaussian to the summed signature without its energy coordinate.
    Predicted count(EL) = (pairs or multiplicity mass of the left class) x (same, right) x 2^-bits.
    Classes holding less than min_frac of the total left x right mass are skipped (their mass is reported)."""
    rng = np.random.default_rng(seed)
    T = X.TS; E0 = int(T[0]); a, b, c, d = order
    def buckets(B):
        e = B.Ssig[:, 0]; wgt = B.mult.astype(np.float64) if weighted else np.ones(B.N)
        return e, wgt, np.bincount(e, weights=wgt, minlength=E0 + 1)[:E0 + 1]
    bk = [buckets(B) for B in blocks]
    def sample_side(i, j, Etot, n):
        (ei, wi, mi), (ej, wj, mj) = bk[i], bk[j]
        ea = np.arange(E0 + 1); eb = Etot - ea
        ok = (eb >= 0) & (eb <= E0)
        mass = np.where(ok, mi * np.where(ok, mj[np.clip(eb, 0, E0)], 0), 0)
        if mass.sum() == 0: return None, 0.0
        pe = rng.choice(E0 + 1, size=n, p=mass / mass.sum())
        out = np.zeros((n, len(T)), np.int64)
        for e in np.unique(pe):
            k_ = int((pe == e).sum())
            ia = np.nonzero(ei == e)[0]; ib = np.nonzero(ej == Etot - e)[0]
            sa = rng.choice(ia, size=k_, p=wi[ia] / wi[ia].sum()); sb = rng.choice(ib, size=k_, p=wj[ib] / wj[ib].sum())
            out[pe == e] = blocks[i].Ssig[sa] + blocks[j].Ssig[sb]
        return out, float(mass.sum())
    # class masses
    tot = []
    for EL in range(E0 + 1):
        (ea_, wa_, ma), (eb_, wb_, mb) = bk[a], bk[b]
        (ec_, wc_, mc), (ed_, wd_, md) = bk[c], bk[d]
        mL = float(np.convolve(ma, mb)[EL]) if EL <= 2 * E0 else 0.0
        mR = float(np.convolve(mc, md)[E0 - EL])
        tot.append((EL, mL, mR))
    allmass = sum(mL * mR for _, mL, mR in tot)
    res = []; skipped = 0.0
    for EL, mL, mR in tot:
        if mL * mR == 0: continue
        if mL * mR < min_frac * allmass: skipped += mL * mR; continue
        Lv, _ = sample_side(a, b, EL, N); Rv, _ = sample_side(c, d, E0 - EL, N)
        V = Lv + Rv
        g = gauss_bits(V[:, 1:], T[1:])
        res.append(dict(EL=EL, mass_log2=round(math.log2(mL * mR), 3), bits=g.get('bits'), rank=g.get('rank'),
                        mahal=g.get('mahal'), pred_log2=round(math.log2(mL * mR) - g['bits'], 3) if g.get('in_lattice') and math.isfinite(g['bits']) else None))
    pred = sum(2.0 ** r['pred_log2'] for r in res if r['pred_log2'] is not None)
    return dict(weighted=weighted, classes=res, pred_total_log2=round(math.log2(pred), 3) if pred else None,
                skipped_mass_frac=round(skipped / allmass, 6))
