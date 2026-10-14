r"""Exact evaluation of list restrictions for the 647 x353 S-join, on a complete S-solution set [M].

Given the complete S-solutions of a joint solution (theory/p647_calib.py), any restriction R of the four
distinct-sigma lists (keep a subset per block) has, exactly,
    work(R)  = the energy-compatible emitted pairs of the join on the kept lists (theory/p647_join.pair_counts,
               the same count join4p visits), for the pairing of the full run;
    S(R)     = the S-solutions of the full lists all of whose sigma are kept (the S-solutions of the restricted
               join are exactly these);
    Y(R)     = sum over S(R) of prod_k mult(sigma_k) = the D-candidates, i.e. pair-level S-solutions.
The cores expected from R are Y(R) 2^-C(D|S) (the D-check success per D-candidate, calibrated separately), so the cost
per core of a rule is proportional to work(R) / Y(R).  Rules: per-block energy caps e <= tau_k, multiplicity floors
mult >= t_k, and the score rule mult(sigma) >= t * 2^(beta (e - e0_k)) (more multiplicity required at higher energy).
usage: python3 theory/p647_restrict.py J"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from p647_join import SL, pair_counts, hash_w, res_w
from nd_cross import orbit_maps

def load_ssols(j, d=None):
    d = d or os.path.join(PC.LOG, "calib")
    npy = os.path.join(d, f"j{j}_ssols.npy")
    if os.path.exists(npy): return np.load(npy)
    S = []
    for l in open(os.path.join(d, f"j{j}_els.jsonl")):
        S.extend(json.loads(l)['sols'])
    return np.array(S, np.int64).reshape(-1, 4)

def work_of(blocks, keep, order, E=645, P_=61):
    """exact energy-compatible left + right pair counts for the kept lists (residues irrelevant: P_ small)"""
    w = hash_w(20, 11); u = res_w(20, P_, 11)
    Ls = [SL(B.Ssig[k_], E, w, u, P_) for B, k_ in zip(blocks, keep)]
    a, b, c, d = order
    cL = pair_counts(Ls[a], Ls[b], E, P_).sum(1); cR = pair_counts(Ls[c], Ls[d], E, P_).sum(1)
    both = (cL > 0) & (cR[::-1] > 0)
    return float(cL[both].sum() + cR[::-1][both].sum())

def evaluate(blocks, S, keepmask, order):
    keep = [np.nonzero(m)[0] for m in keepmask]
    ok = np.ones(len(S), bool)
    for k in range(4): ok &= keepmask[k][S[:, k]]
    Y = float(np.prod([blocks[k].mult[S[ok, k]].astype(np.float64) for k in range(4)], axis=0).sum()) if ok.any() else 0.0
    W = work_of(blocks, keep, order)
    return dict(kept=[int(m.sum()) for m in keepmask], S=int(ok.sum()), Y=Y, W=W, W_over_Y=(W / Y if Y else math.inf))

def main(j, log=print):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    blocks = P.build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
    S = load_ssols(j)
    order = None
    for l in open(os.path.join(PC.LOG, "calib", "calib.jsonl")):
        r = json.loads(l)
        if r['j'] == j: order = r['join_order']
    order = order or [0, 1, 2, 3]
    full = evaluate(blocks, S, [np.ones(B.N, bool) for B in blocks], order)
    log(json.dumps(dict(j=j, rule="full", **full)))
    E = [B.Ssig[:, 0] for B in blocks]; Mu = [B.mult for B in blocks]
    # energies and multiplicities of the S-solutions, per block
    for k in range(4):
        e = E[k][S[:, k]]; m = Mu[k][S[:, k]]
        log(json.dumps(dict(block=k, S_energy_pct=np.percentile(e, [1, 10, 50, 90, 99]).astype(int).tolist(),
                            list_energy_pct=np.percentile(E[k], [1, 10, 50, 90, 99]).astype(int).tolist(),
                            S_mult_mean=round(float(m.mean()), 2), list_mult_mean=round(float(Mu[k].mean()), 2))))
    best = (full['W_over_Y'], 'full')
    # multiplicity floors (same t for all blocks)
    for t in [2, 3, 4, 6, 8, 12, 16]:
        r = evaluate(blocks, S, [Mu[k] >= t for k in range(4)], order)
        log(json.dumps(dict(rule=f"mult>={t}", gain=round(full['W_over_Y'] / r['W_over_Y'], 3), Yfrac=round(r['Y'] / full['Y'], 4), **r)))
        if r['W_over_Y'] < best[0]: best = (r['W_over_Y'], f"mult>={t}")
    # energy caps relative to each block's pair-weighted mean energy
    pm = [float((E[k] * Mu[k]).sum() / Mu[k].sum()) for k in range(4)]
    for dlt in [-12, 0, 12, 24, 36, 48]:
        r = evaluate(blocks, S, [E[k] <= pm[k] + dlt for k in range(4)], order)
        log(json.dumps(dict(rule=f"e<=pm+{dlt}", gain=round(full['W_over_Y'] / r['W_over_Y'], 3), Yfrac=round(r['Y'] / full['Y'], 4), **r)))
        if r['W_over_Y'] < best[0]: best = (r['W_over_Y'], f"e<=pm+{dlt}")
    # score rule
    for t in [1, 2, 4]:
        for beta in [0.02, 0.04, 0.06, 0.08]:
            km = [Mu[k] >= t * np.exp2(beta * (E[k] - pm[k])) for k in range(4)]
            r = evaluate(blocks, S, km, order)
            log(json.dumps(dict(rule=f"score t={t} beta={beta}", gain=round(full['W_over_Y'] / r['W_over_Y'], 3), Yfrac=round(r['Y'] / full['Y'], 4), **r)))
            if r['W_over_Y'] < best[0]: best = (r['W_over_Y'], f"score t={t} beta={beta}")
    log(json.dumps(dict(best_rule=best[1], gain=round(full['W_over_Y'] / best[0], 3))))
    # energy-class selection on top of a rule: rank the left-energy classes by exact yield / exact work
    for rule, km in [("full", [np.ones(B.N, bool) for B in blocks]),
                     (best[1], None)]:
        if km is None:
            if best[1].startswith("mult>="):
                t = int(best[1][6:]); km = [Mu[k] >= t for k in range(4)]
            elif best[1].startswith("e<=pm+"):
                dl = int(best[1][6:]); km = [E[k] <= pm[k] + dl for k in range(4)]
            elif best[1].startswith("score"):
                t = float(best[1].split()[1][2:]); beta = float(best[1].split()[2][5:])
                km = [Mu[k] >= t * np.exp2(beta * (E[k] - pm[k])) for k in range(4)]
            else:
                continue
        el_curve(blocks, S, km, order, rule, log)

def el_curve(blocks, S, km, order, rule, log, E0=645, P_=61):
    """per left-energy class: exact work and exact yield under keep masks km; cumulative curve by yield/work"""
    w = hash_w(20, 11); u = res_w(20, P_, 11)
    keep = [np.nonzero(m)[0] for m in km]
    Ls = [SL(B.Ssig[k_], E0, w, u, P_) for B, k_ in zip(blocks, keep)]
    a, b, c, d = order
    cL = pair_counts(Ls[a], Ls[b], E0, P_).sum(1); cR = pair_counts(Ls[c], Ls[d], E0, P_).sum(1)[::-1]
    both = (cL > 0) & (cR > 0); work = np.where(both, cL + cR, 0).astype(np.float64)
    ok = np.ones(len(S), bool)
    for k in range(4): ok &= km[k][S[:, k]]
    Sk = S[ok]
    EL = blocks[a].Ssig[Sk[:, a], 0] + blocks[b].Ssig[Sk[:, b], 0]
    y = np.prod([blocks[k].mult[Sk[:, k]].astype(np.float64) for k in range(4)], axis=0)
    yEL = np.bincount(EL, weights=y, minlength=E0 + 1)[:E0 + 1]
    nz = np.nonzero(work)[0]
    ratio = np.where(work > 0, yEL / np.maximum(work, 1), 0)
    o = nz[np.argsort(-ratio[nz])]
    cw = np.cumsum(work[o]); cy = np.cumsum(yEL[o]); W = work.sum(); Y = yEL.sum()
    pts = []
    for f in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]:
        i = int(np.searchsorted(cw, f * W)); i = min(i, len(o) - 1)
        pts.append(dict(work_frac=round(float(cw[i] / W), 3), yield_frac=round(float(cy[i] / Y), 4),
                        gain=round(float((cy[i] / cw[i]) / (Y / W)), 3)))
    log(json.dumps(dict(rule=rule, el_classes=len(nz), W=W, Y=Y, curve=pts,
                        top_EL=[int(x) for x in o[:8]])))

if __name__ == "__main__":
    main(int(sys.argv[1]), log=lambda s: print(s, flush=True))
