r"""Same-shape positive controls for the 647 x353 pipeline [M, with E predictions alongside].

Slice: G = Z_m, m = 2 p q (p = qA, q = qB odd primes), M = <mu> with mu mod p generating a proper subgroup M_p of
Z_p^*, mu mod q generating a proper subgroup M_q of Z_q^*, gcd(|M_p|, |M_q|) = 1, so M = M_p x M_q is a product and
fixes Z_2: the shape of 647 x353 (p = 17, |M_p| = 4; q = 19, |M_q| = 3).  The pipeline is the 647 one, unchanged:
  1. stage 1B: every y_B up to the slice's symmetry group (theory/nd_joint_exhaust.enumerate_yB, restart-safe);
  2. joint stage 1 per y_B orbit (JointTest.test): every y_A consistent with y_B solving stage 1A;
  3. cross level per joint solution (theory/p647_cross.py: layer lists, distinct sigma, exact S-join
     theory/p647_join.join4p, D-check, full cyclic PAF check), with the Gaussian predictions of
     theory/p647_calib.gauss_preds next to every exact count;
  4. completion of any core by bordered Goethals-Seidel (theory/nd_doubly.bordered_gs_cyclic), H H^T = N I.
usage: python3 theory/p647_control.py m mu qA qB [--maxyb N] [--maxjoint N] [--stop-at-core] [--dir DIR]
       python3 theory/p647_control.py --list VMAX      (the product-shape cyclic slices with v = 2pq + 1 <= VMAX)"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross, orbit_maps
from nd_joint_exhaust import enumerate_yB
from nd_doubly import bordered_gs_cyclic
import p647_cross as PX
from p647_calib import gauss_preds

def order(a, p):
    k = 1; x = a % p
    while x != 1: x = x * a % p; k += 1
    return k

def product_slices(vmax):
    """(v, p, q, |M_p|, |M_q|, mu) for the product-shape cyclic slices, smallest mu per (p, q, |M_p|, |M_q|)"""
    def isp(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))
    out = []
    P = [p for p in range(5, vmax) if isp(p)]
    for p in P:
        for q in P:
            if q == p or 2 * p * q + 1 > vmax: continue
            m = 2 * p * q
            for dp in [d for d in range(2, p - 1) if (p - 1) % d == 0]:
                for dq in [d for d in range(2, q - 1) if (q - 1) % d == 0]:
                    if math.gcd(dp, dq) != 1: continue
                    mu = next(u for u in range(3, m, 2) if math.gcd(u, m) == 1 and order(u, p) == dp and order(u, q) == dq)
                    out.append((m + 1, p, q, dp, dq, mu))
    return sorted(out)

def run(m, mu, qA, qB, d, maxyb=10 ** 9, maxjoint=10 ** 9, stop_at_core=False, log=print, workers=1):
    os.makedirs(d, exist_ok=True)
    tag = f"{m}_{mu}_{qA}_{qB}"
    X = Cross(m, mu, qA, qB)
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), workers, log=log)
    log(json.dumps(dict(stage="yB", m=m, mu=mu, qA=qA, qB=qB, orbits=len(ys))))
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    outp = os.path.join(d, f"control_{tag}.jsonl")
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            r = json.loads(line)
            if 'jsol' in r: done.add((r['yb'], r['jsol']))
    tot = dict(joint=0, pred_cores=0.0, cores=0, S=0, pred_S=0.0, Spairs=0, pred_Spairs=0.0)
    found = None
    with open(outp, "a") as fo:
        for ib, yb in enumerate(ys[:maxyb]):
            t0 = time.time()
            rec, sols = X.J.test(yb)
            fo.write(json.dumps(dict(yb=ib, joint=rec.get('joint'), S1=rec.get('S_solutions'), secs=round(time.time() - t0, 2))) + "\n"); fo.flush()
            for js, s_ in enumerate(sols[:maxjoint]):
                if (ib, js) in done: continue
                ya = X.J.ya_full(s_['layers'])
                blocks = PX.build(X, yb, ya, maps, NH, oi)
                if blocks is None: continue
                t1 = time.time()
                g = gauss_preds(X, blocks, N=20000)
                r, Ssols, hits, cand = PX.solve(X, blocks)
                cores = []
                for h in hits:
                    xs = PX.assemble(X, blocks, h[1:])
                    if PX.is_core(xs): cores.append(xs)
                out = dict(yb=ib, jsol=js, H=g['H'], Hs=g['Hs'], CSu=g['CSu'], CSp=g['CSp'], C=g['C'], C_rank=g['C_rank'],
                           pred_S_log2=g['pred_S_log2'], S=r['S_solutions'], pred_Spairs_log2=g['pred_Spairs_log2'],
                           Spairs=r['S_pairs'], pred_cores_log2=g['pred_cores_log2'], D_hits=r['D_hits'], cores=len(cores),
                           secs=round(time.time() - t1, 2))
                fo.write(json.dumps(out) + "\n"); fo.flush()
                tot['joint'] += 1; tot['cores'] += len(cores); tot['S'] += r['S_solutions']; tot['Spairs'] += r['S_pairs']
                for key, k2 in (('pred_cores', 'pred_cores_log2'), ('pred_S', 'pred_S_log2'), ('pred_Spairs', 'pred_Spairs_log2')):
                    if g[k2] is not None and math.isfinite(g[k2]): tot[key] += 2.0 ** g[k2]
                if cores and found is None:
                    xs = cores[0]
                    Hm = bordered_gs_cyclic(xs)
                    ok = Hm is not None and np.array_equal(Hm @ Hm.T, Hm.shape[0] * np.eye(Hm.shape[0], dtype=np.int64))
                    found = dict(v=m + 1, N=4 * (m + 1), yb=ib, jsol=js, hadamard=bool(ok))
                    json.dump(dict(**found, mu=mu, qA=qA, qB=qB, core=xs.tolist()), open(os.path.join(d, f"core_{tag}.json"), "w"))
                    if Hm is not None: np.save(os.path.join(d, f"H{4 * (m + 1)}_{tag}.npy"), Hm.astype(np.int8))
                    log(json.dumps(dict(stage="CORE", **found)))
                    if stop_at_core: break
            log(json.dumps(dict(stage="yb_done", yb=ib, joint=rec.get('joint'), **{k: (round(v, 2) if isinstance(v, float) else v) for k, v in tot.items()})))
            if found and stop_at_core: break
    summ = dict(stage="summary", m=m, mu=mu, qA=qA, qB=qB, yB_orbits=len(ys), found=found,
                **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in tot.items()})
    log(json.dumps(summ))
    return summ

if __name__ == "__main__":
    a = sys.argv[1:]
    if "--list" in a:
        for r in product_slices(int(a[a.index("--list") + 1])): print(*r)
        sys.exit()
    opt = {"--maxyb": str(10 ** 9), "--maxjoint": str(10 ** 9), "--dir": "logs/p647/control", "--workers": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith("--"): i += 1
        else: pos.append(int(a[i])); i += 1
    m, mu, qA, qB = pos
    run(m, mu, qA, qB, opt["--dir"], int(opt["--maxyb"]), int(opt["--maxjoint"]), "--stop-at-core" in a,
        log=lambda s: print(s, flush=True), workers=int(opt["--workers"]))
