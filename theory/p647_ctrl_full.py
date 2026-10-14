r"""Full positive-control run on a same-shape product slice (first: 495 x83, m = 494 = 2.13.19) [M, with E alongside].

Same pipeline as 647 x353, unchanged: every y_B orbit (theory/nd_joint_exhaust.enumerate_yB, restart-safe), the exact
joint test of each, the joint solutions reduced to orbits of the stabiliser of their y_B (theory/p647_common), and for
every orbit representative the exact cross level (theory/p647_cross: lists, distinct sigma, S-join p647_join.join4p,
D-check, full PAF check).  Next to every exact count the Gaussian predictions of theory/p647_calib.gauss_preds
(N = 4000): predicted S-solutions, D-candidates and cores.  Every core is completed by bordered Goethals-Seidel and
verified (theory/nd_doubly.bordered_gs_cyclic, then scripts/verify_h2588.py, which works for any such core).
Restart-safe: one JSON line per processed joint solution in logs/p647/control/full_<tag>.jsonl.
The Gaussian predictions are computed for every 10th joint solution (--gauss-every); the exact counts for all.
--classes (from 2026-10-14 10:50): dedupe by the full stabiliser of y_B, per-block reversal included
(theory/p647_classes.py; the old stabiliser misses it, so its orbits come in groups of copies, 8 for 495 x83); the
class maps go to logs/p647/control/classes_<tag>.json.  --maxpairs caps the joint test's pair lists per block (a y_B
over the cap is recorded as undecided); an 'attempt' line before each joint test makes a run killed inside one
(memory) retry that y_B once on restart and then skip it, recorded as undecided.  The joint test itself runs in a
child under an address-space cap (--mem GB, theory/p647_captest.py), so a y_B too big for memory is recorded as
undecided instead of killing the run.
usage: python3 theory/p647_ctrl_full.py m mu qA qB [--yb i,j,...] [--max-cores N] [--gauss-every K] [--classes] [--maxpairs N]
       [--mem GB]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross, orbit_maps
from nd_joint_exhaust import enumerate_yB
from nd_doubly import bordered_gs_cyclic
import p647_common as PC
import p647_cross as PX
from p647_calib import gauss_preds
from p647_captest import capped_test

def main(m, mu, qA, qB, yb_list=None, max_cores=10 ** 9, log=print, N=4000, gauss_every=10, full_stab=False, maxpairs=6 * 10 ** 7,
         mem_gb=6.0):
    d = os.path.join(PC.LOG, "control"); os.makedirs(d, exist_ok=True)
    tag = f"{m}_{mu}_{qA}_{qB}"; outp = os.path.join(d, f"full_{tag}.jsonl")
    X = Cross(m, mu, qA, qB); FA, FB = X.D.FA, X.D.FB
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), 1, log=lambda s: None)
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    done = set(); ncores = 0; attempted = {}
    if os.path.exists(outp):
        for line in open(outp):
            r = json.loads(line)
            if 'jsol' in r: done.add((r['yb'], r['jsol'])); ncores += r.get('cores', 0)
            if 'joint_total' in r: done.add((r['yb'], 'test'))
            if 'attempt' in r: attempted[r['yb']] = attempted.get(r['yb'], 0) + 1
    cpath = os.path.join(d, f"classes_{tag}.json")
    cmaps = json.load(open(cpath)) if os.path.exists(cpath) else {}
    log(json.dumps(dict(stage="start", m=m, mu=mu, qA=qA, qB=qB, yB_orbits=len(ys), done=len(done), cores_so_far=ncores)))
    order = range(len(ys)) if yb_list is None else yb_list
    T = dict(joint=0, S=0, pred_S=0.0, Spairs=0, pred_Spairs=0.0, pred_cores=0.0, pred_cores_cal=0.0, cores=0)
    with open(outp, "a") as fo:
        for ib in order:
            yb = ys[ib]
            if attempted.get(ib, 0) >= 2 and (ib, 'test') not in done:     # killed twice inside the test: give up
                fo.write(json.dumps(dict(yb=ib, joint_total=None, note="run killed inside the joint test (memory); undecided")) + "\n")
                fo.flush(); done.add((ib, 'test')); continue
            if (ib, 'test') not in done:
                fo.write(json.dumps(dict(yb=ib, attempt=1)) + "\n"); fo.flush()
            t0 = time.time(); status, rec, sols = capped_test(X, yb, maxsol=10 ** 7, maxpairs=maxpairs, mem_gb=mem_gb)
            if not sols:
                if (ib, 'test') not in done:
                    extra = {} if rec.get('joint', 0) == 0 else dict(joint_total=None, note=rec.get('note'), pairs_needed=rec.get('pairs_needed'))
                    fo.write(json.dumps(dict(dict(yb=ib, joint_total=0, orbits=0, secs=round(time.time() - t0, 2)), **extra)) + "\n"); fo.flush()
                continue
            if full_stab:
                import p647_classes as PCL
                stab = PCL.stabiliser_full(yb, FB, FA, m)
                cls, rep_idx = PCL.classes(X, sols, stab, FA)
                if str(ib) not in cmaps:
                    cmaps[str(ib)] = dict(class_of=cls, reps=rep_idx, stabiliser=len(stab))
                    json.dump(cmaps, open(cpath + ".tmp", "w")); os.replace(cpath + ".tmp", cpath)
            else:
                stab = PC.stabiliser_on(yb, m, FB, FA)
                keys, rep_idx = PC.canonical_joint(X, sols, stab, FA)
            if (ib, 'test') not in done:
                fo.write(json.dumps(dict(yb=ib, joint_total=len(sols), stabiliser=len(stab), orbits=len(rep_idx),
                                         secs=round(time.time() - t0, 2))) + "\n"); fo.flush()
            log(json.dumps(dict(stage="yb", yb=ib, joint=len(sols), stabiliser=len(stab), orbits=len(rep_idx))))
            for js in rep_idx:
                if (ib, js) in done: continue
                blocks = PX.build(X, yb, X.J.ya_full(sols[js]['layers']), maps, NH, oi)
                if blocks is None: continue
                t1 = time.time()
                nproc = sum(1 for k_ in done if isinstance(k_[1], int)) + T['joint']
                g = gauss_preds(X, blocks, N=N) if nproc % gauss_every == 0 else dict(
                    H=round(sum(math.log2(B.n0 * B.n1) for B in blocks), 3), Hs=round(sum(math.log2(B.N) for B in blocks), 3),
                    CSu=None, CSp=None, C=None, CDgS=None, pred_S_log2=None, pred_Spairs_log2=None, pred_cores_log2=None)
                r, Ssols, hits, cand = PX.solve(X, blocks, P=61)
                cores = []
                for h in hits:
                    xs = PX.assemble(X, blocks, h[1:])
                    if PX.is_core(xs): cores.append(xs)
                out = dict(yb=ib, jsol=int(js), H=g['H'], Hs=g['Hs'], CSu=g['CSu'], CSp=g['CSp'], C=g['C'], CDgS=g['CDgS'],
                           pred_S_log2=g['pred_S_log2'], S=r['S_solutions'], pred_Spairs_log2=g['pred_Spairs_log2'],
                           Spairs=r['S_pairs'], pred_cores_log2=g['pred_cores_log2'], D_hits=r['D_hits'], cores=len(cores),
                           secs=round(time.time() - t1, 2))
                fo.write(json.dumps(out) + "\n"); fo.flush()
                T['joint'] += 1; T['S'] += r['S_solutions']; T['Spairs'] += r['S_pairs']; T['cores'] += len(cores)
                for key, k2 in (('pred_S', 'pred_S_log2'), ('pred_Spairs', 'pred_Spairs_log2'), ('pred_cores', 'pred_cores_log2')) if g['C'] is not None else ():
                    if g[k2] is not None and math.isfinite(g[k2]): T[key] += 2.0 ** g[k2]
                if g['CDgS'] is not None and math.isfinite(g['CDgS']): T['pred_cores_cal'] += r['S_pairs'] * 2.0 ** -g['CDgS']
                for n_, xs in enumerate(cores):
                    ncores += 1
                    Hm = bordered_gs_cyclic(xs)
                    ok = Hm is not None and np.array_equal(Hm @ Hm.T, Hm.shape[0] * np.eye(Hm.shape[0], dtype=np.int64))
                    path = os.path.join(d, f"core_{tag}_yb{ib}_j{js}_{n_}.json")
                    json.dump(dict(m=m, mu=mu, qA=qA, qB=qB, yb_index=ib, jsol=int(js), yb=[list(map(int, y)) for y in yb],
                                   hadamard_internal=bool(ok), core=xs.tolist()), open(path, "w"))
                    log(json.dumps(dict(stage="CORE", v=m + 1, N=4 * (m + 1), yb=ib, jsol=int(js), hadamard_internal=bool(ok), path=path)))
                if T['joint'] % 25 == 0:
                    log(json.dumps(dict(stage="progress", **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in T.items()})))
                if ncores >= max_cores: break
            if ncores >= max_cores: break
    log(json.dumps(dict(stage="summary", **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in T.items()}, cores_total=ncores)))

if __name__ == "__main__":
    a = sys.argv[1:]
    pos = [int(x) for x in a[:4]]
    ybl = [int(x) for x in a[a.index("--yb") + 1].split(",")] if "--yb" in a else None
    mc = int(a[a.index("--max-cores") + 1]) if "--max-cores" in a else 10 ** 9
    ge = int(a[a.index("--gauss-every") + 1]) if "--gauss-every" in a else 10
    mp = int(a[a.index("--maxpairs") + 1]) if "--maxpairs" in a else 6 * 10 ** 7
    mg = float(a[a.index("--mem") + 1]) if "--mem" in a else 6.0
    main(*pos, yb_list=ybl, max_cores=mc, log=lambda s: print(s, flush=True), gauss_every=ge, full_stab="--classes" in a, maxpairs=mp,
         mem_gb=mg)
