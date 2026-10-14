r"""Re-derive the H(2588) core from scratch: one joint solution, one left-energy class [M].

Rebuilds the cross level of joint solution J of y_B #K (theory/p647_cross.build: per-layer lists, distinct sigma,
S-signatures, layer PAFs), runs the exact S-join (theory/p647_join.join4p) for the single left-energy class EL, D-checks
every S-solution (theory/p647_cross.dcheck), assembles every D-hit on Z_646, keeps those whose full cyclic PAF sum is
(2584, -4, ..., -4), and compares them with the recorded core (logs/p647/CORE_*.json), up to the block order and
sign normalisation of theory/p647_finish.py.
usage: python3 theory/p647_core_redo.py [--yb 0] [--j 1508] [--el 317] [--core logs/p647/CORE_c2_j1508_0.json]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from p647_join import join4p
from nd_cross import orbit_maps

def norm(X):
    X = np.asarray(X, np.int64)
    k2 = int(np.nonzero(np.abs(X.sum(1)) == 2)[0][0]); X = X[[k2] + [k for k in range(4) if k != k2]]
    return -X if X[0].sum() < 0 else X

def main(k=0, j=1508, EL=317, core=None, log=print):
    t0 = time.time()
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, k)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    blocks = P.build(X, ys[k], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
    Ssig = [B.Ssig for B in blocks]
    _, st0 = join4p(Ssig, X.TS, stats_only=True)
    found = []
    def done_cb(e, s_): found.extend(s_)
    join4p(Ssig, X.TS, order=st0['order'], els=[EL], done_cb=done_cb)
    S = np.array([list(t) for t in found], np.int64)
    TE = 2 * (X.TS + X.TD)
    hits, cand = P.dcheck(blocks, S, TE)
    cores = [xs for xs in (P.assemble(X, blocks, h[1:]) for h in hits) if P.is_core(xs)]
    out = dict(yb=k, j=j, EL=EL, S_solutions=len(S), D_candidates=int(cand.sum()), D_hits=len(hits), cores=len(cores),
               secs=round(time.time() - t0, 1))
    if core:
        ref = norm(json.load(open(core))['core'])
        out['matches_recorded_core'] = any(np.array_equal(norm(xs), ref) for xs in cores)
    log(json.dumps(out))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    g = lambda f, d: a[a.index(f) + 1] if f in a else d
    main(int(g("--yb", 0)), int(g("--j", 1508)), int(g("--el", 317)),
         g("--core", os.path.join(PC.LOG, "CORE_c2_j1508_0.json")), log=lambda s: print(s, flush=True))
