r"""Re-derive, export and verify the cores of a control slice found by theory/p647_ctrl_full.py [M].

For each (y_B index, joint solution) given: rebuild the cross level of that joint solution (theory/p647_cross),
solve it exactly (S-join, D-check), assemble every D-hit on Z_m and keep those whose full cyclic PAF sum is
(4m, -4, ..., -4); write the core JSON logs/p647/control/core_<tag>_yb<i>_j<j>_<n>.json if it is missing (a run killed
between its record and its core file leaves it missing), and the text file core_<tag>_yb<i>_j<j>_<n>.txt (A=, B=,
C=, D= lines; the block with row sum +-2 first, made +2), which scripts/verify_h2588.py checks independently
(no code shared with theory/; it builds H(4(m+1)) by bordered Goethals-Seidel and checks H H^T = N I exactly).
usage: python3 theory/p647_ctrl_core.py m mu qA qB yb:jsol [yb:jsol ...]"""
import sys, os, json, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross, orbit_maps
from nd_joint_exhaust import enumerate_yB
from nd_doubly import bordered_gs_cyclic
import p647_common as PC
import p647_cross as PX
from p647_captest import capped_test
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main(m, mu, qA, qB, pairs, log=print):
    d = os.path.join(PC.LOG, "control"); tag = f"{m}_{mu}_{qA}_{qB}"
    X = Cross(m, mu, qA, qB)
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), 1, log=lambda s: None)
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    cache = {}; out = []
    for ib, js in pairs:
        if ib not in cache:
            st, rec, sols = capped_test(X, ys[ib], maxsol=10 ** 7, maxpairs=4 * 10 ** 6, mem_gb=5.0); assert st == 'ok', rec
            cache[ib] = sols
        sols = cache[ib]
        blocks = PX.build(X, ys[ib], X.J.ya_full(sols[js]['layers']), maps, NH, oi)
        r, Ss, hits, cand = PX.solve(X, blocks, P=61)
        cores = [xs for xs in (PX.assemble(X, blocks, h[1:]) for h in hits) if PX.is_core(xs)]
        log(json.dumps(dict(yb=ib, jsol=js, S=r['S_solutions'], Spairs=r['S_pairs'], D_hits=r['D_hits'], cores=len(cores))))
        for n_, xs in enumerate(cores):
            jp = os.path.join(d, f"core_{tag}_yb{ib}_j{js}_{n_}.json")
            if not os.path.exists(jp):
                Hm = bordered_gs_cyclic(xs)
                ok = Hm is not None and np.array_equal(Hm @ Hm.T, Hm.shape[0] * np.eye(Hm.shape[0], dtype=np.int64))
                json.dump(dict(m=m, mu=mu, qA=qA, qB=qB, yb_index=ib, jsol=int(js), yb=[list(map(int, y)) for y in ys[ib]],
                               hadamard_internal=bool(ok), core=xs.tolist(), note="re-derived by theory/p647_ctrl_core.py"), open(jp, "w"))
            X4 = np.array(json.load(open(jp))['core'], np.int64)
            k2 = int(np.nonzero(np.abs(X4.sum(1)) == 2)[0][0]); X4 = X4[[k2] + [k for k in range(4) if k != k2]]
            if X4[0].sum() < 0: X4 = -X4
            tp = jp[:-5] + ".txt"
            with open(tp, "w") as f:
                f.write(f"# b = 1 core for H({4 * (m + 1)}): four +-1 sequences on cyclic Z_{m}, control slice {m + 1} x{mu} "
                        f"(y_B #{ib}, joint solution {js}); row sums {X4.sum(1).tolist()}\n")
                for name, x in zip("ABCD", X4): f.write(name + "=" + "".join("+" if v > 0 else "-" for v in x) + "\n")
            npz = os.path.join(d, f"H{4 * (m + 1)}_yb{ib}_j{js}_{n_}_packed.npz")
            res = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "verify_h2588.py"), tp, npz], capture_output=True, text=True)
            last = res.stdout.strip().splitlines()[-1] if res.stdout.strip() else res.stderr[-300:]
            v = json.loads(last) if last.startswith("{") else dict(hadamard=False, err=last)
            log(json.dumps(dict(stage="verified", core=os.path.relpath(tp, ROOT), independent_hadamard=v.get('hadamard'),
                                N=v.get('N'), sha256_of_packed_rows=v.get('sha256_of_packed_rows'), border=v.get('border'))))
            out.append(v)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]
    pr = [tuple(int(t) for t in x.split(":")) for x in a[4:]]
    main(*(int(x) for x in a[:4]), pr, log=lambda s: print(s, flush=True))
