"""Stage B: screen the twisted b = 1 cyclic slices of the open orders (report/prompt_b1_twist_design.md §3).

For each v, every cyclic subgroup U = <u> of (Z/m)* (m = v - 1) and every twist ('BC' if d even, 'BCD' if 3 | d;
'id' with --id for reference) with RMIN <= orbits <= RMAX:
  1. exact row-sum feasibility (DP);
  2. exact quotient screens [P] (scripts/b1tw_lag0.py): for every subgroup H of order K, the all-lags quotient test when
     its enumeration is small, else the lag-0 energy test (energy only when the row-sum table would be large);
  3. for survivors, the Gaussian log2 E and conditional constraint bits (scripts/b1tw_lib.py gaussian).
Appends JSON lines to OUT, skipping (m, u, tw) already there (restart-safe).
usage: b1tw_screen.py OUT v [v ...] [--rmin R] [--rmax R] [--id] [--ns N]"""
import sys, os, json, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice, gaussian, cyclic_subgroups, log2_rowsum_configs
from b1tw_lag0 import quotient_screen, divisors, cells_of_group, letters


def energy_only(sl, K):
    """lag-0 energy test without row sums (1-D DP): can sum_groups mult * sum_cells w y^2 equal T0?"""
    F = sl.m // K; T0 = 4 * sl.m - 4 * (K - 1)
    acc = np.zeros(T0 + 1, bool); acc[0] = True
    for g in sl.groups:
        rep, mult, cells = cells_of_group(sl, g, F)
        for w, ns in cells:
            es = sorted({mult * w * y * y for y in letters(ns)})
            nxt = np.zeros_like(acc)
            for e in es:
                if e <= T0: nxt[e:] |= acc[:T0 + 1 - e]
            acc = nxt
    return bool(acc[T0])


def quotient_ok(sl, K):
    """all-lags test if small enough, else lag-0 energy; returns (verdict, which)"""
    q = quotient_screen(sl, K, maxvec=300_000)
    if q is not None: return q, "all-lags"
    return energy_only(sl, K), "lag0-energy"


def main():
    args = sys.argv[1:]
    OUT = args[0]
    def opt(f, d, ty=int):
        return ty(args[args.index(f) + 1]) if f in args else d
    rmin, rmax, ns = opt("--rmin", 8), opt("--rmax", 150), opt("--ns", 20000)
    tws = ("id", "BC", "BCD") if "--id" in args else ("BC", "BCD")
    vs = []
    for a in args[1:]:
        if a.startswith("--"): break
        vs.append(int(a))
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            z = json.loads(l); done.add((z["m"], z["u"], z["tw"]))
    out = open(OUT, "a")
    for v in vs:
        m = v - 1
        for d, u, U in cyclic_subgroups(m):
            if d == 1: continue
            for tw in tws:
                if (tw == "BC" and d % 2) or (tw == "BCD" and d % 3) or (m, u, tw) in done: continue
                sl = Slice(m, u, tw)
                if not (rmin <= sl.r <= rmax): continue
                t0 = time.time()
                rec = dict(v=v, m=m, u=u, d=d, tw=tw, orbits=sl.r, classes=sl.L, minus1_in_U=sl.minus1_in_U,
                           group_bits={g: len(b) for g, b in sl.groups.items()})
                lrs = log2_rowsum_configs(sl); rec["log2_rowsum_configs"] = round(lrs, 2)
                if lrs == float("-inf"):
                    rec.update(alive=False, killed="row sums")
                else:
                    killed = None; tested = []
                    for K in sorted([k for k in divisors(m) if 1 < k < m], reverse=True):   # small quotients first
                        ok, how = quotient_ok(sl, K); tested.append([K, how])
                        if not ok:
                            killed = dict(K=K, F=m // K, test=how); break
                    rec.update(alive=killed is None, killed=killed, tested=tested)
                    if killed is None:
                        ga = gaussian(sl, N=ns)
                        rec.update({k: ga.get(k) for k in ("feasible", "why", "n_cond", "log2E", "steps", "mahal")})
                rec["secs"] = round(time.time() - t0, 1)
                out.write(json.dumps(rec) + "\n"); out.flush()
                print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
