#!/usr/bin/env python3
"""Bounded exact tests of the cubic-root-switch proposal (no heuristic estimates).

Seed gate enumerates every lambda, twist, root completion, and nonzero dilation.
It retains distinct PAF vectors, using full tuples (no probabilistic hash test).
A zero kills this seed class for BOTH quotients. A positive returns one seed
quadruple only. Phase failure is therefore scoped to that quadruple.
Coverage mode tests compression only, not correlations. Caps mean UNKNOWN.
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path
from openai_4p3_build import build, elliptic_seed, prime, quadratic_table
from openai_4p3_compression import TARGETS
from openai_4p3_verify import check_core


class Cap(Exception):
    pass


class Budget:
    def __init__(self, seconds, pairs, vectors):
        self.deadline = time.monotonic()+seconds
        self.pairs, self.vectors = pairs, vectors

    def check(self):
        if time.monotonic() >= self.deadline:
            raise Cap("time_limit")


def paf(a):
    p = len(a)
    return tuple(sum(a[k]*a[(k+t) % p] for k in range(p)) for t in range(1,(p+1)//2))


def seed_gate(p, budget):
    bases, defects = {}, {}
    for lam in range(2,p):
        budget.check()
        for sign,mask,scale in itertools.product((-1,1),range(8),range(1,p)):
            budget.check()
            par = dict(lam=lam,sign=sign,mask=mask,scale=scale)
            a = elliptic_seed(p,par)
            b = paf(a)
            bases.setdefault(b,par)
            for r in range(p):
                if a[r] != -1 or r*scale % p not in (0,1,lam):
                    continue
                d = a.copy()
                d[r] = 1
                v = tuple((3*x+y)//4 for x,y in zip(b,paf(d)))
                defects.setdefault(v,dict(par,root=r))
            if len(bases)+len(defects) > budget.vectors:
                raise Cap("vector_limit")
    bs, ds = list(bases.items()), list(defects.items())
    pairs = len(ds)*(len(ds)+1)//2
    out = dict(p=p,base_vectors=len(bs),defect_vectors=len(ds),pairs=pairs)
    if pairs > budget.pairs or len(bs)*len(ds) > budget.pairs:
        return dict(out,status="capped",reason="pair_limit")
    table = {}
    for j,(a,_) in enumerate(ds):
        budget.check()
        for k in range(j,len(ds)):
            table.setdefault(tuple(x+y for x,y in zip(a,ds[k][0])),(j,k))
    for a,pa in bs:
        budget.check()
        for b,pb in ds:
            key = tuple(-3-x-y for x,y in zip(a,b))
            if key in table:
                j,k = table[key]
                return dict(out,status="pass_seed_gate",parameters=[pa,pb,ds[j][1],ds[k][1]])
    return dict(out,status="no_seed_quadruple")


def phase_candidates(p,q,par,budget):
    b = elliptic_seed(p,par)
    special = b.copy()
    if "root" in par:
        special[par["root"]] = 1
    candidates = {}
    for h1,h2,h3 in itertools.product(range(p),repeat=3):
        budget.check()
        hs = [0,h1,h2,h3]
        a = [special]+[[b[(k-h) % p] for k in range(p)] for h in hs[1:]]
        vals = []
        for d in range(1,4):
            for t in range(p):
                vals.append(sum(a[u][k]*a[u ^ d if q == "v4" else (u+d) % 4][(k+t) % p]
                                for u in range(4) for k in range(p)))
        candidates.setdefault(tuple(vals),hs)
        if len(candidates) > budget.vectors:
            raise Cap("vector_limit")
    return list(candidates.items())


def phase_join(rec,q,budget):
    p,pars = rec["p"],rec["parameters"]
    cs = [phase_candidates(p,q,par,budget) for par in pars]
    out = dict(p=p,q=q,lists=[len(c) for c in cs],scope="first seed-gate witness only")
    if len(cs[0])*len(cs[1]) > budget.pairs or len(cs[2])*len(cs[3]) > budget.pairs:
        return dict(out,status="capped",reason="pair_limit")
    pairs = {}
    for a,ha in cs[0]:
        budget.check()
        for b,hb in cs[1]:
            pairs.setdefault(tuple(x+y for x,y in zip(a,b)),(ha,hb))
    for c,hc in cs[2]:
        budget.check()
        for d,hd in cs[3]:
            key = tuple(4-x-y for x,y in zip(c,d))
            if key not in pairs:
                continue
            ha,hb = pairs[key]
            params = [dict(par,shifts=h) for par,h in zip(pars,(ha,hb,hc,hd))]
            core = build(p,q,"elliptic",params)
            verified = check_core(core)
            if not verified["ok"]:
                raise RuntimeError("Independent PAF check rejected phase join")
            return dict(out,status="core",parameters=params,core=core,verification=verified)
    return dict(out,status="no_shifts_for_selected_seeds")


def trace_gate(p,budget):
    L = quadratic_table(p)
    traces = {}
    for lam in range(2,p):
        budget.check()
        T = sum(L[u*(u-1)*(u-lam) % p] for u in range(p))
        for sign in (-1,1):
            traces.setdefault(sign*T,dict(lam=lam,sign=sign))
    const,spike = {},{}
    for T,par in traces.items():
        for mask in range(8):
            rs = [1-2*((mask >> j) & 1) for j in range(3)]
            z = T+sum(rs)
            params = dict(par,mask=mask,scale=1)
            const.setdefault(z,params)
            if mask:
                root = [0,1,par["lam"]][rs.index(-1)]
                spike.setdefault(z,dict(params,root=root))
    ds = sorted(spike)
    pairs = {}
    for b,c in itertools.product(ds,repeat=2):
        pairs.setdefault(b*(b+1)+c*(c+1),(b,c))
    for a,b in itertools.product(sorted(const),ds):
        key = p-a*a-b*(b+1)
        if key in pairs:
            c,d = pairs[key]
            return dict(p=p,N=16*p+12,status="pass_compression_only",traces=len(traces),
                        witness=[a,b,c,d],parameters=[const[a],spike[b],spike[c],spike[d]],
                        full_core_status="unknown")
    return dict(p=p,N=16*p+12,status="fail_compression",traces=len(traces),full_core_status="fails")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coverage",action="store_true")
    ap.add_argument("--primes",default="3,5,7,11,13,17,19")
    ap.add_argument("--phase-max-p",type=int,default=13)
    ap.add_argument("--max-pairs",type=int,default=2000000)
    ap.add_argument("--max-vectors",type=int,default=20000)
    ap.add_argument("--seconds",type=float,default=540)
    ap.add_argument("--out",required=True)
    ap.add_argument("--artifact-dir",default="logs/openai_4p3")
    args = ap.parse_args()
    if not 0 < args.seconds <= 600:
        ap.error("Each invocation must have a time budget of at most 600 seconds")
    ps = TARGETS if args.coverage else list(map(int,args.primes.split(",")))
    if not all(prime(p) for p in ps):
        ap.error("All p must be odd primes")
    budget = Budget(args.seconds,args.max_pairs,args.max_vectors)
    result = dict(command=" ".join(sys.argv),seed=None,arithmetic="Python integers",
                  mode="compression_only" if args.coverage else "seed_then_selected_phases",records=[])
    out,adir = Path(args.out),Path(args.artifact_dir)
    out.parent.mkdir(parents=True,exist_ok=True)
    adir.mkdir(parents=True,exist_ok=True)
    start = time.monotonic()
    for p in ps:
        try:
            rec = trace_gate(p,budget) if args.coverage else seed_gate(p,budget)
            if rec["status"] == "pass_seed_gate" and p <= args.phase_max_p:
                rec["phases"] = []
                for q in ("v4","z4"):
                    phase = phase_join(rec,q,budget)
                    if "core" in phase:
                        path = adir / f"core_p{p}_{q}.json"
                        path.write_text(json.dumps(phase.pop("core"),indent=2)+"\n")
                        phase["core_file"] = str(path)
                    rec["phases"].append(phase)
        except Cap as ex:
            # Capped means undecided; never turn an incomplete enumeration into a zero.
            rec = dict(p=p,status="capped",reason=str(ex))
        result["records"].append(rec)
        result["seconds"] = round(time.monotonic()-start,3)
        out.write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps({k:v for k,v in rec.items() if k not in ("parameters","phases")}),flush=True)
        if rec.get("reason") == "time_limit":
            result["unfinished_primes"] = ps[len(result["records"]):]
            out.write_text(json.dumps(result,indent=2)+"\n")
            break


if __name__ == "__main__":
    main()
