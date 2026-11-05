#!/usr/bin/env python3
"""Exact compression gates for the k4 family. See report/openai_4p3_answer.md.
No assertion here is an existence assertion for full cores.
The embedded HOST_JS was executed in this session's V8 runtime; the independent
Python implementation below is provided for reproduction and cross-checking.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

TARGETS = [149,179,257,271,331,359,367,389,457,461,467,509,521,523,571,599]

def normal_forms(p):
    # Each nonconstant row has three t's and one t+2, t odd.
    bound = math.isqrt(p + 3) + 2
    odds = [x for x in range(-bound, bound + 1) if x % 2]
    for a in odds:
        if a <= 0 or a*a > p:
            continue
        for b,c,d in itertools.product(odds, repeat=3):
            if a*a + b*(b+1) + c*(c+1) + d*(d+1) == p:
                yield [[a]*4] + [[t+2,t,t,t] for t in (b,c,d)]

def balanced_minima(p):
    best_total = best_cell = None
    tw = cw = None
    for z in normal_forms(p):
        costs = [(abs(x)-1)//2 for row in z for x in row]
        total, cell = sum(costs), max(costs)
        if best_total is None or total < best_total:
            best_total, tw = total, z
        if best_cell is None or cell < best_cell:
            best_cell, cw = cell, z
    if tw is None:
        raise RuntimeError("No compression representation found")
    return dict(p=p, N=16*p+12, min_total=best_total, min_per_cell=best_cell,
                total_witness=tw, cell_witness=cw)

def cyclotomic_sums(p, e, budget):
    if (p-1) % e:
        raise ValueError("e must divide p-1")
    s = (p-1)//e
    out = {}
    for negatives in range(e+1):
        for zero in (-1,1):
            base = zero + s*(e-2*negatives)
            plus, minus = (p+base)//2, (p-base)//2
            for up in range(min(budget,minus)+1):
                for down in range(min(budget-up,plus)+1):
                    value, cost = base+2*(up-down), up+down
                    out[value] = min(cost, out.get(value, cost))
    return out

def cyclotomic_gate(p, e, budget):
    sums = cyclotomic_sums(p, e, budget)
    for z in normal_forms(p):
        if all(v in sums for row in z for v in row):
            return dict(status="pass_compression_only", witness=z)
    return dict(status="fail_compression")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", default=",".join(map(str,TARGETS)))
    ap.add_argument("--indices", default="2,4,6,8,12,16")
    ap.add_argument("--max-c", type=int, default=8)
    ap.add_argument("--out")
    args = ap.parse_args()
    records = []
    for p in map(int,args.primes.split(",")):
        rec = balanced_minima(p)
        rec["cyclotomic"] = {}
        forms = list(normal_forms(p))
        for e in map(int,args.indices.split(",")):
            if (p-1) % e:
                rec["cyclotomic"][str(e)] = dict(status="not_defined")
                continue
            for c in range(args.max_c+1):
                sums = cyclotomic_sums(p,e,c)
                witness = next((z for z in forms if all(v in sums for row in z for v in row)),None)
                if witness is not None:
                    rec["cyclotomic"][str(e)] = dict(min_per_cell=c,witness=witness,
                                                       status="pass_compression_only")
                    break
            else:
                rec["cyclotomic"][str(e)] = dict(status="fails_through_budget",budget=args.max_c)
        records.append(rec)
    result = dict(command="python3 theory/openai_4p3_compression.py",
                  seed=None, arithmetic="Python integers", records=records)
    text = json.dumps(result,indent=2)+"\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True,exist_ok=True)
        Path(args.out).write_text(text)
    print(text,end="")

HOST_JS = r"""

function gcd(a,b){while(b){[a,b]=[b,a%b]}return a;}
function cost(z){return Math.max(0,(Math.abs(z)-1)/2);}
const ps=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,149,179,257,271,331,359,367,389,457,461,467,509,521,523,571,599];
const out=[];
for (const p of ps){
 const lim=Math.floor(Math.sqrt(p+3))+2;const odds=[];for(let x=-lim;x<=lim;x++)if(Math.abs(x)%2)odds.push(x);
 let best=Infinity,witness=null,cbest=Infinity,cwitness=null;
 for(const a of odds)if(a>0&&a*a<=p)for(const b of odds)for(const c of odds)for(const d of odds){
  if(a*a+b*(b+1)+c*(c+1)+d*(d+1)!==p)continue;
  const z=[[a,a,a,a],...[b,c,d].map(x=>[x+2,x,x,x])];
  const F=z.flat().reduce((s,x)=>s+cost(x),0),C=Math.max(...z.flat().map(cost));
  if(F<best){best=F;witness=z;}if(C<cbest){cbest=C;cwitness=z;}
 }
 out.push({p,N:16*p+12,min_total:best,min_per_cell:cbest,total_witness:witness,cell_witness:cwitness});
}
return out;
"""

if __name__ == "__main__":
    main()
