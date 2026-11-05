#!/usr/bin/env python3
"""Reference builders for report/openai_4p3_answer.md. Outputs candidates, not certificates.

Coordinates are q*p+k, q=0,1,2,3; q addition is XOR (v4) or modulo 4 (z4).
Elliptic parameters: four rows with lam, sign, mask, scale, shifts[4]; rows 1..3
also have root, a position where the seed is -1 and the cubic vanishes.
Cyclotomic parameters: cells[4][4], each with e, signs[e], zero, shift, scale,
and flips (distinct positions in 0..p-1). All parameters are literal integers.
Credit: the completion is EKLO Theorem 7 in Kumar's form; border tooling is
Claude Code's theory/bt_border.py. This file only builds proposed core sequences.
"""
import argparse
import json
import math
from pathlib import Path


def prime(p):
    return p >= 3 and p % 2 and all(p % d for d in range(3, math.isqrt(p)+1, 2))


def quadratic_table(p):
    a = [-1]*p
    a[0] = 0
    for x in range(1, p):
        a[x*x % p] = 1
    return a


def elliptic_seed(p, par):
    lam, sign, mask, scale = (par[k] for k in ("lam", "sign", "mask", "scale"))
    if not (2 <= lam < p and sign in (-1, 1) and 0 <= mask < 8 and 1 <= scale < p):
        raise ValueError("Invalid cubic parameters")
    roots = [0, 1, lam]
    L = quadratic_table(p)
    out = []
    for k in range(p):
        u = k*scale % p
        out.append(1-2*((mask >> roots.index(u)) & 1) if u in roots
                   else sign*L[u*(u-1)*(u-lam) % p])
    return out


def elliptic(p, params):
    if len(params) != 4:
        raise ValueError("Exactly four rows required")
    xs = []
    for i, par in enumerate(params):
        b = elliptic_seed(p, par)
        hs = par.get("shifts", [0]*4)
        if len(hs) != 4 or any(not isinstance(h, int) for h in hs):
            raise ValueError("Four integer shifts required")
        special = b.copy()
        if i:
            r = par["root"]
            if not (0 <= r < p and b[r] == -1 and r*par["scale"] % p in (0, 1, par["lam"])):
                raise ValueError("The changed position must be a negative root")
            special[r] = 1
        xs.append([((special if q == 0 else b)[(k-hs[q]) % p])
                   for q in range(4) for k in range(p)])
    return xs


def primitive_root(p):
    factors = [d for d in range(2, p) if (p-1) % d == 0
               and all(d % j for j in range(2, math.isqrt(d)+1))]
    return next(g for g in range(2, p) if all(pow(g, (p-1)//d, p) != 1 for d in factors))


def cyclotomic(p, params):
    cells = params["cells"]
    if len(cells) != 4 or any(len(row) != 4 for row in cells):
        raise ValueError("cells must be 4 by 4")
    g = primitive_root(p)
    log = {}
    u = 1
    for j in range(p-1):
        log[u] = j
        u = u*g % p
    xs = []
    for row in cells:
        x = []
        for par in row:
            e, signs = par["e"], par["signs"]
            eta, h, scale = par.get("zero", 1), par.get("shift", 0), par.get("scale", 1)
            flips = par.get("flips", [])
            if not (e > 0 and (p-1) % e == 0 and len(signs) == e
                    and all(s in (-1, 1) for s in signs) and eta in (-1, 1) and 1 <= scale < p):
                raise ValueError("Invalid cyclotomic parameters")
            if len(set(flips)) != len(flips) or any(not 0 <= k < p for k in flips):
                raise ValueError("flips must be distinct positions")
            fs = set(flips)
            for k in range(p):
                u = scale*(k-h) % p
                value = eta if not u else signs[log[u] % e]
                x.append(-value if k in fs else value)
        xs.append(x)
    return xs


EXAMPLE_SEEDS = {
    3: [(2,-1,1,1), (2,-1,1,1,0), (2,-1,3,1,0), (2,-1,3,1,0)],
    5: [(2,-1,0,1), (2,-1,1,2,0), (2,1,3,1,0), (2,1,3,2,0)],
    7: [(2,-1,1,1), (2,-1,1,2,0), (2,-1,1,3,0), (2,-1,2,2,4)],
}
EXAMPLE_SHIFTS = {
    (3,"v4"): [[0,1,1,2],[0,0,1,2],[0,0,0,0],[0,0,1,2]],
    (3,"z4"): [[0,0,2,1],[0,0,1,2],[0,0,0,0],[0,0,1,2]],
    (5,"v4"): [[0,0,0,2],[0,0,4,3],[0,0,1,3],[0,0,0,1]],
    (5,"z4"): [[0,1,0,3],[0,0,0,1],[0,0,0,2],[0,0,0,4]],
    (7,"v4"): [[0,2,0,5],[0,0,1,3],[0,0,0,1],[0,0,6,5]],
}


def example_params(p, q):
    shifts = EXAMPLE_SHIFTS[(p, q)]
    return [dict(zip(("lam","sign","mask","scale","root"), seed), shifts=h)
            for seed,h in zip(EXAMPLE_SEEDS[p], shifts)]


def build(p, q, ansatz, params):
    if not prime(p) or q not in ("v4", "z4"):
        raise ValueError("An odd prime and q=v4 or z4 are required")
    xs = elliptic(p, params) if ansatz == "elliptic" else cyclotomic(p, params)
    return dict(p=p, q=q, ansatz=ansatz, parameters=params, sequences=xs,
                compressions=[[sum(x[j*p:(j+1)*p]) for j in range(4)] for x in xs],
                status="candidate; run independent verifier")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--p", type=int)
    ap.add_argument("--q", choices=("v4","z4"), default="v4")
    ap.add_argument("--ansatz", choices=("elliptic","cyclotomic"), default="elliptic")
    ap.add_argument("--params", help="JSON parameters (list for elliptic, cells object for cyclotomic)")
    ap.add_argument("--example", type=int, choices=(3,5,7))
    ap.add_argument("--out")
    args = ap.parse_args()
    if args.example:
        if (args.example,args.q) not in EXAMPLE_SHIFTS:
            ap.error("No stored example for that (p,q)")
        p, ansatz, params = args.example, "elliptic", example_params(args.example,args.q)
    else:
        if not args.p or not args.params:
            ap.error("Use --example, or --p and --params")
        p, ansatz, params = args.p, args.ansatz, json.loads(Path(args.params).read_text())
    text = json.dumps(build(p,args.q,ansatz,params),indent=2)+"\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True,exist_ok=True)
        Path(args.out).write_text(text)
    else:
        print(text,end="")


if __name__ == "__main__":
    main()
