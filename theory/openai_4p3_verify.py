#!/usr/bin/env python3
"""Independent integer verification. Imports NO construction/search/border code.

Accepts core JSON with p,q,sequences, or matrix JSON with H. Every product and
sum below uses Python integers. A successful matrix check checks all N^2 entries
of H H^T, including the diagonal. No FFT or rounded floating point is used.
"""
import argparse
import json
import sys
from pathlib import Path


def check_core(obj):
    p, q, xs = obj["p"], obj["q"], obj["sequences"]
    if type(p) is not int or p < 3 or q not in ("v4","z4"):
        raise ValueError("Invalid p or q")
    if len(xs) != 4 or any(len(x) != 4*p or any(type(v) is not int or v not in (-1,1) for v in x) for x in xs):
        raise ValueError("Invalid core dimensions/entries")
    totals, bad = [], []
    for d in range(4):
        for t in range(p):
            value = 0
            for x in xs:
                for u in range(4):
                    v = u ^ d if q == "v4" else (u+d) % 4
                    for k in range(p):
                        value += x[u*p+k]*x[v*p+(k+t) % p]
            wanted = (16*p if t == 0 else -12) if d == 0 else 4
            totals.append(value)
            if value != wanted:
                bad.append(dict(d=d,t=t,actual=value,wanted=wanted))
    return dict(ok=not bad,p=p,q=q,checked_lags=4*p,totals=totals,bad=bad)


def check_matrix(H):
    n = len(H)
    if not n or any(len(r) != n or any(type(v) is not int or v not in (-1,1) for v in r) for r in H):
        raise ValueError("Invalid square binary matrix")
    for i in range(n):
        for j in range(n):
            value = sum(H[i][k]*H[j][k] for k in range(n))
            if value != (n if i == j else 0):
                return dict(ok=False,n=n,i=i,j=j,dot=value)
    return dict(ok=True,n=n,checked_entries=n*n,arithmetic="Python integer row dot products")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("--out")
    args = ap.parse_args()
    obj = json.loads(Path(args.input).read_text())
    result = check_matrix(obj["H"]) if "H" in obj else check_core(obj)
    result.update(command=" ".join(sys.argv), seed=None)
    text = json.dumps(result,indent=2)+"\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True,exist_ok=True)
        Path(args.out).write_text(text)
    print(text,end="")
    if not result["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
