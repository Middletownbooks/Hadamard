#!/usr/bin/env python3
"""Small exact border controls, using Claude Code's bt_border/bt_engine tooling.

Credit: EKLO Theorem 7 in Kumar's form. This wrapper normalizes the compression
and uses a universal finite certificate for (0,2 delta_0,2 delta_0,2 delta_0).
Matrix verification belongs to the separate openai_4p3_verify.py, never to this
construction script. See the report for the proof that block constants vanish.
"""
import argparse
from collections import Counter
import json
import sys
from pathlib import Path
import numpy as np
from bt_engine import BStruct
from bt_border import assemble, border_for_compressions, gs_array, QData


def structure(p,q):
    return BStruct([2,2,p],[(0,0,1)],name="openai_4p3") if q == "v4" else BStruct([4,p],[(0,1)],name="openai_4p3")


def normalize(core):
    p,q,xs = core["p"],core["q"],core["sequences"]
    constant,spikes = [],[]
    for i,x in enumerate(xs):
        z = [sum(x[j*p:(j+1)*p]) for j in range(4)]
        c = Counter(z)
        if len(c) == 1:
            constant.append((i,0,1))
        elif sorted(c.values()) == [1,3]:
            b = next(v for v,n in c.items() if n == 3)
            single = next(v for v,n in c.items() if n == 1)
            if abs(single-b) != 2:
                raise ValueError("Compression is not in the proved normal form")
            spikes.append((i,z.index(single),(single-b)//2))
        else:
            raise ValueError("Compression is not in the proved normal form")
    if len(constant) != 1 or len(spikes) != 3:
        raise ValueError("Need one constant row and three singleton rows")
    normalized = []
    operations = constant+spikes
    for i,s,sign in operations:
        normalized.append([sign*xs[i][(j ^ s if q == "v4" else (j+s) % 4)*p+k]
                           for j in range(4) for k in range(p)])
    z = [[sum(x[j*p:(j+1)*p]) for j in range(4)] for x in normalized]
    return dict(p=p,q=q,sequences=normalized,compressions=z,normalization=operations)


def certificate(q):
    st = structure(3,q)
    z = np.array([[0,0,0,0],[2,0,0,0],[2,0,0,0],[2,0,0,0]],dtype=np.int64)
    found = border_for_compressions(st,z,max_frames=76,tmax=2.0)
    if found is None:
        raise RuntimeError("Expected normal-form certificate was not found")
    X,Rh,Ch,frame = found
    Q = QData(st)
    Mc = gs_array(z,Q.add,Q.neg)
    assert np.all(X@X.T == 12*np.eye(12,dtype=np.int64))
    assert np.all(Rh@Rh.T == 16*np.eye(12,dtype=np.int64))
    assert all(np.all(Rh[:,4*j:4*j+4].sum(axis=1) == 0) for j in range(4))
    assert np.all(Ch@Ch.T == np.kron(np.eye(4,dtype=np.int64),Q.Gmat))
    assert np.all(X@Ch.T+Rh@Mc.T == 0)
    assert np.all(np.abs(X) == 1) and np.all(np.abs(Rh) == 1) and np.all(np.abs(Ch) == 1)
    return dict(q=q,X=X.tolist(),Rh=Rh.tolist(),Ch=Ch.tolist(),Mc=Mc.tolist(),
                frame_index_zero_based=int(frame),integer_identities_pass=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("core",nargs="?")
    ap.add_argument("--certificates",action="store_true")
    ap.add_argument("--max-order",type=int,default=512)
    ap.add_argument("--out",required=True)
    args = ap.parse_args()
    if args.certificates:
        result = dict(certificates=[certificate(q) for q in ("v4","z4")])
    else:
        if not args.core:
            ap.error("Supply a core JSON, or --certificates")
        core = normalize(json.loads(Path(args.core).read_text()))
        p,q = core["p"],core["q"]
        if 16*p+12 > args.max_order:
            ap.error("Order exceeds --max-order; this wrapper defaults to small controls")
        cert = certificate(q)
        H = assemble(structure(p,q),core["sequences"],
                     *[np.array(cert[k],dtype=np.int64) for k in ("X","Rh","Ch")])
        result = dict(p=p,q=q,normalized_core=core,H=H.tolist(),
                      status="assembled; independent integer matrix verification required")
    result.update(command=" ".join(sys.argv),seed=None)
    out = Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(output=str(out),mode="certificates" if args.certificates else "matrix")))


if __name__ == "__main__":
    main()
