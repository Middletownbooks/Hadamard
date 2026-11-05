#!/usr/bin/env python3
"""Reproduce short controls and independently check all small matrices.
Construction imports are confined to assembly; the verifier module imports none.
Run the compression and coverage commands in the report first.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from bt_border import assemble
from openai_4p3_border import certificate, normalize, structure
from openai_4p3_build import EXAMPLE_SHIFTS, build, example_params, cyclotomic
from openai_4p3_compression import cyclotomic_sums
from openai_4p3_verify import check_core, check_matrix


def compression_check(z,p,q):
    # Direct independent quotient correlations, no normal-form formula.
    totals = [sum(row[u]*row[u ^ d if q == "v4" else (u+d) % 4]
                  for row in z for u in range(4)) for d in range(4)]
    assert totals == [4*p+12,4*p,4*p,4*p],(z,p,q,totals)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir",default="logs/openai_4p3")
    args = ap.parse_args()
    base = Path(args.dir)
    start = time.monotonic()
    host = json.loads((base/"compression_host.json").read_text())
    py = json.loads((base/"compression_python.json").read_text())
    hm = {r["p"]:r for r in host["records"]}
    for r in py["records"]:
        for key in ("min_total","min_per_cell"):
            assert r[key] == hm[r["p"]][key],(r["p"],key)
        for q in ("v4","z4"):
            compression_check(r["total_witness"],r["p"],q)
            compression_check(r["cell_witness"],r["p"],q)
            for g in r["cyclotomic"].values():
                if "witness" in g:
                    compression_check(g["witness"],r["p"],q)
    coverage = json.loads((base/"coverage.json").read_text())
    for r in coverage["records"]:
        assert r["status"] == "pass_compression_only"
        for q in ("v4","z4"):
            core = build(r["p"],q,"elliptic",r["parameters"])
            compression_check(core["compressions"],r["p"],q)
    # Independently enumerate all class signs and zero/one-position repairs in
    # small cyclotomic cells; compare the actually built sums with the gate.
    cyclotomic_controls = []
    for p in (5,7,13):
        for e in (2,4,6):
            if (p-1) % e:
                continue
            actual = set()
            for mask in range(1 << e):
                signs = [1-2*((mask >> j) & 1) for j in range(e)]
                for eta in (-1,1):
                    for flips in [[]]+[[k] for k in range(p)]:
                        cell = dict(e=e,signs=signs,zero=eta,scale=2,shift=1,flips=flips)
                        xs = cyclotomic(p,dict(cells=[[cell]*4 for _ in range(4)]))
                        actual.add(sum(xs[0][:p]))
            assert actual == set(cyclotomic_sums(p,e,1)),(p,e)
            cyclotomic_controls.append([p,e,len(actual)])
    certs = {q:certificate(q) for q in ("v4","z4")}
    records = []
    for (p,q) in EXAMPLE_SHIFTS:
        core = build(p,q,"elliptic",example_params(p,q))
        cv = check_core(core)
        assert cv["ok"]
        normal = normalize(core)
        assert check_core(normal)["ok"]
        cert = certs[q]
        H = assemble(structure(p,q),normal["sequences"],
                     *[np.array(cert[k],dtype=np.int64) for k in ("X","Rh","Ch")]).tolist()
        hv = check_matrix(H)
        assert hv["ok"]
        (base/f"core_p{p}_{q}.json").write_text(json.dumps(core,indent=2)+"\n")
        if (p,q) == (5,"v4"):
            (base/"matrix_p5_v4.json").write_text(json.dumps(dict(H=H,p=p,q=q),indent=2)+"\n")
        # A deliberately corrupted core and matrix must each fail.
        badx = dict(core,sequences=[row.copy() for row in core["sequences"]])
        badx["sequences"][0][0] *= -1
        assert not check_core(badx)["ok"]
        H[0][0] *= -1
        assert not check_matrix(H)["ok"]
        records.append(dict(p=p,q=q,core=cv,matrix=hv,negative_controls_rejected=True))
    result = dict(command=" ".join(sys.argv),seed=None,records=records,
                  compression_host_python_agree=True,all_compression_witnesses_checked=True,
                  cyclotomic_builder_controls=cyclotomic_controls,
                  target_cubic_compressions_checked=32,
                  seconds=round(time.monotonic()-start,3))
    (base/"checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(ok=True,matrices=[(r["p"],r["q"],r["matrix"]["n"]) for r in records],
                          seconds=result["seconds"])))


if __name__ == "__main__":
    main()
