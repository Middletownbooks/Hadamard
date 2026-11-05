#!/usr/bin/env python3
"""Exact structural tests of the 1916/3964 controls; no seed search.

For every Q-cell and every CRT translation, measure the distance from a quadratic
product on the nonzero rectangle (axes arbitrary), and from an arbitrary 3x3
quadratic-class sign table. Also test arbitrary balanced rank-one interiors.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from openai_identity_verify import paf_check, primitive_energy

ROOT = Path(__file__).resolve().parents[1]


def legendre(p):
    return np.array([0] + [1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1,p)], np.int64)


def crt_cells(rows, p, r):
    ix = np.array([[next(k for k in range(p*r) if k % p == a and k % r == b)
                    for b in range(r)] for a in range(p)])
    return np.asarray(rows, np.int64).reshape(4,4,p*r)[:,:,ix]


def cell_distances(cell, p, r):
    chp, chr = legendre(p), legendre(r)
    product = chp[:,None] * chr[None,:]
    interior = product != 0
    min_product = int(interior.sum())
    min_table = p*r
    balanced_rank1_centres = []
    for a,b in itertools.product(range(p),range(r)):
        f = cell[np.ix_((np.arange(p)+a)%p, (np.arange(r)+b)%r)]
        mismatch = int((f[interior] != product[interior]).sum())
        min_product = min(min_product, mismatch, int(interior.sum())-mismatch)
        error = 0
        for cp,cr in itertools.product((-1,0,1),repeat=2):
            vals = f[np.ix_(chp==cp,chr==cr)]
            plus = int((vals==1).sum())
            error += min(plus, vals.size-plus)
        min_table = min(min_table,error)
        f0 = f[1:,1:]
        u,v = f0[:,0], f0[0,:]*f0[0,0]
        if np.array_equal(f0,u[:,None]*v[None,:]) and int(u.sum())==0 and int(v.sum())==0:
            balanced_rank1_centres.append([a,b])
    return dict(min_product_interior_flips=min_product, min_quadratic_table_flips=min_table,
                balanced_rank1_centres=balanced_rank1_centres)


def coupled_group(modulus, generators):
    identity = (tuple(range(4)), 1)
    seen = {identity}
    todo = [identity]
    while todo:
        perm,u = todo.pop()
        for g,v in generators:
            nxt = (tuple(g[j] for j in perm), u*v % modulus)
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return sorted(seen)


def analyse(name, rows, p, r, generators):
    m=p*r
    cells=crt_cells(rows,p,r)
    rec=dict(name=name, q='v4', factors=[p,r],m=m,sequences=cells.reshape(4,4*m).tolist())
    paf=paf_check(rec)
    assert paf['ok'],paf
    group=coupled_group(m,generators)
    align=[]
    raw=np.array(rows,np.int64).reshape(4,4,m)
    for row in raw:
        found=[]
        for shift in range(4):
            y=row[[q^shift for q in range(4)]]
            if all(np.array_equal(y, y[np.ix_(perm,(u*np.arange(m))%m)]) for perm,u in group):
                found.append(shift)
        align.append(found)
    projected={u for perm,u in group}
    kernel={u for perm,u in group if perm==tuple(range(4))}
    distances=[[cell_distances(cells[i,q],p,r) for q in range(4)] for i in range(4)]
    return dict(name=name,p=p,r=r,m=m,exact_core_check=paf,
                compressions=raw.sum(axis=2).tolist(), primitive_energy=primitive_energy(rec),
                symmetry=dict(order=len(group),K_projection_order=len(projected),
                              K_projection_index=math.prod([p-1,r-1])//len(projected),
                              Q_trivial_kernel_order=len(kernel),valid_Q_translations_per_row=align,
                              generators=[dict(Q_permutation=list(a),K_multiplier=b) for a,b in generators]),
                cells=distances)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--kumar',required=True,help='local hadamard/3964/certificate.json')
    ap.add_argument('--out',default='logs/openai_identity/known_cores.json')
    a=ap.parse_args()
    lines=(ROOT/'data/core_v479_b3.txt').read_text().splitlines()
    seqs={k.strip():v.strip() for line in lines[1:] if '=' in line for k,v in [line.split('=',1)]}
    rows1916=[[1 if ch=='+' else -1 for ch in seqs[k]] for k in 'ABCD']
    c=json.loads(Path(a.kumar).read_text())
    assert c['group']==[2,2,247]
    rows3964=[[1 if ch=='+' else -1 for ch in s] for s in c['sequences']]
    results=[analyse('1916',rows1916,7,17,[([0,2,3,1],16),([0,1,2,3],13)]),
             analyse('3964',rows3964,13,19,[([0,3,1,2],66),([0,1,2,3],35)])]
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(dict(script='theory/openai_identity_known.py',results=results),indent=2)+'\n')
    for rec in results:
        print(json.dumps({k:rec[k] for k in ('name','compressions','symmetry')}))
        print(json.dumps(dict(name=rec['name'],product_distances=[[c['min_product_interior_flips'] for c in row] for row in rec['cells']],
                              quadratic_table_distances=[[c['min_quadratic_table_flips'] for c in row] for row in rec['cells']],
                              balanced_rank1_cells=sum(bool(c['balanced_rank1_centres']) for row in rec['cells'] for c in row))))


if __name__=='__main__':
    main()
